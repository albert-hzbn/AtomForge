#include "util/PathUtils.h"
#include "ElementData.h"
#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

void pushDirectoryHistory(std::vector<std::string>& history, int& historyIndex, const std::string& dir)
{
    if (historyIndex + 1 < (int)history.size())
        history.erase(history.begin() + historyIndex + 1, history.end());

    history.push_back(dir);
    historyIndex = (int)history.size() - 1;
}

bool hasExtension(const std::string& name, const std::string& extension)
{
    if (extension.empty())
        return true;
    if (name.size() < extension.size())
        return false;
    return name.compare(name.size() - extension.size(), extension.size(), extension) == 0;
}

std::string replaceFileExtension(const std::string& filename,
                                 const std::string& extension,
                                 const std::string& fallbackBase)
{
    const std::size_t dot = filename.find_last_of('.');
    std::string base = (dot != std::string::npos) ? filename.substr(0, dot) : filename;
    if (base.empty())
        base = fallbackBase;
    return base + extension;
}

void updateFilenameWithExtension(char* filenameBuffer,
                                 std::size_t bufferSize,
                                 const std::string& extension,
                                 const std::string& fallbackBase)
{
    const std::string currentName(filenameBuffer);
    const std::string updatedName = replaceFileExtension(currentName, extension, fallbackBase);
    std::snprintf(filenameBuffer, bufferSize, "%s", updatedName.c_str());
}

int atomicNumberFromSymbol(const std::string& symbol)
{
    if (symbol.empty())
        return -1;

    std::string lowered = symbol;
    for (char& ch : lowered)
        ch = (char)std::tolower((unsigned char)ch);

    for (int z = 1; z <= 118; ++z)
    {
        std::string candidate = elementSymbol(z);
        for (char& ch : candidate)
            ch = (char)std::tolower((unsigned char)ch);
        if (candidate == lowered)
            return z;
    }

    return -1;
}

std::string normalizePathSeparators(const std::string& path)
{
    std::string out = path;
    std::replace(out.begin(), out.end(), '\\', '/');
    return out;
}

std::string joinPath(const std::string& base, const std::string& name)
{
    if (base.empty() || base == ".")
        return name;
    if (base.back() == '/' || base.back() == '\\')
        return base + name;
    return base + "/" + name;
}

bool isDriveRootPath(const std::string& path)
{
    if (path.size() < 2)
        return false;
    if (!std::isalpha((unsigned char)path[0]) || path[1] != ':')
        return false;
    return (path.size() == 2) || (path.size() == 3 && (path[2] == '/' || path[2] == '\\'));
}

std::string parentPath(const std::string& path)
{
    if (path.empty() || path == ".")
        return ".";

    std::string out = path;
    while (out.size() > 1 && (out.back() == '/' || out.back() == '\\'))
    {
        if (isDriveRootPath(out))
            return normalizePathSeparators(out);
        out.pop_back();
    }

    if (out == "/" || out == "\\" || isDriveRootPath(out))
        return normalizePathSeparators(out);

    std::size_t pos = out.find_last_of("/\\");
    if (pos == std::string::npos)
        return ".";
    if (pos == 0)
        return out.substr(0, 1);
    if (pos == 2 && std::isalpha((unsigned char)out[0]) && out[1] == ':')
        return normalizePathSeparators(out.substr(0, 3));
    return out.substr(0, pos);
}

std::string detectHomePath()
{
#ifdef _WIN32
    if (const char* userProfile = std::getenv("USERPROFILE"))
        return normalizePathSeparators(userProfile);

    const char* homeDrive = std::getenv("HOMEDRIVE");
    const char* homePath = std::getenv("HOMEPATH");
    if (homeDrive && homePath)
        return normalizePathSeparators(std::string(homeDrive) + homePath);
#endif

    if (const char* home = std::getenv("HOME"))
        return normalizePathSeparators(home);

    return ".";
}

void appendUniquePath(std::vector<std::string>& paths, const std::string& value)
{
    if (value.empty())
        return;

    std::string normalized = normalizePathSeparators(value);
    if (isDriveRootPath(normalized) && normalized.size() == 2)
        normalized += "/";

    if (std::find(paths.begin(), paths.end(), normalized) == paths.end())
        paths.push_back(normalized);
}

std::string toLowerStr(const std::string& s)
{
    std::string out = s;
    for (auto& c : out)
        c = (char)std::tolower((unsigned char)c);
    return out;
}

#ifdef _WIN32
static std::string toNativePath(const std::string& path)
{
    std::string native = path;
    std::replace(native.begin(), native.end(), '/', '\\');
    return native;
}
#endif

bool loadDirectoryEntries(const std::string& directory,
                          bool filterFiles,
                          const std::function<bool(const std::string&)>& includeFile,
                          std::vector<DirectoryEntry>& entries)
{
#ifdef _WIN32
    std::string nativeDir = toNativePath(directory.empty() ? "." : directory);
    if (!nativeDir.empty() && nativeDir.back() != '\\' && nativeDir.back() != '/')
        nativeDir.push_back('\\');
    std::string searchPattern = nativeDir + "*";

    WIN32_FIND_DATAA findData;
    HANDLE handle = FindFirstFileA(searchPattern.c_str(), &findData);
    if (handle == INVALID_HANDLE_VALUE)
        return false;

    do
    {
        std::string name(findData.cFileName);
        if (name == "." || name == "..")
            continue;

        bool isDir = (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        if (!isDir && filterFiles && !includeFile(name))
            continue;

        entries.emplace_back(name, isDir);
    } while (FindNextFileA(handle, &findData));

    FindClose(handle);
#else
    DIR* d = opendir(directory.empty() ? "." : directory.c_str());
    if (!d)
        return false;

    struct dirent* de;
    while ((de = readdir(d)) != NULL)
    {
        std::string name(de->d_name);
        if (name == "." || name == "..")
            continue;

        std::string full = joinPath(directory, name);
        struct stat st;
        bool isDir = (stat(full.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
        if (!isDir && filterFiles && !includeFile(name))
            continue;

        entries.emplace_back(name, isDir);
    }
    closedir(d);
#endif

    std::sort(entries.begin(), entries.end(), [](const DirectoryEntry& a, const DirectoryEntry& b) {
        const auto& [nameA, isDirA] = a;
        const auto& [nameB, isDirB] = b;
        if (isDirA != isDirB)
            return isDirA > isDirB;
        return nameA < nameB;
    });

    return true;
}

void drawDirectoryEntries(const std::vector<DirectoryEntry>& entries,
                          char* selectedFilename,
                          int idBase,
                          const std::function<void(const std::string&)>& onEnterDirectory,
                          bool* fileDoubleClicked)
{
    for (size_t i = 0; i < entries.size(); ++i)
    {
        const auto& [name, isDir] = entries[i];

        ImGui::PushID(static_cast<int>(i) + idBase);
        if (isDir)
        {
            // Leave room for the hand-drawn folder icon (5 spaces ≈ icon width)
            const ImVec4 windowBg = ImGui::GetStyleColorVec4(ImGuiCol_WindowBg);
            const float bgLuma = 0.2126f * windowBg.x + 0.7152f * windowBg.y + 0.0722f * windowBg.z;
            const ImVec4 folderTextColor = (bgLuma > 0.55f)
                ? ImVec4(0.16f, 0.22f, 0.30f, 1.0f)
                : ImVec4(0.55f, 0.78f, 1.00f, 1.0f);

            ImGui::PushStyleColor(ImGuiCol_Text, folderTextColor);
            std::string label = std::string("     ") + name + "##dir";
            const bool clicked = ImGui::Selectable(label.c_str());
            ImGui::PopStyleColor();

            // Draw a two-rectangle folder icon over the reserved leading space
            const ImVec2 itemMin  = ImGui::GetItemRectMin();
            const ImVec2 itemSize = ImGui::GetItemRectSize();
            const float  th  = ImGui::GetTextLineHeight();
            const float  iw  = th * 1.10f;
            const float  x0  = itemMin.x + 2.0f;
            const float  y0  = itemMin.y + (itemSize.y - th) * 0.5f;
            ImDrawList*  dl  = ImGui::GetWindowDrawList();
            const ImU32  col = IM_COL32(100, 175, 255, 225);
            // tab nub (top-left corner)
            dl->AddRectFilled(ImVec2(x0,              y0 + th * 0.14f),
                              ImVec2(x0 + iw * 0.48f, y0 + th * 0.36f), col, 2.0f);
            // folder body
            dl->AddRectFilled(ImVec2(x0,        y0 + th * 0.28f),
                              ImVec2(x0 + iw,   y0 + th * 0.92f), col, 2.0f);

            if (clicked)
                onEnterDirectory(name);
        }
        else
        {
            std::string label = "  " + name + "##file";
            bool selected = (std::string(selectedFilename) == name);
            if (ImGui::Selectable(label.c_str(), selected))
                std::snprintf(selectedFilename, 1024, "%s", name.c_str());
            if (fileDoubleClicked && ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(0))
                *fileDoubleClicked = true;
        }
        ImGui::PopID();
    }
}
