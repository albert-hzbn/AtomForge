#include "ui/ResponsiveLayout.h"
#include "ui/PathPicker.h"
#include "util/PathUtils.h"
#include "imgui.h"
#include "ui/DialogLayout.h"

#include <cstdio>
#include <algorithm>
#include <filesystem>

void PathPicker::navigate(const std::string& path)
{
    m_directory = normalizePathSeparators(path);
    pushDirectoryHistory(m_history,m_historyIndex,m_directory);
    std::snprintf(m_path,sizeof(m_path),"%s",m_directory.c_str());
    m_error.clear();
    m_overwrite.clear();
    if (!m_save) m_filename[0] = '\0';
}

void PathPicker::open(const std::string& title, bool save, const std::string& initialPath)
{
    m_title = title;
    m_save = save;
    m_requested = true;
    m_filename[0] = '\0';
    m_error.clear();
    m_overwrite.clear();
    const auto path = std::filesystem::u8path(initialPath);
    m_extension = save ? path.extension().u8string() : "";
    if (!initialPath.empty())
    {
        std::error_code error;
        if (std::filesystem::is_directory(path,error)) navigate(initialPath);
        else
        {
            if (path.has_parent_path()) navigate(path.parent_path().u8string());
            std::snprintf(m_filename,sizeof(m_filename),"%s",path.filename().u8string().c_str());
        }
    }
    if (m_directory.empty()) navigate(detectHomePath());
}

std::optional<std::string> PathPicker::draw()
{
    std::optional<std::string> chosen;
    const std::string title = m_title + "###ElectronicPathPicker";
    if (m_requested) { ImGui::OpenPopup(title.c_str()); m_requested = false; }
    responsive::windowSize(ImVec2(800,560),ImGuiCond_Appearing);
    responsive::windowConstraints(ImVec2(640,480),ImVec2(FLT_MAX,FLT_MAX));
    bool open = true;
    if (responsive::beginModal(title.c_str(),&open,ImGuiWindowFlags_NoCollapse))
    {
        if (responsive::button("<") && m_historyIndex > 0)
        {
            m_directory = m_history[--m_historyIndex];
            std::snprintf(m_path,sizeof(m_path),"%s",m_directory.c_str());
            if (!m_save) m_filename[0] = '\0';
            m_error.clear();
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine();
        if (responsive::button(">") && m_historyIndex+1 < static_cast<int>(m_history.size()))
        {
            m_directory = m_history[++m_historyIndex];
            std::snprintf(m_path,sizeof(m_path),"%s",m_directory.c_str());
            if (!m_save) m_filename[0] = '\0';
            m_error.clear();
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine();
        if (responsive::button("Up")) navigate(parentPath(m_directory));
        ImGui::SameLine();
        ImGui::SetNextItemWidth(responsive::dp(-1));
        if (ImGui::InputText("##folder",m_path,sizeof(m_path),ImGuiInputTextFlags_EnterReturnsTrue)) navigate(m_path);
        ImGui::Separator();
        const float height = std::max(100.0f,ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing()*3 - ImGui::GetTextLineHeightWithSpacing()*2);
        responsive::beginChild("Locations",ImVec2(responsive::dp(150),height),true);
        ImGui::TextDisabled("Locations");
        if (ImGui::Selectable("Home")) navigate(detectHomePath());
        for (const char* name : {"Desktop","Documents","Downloads"})
        {
            const auto location = std::filesystem::u8path(detectHomePath())/name;
            std::error_code error;
            if (std::filesystem::is_directory(location,error) && ImGui::Selectable(name))
                navigate(location.u8string());
        }
        std::error_code currentError;
        const auto current = std::filesystem::current_path(currentError);
        if (!currentError && ImGui::Selectable("Working folder")) navigate(current.u8string());
        ImGui::Spacing();
        ImGui::TextDisabled("Drives");
        ImGui::Separator();
        for (const auto& root : getDriveRoots()) if (ImGui::Selectable(root.c_str())) navigate(root);
        ImGui::EndChild();
        ImGui::SameLine();
        responsive::beginChild("Files",ImVec2(0,height),true);
        std::vector<DirectoryEntry> entries;
        bool accepted = false;
        // Include extensionless VASP outputs and arbitrarily named calculations.
        if (loadDirectoryEntries(m_directory,true,[](const std::string&) { return true; },entries))
            drawDirectoryEntries(entries,m_filename,0,[&](const std::string& name) { navigate(joinPath(m_directory,name)); },&accepted);
        else ImGui::TextDisabled("Unable to open folder");
        ImGui::EndChild();
        ImGui::Spacing();
        ImGui::TextWrapped("%s",m_save ? "Choose a folder and filename." : "Select a file, then Open. Double-click a file to open it directly.");
        ImGui::Separator();
        ImGui::SetNextItemWidth(responsive::dp(-85));
        if (m_save)
        {
            if (ImGui::InputText("Filename",m_filename,sizeof(m_filename),ImGuiInputTextFlags_EnterReturnsTrue)) accepted = true;
        }
        else ImGui::TextWrapped("Selected: %s",m_filename[0] ? m_filename : "No file selected");
        if (!m_error.empty()) ImGui::TextWrapped("%s",m_error.c_str());
        ImGui::SetCursorPosY(std::max(ImGui::GetCursorPosY(),ImGui::GetWindowHeight()-ImGui::GetStyle().WindowPadding.y-ImGui::GetFrameHeight()));
        ImGui::SetCursorPosX(std::max(ImGui::GetStyle().WindowPadding.x,ImGui::GetWindowWidth()-ImGui::GetStyle().WindowPadding.x-200-ImGui::GetStyle().ItemSpacing.x));
        if (dialogLayout::primaryButton(m_save ? "Save" : "Open",ImVec2(100,0))) accepted = true;
        ImGui::SameLine();
        if (responsive::button("Cancel",responsive::size(100,0))) ImGui::CloseCurrentPopup();
        if (accepted)
        {
            auto path = std::filesystem::u8path(m_directory) / std::filesystem::u8path(m_filename);
            std::error_code error;
            const bool directory = std::filesystem::is_directory(path,error);
            if (m_save && !directory && !m_extension.empty()) path.replace_extension(std::filesystem::u8path(m_extension));
            if (!m_filename[0]) m_error = "Choose a filename.";
            else if (directory) navigate(path.u8string());
            else if (!m_save && !std::filesystem::is_regular_file(path,error)) m_error = "Choose an existing file.";
            else if (m_save && std::filesystem::exists(path,error))
            { m_overwrite = path.u8string(); ImGui::OpenPopup("Replace existing file?"); }
            else { chosen = path.u8string(); ImGui::CloseCurrentPopup(); }
        }
        if (responsive::beginModal("Replace existing file?",nullptr,ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::TextWrapped("Replace %s?",m_overwrite.c_str());
            if (responsive::button("Replace")) { chosen = m_overwrite; ImGui::CloseCurrentPopup(); }
            ImGui::SameLine();
            if (responsive::button("Cancel##replace")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
            if (chosen) ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    if (chosen)
    {
        // Keep subsequent browsing in the selected folder, including absolute filenames.
        navigate(std::filesystem::u8path(*chosen).parent_path().u8string());
    }
    return chosen;
}
