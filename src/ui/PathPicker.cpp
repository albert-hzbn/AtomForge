#include "ui/PathPicker.h"
#include "util/PathUtils.h"
#include "imgui.h"

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
    ImGui::SetNextWindowSize(ImVec2(720,500),ImGuiCond_Appearing);
    bool open = true;
    if (ImGui::BeginPopupModal(title.c_str(),&open))
    {
        if (ImGui::Button("<") && m_historyIndex > 0)
        {
            m_directory = m_history[--m_historyIndex];
            std::snprintf(m_path,sizeof(m_path),"%s",m_directory.c_str());
        }
        ImGui::SameLine();
        if (ImGui::Button(">") && m_historyIndex+1 < static_cast<int>(m_history.size()))
        {
            m_directory = m_history[++m_historyIndex];
            std::snprintf(m_path,sizeof(m_path),"%s",m_directory.c_str());
        }
        ImGui::SameLine();
        if (ImGui::Button("Up")) navigate(parentPath(m_directory));
        ImGui::SameLine();
        ImGui::SetNextItemWidth(-1);
        if (ImGui::InputText("##folder",m_path,sizeof(m_path),ImGuiInputTextFlags_EnterReturnsTrue)) navigate(m_path);
        const float height = std::max(100.0f,ImGui::GetContentRegionAvail().y - 115);
        ImGui::BeginChild("Locations",ImVec2(125,height),true);
        ImGui::TextDisabled("Locations");
        if (ImGui::Selectable("Home")) navigate(detectHomePath());
        ImGui::Separator();
        for (const auto& root : getDriveRoots()) if (ImGui::Selectable(root.c_str())) navigate(root);
        ImGui::EndChild();
        ImGui::SameLine();
        ImGui::BeginChild("Files",ImVec2(0,height),true);
        std::vector<DirectoryEntry> entries;
        bool accepted = false;
        // Include extensionless VASP outputs and arbitrarily named calculations.
        if (loadDirectoryEntries(m_directory,true,[](const std::string&) { return true; },entries))
            drawDirectoryEntries(entries,m_filename,0,[&](const std::string& name) { navigate(joinPath(m_directory,name)); },&accepted);
        else ImGui::TextDisabled("Unable to open folder");
        ImGui::EndChild();
        ImGui::SetNextItemWidth(-85);
        if (ImGui::InputText("Filename",m_filename,sizeof(m_filename),ImGuiInputTextFlags_EnterReturnsTrue)) accepted = true;
        if (!m_error.empty()) ImGui::TextWrapped("%s",m_error.c_str());
        if (ImGui::Button(m_save ? "Save" : "Open",ImVec2(100,0))) accepted = true;
        ImGui::SameLine();
        if (ImGui::Button("Cancel",ImVec2(100,0))) ImGui::CloseCurrentPopup();
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
        if (ImGui::BeginPopupModal("Replace existing file?",nullptr,ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::TextWrapped("Replace %s?",m_overwrite.c_str());
            if (ImGui::Button("Replace")) { chosen = m_overwrite; ImGui::CloseCurrentPopup(); }
            ImGui::SameLine();
            if (ImGui::Button("Cancel##replace")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
            if (chosen) ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    return chosen;
}
