#pragma once

#include <optional>
#include <string>
#include <vector>

// Modal adapter around the directory-list and navigation helpers shared with
// FileBrowser's structure Open/Save dialogs. State belongs to each caller.
class PathPicker
{
public:
    void open(const std::string& title, bool save, const std::string& initialPath);
    std::optional<std::string> draw();

private:
    void navigate(const std::string& path);
    std::string m_title;
    std::string m_directory;
    std::string m_error;
    std::string m_overwrite;
    std::string m_extension;
    std::vector<std::string> m_history;
    int m_historyIndex = -1;
    char m_path[2048]{};
    char m_filename[1024]{};
    bool m_requested = false;
    bool m_save = false;
};
