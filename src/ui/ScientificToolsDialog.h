#pragma once

#include "model/Structure.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include <array>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

class ScientificToolsDialog
{
public:
    void drawMenuItems(const char* category);
    void draw(const Structure& structure, const std::function<void(Structure&)>& loadResult);
private:
    struct Field {
        std::array<char, 8192> value{};
        bool enabled = true;
        std::string path;
        bool useFile = false;
        bool selectColumn = false;
        int column = 0;
        std::array<char, 256> field{};
    };
    struct Result { std::filesystem::path output; std::string report; std::filesystem::path structures; };
    void selectTool(int index);
    bool m_open = false;
    bool m_showResults = false;
    bool m_showInputs = false;
    int m_tool = -1;
    int m_pickerTarget = -1;
    std::array<char, 2048> m_python{};
    std::vector<Field> m_fields;
    PathPicker m_picker;
    atomforge::BackgroundTask<Result> m_task;
    std::string m_error;
    Result m_result;
};
