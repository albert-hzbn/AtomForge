#include "ui/ScientificToolsDialog.h"
#include "ui/ScienceCatalog.h"
#include "ui/ResponsiveLayout.h"
#include "util/ApplicationPaths.h"
#include "util/ScientificProcess.h"
#include "io/Trajectory.h"
#include "imgui.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace
{
std::string quote(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c);
        else out << static_cast<char>(c);
    }
    out << '"';
    return out.str();
}

std::string structureJson(const Structure& structure)
{
    if (structure.atoms.empty()) throw std::runtime_error("The active structure is empty");
    std::ostringstream out;
    out << std::setprecision(17) << "{\"symbols\":[";
    for (std::size_t i = 0; i < structure.atoms.size(); ++i) {
        if (i) out << ',';
        out << quote(structure.atoms[i].symbol);
    }
    out << "],\"positions\":[";
    for (std::size_t i = 0; i < structure.atoms.size(); ++i) {
        if (i) out << ',';
        const auto& atom = structure.atoms[i];
        out << '[' << atom.x << ',' << atom.y << ',' << atom.z << ']';
    }
    out << ']';
    if (structure.hasUnitCell) {
        out << ",\"cell\":[";
        for (int i = 0; i < 3; ++i) {
            if (i) out << ',';
            const auto& vector = structure.cellVectors[static_cast<std::size_t>(i)];
            out << '[' << vector[0] << ',' << vector[1] << ',' << vector[2] << ']';
        }
        out << ']';
    }
    out << '}';
    return out.str();
}

std::filesystem::path entryPoint()
{
    const auto app = applicationDirectory();
    for (const auto& root : {app / "python", app.parent_path() / "python",
                            app.parent_path() / "share" / "atomforge" / "python",
                            std::filesystem::path(ATOMFORGE_SOURCE_PYTHON)}) {
        const auto entry = root / "atomforge" / "science" / "_desktop_entry.py";
        if (std::filesystem::is_regular_file(entry)) return entry;
    }
    throw std::runtime_error("Scientific Python package is missing from this installation");
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    std::string text(65536, '\0');
    input.read(text.data(), static_cast<std::streamsize>(text.size()));
    text.resize(static_cast<std::size_t>(input.gcount()));
    return text;
}
}

void ScientificToolsDialog::selectTool(int index)
{
    m_tool = index;
    m_fields.clear();
    for (const auto& parameter : scienceToolCatalog()[static_cast<std::size_t>(index)].parameters) {
        Field field;
        std::snprintf(field.value.data(), field.value.size(), "%s", parameter.value);
        field.enabled = parameter.required || parameter.value[0] != '\0';
        field.useFile = (std::strcmp(parameter.kind, "data") == 0 ||
                         std::strcmp(parameter.kind, "structure") == 0) && parameter.value[0] == '\0';
        m_fields.push_back(field);
    }
    m_error.clear();
    m_result = {};
    m_showInputs = true;
    m_showResults = false;
}

void ScientificToolsDialog::drawMenuItems(const char* category)
{
    const auto& catalog = scienceToolCatalog();
    for (std::size_t i = 0; i < catalog.size(); ++i)
        if (std::strcmp(catalog[i].category, category) == 0 &&
            ImGui::MenuItem(catalog[i].title, nullptr, m_open && m_tool == static_cast<int>(i), !m_task.running())) {
            selectTool(static_cast<int>(i));
            m_open = true;
        }
}

void ScientificToolsDialog::draw(const Structure& structure, const std::function<void(Structure&)>& loadResult)
{
    if (m_task.poll()) {
        m_error = m_task.error();
        if (m_task.result()) {
            m_result = std::move(*m_task.result());
            m_showResults = true;
        }
        m_task.clearResult();
    }
    if (!m_open) return;
    if (m_tool < 0) selectTool(0);
    if (!m_python[0]) {
#ifdef _WIN32
        std::snprintf(m_python.data(), m_python.size(), "python");
#else
        std::snprintf(m_python.data(), m_python.size(), "python3");
#endif
    }
    const auto& catalog = scienceToolCatalog();
    const auto& tool = catalog[static_cast<std::size_t>(m_tool)];
    const std::string title = std::string(tool.title) + "###Scientific analysis";
    responsive::windowSize(ImVec2(760, 700), ImGuiCond_FirstUseEver);
    responsive::windowConstraints(ImVec2(420, 320), ImVec2(1400, 1200));
    if (responsive::begin(title.c_str(), &m_open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, responsive::size(10, 6));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, responsive::size(10, 8));
        ImGui::TextDisabled("Analysis / %s", tool.category);
        const std::string summary = std::string(tool.help).substr(0, std::string(tool.help).find('\n'));
        ImGui::TextWrapped("%s", summary.c_str());
        if (m_task.running()) {
            ImGui::TextUnformatted("Calculation running...");
            ImGui::SameLine();
            if (responsive::button("Cancel")) m_task.cancel();
        }
        if (!m_error.empty() && ImGui::TreeNodeEx("Calculation error", ImGuiTreeNodeFlags_DefaultOpen)) {
            responsive::beginChild("Error details", ImVec2(0, responsive::dp(90)), true);
            ImGui::TextWrapped("%s", m_error.c_str());
            ImGui::EndChild();
            ImGui::TreePop();
        }
        ImGui::Spacing();
        if (ImGui::BeginTabBar("Scientific workflow")) {
        if (ImGui::BeginTabItem("Inputs", nullptr, m_showInputs ? ImGuiTabItemFlags_SetSelected : 0)) {
        m_showInputs = false;
        ImGui::BeginDisabled(m_task.running());
        const float parameterHeight = std::max(responsive::dp(80), ImGui::GetContentRegionAvail().y -
            2 * ImGui::GetFrameHeightWithSpacing());
        responsive::beginChild("Scientific parameters", ImVec2(0, parameterHeight), true);
        const char* sections[] = {"Input data", "Calculation settings", "Optional inputs"};
        for (int section = 0; section < 3; ++section) {
        bool headingShown = false;
        for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
            const auto& parameter = tool.parameters[i];
            auto& field = m_fields[i];
            const std::string kind = parameter.kind;
            const bool optional = !parameter.required && parameter.value[0] == '\0';
            const int fieldSection = optional ? 2 :
                ((kind == "data" || kind == "structure" || kind == "calculator") ? 0 : 1);
            if (fieldSection != section) continue;
            if (!headingShown) {
                ImGui::SeparatorText(sections[section]);
                headingShown = true;
            }
            ImGui::PushID(static_cast<int>(i));
            if (optional) ImGui::Checkbox(parameter.label, &field.enabled);
            else if (kind != "bool") ImGui::TextWrapped("%s", parameter.label);
            if (optional && !field.enabled) {
                ImGui::Spacing();
                ImGui::PopID();
                continue;
            }
            ImGui::BeginDisabled(!field.enabled);
            if (kind == "bool") {
                bool value = std::string(field.value.data()) == "true";
                if (ImGui::Checkbox(parameter.label, &value)) std::snprintf(field.value.data(), field.value.size(), "%s", value ? "true" : "false");
            } else if (kind == "axes") {
                std::istringstream input(field.value.data());
                char bracket;
                input >> bracket;
                bool values[3]{};
                for (int axis = 0; axis < 3; ++axis) {
                    std::string token;
                    std::getline(input, token, axis == 2 ? ']' : ',');
                    values[axis] = token.find("true") != std::string::npos;
                }
                const char* labels[] = {"a", "b", "c"};
                bool changed = false;
                for (int axis = 0; axis < 3; ++axis) {
                    if (axis) ImGui::SameLine();
                    changed |= ImGui::Checkbox(labels[axis], &values[axis]);
                }
                if (changed) std::snprintf(field.value.data(), field.value.size(), "[%s,%s,%s]",
                    values[0] ? "true" : "false", values[1] ? "true" : "false", values[2] ? "true" : "false");
            } else if (kind == "range" || kind == "vector") {
                const int count = kind == "range" ? 2 : 3;
                double values[3]{};
                std::string text = field.value.data();
                for (char& c : text) if (c == '[' || c == ']' || c == ',') c = ' ';
                std::istringstream input(text);
                for (int component = 0; component < count; ++component) input >> values[component];
                ImGui::SetNextItemWidth(-1);
                if (ImGui::InputScalarN("##components", ImGuiDataType_Double, values, count, nullptr, nullptr, "%.8g")) {
                    std::ostringstream output;
                    output << std::setprecision(17) << '[';
                    for (int component = 0; component < count; ++component) {
                        if (component) output << ',';
                        output << values[component];
                    }
                    output << ']';
                    std::snprintf(field.value.data(), field.value.size(), "%s", output.str().c_str());
                }
                ImGui::TextDisabled(kind == "range" ? "Start / End" : "kx / ky / kz");
            } else if (kind == "window") {
                int selected = std::string(field.value.data()) == "\"hann\"" ? 0 : 1;
                ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340)));
                if (ImGui::Combo("##window", &selected, "Hann\0Rectangular\0"))
                    std::snprintf(field.value.data(), field.value.size(), "%s", selected == 0 ? "\"hann\"" : "\"none\"");
            } else if (kind == "float" || kind == "int") {
                ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340)));
                ImGui::InputText("##number", field.value.data(), field.value.size(), ImGuiInputTextFlags_CharsScientific);
            } else {
                if (kind == "data" || kind == "structure") {
                    int source = field.useFile ? 0 : 1;
                    ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(170)));
                    if (ImGui::Combo("##source", &source, "Load from file\0Enter values\0")) field.useFile = source == 0;
                    if (field.useFile) {
                        if (ImGui::GetContentRegionAvail().x > responsive::dp(320)) ImGui::SameLine();
                        if (responsive::button("Browse...")) { m_pickerTarget = static_cast<int>(i); m_picker.open(parameter.label, false, field.path); }
                        ImGui::TextWrapped("%s", field.path.empty() ? "Choose a file" : field.path.c_str());
                        if (kind == "data" && ImGui::TreeNode("File options")) {
                            ImGui::SetNextItemWidth(-1);
                            ImGui::InputTextWithHint("##field", "JSON field, e.g. result.msd_A2 (optional)", field.field.data(), field.field.size());
                            ImGui::Checkbox("Select numeric column", &field.selectColumn);
                            if (field.selectColumn) { ImGui::SetNextItemWidth(-1); ImGui::InputInt("##column", &field.column); ImGui::TextDisabled("Column index starts at zero"); }
                            ImGui::TreePop();
                        }
                    }
                    if (kind == "structure" && !structure.atoms.empty() && responsive::button("Use active structure")) {
                        try {
                        // Large structures go through a file, avoiding a fixed-size text buffer.
                        const auto root = std::filesystem::temp_directory_path() / "AtomForge-science";
                        std::filesystem::create_directories(root);
                        const auto file = root / ("structure-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
                        std::ofstream output(file); output << structureJson(structure);
                        if (!output) throw std::runtime_error("Cannot save active structure for analysis");
                        field.path = file.u8string(); field.useFile = true;
                        } catch (const std::exception& error) { m_error = error.what(); }
                    }
                }
                if (!field.useFile || kind == "calculator" || kind == "string") {
                    ImGui::SetNextItemWidth(-1);
                    ImGui::InputText("##value", field.value.data(), field.value.size());
                    if (kind == "calculator") ImGui::TextWrapped("Explicit ASE calculator configuration. EMT is an example for supported metals; select a calculator appropriate to your system.");
                    else if (kind == "data") ImGui::TextDisabled("Inline JSON array, or browse numeric CSV/NPY/JSON/structure data.");
                }
            }
            ImGui::EndDisabled();
            ImGui::Spacing();
            ImGui::PopID();
        }
        }
        ImGui::EndChild();
        if (responsive::button("Run calculation")) {
            try {
                std::ostringstream request; request << '{';
                bool first = true;
                for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
                    const auto& field = m_fields[i];
                    if (!field.enabled) continue;
                    const auto& parameter = tool.parameters[i];
                    std::string value = field.value.data();
                    if (field.useFile) {
                        if (field.path.empty()) throw std::runtime_error("Choose a file for " + std::string(parameter.name));
                        {
                            value = "{\"file\":" + quote(field.path);
                            if (field.field[0]) value += ",\"field\":" + quote(field.field.data());
                            if (field.selectColumn) value += ",\"column\":" + std::to_string(field.column);
                            value += '}';
                        }
                    }
                    if (value.empty()) throw std::runtime_error("Provide " + std::string(parameter.name));
                    if (!first) request << ',';
                    first = false;
                    request << quote(parameter.name) << ':' << value;
                }
                request << '}';
                const auto entry = entryPoint();
                const auto directory = std::filesystem::temp_directory_path() / "AtomForge-science" /
                    ("run-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
                std::filesystem::create_directories(directory);
                { std::ofstream output(directory / "request.json"); output << request.str(); if (!output) throw std::runtime_error("Cannot save calculation parameters"); }
                const std::string python = m_python.data(), id = tool.id;
                m_error.clear();
                m_result = {};
                m_task.start([directory, entry, python, id] {
                    const auto output = directory / "result.json", report = directory / "report.txt", log = directory / "process.log";
                    const auto structures = directory / "frames.extxyz";
                    const int code = runScientificProcess({python, entry.u8string(), id, "--input", (directory / "request.json").u8string(),
                        "--output", output.u8string(), "--report", report.u8string(), "--structures", structures.u8string()}, log);
                    if (code != 0) throw std::runtime_error("Calculation failed (" + std::to_string(code) + ").\n" + readText(log));
                    return Result{output, readText(report), std::filesystem::exists(structures) ? structures : std::filesystem::path{}};
                });
            } catch (const std::exception& error) { m_error = error.what(); }
        }
        ImGui::EndDisabled();
        ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Results", nullptr, m_showResults ? ImGuiTabItemFlags_SetSelected : 0)) {
        m_showResults = false;
        if (!m_result.output.empty()) {
            if (responsive::button("Save results...")) { m_pickerTarget = -3; m_picker.open("Save scientific results", true, "analysis.json"); }
            if (!m_result.structures.empty()) {
                if (ImGui::GetContentRegionAvail().x > responsive::dp(520)) ImGui::SameLine();
                if (responsive::button("Save structures...")) { m_pickerTarget = -4; m_picker.open("Save result structures or trajectory", true, "frames.extxyz"); }
                if (responsive::button("Open final structure in a new tab")) {
                    try {
                        auto frames = loadXyzTrajectory(m_result.structures.u8string());
                        if (!frames.empty()) loadResult(frames.back());
                    } catch (const std::exception& error) { m_error = error.what(); }
                }
            }
            responsive::beginChild("Scientific results", ImVec2(0, 0), true);
            ImGui::TextWrapped("%s", m_result.report.c_str());
            ImGui::EndChild();
        } else {
            ImGui::TextWrapped(m_task.running() ? "Calculation running. Results will appear here when it finishes." :
                "No results yet. Load your data on the Inputs tab and run the calculation.");
        }
        ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Method")) {
            responsive::beginChild("Method details", ImVec2(0, 0));
            ImGui::TextWrapped("%s", tool.help);
            ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Environment")) {
            ImGui::BeginDisabled(m_task.running());
            ImGui::TextWrapped("Choose the Python environment used by all scientific analysis tools. It must have the atomforge-py science dependencies installed.");
            ImGui::Spacing();
            ImGui::TextUnformatted("Python interpreter");
            ImGui::SetNextItemWidth(-1);
            ImGui::InputText("##python", m_python.data(), m_python.size());
            if (responsive::button("Browse interpreter...")) { m_pickerTarget = -2; m_picker.open("Select Python interpreter", false, m_python.data()); }
            ImGui::Spacing();
            ImGui::TextWrapped("Calculations run in a separate process. You can continue working with your structure while an analysis runs.");
            ImGui::EndDisabled();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        }
        if (auto path = m_picker.draw()) {
            try {
                if (m_pickerTarget == -2) std::snprintf(m_python.data(), m_python.size(), "%s", path->c_str());
                else if (m_pickerTarget == -3) std::filesystem::copy_file(m_result.output, std::filesystem::u8path(*path), std::filesystem::copy_options::overwrite_existing);
                else if (m_pickerTarget == -4) std::filesystem::copy_file(m_result.structures, std::filesystem::u8path(*path), std::filesystem::copy_options::overwrite_existing);
                else if (m_pickerTarget >= 0 && static_cast<std::size_t>(m_pickerTarget) < m_fields.size()) {
                    auto& field = m_fields[static_cast<std::size_t>(m_pickerTarget)]; field.path = *path; field.useFile = true;
                }
            } catch (const std::exception& error) { m_error = error.what(); }
        }
        ImGui::PopStyleVar(2);
    }
    ImGui::End();
}
