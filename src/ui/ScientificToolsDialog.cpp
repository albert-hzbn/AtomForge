#include "ui/ScientificToolsDialog.h"
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "ui/ResponsiveLayout.h"
#include "io/StructureLoader.h"
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
using atomforge::science::Json;

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream output(path, std::ios::binary);
    output << text;
    if (!output) throw std::runtime_error("Cannot save " + path.u8string());
}

Json potentialJson(int potential, double epsilon, double sigma, double cutoff)
{
    Json result = Json::object();
    if (potential == 0) {
        result["potential"] = "EMT";
        return result;
    }
    result["potential"] = "LennardJones";
    result["epsilon"] = epsilon;
    result["sigma"] = sigma;
    result["cutoff"] = cutoff;
    return result;
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

bool ScientificToolsDialog::open(const std::string& toolId)
{
    const auto& catalog = scienceToolCatalog();
    for (std::size_t i = 0; i < catalog.size(); ++i)
        if (toolId == catalog[i].id && !m_task.running()) {
            selectTool(static_cast<int>(i));
            m_open = true;
            return true;
        }
    return false;
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
            const float progress = m_task.progress();
            if (progress >= 0.0f) {
                ImGui::ProgressBar(progress, ImVec2(std::max(responsive::dp(120), ImGui::GetContentRegionAvail().x - responsive::dp(110)), 0), "Calculating...");
            } else ImGui::TextUnformatted("Calculation running...");
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
                        std::ofstream output(file); output << atomforge::science::structureJson(structure).dump();
                        if (!output) throw std::runtime_error("Cannot save active structure for analysis");
                        field.path = file.u8string(); field.useFile = true;
                        } catch (const std::exception& error) { m_error = error.what(); }
                    }
                }
                if (kind == "calculator") {
                    ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340)));
                    ImGui::Combo("##potential", &field.potential, "EMT (effective-medium theory)\0Lennard-Jones 12-6\0");
                    if (field.potential == 0) {
                        ImGui::TextWrapped("Built-in EMT for Al, Cu, Ag, Au, Ni, Pd and Pt (H, C, N and O are illustrative only).");
                    } else {
                        const float width = std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340));
                        ImGui::SetNextItemWidth(width); ImGui::InputDouble("Well depth epsilon (eV)", &field.epsilon, 0, 0, "%.6g");
                        ImGui::SetNextItemWidth(width); ImGui::InputDouble("Zero-crossing sigma (Angstrom)", &field.sigma, 0, 0, "%.6g");
                        ImGui::SetNextItemWidth(width); ImGui::InputDouble("Cutoff (Angstrom)", &field.cutoff, 0, 0, "%.6g");
                        ImGui::TextDisabled("Energy is shifted to zero at the cutoff.");
                    }
                } else if (!field.useFile || kind == "string") {
                    ImGui::SetNextItemWidth(-1);
                    ImGui::InputText("##value", field.value.data(), field.value.size());
                    if (kind == "data") ImGui::TextDisabled("Inline JSON array, or browse numeric CSV/NPY/JSON/structure data.");
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
                Json request = Json::object();
                for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
                    const auto& field = m_fields[i];
                    if (!field.enabled) continue;
                    const auto& parameter = tool.parameters[i];
                    const std::string kind = parameter.kind;
                    if (kind == "calculator") {
                        request[parameter.name] = potentialJson(field.potential, field.epsilon, field.sigma, field.cutoff);
                    } else if (field.useFile) {
                        if (field.path.empty()) throw std::runtime_error("Choose a file for " + std::string(parameter.label));
                        Json reference = Json::object();
                        reference["file"] = field.path;
                        if (kind == "data" && field.field[0]) reference["field"] = std::string(field.field.data());
                        if (kind == "data" && field.selectColumn) reference["column"] = field.column;
                        request[parameter.name] = reference;
                    } else {
                        const std::string value = field.value.data();
                        if (value.find_first_not_of(" \t\r\n") == std::string::npos) throw std::runtime_error("Provide " + std::string(parameter.label));
                        try { request[parameter.name] = Json::parse(value); }
                        catch (const std::exception& error) { throw std::runtime_error(std::string(parameter.label) + ": " + error.what()); }
                    }
                }
                const auto directory = std::filesystem::temp_directory_path() / "AtomForge-science" /
                    ("run-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
                std::filesystem::create_directories(directory);
                const std::string id = tool.id;
                m_error.clear();
                m_result = {};
                m_task.start([directory, request, id] {
                    const auto reader = [](const std::filesystem::path& path) {
                        Structure loaded;
                        std::string error;
                        if (!loadStructureFromFile(path.u8string(), loaded, error)) throw std::runtime_error(error);
                        return loaded;
                    };
                    const auto output = atomforge::science::runTool(id, request, directory, reader);
                    const auto resultPath = directory / "result.json", structures = directory / "frames.extxyz";
                    writeText(resultPath, atomforge::science::resultDocument(id, request, output).dump(2) + "\n");
                    if (!output.frames.empty()) atomforge::science::writeExtxyz(structures, output.frames, output.velocities, output.times);
                    return Result{resultPath, atomforge::science::resultReport(id, output.result),
                                  output.frames.empty() ? std::filesystem::path{} : structures};
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
                        auto frames = atomforge::science::readFrames(m_result.structures);
                        if (!frames.empty()) loadResult(frames.back().structure);
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
        ImGui::EndTabBar();
        }
        if (auto path = m_picker.draw()) {
            try {
                if (m_pickerTarget == -3) std::filesystem::copy_file(m_result.output, std::filesystem::u8path(*path), std::filesystem::copy_options::overwrite_existing);
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
