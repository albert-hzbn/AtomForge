#include "ui/ScientificToolsDialog.h"
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "ui/ResponsiveLayout.h"
#include "ui/SciencePlot.h"
#include "io/StructureLoader.h"
#include "io/Trajectory.h"
#include "imgui.h"
#include <algorithm>
#include <cctype>
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

std::string readText(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read " + path.u8string());
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

std::filesystem::path newRunDirectory()
{
    const auto directory = std::filesystem::temp_directory_path() / "AtomForge-science" /
        ("run-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    return directory;
}

std::vector<std::pair<std::string, std::string>> generatedFiles(const Json& result)
{
    std::vector<std::pair<std::string, std::string>> files;
    if (const Json* generated = result.find("files"); generated && generated->isObject())
        for (const auto& [name, text] : generated->members()) files.push_back({name, text.string()});
    return files;
}

// Large structures go through a file, avoiding a fixed-size text buffer.
std::string saveActiveStructure(const Structure& structure)
{
    const auto root = std::filesystem::temp_directory_path() / "AtomForge-science";
    std::filesystem::create_directories(root);
    const auto file = root / ("structure-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
    std::ofstream output(file);
    output << atomforge::science::structureJson(structure).dump();
    if (!output) throw std::runtime_error("Cannot save active structure for analysis");
    return file.u8string();
}

std::string cellJson(const Structure& structure)
{
    Json cell = Json::array();
    for (const auto& row : structure.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
    return cell.dump();
}

Json potentialJson(int potential, double epsilon, double sigma, double cutoff, const std::string& file, int format)
{
    Json result = Json::object();
    if (potential == 0) {
        result["potential"] = "EMT";
        return result;
    }
    if (potential == 2) {
        static const char* formats[] = {"auto", "setfl", "fs", "funcfl"};
        if (file.empty()) throw std::runtime_error("Choose an EAM potential file");
        result["potential"] = "EAM";
        result["file"] = file;
        result["format"] = formats[std::clamp(format, 0, 3)];
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
        field.useFile = ((std::strcmp(parameter.kind, "data") == 0 ||
                          std::strcmp(parameter.kind, "structure") == 0) && parameter.value[0] == '\0') ||
                        std::strcmp(parameter.kind, "file") == 0;
        m_fields.push_back(field);
    }
    m_error.clear();
    m_result = {};
}

std::string ScientificToolsDialog::snapshot() const
{
    Json state = Json::object();
    state["open"] = m_open;
    if (m_tool < 0) return state.dump();
    state["tool"] = scienceToolCatalog()[static_cast<std::size_t>(m_tool)].id;
    Json fields = Json::array();
    for (const auto& field : m_fields) {
        Json item = Json::object();
        item["value"] = std::string(field.value.data());
        item["enabled"] = field.enabled;
        item["path"] = field.path;
        item["use_file"] = field.useFile;
        item["select_column"] = field.selectColumn;
        item["column"] = field.column;
        item["field"] = std::string(field.field.data());
        item["potential"] = field.potential;
        item["epsilon"] = field.epsilon;
        item["sigma"] = field.sigma;
        item["cutoff"] = field.cutoff;
        item["potential_file"] = field.potentialFile;
        item["eam_format"] = field.eamFormat;
        fields.push(item);
    }
    state["fields"] = fields;
    // Results live in temporary files; they are embedded so the project is self-contained.
    try {
        if (!m_result.output.empty()) {
            Json document = Json::parse(readText(m_result.output));
            Json frames;
            if (!m_result.structures.empty()) frames = readText(m_result.structures);
            state["result"] = document;
            state["frames"] = frames;
            state["property"] = m_property;
        }
    } catch (const std::exception&) {
        // Temporary result files removed outside AtomForge: save the inputs only.
    }
    return state.dump();
}

void ScientificToolsDialog::restore(const std::string& text)
{
    if (text.empty() || m_task.running()) return;
    try {
        const Json state = Json::parse(text);
        const Json* tool = state.find("tool");
        const auto& catalog = scienceToolCatalog();
        int index = -1;
        for (std::size_t i = 0; tool && tool->isString() && i < catalog.size(); ++i)
            if (tool->string() == catalog[i].id) index = static_cast<int>(i);
        if (index < 0) { m_open = false; m_tool = -1; m_result = {}; return; }
        selectTool(index);
        const auto copyText = [](auto& buffer, const Json* value) {
            if (value && value->isString()) std::snprintf(buffer.data(), buffer.size(), "%s", value->string().c_str());
        };
        const auto number = [](const Json& item, const char* key, double fallback) {
            const Json* value = item.find(key);
            return value && value->isNumber() ? value->number() : fallback;
        };
        const auto flag = [](const Json& item, const char* key, bool fallback) {
            const Json* value = item.find(key);
            return value && value->isBool() ? value->boolean() : fallback;
        };
        const auto string = [](const Json& item, const char* key) {
            const Json* value = item.find(key);
            return value && value->isString() ? value->string() : std::string();
        };
        // Inputs are restored only when the tool still has the same parameters.
        if (const Json* fields = state.find("fields"); fields && fields->isArray() && fields->size() == m_fields.size()) {
            for (std::size_t i = 0; i < m_fields.size(); ++i) {
                const Json& item = fields->items()[i];
                if (!item.isObject()) continue;
                auto& field = m_fields[i];
                copyText(field.value, item.find("value"));
                field.enabled = flag(item, "enabled", field.enabled);
                field.path = string(item, "path");
                field.useFile = flag(item, "use_file", field.useFile);
                field.selectColumn = flag(item, "select_column", false);
                field.column = static_cast<int>(number(item, "column", 0));
                copyText(field.field, item.find("field"));
                field.potential = static_cast<int>(number(item, "potential", 0));
                field.epsilon = number(item, "epsilon", field.epsilon);
                field.sigma = number(item, "sigma", field.sigma);
                field.cutoff = number(item, "cutoff", field.cutoff);
                field.potentialFile = string(item, "potential_file");
                field.eamFormat = static_cast<int>(number(item, "eam_format", 0));
            }
        }
        if (const Json* document = state.find("result"); document && document->isObject() && document->contains("result")) {
            const auto directory = newRunDirectory();
            const auto resultPath = directory / "result.json";
            writeText(resultPath, document->dump(2) + "\n");
            std::filesystem::path structures;
            if (const Json* frames = state.find("frames"); frames && frames->isString() && !frames->string().empty()) {
                structures = directory / "frames.extxyz";
                writeText(structures, frames->string());
            }
            const std::string id = catalog[static_cast<std::size_t>(index)].id;
            const Json& result = document->at("result");
            m_result = Result{resultPath, atomforge::science::resultReport(id, result), structures,
                              atomforge::science::resultPlots(id, result),
                              atomforge::science::perAtomProperties(id, result),
                              generatedFiles(result), summaryOf(result), matricesOf(result)};
            m_property = static_cast<int>(number(state, "property", 0));
        }
        m_open = flag(state, "open", false);
    } catch (const std::exception& error) {
        m_error = std::string("Saved scientific tool state could not be restored: ") + error.what();
    }
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

const char* ScientificToolsDialog::companionCell(const std::string& parameter)
{
    if (parameter == "positions" || parameter == "reference_sites") return "cell";
    if (parameter == "reference") return "reference_cell";
    if (parameter == "current") return "current_cell";
    return nullptr;
}

bool ScientificToolsDialog::acceptsTrajectory(const std::string& tool, const std::string& parameter)
{
    if (tool == "msd" || tool == "structure-factor") return parameter == "positions";
    if (tool == "vacf" || tool == "vibrational-spectrum") return parameter == "velocities";
    return false;
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

ScientificToolsDialog::Layout ScientificToolsDialog::layoutFor(const std::string& tool)
{
    static const char* atoms[] = {"structure-type", "dislocation-lines", "cluster-analysis", "void-analysis",
                                  "centrosymmetry", "bond-order", "wigner-seitz", "local-strain"};
    static const char* simulations[] = {"relax", "nvt", "npt", "neb", "phonons"};
    static const char* generators[] = {"dft-inputs", "lammps-export"};
    static const char* trajectories[] = {"msd", "diffusion", "vacf", "vibrational-spectrum", "structure-factor"};
    static const char* properties[] = {"band-gap", "effective-mass", "work-function", "equation-of-state",
                                       "elastic-tensor", "harmonic-thermodynamics"};
    for (const char* id : atoms) if (tool == id) return Layout::Atoms;
    for (const char* id : trajectories) if (tool == id) return Layout::Trajectory;
    for (const char* id : properties) if (tool == id) return Layout::Properties;
    for (const char* id : simulations) if (tool == id) return Layout::Simulation;
    for (const char* id : generators) if (tool == id) return Layout::Generator;
    return Layout::Plot;
}

const char* ScientificToolsDialog::windowId(Layout layout)
{
    switch (layout) {
    case Layout::Atoms: return "###Atom analysis";
    case Layout::Simulation: return "###Simulation run";
    case Layout::Generator: return "###Input generator";
    case Layout::Trajectory: return "###Trajectory analysis";
    case Layout::Properties: return "###Property calculator";
    case Layout::Plot: break;
    }
    return "###Plot analysis";
}

std::string ScientificToolsDialog::prettyLabel(const std::string& key)
{
    // Trailing unit tokens of a result field name become a suffix: max_force_eV_per_A -> Max force (eV/A).
    static const std::pair<const char*, const char*> units[] = {
        {"eV", "eV"}, {"A", "A"}, {"A2", "A^2"}, {"A3", "A^3"}, {"K", "K"}, {"GPa", "GPa"}, {"fs", "fs"},
        {"ps", "ps"}, {"deg", "deg"}, {"THz", "THz"}, {"amu", "amu"}, {"kV", "kV"}, {"mrad", "mrad"},
        {"atom", "atom"}, {"inv", "1/"}, {"cm", "cm"}, {"meV", "meV"}};
    std::vector<std::string> tokens;
    std::string token;
    for (char c : key) {
        if (c == '_') { if (!token.empty()) tokens.push_back(token); token.clear(); }
        else token += c;
    }
    if (!token.empty()) tokens.push_back(token);
    const auto unitOf = [&](const std::string& word) -> const char* {
        for (const auto& [suffix, shown] : units) if (word == suffix) return shown;
        return nullptr;
    };
    // The unit is the longest tail of unit words joined by "per"; at least one word stays in the label.
    std::size_t first = tokens.size();
    while (first > 1 && (unitOf(tokens[first - 1]) || tokens[first - 1] == "per")) --first;
    while (first < tokens.size() && tokens[first] == "per") ++first;
    std::string label, unit;
    for (std::size_t k = 0; k < tokens.size(); ++k) {
        if (k < first) label += (label.empty() ? "" : " ") + tokens[k];
        else unit += tokens[k] == "per" ? "/" : unitOf(tokens[k]);
    }
    if (!label.empty()) label[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(label[0])));
    return unit.empty() ? label : label + " (" + unit + ")";
}

std::vector<std::pair<std::string, std::string>> ScientificToolsDialog::summaryOf(const Json& result)
{
    std::vector<std::pair<std::string, std::string>> rows;
    if (!result.isObject()) return rows;
    for (const auto& [key, value] : result.members()) {
        std::string text;
        if (value.isNumber()) {
            char buffer[64];
            std::snprintf(buffer, sizeof(buffer), "%.6g", value.number());
            text = buffer;
        } else if (value.isBool()) text = value.boolean() ? "yes" : "no";
        else if (value.isString() && value.string().size() <= 80 && value.string().find('\n') == std::string::npos) text = value.string();
        else continue;
        rows.push_back({prettyLabel(key), text});
    }
    return rows;
}

std::vector<std::pair<std::string, std::vector<std::vector<double>>>> ScientificToolsDialog::matricesOf(const Json& result)
{
    // Small square or rectangular numeric tables (tensors, cells), shown as grids.
    std::vector<std::pair<std::string, std::vector<std::vector<double>>>> matrices;
    if (!result.isObject()) return matrices;
    for (const auto& [key, value] : result.members()) {
        if (!value.isArray() || value.size() < 2 || value.size() > 9) continue;
        std::vector<std::vector<double>> rows;
        bool numeric = true;
        for (const auto& row : value.items()) {
            if (!row.isArray() || row.size() < 2 || row.size() > 9 || (!rows.empty() && row.size() != rows.front().size())) { numeric = false; break; }
            std::vector<double> numbers;
            for (const auto& item : row.items()) {
                if (!item.isNumber()) { numeric = false; break; }
                numbers.push_back(item.number());
            }
            if (!numeric) break;
            rows.push_back(numbers);
        }
        if (!numeric) continue;
        // Round-off below 1e-9 of the largest entry is shown as zero.
        double largest = 0;
        for (const auto& row : rows) for (double v : row) largest = std::max(largest, std::abs(v));
        for (auto& row : rows) for (double& v : row) if (std::abs(v) <= 1e-9 * largest) v = 0;
        matrices.push_back({prettyLabel(key), rows});
    }
    return matrices;
}

bool ScientificToolsDialog::isSource(std::size_t index) const
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    const auto& parameter = tool.parameters[index];
    const std::string kind = parameter.kind;
    return kind == "data" || kind == "structure" || kind == "file";
}

bool ScientificToolsDialog::suppliedByActive(std::size_t index) const
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    const int input = activeInput();
    if (!m_useActive || input < 0 || layoutFor(tool.id) != Layout::Atoms) return false;
    if (static_cast<int>(index) == input) return true;
    const char* companion = companionCell(tool.parameters[static_cast<std::size_t>(input)].name);
    const std::string name = tool.parameters[index].name;
    return companion && (name == companion || name == "pbc");
}

// First structure or position input: filled from the active structure in the atom layout.
int ScientificToolsDialog::activeInput() const
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        const std::string kind = tool.parameters[i].kind, name = tool.parameters[i].name;
        if (kind == "structure" || (kind == "data" && name == "positions")) return static_cast<int>(i);
    }
    return -1;
}

void ScientificToolsDialog::useActiveStructure(std::size_t index, const Structure& structure)
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    auto& field = m_fields[index];
    field.path = saveActiveStructure(structure);
    field.useFile = true;
    field.enabled = true;
    field.field[0] = '\0';
    field.selectColumn = false;
    // Periodic local analyses also need the matching cell and axes.
    const char* companion = std::string(tool.parameters[index].kind) == "data" ? companionCell(tool.parameters[index].name) : nullptr;
    for (std::size_t j = 0; companion && structure.hasUnitCell && j < tool.parameters.size(); ++j) {
        auto& other = m_fields[j];
        if (std::string(tool.parameters[j].name) == companion) {
            std::snprintf(other.value.data(), other.value.size(), "%s", cellJson(structure).c_str());
            other.useFile = false;
            other.enabled = true;
        } else if (std::string(tool.parameters[j].name) == "pbc") {
            std::snprintf(other.value.data(), other.value.size(), "%s", "[true,true,true]");
            other.enabled = true;
        }
    }
}

void ScientificToolsDialog::drawField(std::size_t i, const Structure& structure)
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    const auto& parameter = tool.parameters[i];
    auto& field = m_fields[i];
    const std::string kind = parameter.kind;
    const bool optional = !parameter.required && parameter.value[0] == '\0';
    ImGui::PushID(static_cast<int>(i));
    if (optional) ImGui::Checkbox(parameter.label, &field.enabled);
    else if (kind != "bool") ImGui::TextWrapped("%s", parameter.label);
    if (optional && !field.enabled) {
        ImGui::Spacing();
        ImGui::PopID();
        return;
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
        if (kind == "data" || kind == "structure" || kind == "file") {
            int source = field.useFile ? 0 : 1;
            if (kind != "file") {
                ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(170)));
                if (ImGui::Combo("##source", &source, "Load from file\0Enter values\0")) field.useFile = source == 0;
            }
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
            const std::string name = parameter.name;
            const bool cellInput = name == "cell" || name == "reference_cell" || name == "current_cell";
            if (kind == "data" && cellInput && structure.hasUnitCell && responsive::button("Use active cell")) {
                std::snprintf(field.value.data(), field.value.size(), "%s", cellJson(structure).c_str());
                field.useFile = false;
            }
            if (!structure.atoms.empty() && (kind == "structure" || companionCell(name)) && responsive::button("Use active structure")) {
                try { useActiveStructure(i, structure); }
                catch (const std::exception& error) { m_error = error.what(); }
            }
            if (kind == "data" && !m_trajectory.empty() && acceptsTrajectory(tool.id, name)) {
                if (responsive::button("Use loaded trajectory")) {
                    field.path = m_trajectory;
                    field.useFile = true;
                    field.field[0] = '\0';
                    field.selectColumn = false;
                }
                ImGui::TextDisabled("Loaded: %s", std::filesystem::u8path(m_trajectory).filename().u8string().c_str());
            }
        }
        if (kind == "calculator") {
            ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340)));
            ImGui::Combo("##potential", &field.potential, "EMT (effective-medium theory)\0Lennard-Jones 12-6\0EAM / Finnis-Sinclair file\0");
            if (field.potential == 2) {
                if (responsive::button("Browse potential...")) { m_pickerTarget = -100 - static_cast<int>(i); m_picker.open("EAM potential file", false, field.potentialFile); }
                ImGui::TextWrapped("%s", field.potentialFile.empty() ? "Choose a LAMMPS-format .eam.alloy, .eam.fs or .eam file" : field.potentialFile.c_str());
                ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(340)));
                ImGui::Combo("Table format", &field.eamFormat, "Detect from file name\0setfl (eam/alloy)\0Finnis-Sinclair (eam/fs)\0funcfl (single element eam)\0");
                ImGui::TextDisabled("Elements are matched by symbol; the file's cutoff is used.");
            } else if (field.potential == 0) {
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

void ScientificToolsDialog::drawFields(const Structure& structure, int group, int skip)
{
    // group: 0 sources (data, structures, files), 1 settings, 2 optional, 3 potential.
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        if (static_cast<int>(i) == skip || suppliedByActive(i)) continue;
        const auto& parameter = tool.parameters[i];
        const std::string kind = parameter.kind;
        const bool optional = !parameter.required && parameter.value[0] == '\0';
        const int fieldGroup = kind == "calculator" ? 3 : optional ? 2 : isSource(i) ? 0 : 1;
        if (fieldGroup == group) drawField(i, structure);
    }
}

bool ScientificToolsDialog::hasFields(int group, int skip) const
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        if (static_cast<int>(i) == skip || suppliedByActive(i)) continue;
        const auto& parameter = tool.parameters[i];
        const std::string kind = parameter.kind;
        const bool optional = !parameter.required && parameter.value[0] == '\0';
        if ((kind == "calculator" ? 3 : optional ? 2 : isSource(i) ? 0 : 1) == group) return true;
    }
    return false;
}

void ScientificToolsDialog::drawFieldGrid(const Structure& structure, int group, int skip, int columns)
{
    // Lays a group of fields across columns, keeping short settings on one row.
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    std::vector<std::size_t> members;
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        if (static_cast<int>(i) == skip) continue;
        const auto& parameter = tool.parameters[i];
        const std::string kind = parameter.kind;
        const bool optional = !parameter.required && parameter.value[0] == '\0';
        if ((kind == "calculator" ? 3 : optional ? 2 : isSource(i) ? 0 : 1) == group) members.push_back(i);
    }
    if (members.empty()) return;
    columns = std::max(1, std::min(columns, static_cast<int>(members.size())));
    if (ImGui::GetContentRegionAvail().x < responsive::dp(300) * columns) columns = std::max(1, static_cast<int>(ImGui::GetContentRegionAvail().x / responsive::dp(300)));
    const std::string id = "Field grid " + std::to_string(group);
    if (ImGui::BeginTable(id.c_str(), columns, ImGuiTableFlags_SizingStretchSame)) {
        for (std::size_t n = 0; n < members.size(); ++n) {
            if (n % static_cast<std::size_t>(columns) == 0) ImGui::TableNextRow();
            ImGui::TableNextColumn();
            drawField(members[n], structure);
        }
        ImGui::EndTable();
    }
}

void ScientificToolsDialog::startCalculation(const Structure& structure)
{
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    try {
        if (layoutFor(tool.id) == Layout::Atoms && m_useActive) {
            const int input = activeInput();
            if (input >= 0) {
                if (structure.atoms.empty()) throw std::runtime_error("The active structure has no atoms; load a structure or choose a file instead");
                useActiveStructure(static_cast<std::size_t>(input), structure);
            }
        }
        Json request = Json::object();
        for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
            const auto& field = m_fields[i];
            if (!field.enabled) continue;
            const auto& parameter = tool.parameters[i];
            const std::string kind = parameter.kind;
            if (kind == "calculator") {
                request[parameter.name] = potentialJson(field.potential, field.epsilon, field.sigma, field.cutoff, field.potentialFile, field.eamFormat);
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
                // Plain text (e.g. an element symbol) is accepted without JSON quotes.
                if (kind == "string" && value.find_first_not_of(" \t") != std::string::npos && value[value.find_first_not_of(" \t")] != '"')
                    request[parameter.name] = value.substr(value.find_first_not_of(" \t"), value.find_last_not_of(" \t\r\n") - value.find_first_not_of(" \t") + 1);
                else try { request[parameter.name] = Json::parse(value); }
                catch (const std::exception& error) { throw std::runtime_error(std::string(parameter.label) + ": " + error.what()); }
            }
        }
        const auto directory = newRunDirectory();
        const std::string id = tool.id;
        m_error.clear();
        m_result = {};
        m_plot = 0;
        m_file = 0;
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
                          output.frames.empty() ? std::filesystem::path{} : structures,
                          atomforge::science::resultPlots(id, output.result),
                          atomforge::science::perAtomProperties(id, output.result),
                          generatedFiles(output.result), summaryOf(output.result), matricesOf(output.result)};
        });
    } catch (const std::exception& error) { m_error = error.what(); }
}

void ScientificToolsDialog::drawStatus()
{
    if (m_task.running()) {
        const float progress = m_task.progress();
        const float width = std::max(responsive::dp(120), ImGui::GetContentRegionAvail().x - responsive::dp(110));
        if (progress >= 0.0f) ImGui::ProgressBar(progress, ImVec2(width, 0), "Calculating...");
        else ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(width, 0), "Calculating...");
        ImGui::SameLine();
        if (responsive::button("Cancel")) m_task.cancel();
    }
    if (!m_error.empty()) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.95f, 0.45f, 0.40f, 1.0f));
        ImGui::TextWrapped("%s", m_error.c_str());
        ImGui::PopStyleColor();
    }
}

bool ScientificToolsDialog::drawRunButton(const char* label, const Structure& structure, float width)
{
    ImGui::BeginDisabled(m_task.running());
    const bool pressed = responsive::button(label, ImVec2(width, responsive::dp(34)));
    ImGui::EndDisabled();
    if (pressed) startCalculation(structure);
    return pressed;
}

void ScientificToolsDialog::drawSummary(const char* id, int pairsPerRow)
{
    if (m_result.summary.empty()) return;
    pairsPerRow = std::max(1, std::min(pairsPerRow, static_cast<int>(m_result.summary.size())));
    if (ImGui::BeginTable(id, 2 * pairsPerRow, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp)) {
        for (int pair = 0; pair < pairsPerRow; ++pair) {
            ImGui::TableSetupColumn("Quantity", ImGuiTableColumnFlags_WidthStretch, 0.6f);
            ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch, 0.4f);
        }
        for (std::size_t n = 0; n < m_result.summary.size(); ++n) {
            if (n % static_cast<std::size_t>(pairsPerRow) == 0) ImGui::TableNextRow();
            const auto& [label, value] = m_result.summary[n];
            ImGui::TableNextColumn(); ImGui::TextWrapped("%s", label.c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(value.c_str());
        }
        ImGui::EndTable();
    }
}

void ScientificToolsDialog::drawPropertyControls(const Structure& structure, const ColourAtoms& colourAtoms)
{
    if (m_result.properties.empty() || !colourAtoms) return;
    m_property = std::clamp(m_property, 0, static_cast<int>(m_result.properties.size()) - 1);
    for (std::size_t p = 0; p < m_result.properties.size(); ++p) {
        const auto& property = m_result.properties[p];
        const bool matches = property.values.size() == structure.atoms.size();
        ImGui::BeginDisabled(!matches);
        if (ImGui::RadioButton(property.name.c_str(), static_cast<int>(p) == m_property)) {
            m_property = static_cast<int>(p);
            colourAtoms(property.name, property.values);
        }
        ImGui::EndDisabled();
    }
    const auto& selected = m_result.properties[static_cast<std::size_t>(m_property)];
    if (selected.values.size() != structure.atoms.size())
        ImGui::TextDisabled("The active structure has %zu atoms; this result has %zu.", structure.atoms.size(), selected.values.size());
    else if (responsive::button("Colour atoms in the view")) colourAtoms(selected.name, selected.values);
}

void ScientificToolsDialog::drawPlots(float height)
{
    if (m_result.plots.empty()) return;
    m_plot = std::clamp(m_plot, 0, static_cast<int>(m_result.plots.size()) - 1);
    if (m_result.plots.size() > 1 && ImGui::BeginTabBar("Plots")) {
        for (std::size_t i = 0; i < m_result.plots.size(); ++i) {
            const std::string name = (m_result.plots[i].title.empty() ? "Plot " + std::to_string(i + 1) : m_result.plots[i].title) + "##plot" + std::to_string(i);
            if (ImGui::BeginTabItem(name.c_str())) { m_plot = static_cast<int>(i); ImGui::EndTabItem(); }
        }
        ImGui::EndTabBar();
    }
    const float plotHeight = std::max(responsive::dp(200), height - ImGui::GetFrameHeightWithSpacing());
    uiPlot::drawSciencePlot("##plot", m_result.plots[static_cast<std::size_t>(m_plot)], plotHeight / responsive::scale());
    if (responsive::button("Export plot data...")) {
        m_pickerTarget = -10 - m_plot;
        m_picker.open("Save plot data", true, "plot.csv");
    }
}

void ScientificToolsDialog::drawSaveButtons(const std::function<void(Structure&)>& loadResult)
{
    if (responsive::button("Save results...")) { m_pickerTarget = -3; m_picker.open("Save scientific results", true, "analysis.json"); }
    if (!m_result.structures.empty()) {
        if (ImGui::GetContentRegionAvail().x > responsive::dp(380)) ImGui::SameLine();
        if (responsive::button("Save structures...")) { m_pickerTarget = -4; m_picker.open("Save result structures or trajectory", true, "frames.extxyz"); }
        if (ImGui::GetContentRegionAvail().x > responsive::dp(380)) ImGui::SameLine();
        if (responsive::button("Open final structure in a new tab")) {
            try {
                auto frames = atomforge::science::readFrames(m_result.structures);
                if (!frames.empty()) loadResult(frames.back().structure);
            } catch (const std::exception& error) { m_error = error.what(); }
        }
    }
}

void ScientificToolsDialog::drawDetails()
{
    if (m_result.report.empty()) return;
    if (ImGui::CollapsingHeader("Full report")) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(m_result.report.c_str());
        ImGui::PopTextWrapPos();
    }
}

// Per-atom analyses: one panel next to the viewport, fed by the active structure.
void ScientificToolsDialog::drawAtomLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult, const ColourAtoms& colourAtoms)
{
    const int input = activeInput();
    ImGui::BeginDisabled(m_task.running());
    ImGui::SeparatorText("Atoms");
    if (input >= 0) {
        ImGui::Checkbox("Analyse the active structure", &m_useActive);
        if (m_useActive) {
            if (structure.atoms.empty()) ImGui::TextDisabled("No structure is loaded.");
            else ImGui::TextDisabled("%zu atoms, %s", structure.atoms.size(), structure.hasUnitCell ? "periodic cell" : "no cell (open boundaries)");
        } else drawField(static_cast<std::size_t>(input), structure);
    }
    const int skip = m_useActive ? input : -1;
    if (hasFields(0, skip)) drawFields(structure, 0, skip);
    if (hasFields(1, skip)) { ImGui::SeparatorText("Settings"); drawFields(structure, 1, skip); }
    if (hasFields(2, skip) && ImGui::CollapsingHeader("More options")) drawFields(structure, 2, skip);
    ImGui::EndDisabled();
    ImGui::Spacing();
    drawRunButton("Analyse", structure, -1);
    drawStatus();
    if (!m_result.output.empty()) {
        ImGui::SeparatorText("Result");
        drawSummary("Atom summary");
        if (!m_result.properties.empty() && colourAtoms) {
            ImGui::SeparatorText("Colour atoms by");
            drawPropertyControls(structure, colourAtoms);
            ImGui::Checkbox("Colour automatically after each analysis", &m_autoColour);
        }
        if (!m_result.plots.empty() && ImGui::CollapsingHeader("Plots", ImGuiTreeNodeFlags_DefaultOpen)) drawPlots(responsive::dp(260));
        drawSaveButtons(loadResult);
        drawDetails();
    }
}

// Curves and spectra: settings sidebar and a large plot area.
void ScientificToolsDialog::drawPlotLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult)
{
    const bool stack = responsive::stacked(900.0f);
    const float sidebar = stack ? 0.0f : std::min(responsive::dp(380), ImGui::GetContentRegionAvail().x * 0.42f);
    const float height = stack ? responsive::dp(360) : ImGui::GetContentRegionAvail().y;
    responsive::beginChild("Plot settings", ImVec2(sidebar, height), true);
    ImGui::BeginDisabled(m_task.running());
    if (hasFields(0, -1)) { ImGui::SeparatorText("Data"); drawFields(structure, 0, -1); }
    if (hasFields(1, -1)) { ImGui::SeparatorText("Settings"); drawFields(structure, 1, -1); }
    if (hasFields(2, -1) && ImGui::CollapsingHeader("Optional inputs")) drawFields(structure, 2, -1);
    ImGui::EndDisabled();
    ImGui::Spacing();
    drawRunButton("Compute", structure, -1);
    drawStatus();
    ImGui::EndChild();
    if (!stack) ImGui::SameLine();
    responsive::beginChild("Plot area", ImVec2(0, stack ? 0 : height), false);
    if (m_result.output.empty()) {
        ImGui::TextDisabled(m_task.running() ? "Calculating..." : "Results appear here.");
    } else {
        const float plotHeight = m_result.plots.empty() ? 0.0f :
            std::max(responsive::dp(240), ImGui::GetContentRegionAvail().y * 0.62f);
        drawPlots(plotHeight);
        ImGui::Spacing();
        drawSummary("Plot summary");
        drawSaveButtons(loadResult);
        drawDetails();
    }
    ImGui::EndChild();
}

// Trajectory analyses: the trajectory source across the top, settings in a row, a wide plot.
void ScientificToolsDialog::drawTrajectoryLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult)
{
    ImGui::BeginDisabled(m_task.running());
    responsive::beginChild("Trajectory source", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
    ImGui::SeparatorText("Trajectory");
    if (m_trajectory.empty()) ImGui::TextDisabled("Load frames in Trajectory playback to use them here, or choose files below.");
    else ImGui::Text("Playback trajectory: %s", std::filesystem::u8path(m_trajectory).filename().u8string().c_str());
    drawFieldGrid(structure, 0, -1, 2);
    ImGui::EndChild();
    if (hasFields(1, -1)) { ImGui::SeparatorText("Settings"); drawFieldGrid(structure, 1, -1, 3); }
    if (hasFields(2, -1) && ImGui::TreeNode("Optional inputs")) { drawFieldGrid(structure, 2, -1, 3); ImGui::TreePop(); }
    ImGui::EndDisabled();
    drawRunButton("Analyse trajectory", structure, responsive::dp(220));
    drawStatus();
    if (m_result.output.empty()) return;
    ImGui::Separator();
    responsive::beginChild("Trajectory result", ImVec2(0, 0), false);
    drawPlots(std::max(responsive::dp(260), ImGui::GetContentRegionAvail().y * 0.6f));
    drawSummary("Trajectory summary");
    drawSaveButtons(loadResult);
    drawDetails();
    ImGui::EndChild();
}

// Material properties: inputs on one side, the computed quantities as large values and tensors.
void ScientificToolsDialog::drawPropertiesLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult)
{
    const bool stack = responsive::stacked(860.0f);
    const float formWidth = stack ? 0.0f : ImGui::GetContentRegionAvail().x * 0.42f;
    const float height = stack ? responsive::dp(330) : ImGui::GetContentRegionAvail().y;
    responsive::beginChild("Property inputs", ImVec2(formWidth, height), true);
    ImGui::BeginDisabled(m_task.running());
    if (hasFields(0, -1)) { ImGui::SeparatorText("Data"); drawFields(structure, 0, -1); }
    if (hasFields(1, -1)) { ImGui::SeparatorText("Settings"); drawFields(structure, 1, -1); }
    if (hasFields(2, -1) && ImGui::TreeNode("Optional inputs")) { drawFields(structure, 2, -1); ImGui::TreePop(); }
    ImGui::EndDisabled();
    drawRunButton("Calculate", structure, -1);
    drawStatus();
    ImGui::EndChild();
    if (!stack) ImGui::SameLine();
    responsive::beginChild("Property values", ImVec2(0, stack ? 0 : height), true);
    if (m_result.output.empty()) {
        ImGui::TextDisabled(m_task.running() ? "Calculating..." : "Computed properties appear here.");
        ImGui::EndChild();
        return;
    }
    // Headline quantities as cards: label above a large value.
    const int columns = std::max(1, static_cast<int>(ImGui::GetContentRegionAvail().x / responsive::dp(190)));
    if (!m_result.summary.empty() && ImGui::BeginTable("Property cards", columns, ImGuiTableFlags_SizingStretchSame)) {
        for (const auto& [label, value] : m_result.summary) {
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label.c_str());
            ImGui::SetWindowFontScale(1.45f);
            ImGui::TextUnformatted(value.c_str());
            ImGui::SetWindowFontScale(1.0f);
            ImGui::Spacing();
        }
        ImGui::EndTable();
    }
    for (std::size_t m = 0; m < m_result.matrices.size(); ++m) {
        const auto& [label, rows] = m_result.matrices[m];
        ImGui::SeparatorText(label.c_str());
        const std::string id = "Matrix " + std::to_string(m);
        if (ImGui::BeginTable(id.c_str(), static_cast<int>(rows.front().size()), ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchSame)) {
            for (const auto& row : rows) {
                ImGui::TableNextRow();
                for (double value : row) { ImGui::TableNextColumn(); ImGui::Text("%.5g", value); }
            }
            ImGui::EndTable();
        }
    }
    if (!m_result.plots.empty() && ImGui::CollapsingHeader("Plots", ImGuiTreeNodeFlags_DefaultOpen)) drawPlots(responsive::dp(280));
    drawSaveButtons(loadResult);
    drawDetails();
    ImGui::EndChild();
}

// Simulations: numbered setup steps, an inline run bar and a run summary.
void ScientificToolsDialog::drawSimulationLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult)
{
    responsive::beginChild("Simulation steps", ImVec2(0, m_result.output.empty() ? -responsive::dp(96) : responsive::dp(220)), true);
    ImGui::BeginDisabled(m_task.running());
    int step = 1;
    const auto heading = [&step](const char* text) { ImGui::SeparatorText((std::to_string(step++) + "  " + text).c_str()); };
    if (hasFields(0, -1)) { heading("System"); drawFields(structure, 0, -1); }
    if (hasFields(3, -1)) { heading("Interatomic potential"); drawFields(structure, 3, -1); }
    if (hasFields(1, -1)) { heading("Protocol"); drawFields(structure, 1, -1); }
    if (hasFields(2, -1) && ImGui::TreeNode("Advanced options")) { drawFields(structure, 2, -1); ImGui::TreePop(); }
    ImGui::EndDisabled();
    ImGui::EndChild();
    drawRunButton(m_task.running() ? "Running..." : "Run simulation", structure, responsive::dp(200));
    drawStatus();
    if (!m_result.output.empty()) {
        responsive::beginChild("Simulation result", ImVec2(0, 0), false);
        ImGui::SeparatorText("Run summary");
        const bool wide = ImGui::GetContentRegionAvail().x > responsive::dp(760);
        if (wide && !m_result.plots.empty()) {
            responsive::beginChild("Run numbers", ImVec2(ImGui::GetContentRegionAvail().x * 0.4f, responsive::dp(320)), false);
            drawSummary("Simulation summary");
            ImGui::EndChild();
            ImGui::SameLine();
            responsive::beginChild("Run plots", ImVec2(0, responsive::dp(320)), false);
            drawPlots(responsive::dp(280));
            ImGui::EndChild();
        } else {
            drawSummary("Simulation summary", wide ? 2 : 1);
            drawPlots(responsive::dp(260));
        }
        drawSaveButtons(loadResult);
        drawDetails();
        ImGui::EndChild();
    }
}

// Input generators: a short form and a browser for the generated files.
void ScientificToolsDialog::drawGeneratorLayout(const Structure& structure)
{
    responsive::beginChild("Generator form", ImVec2(0, std::max(responsive::dp(220), ImGui::GetContentRegionAvail().y * 0.42f)), true);
    ImGui::BeginDisabled(m_task.running());
    if (ImGui::BeginTable("Generator columns", 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("Source", ImGuiTableColumnFlags_WidthStretch, 0.36f);
        ImGui::TableSetupColumn("Options", ImGuiTableColumnFlags_WidthStretch, 0.64f);
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SeparatorText("Structure and data");
        drawFields(structure, 0, -1);
        if (hasFields(3, -1)) { ImGui::SeparatorText("Potential"); drawFields(structure, 3, -1); }
        ImGui::TableNextColumn();
        ImGui::SeparatorText("Options");
        drawFieldGrid(structure, 1, -1, 2);
        if (hasFields(2, -1) && ImGui::TreeNode("Optional")) { drawFieldGrid(structure, 2, -1, 2); ImGui::TreePop(); }
        ImGui::EndTable();
    }
    ImGui::EndDisabled();
    ImGui::EndChild();
    drawRunButton("Generate files", structure, responsive::dp(200));
    drawStatus();
    ImGui::SeparatorText("Generated files");
    if (m_result.files.empty()) {
        ImGui::TextDisabled(m_result.output.empty() ? "Generate to preview the files here." : "This run produced no files.");
        return;
    }
    m_file = std::clamp(m_file, 0, static_cast<int>(m_result.files.size()) - 1);
    const float listWidth = std::min(responsive::dp(220), ImGui::GetContentRegionAvail().x * 0.3f);
    const float height = std::max(responsive::dp(200), ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 2);
    responsive::beginChild("File list", ImVec2(listWidth, height), true);
    for (std::size_t i = 0; i < m_result.files.size(); ++i) {
        const auto& [name, text] = m_result.files[i];
        if (ImGui::Selectable(name.c_str(), static_cast<int>(i) == m_file)) m_file = static_cast<int>(i);
        ImGui::TextDisabled("%zu lines", static_cast<std::size_t>(std::count(text.begin(), text.end(), '\n')));
    }
    ImGui::EndChild();
    ImGui::SameLine();
    auto& [name, text] = m_result.files[static_cast<std::size_t>(m_file)];
    ImGui::InputTextMultiline("##preview", text.data(), text.size() + 1, ImVec2(-1, height), ImGuiInputTextFlags_ReadOnly);
    if (responsive::button("Save this file...")) { m_pickerTarget = -6; m_picker.open("Save generated file", true, name); }
    ImGui::SameLine();
    if (responsive::button("Save all files...")) { m_pickerTarget = -5; m_picker.open("Choose the folder (save as the first file)", true, m_result.files.front().first); }
    ImGui::SameLine();
    if (responsive::button("Copy to clipboard")) ImGui::SetClipboardText(text.c_str());
    ImGui::SameLine();
    if (responsive::button("Save results...")) { m_pickerTarget = -3; m_picker.open("Save scientific results", true, "analysis.json"); }
}

void ScientificToolsDialog::draw(const Structure& structure, const std::function<void(Structure&)>& loadResult,
                                 const ColourAtoms& colourAtoms)
{
    if (m_task.poll()) {
        m_error = m_task.error();
        if (m_task.result()) {
            m_result = std::move(*m_task.result());
            // Atom analyses colour the view straight away.
            if (m_tool >= 0 && layoutFor(scienceToolCatalog()[static_cast<std::size_t>(m_tool)].id) == Layout::Atoms &&
                m_autoColour && colourAtoms && !m_result.properties.empty()) {
                m_property = std::clamp(m_property, 0, static_cast<int>(m_result.properties.size()) - 1);
                const auto& property = m_result.properties[static_cast<std::size_t>(m_property)];
                if (property.values.size() == structure.atoms.size()) colourAtoms(property.name, property.values);
            }
        }
        m_task.clearResult();
    }
    if (!m_open) return;
    if (m_tool < 0) selectTool(0);
    const auto& tool = scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
    const Layout layout = layoutFor(tool.id);
    const std::string title = std::string(tool.title) + windowId(layout);
    switch (layout) {
    case Layout::Atoms: responsive::windowSize(ImVec2(470, 760), ImGuiCond_FirstUseEver); break;
    case Layout::Plot: responsive::windowSize(ImVec2(1120, 700), ImGuiCond_FirstUseEver); break;
    case Layout::Trajectory: responsive::windowSize(ImVec2(1000, 820), ImGuiCond_FirstUseEver); break;
    case Layout::Properties: responsive::windowSize(ImVec2(1040, 680), ImGuiCond_FirstUseEver); break;
    case Layout::Simulation: responsive::windowSize(ImVec2(900, 800), ImGuiCond_FirstUseEver); break;
    case Layout::Generator: responsive::windowSize(ImVec2(980, 760), ImGuiCond_FirstUseEver); break;
    }
    responsive::windowConstraints(ImVec2(380, 320), ImVec2(1800, 1400));
    if (responsive::begin(title.c_str(), &m_open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, responsive::size(10, 6));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, responsive::size(10, 8));
        ImGui::TextDisabled("%s", tool.category);
        const std::string summary = std::string(tool.help).substr(0, std::string(tool.help).find('\n'));
        ImGui::TextWrapped("%s", summary.c_str());
        ImGui::Spacing();
        switch (layout) {
        case Layout::Atoms: drawAtomLayout(structure, loadResult, colourAtoms); break;
        case Layout::Plot: drawPlotLayout(structure, loadResult); break;
        case Layout::Trajectory: drawTrajectoryLayout(structure, loadResult); break;
        case Layout::Properties: drawPropertiesLayout(structure, loadResult); break;
        case Layout::Simulation: drawSimulationLayout(structure, loadResult); break;
        case Layout::Generator: drawGeneratorLayout(structure); break;
        }
        if (auto path = m_picker.draw()) {
            try {
                if (m_pickerTarget == -5) {
                    const auto folder = std::filesystem::u8path(*path).parent_path();
                    for (const auto& [name, text] : m_result.files) writeText(folder / std::filesystem::u8path(name), text);
                } else if (m_pickerTarget == -6 && static_cast<std::size_t>(m_file) < m_result.files.size())
                    writeText(std::filesystem::u8path(*path), m_result.files[static_cast<std::size_t>(m_file)].second);
                else if (m_pickerTarget == -3) std::filesystem::copy_file(m_result.output, std::filesystem::u8path(*path), std::filesystem::copy_options::overwrite_existing);
                else if (m_pickerTarget <= -10 && m_pickerTarget > -100 && static_cast<std::size_t>(-10 - m_pickerTarget) < m_result.plots.size())
                    writeText(std::filesystem::u8path(*path), atomforge::science::plotCsv(m_result.plots[static_cast<std::size_t>(-10 - m_pickerTarget)]));
                else if (m_pickerTarget == -4) std::filesystem::copy_file(m_result.structures, std::filesystem::u8path(*path), std::filesystem::copy_options::overwrite_existing);
                else if (m_pickerTarget <= -100 && static_cast<std::size_t>(-100 - m_pickerTarget) < m_fields.size())
                    m_fields[static_cast<std::size_t>(-100 - m_pickerTarget)].potentialFile = *path;
                else if (m_pickerTarget >= 0 && static_cast<std::size_t>(m_pickerTarget) < m_fields.size()) {
                    auto& field = m_fields[static_cast<std::size_t>(m_pickerTarget)]; field.path = *path; field.useFile = true;
                }
            } catch (const std::exception& error) { m_error = error.what(); }
        }
        ImGui::PopStyleVar(2);
    }
    ImGui::End();
}
