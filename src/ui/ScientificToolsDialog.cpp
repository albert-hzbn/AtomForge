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

// Index order of the desktop's interatomic-potential selector.
enum PotentialChoice { kEmt, kLennardJones, kEam, kTersoff, kStillingerWeber, kBuckingham };

Json potentialJson(int potential, double epsilon, double sigma, double cutoff, const std::string& file, int format,
                   const std::string& options)
{
    Json result = Json::object();
    if (potential == kTersoff || potential == kStillingerWeber) {
        result["potential"] = potential == kTersoff ? "Tersoff" : "StillingerWeber";
        if (!file.empty()) result["file"] = file;  // empty: built-in Si
        return result;
    }
    if (potential == kBuckingham) {
        try { result = Json::parse(options); }
        catch (const std::exception& error) { throw std::runtime_error(std::string("Buckingham parameters: ") + error.what()); }
        if (!result.isObject()) throw std::runtime_error("Buckingham parameters must be a JSON object");
        result["potential"] = "Buckingham";
        return result;
    }
    if (potential == kEmt) {
        result["potential"] = "EMT";
        return result;
    }
    if (potential == kEam) {
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

ScientificToolsDialog::Field::Field()
{
    std::snprintf(potentialOptions.data(), potentialOptions.size(), "%s",
        "{\"charges\": {\"Mg\": 2, \"O\": -2}, \"cutoff\": 10,\n"
        " \"pairs\": [{\"elements\": [\"Mg\", \"O\"], \"A\": 821.6, \"rho\": 0.3242, \"C\": 0},\n"
        "           {\"elements\": [\"O\", \"O\"], \"A\": 22764, \"rho\": 0.149, \"C\": 27.88}]}");
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
    m_kept.clear();
}

std::string ScientificToolsDialog::snapshot() const
{
    Json state = Json::object();
    state["open"] = m_open;
    if (m_tool < 0) return state.dump();
    state["tool"] = currentTool().id;
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
        item["potential_options"] = std::string(field.potentialOptions.data());
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
        const int index = tool && tool->isString() ? catalogIndex(tool->string()) : -1;
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
                copyText(field.potentialOptions, item.find("potential_options"));
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
            const std::string id = currentTool().id;
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
    if (tool == "trajectory-structure") return parameter == "trajectory_file";
    return false;
}

bool ScientificToolsDialog::open(const std::string& toolId)
{
    const int index = catalogIndex(toolId);
    if (index < 0 || m_task.running()) return false;
    selectTool(index);
    m_open = true;
    return true;
}

ScientificToolsDialog::Layout ScientificToolsDialog::layoutFor(const std::string& tool)
{
    const ScienceToolDef* found = findScienceTool(tool);
    return found ? found->view : Layout::Plot;
}

std::string ScientificToolsDialog::cellJson(const Structure& structure)
{
    Json cell = Json::array();
    for (const auto& row : structure.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
    return cell.dump();
}

int ScientificToolsDialog::catalogIndex(const std::string& id)
{
    const ScienceToolDef* found = findScienceTool(id);
    return found ? static_cast<int>(found - scienceToolCatalog().data()) : -1;
}

const ScienceToolDef& ScientificToolsDialog::currentTool() const
{
    return scienceToolCatalog()[static_cast<std::size_t>(m_tool)];
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
    const auto& tool = currentTool();
    const auto& parameter = tool.parameters[index];
    const std::string kind = parameter.kind;
    return kind == "data" || kind == "structure" || kind == "file";
}

bool ScientificToolsDialog::suppliedByActive(std::size_t index) const
{
    const auto& tool = currentTool();
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
    const auto& tool = currentTool();
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        const std::string kind = tool.parameters[i].kind, name = tool.parameters[i].name;
        if (kind == "structure" || (kind == "data" && name == "positions")) return static_cast<int>(i);
    }
    return -1;
}

void ScientificToolsDialog::useActiveStructure(std::size_t index, const Structure& structure)
{
    const auto& tool = currentTool();
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

bool ScientificToolsDialog::hasFields(int group, int skip) const
{
    const auto& tool = currentTool();
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        if (static_cast<int>(i) == skip || suppliedByActive(i)) continue;
        const auto& parameter = tool.parameters[i];
        const std::string kind = parameter.kind;
        const bool optional = !parameter.required && parameter.value[0] == '\0';
        if ((kind == "calculator" ? 3 : optional ? 2 : isSource(i) ? 0 : 1) == group) return true;
    }
    return false;
}

void ScientificToolsDialog::startCalculation(const Structure& structure)
{
    const auto& tool = currentTool();
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
                request[parameter.name] = potentialJson(field.potential, field.epsilon, field.sigma, field.cutoff, field.potentialFile, field.eamFormat, field.potentialOptions.data());
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
                // Plain text (e.g. an element symbol) is accepted without JSON quotes;
                // text starting with a quote, { or [ is read as JSON.
                const char first = value[value.find_first_not_of(" \t\r\n")];
                if (kind == "string" && first != '"' && first != '{' && first != '[')
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

void ScientificToolsDialog::draw(const Structure& structure, const std::function<void(Structure&)>& loadResult,
                                 const ColourAtoms& colourAtoms)
{
    if (m_task.poll()) {
        m_error = m_task.error();
        if (m_task.result()) {
            m_result = std::move(*m_task.result());
            // Atom analyses colour the view straight away.
            if (m_tool >= 0 && layoutFor(currentTool().id) == Layout::Atoms &&
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
    const auto& tool = currentTool();
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
                    writeText(std::filesystem::u8path(*path), atomforge::science::plotCsv(displayedPlot(static_cast<std::size_t>(-10 - m_pickerTarget))));
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
