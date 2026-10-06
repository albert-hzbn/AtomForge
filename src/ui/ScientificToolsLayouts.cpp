// Drawing of the scientific tool windows: parameter widgets, result views and
// one layout per tool view (see ScienceToolView). State, persistence and
// running live in ScientificToolsDialog.cpp.
#include "ui/ScientificToolsDialog.h"
#include "science/ScienceData.h"
#include "ui/ResponsiveLayout.h"
#include "ui/SciencePlot.h"
#include "imgui.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>

using atomforge::science::Json;

void ScientificToolsDialog::drawField(std::size_t i, const Structure& structure)
{
    const auto& tool = currentTool();
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
            if ((kind == "data" || kind == "file") && !m_trajectory.empty() && acceptsTrajectory(tool.id, name)) {
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
            ImGui::Combo("##potential", &field.potential,
                         "EMT (effective-medium theory)\0Lennard-Jones 12-6\0EAM / Finnis-Sinclair file\0"
                         "Tersoff (covalent)\0Stillinger-Weber (covalent)\0Buckingham + Coulomb (ionic)\0");
            if (field.potential == 3 || field.potential == 4) {
                const bool tersoff = field.potential == 3;
                if (responsive::button("Browse parameter file...")) { m_pickerTarget = -100 - static_cast<int>(i); m_picker.open(tersoff ? "Tersoff parameter file" : "Stillinger-Weber parameter file", false, field.potentialFile); }
                if (!field.potentialFile.empty()) {
                    ImGui::SameLine();
                    if (responsive::button("Use built-in Si")) field.potentialFile.clear();
                }
                ImGui::TextWrapped("%s", field.potentialFile.empty()
                    ? (tersoff ? "Built-in Si (Tersoff 1988). Choose a LAMMPS .tersoff file for other elements, e.g. SiC.tersoff."
                               : "Built-in Si (Stillinger-Weber 1985). Choose a LAMMPS .sw file for other elements.")
                    : field.potentialFile.c_str());
            } else if (field.potential == 5) {
                ImGui::TextWrapped("Pairs A exp(-r/rho) - C/r^6 and fixed charges (e); Coulomb terms are Ewald-summed in periodic cells.");
                ImGui::InputTextMultiline("##buckingham", field.potentialOptions.data(), field.potentialOptions.size(),
                                          ImVec2(-1, ImGui::GetTextLineHeight() * 6));
            } else if (field.potential == 2) {
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
    const auto& tool = currentTool();
    for (std::size_t i = 0; i < tool.parameters.size(); ++i) {
        if (static_cast<int>(i) == skip || suppliedByActive(i)) continue;
        const auto& parameter = tool.parameters[i];
        const std::string kind = parameter.kind;
        const bool optional = !parameter.required && parameter.value[0] == '\0';
        const int fieldGroup = kind == "calculator" ? 3 : optional ? 2 : isSource(i) ? 0 : 1;
        if (fieldGroup == group) drawField(i, structure);
    }
}

void ScientificToolsDialog::drawFieldGrid(const Structure& structure, int group, int skip, int columns)
{
    // Lays a group of fields across columns, keeping short settings on one row.
    const auto& tool = currentTool();
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
    const float plotHeight = std::max(responsive::dp(200), height - 2 * ImGui::GetFrameHeightWithSpacing());
    uiPlot::drawSciencePlot("##plot", displayedPlot(static_cast<std::size_t>(m_plot)), plotHeight / responsive::scale());
    if (responsive::button("Export plot data...")) {
        m_pickerTarget = -10 - m_plot;
        m_picker.open("Save plot data", true, "plot.csv");
    }
    if (!m_kept.empty()) {
        ImGui::SameLine();
        ImGui::Checkbox(("Overlay " + std::to_string(m_kept.size()) + " kept run" + (m_kept.size() > 1 ? "s" : "")).c_str(), &m_overlay);
    }
}

atomforge::science::PlotSpec ScientificToolsDialog::displayedPlot(std::size_t index) const
{
    const auto& current = m_result.plots.at(index);
    if (!m_overlay || m_kept.empty()) return current;
    std::vector<std::pair<std::string, std::vector<atomforge::science::PlotSpec>>> others;
    for (const auto& run : m_kept) others.push_back({run.label, run.plots});
    return atomforge::science::overlayPlots(current, "Current", others);
}

void ScientificToolsDialog::keepForComparison()
{
    if (m_result.output.empty()) return;
    m_kept.push_back({"Run " + std::to_string(m_kept.size() + 1), m_result.plots, m_result.summary});
}

void ScientificToolsDialog::drawComparison()
{
    if (m_kept.empty() || m_result.summary.empty()) return;
    if (!ImGui::CollapsingHeader(("Compare with kept runs (" + std::to_string(m_kept.size()) + ")").c_str(), ImGuiTreeNodeFlags_DefaultOpen)) return;
    const int columns = static_cast<int>(m_kept.size()) + 2;
    if (ImGui::BeginTable("Run comparison", columns, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollX |
                                                     ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("Quantity");
        for (const auto& run : m_kept) ImGui::TableSetupColumn(run.label.c_str());
        ImGui::TableSetupColumn("Current");
        ImGui::TableHeadersRow();
        for (const auto& [label, value] : m_result.summary) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(label.c_str());
            for (const auto& run : m_kept) {
                ImGui::TableNextColumn();
                const auto found = std::find_if(run.summary.begin(), run.summary.end(), [&](const auto& row) { return row.first == label; });
                ImGui::TextUnformatted(found == run.summary.end() ? "-" : found->second.c_str());
            }
            ImGui::TableNextColumn(); ImGui::TextUnformatted(value.c_str());
        }
        ImGui::EndTable();
    }
    if (responsive::button("Clear kept runs")) m_kept.clear();
}

void ScientificToolsDialog::drawSaveButtons(const std::function<void(Structure&)>& loadResult)
{
    drawComparison();
    if (responsive::button("Keep for comparison")) keepForComparison();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Keep this result; change inputs and run again to overlay plots and compare values.");
    ImGui::SameLine();
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
    if (hasFields(3, -1)) { ImGui::SeparatorText("Interatomic potential"); drawFields(structure, 3, -1); }
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
