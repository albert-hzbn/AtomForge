// The command palette (Ctrl+K): finds any menu action by name and runs it.
#include "FileBrowser.h"
#include "ImGuiSetup.h"
#include "ui/DialogLayout.h"
#include "ui/MenuParity.h"
#include "ui/ResponsiveLayout.h"
#include "ui/UiScript.h"
#include "science/ScienceCatalog.h"
#include "pipeline/Pipeline.h"

#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cstdio>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Command
{
    std::string label, path;
    std::function<void()> run;
    bool enabled = true;
};

std::string lower(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

// Every word of the query must occur in the label or the menu path; matches
// at the start of the label, then in the label, rank first. -1: no match.
int score(const Command& command, const std::vector<std::string>& words)
{
    const std::string label = lower(command.label), path = lower(command.path);
    int total = 0;
    for (const auto& word : words) {
        const auto inLabel = label.find(word);
        if (inLabel == 0) total += 30;
        else if (inLabel != std::string::npos) total += label[inLabel - 1] == ' ' ? 20 : 10;
        else if (path.find(word) != std::string::npos) total += 2;
        else return -1;
    }
    return total + (command.enabled ? 1 : 0);
}
}

void FileBrowser::drawCommandPalette(Structure& structure, EditMenuDialogs& editMenuDialogs,
                                     const std::function<void(Structure&)>& updateBuffers)
{
    const ImGuiIO& io = ImGui::GetIO();
    if (!io.WantTextInput && io.KeyCtrl && (ImGui::IsKeyPressed(ImGuiKey_K, false) || (io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_P, false))))
        showCommandPalette = true;
    if (showCommandPalette) {
        ImGui::OpenPopup("##commandPalette");
        showCommandPalette = false;
        commandQuery[0] = '\0';
        commandSelected = 0;
    }

    const float width = std::min(responsive::dp(620.0f), io.DisplaySize.x - responsive::dp(32.0f));
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.12f), ImGuiCond_Always, ImVec2(0.5f, 0.0f));
    ImGui::SetNextWindowSize(ImVec2(width, 0));
    if (!ImGui::BeginPopup("##commandPalette", ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize)) return;

    // The commands, as in the menus.
    const bool atoms = !structure.atoms.empty();
    std::vector<Command> commands;
    const auto add = [&](const char* path, const std::string& label, std::function<void()> run, bool enabled = true) {
        commands.push_back({label, path, std::move(run), enabled});
    };
    add("File", "Open", [this] { openStructurePopup = true; });
    add("File", "Save", [this] { requestSaveFile = true; }, atoms);
    add("File", "Save As", [this] { saveStructurePopup = true; }, atoms);
    add("File", "Export Image", [this] { exportImageDialog(); }, atoms);
    add("File", "Close", [this] { requestCloseStructure = true; }, atoms);
    add("File", "Open Project", [this] { projectAction = 1; projectPicker.open("Open AtomForge project", false, "project.afproject"); });
    add("File", "Save Project As", [this] { projectAction = 2; projectPicker.open("Save AtomForge project", true, "project.afproject"); });
    add("File", "Quit", [this] { requestQuit = true; });
    add("Edit", undoDescription.empty() ? std::string("Undo") : "Undo " + undoDescription, [this] { requestUndo = true; }, !undoDescription.empty());
    add("Edit", redoDescription.empty() ? std::string("Redo") : "Redo " + redoDescription, [this] { requestRedo = true; }, !redoDescription.empty());
    add("Build", "Bulk Crystal", [this] { bulkCrystalDialog.open(); });
#if ATOMFORGE_ENABLE_SSS_BUILDER
    add("Build", "Substitutional Solid Solution", [this] { substitutionalSolidSolutionDialog.open(); });
#endif
#if ATOMFORGE_ENABLE_SFE_BUILDER
    add("Build", "Stacking Faults", [this] { stackingFaultDialog.open(); });
#endif
    add("Build", "CSL Grain Boundary", [this] { cslDialog.open(); });
    add("Build", "Interface Builder", [this] { interfaceBuilderDialog.open(); });
    add("Build", "Nanocrystal", [this] { nanoCrystalDialog.open(); });
    add("Build", "Custom Structure", [this] { customStructureDialog.open(); });
    add("Build", "Polycrystal", [this] { polyCrystalDialog.open(); });
    add("Build", "Amorphous Structure", [this] { amorphousBuilderDialog.open(); });
    add("Edit", "Edit Structure", [&editMenuDialogs] { editMenuDialogs.openEditStructure(); }, atoms);
    add("Edit", "Insert Dislocation", [this] { dislocationDialog.open(); });
    add("Edit", "Add Interstitial Atoms", [this] { interstitialAtomsDialog.open(); }, atoms);
    add("Edit", "Transform Structure", [this] { transformDialog.open(); }, structure.hasUnitCell);
    add("Edit", "Merge Structures", [this] { mergeStructuresDialog.open(); });
    add("Edit", "Cell Sculptor", [this] { cellSculptorDialog.open(); });
    add("Edit", "Structure pipeline", [this] { pipelineDialog.show(); }, atoms);
    add("Edit", "Box Select Mode", [this] { boxSelectMode = !boxSelectMode; if (boxSelectMode) lassoSelectMode = false; });
    add("Edit", "Lasso Select Mode", [this] { lassoSelectMode = !lassoSelectMode; if (lassoSelectMode) boxSelectMode = false; });
    for (const auto& location : menuLocations()) {
        if (location.dialog || !atomforge::pipeline::findModifierType(location.step)) continue;
        const std::string path = std::string(location.menu) + (location.submenu[0] ? std::string(" > ") + location.submenu : "");
        const std::string step = location.step;
        commands.push_back({location.label, path, [this, step] { operationDialog.open(step); },
                            atoms || step == "build-bulk"});
    }
    add("Analysis > Structure and defects", "Common Neighbour Analysis", [this] { cnaDialog.open(); }, atoms);
    add("Analysis > Structure and defects", "Radial Distribution Function", [this] { rdfDialog.open(); }, atoms);
    add("Analysis > Structure and defects", "Angular Distribution Function", [this] { angularDistributionDialog.open(); }, atoms);
    add("Analysis > Structure and defects", "Short-Range Order (SRO)", [this] { shortRangeOrderDialog.openRequested = true; }, atoms);
    add("Analysis > Structure and defects", "Dislocation Post-Processing", [this] { dislocationAnalysisDialog.open(); });
    add("Analysis > Trajectories and transport", "Trajectory playback", [this] { trajectoryDialog.show(); });
    add("Analysis > Electronic structure", "Electronic Post-processing", [this] { electronicDialog.open(); });
    add("Analysis > Electronic structure", "Wannier Analysis", [this] { wannierDialog.open(); });
    add("Analysis > Electronic structure", "LOBSTER Bonding Analysis", [this] { lobsterDialog.open(); });
    for (const auto& tool : scienceToolCatalog()) {
        const std::string id = tool.id;
        commands.push_back({tool.title, std::string("Analysis > ") + tool.category, [this, id] { scientificToolsDialog.open(id); }});
    }
    const auto toggle = [&](const char* label, bool& flag, bool enabled) {
        add("View > Show", label, [&flag] { flag = !flag; }, enabled);
    };
    toggle("Atoms", showAtoms, atoms);
    toggle("Element Labels", showElementLabels, true);
    toggle("Element Legend", showElementLegend, true);
    toggle("Bounding Box", showBoundingBox, structure.hasUnitCell);
    toggle("Lattice Planes", showLatticePlanes, structure.hasUnitCell);
    toggle("Miller Directions", showMillerDirections, structure.hasUnitCell);
    toggle("Voronoi Volumes", showVoronoi, atoms);
    toggle("Dislocation Lines", showDislocationLines, atoms);
    add("View > Projection", "Perspective", [this] { viewMode = ViewMode::Isometric; requestResetDefaultView = true; });
    add("View > Projection", "Orthographic", [this] { viewMode = ViewMode::Orthographic; requestResetDefaultView = true; });
    add("View", "Reset View", [this] { requestResetDefaultView = true; });
    const auto colourBy = [&](const char* label, AtomColorMode mode, bool enabled) {
        add("View > Colour Atoms By", label, [this, mode, &structure, &updateBuffers] {
            if (atomColorMode == mode) return;
            atomColorMode = mode;
            atomColorModeJustChanged = true;
            updateBuffers(structure);
        }, enabled);
    };
    colourBy("Element", AtomColorMode::ElementType, true);
    colourBy("Crystal Orientation", AtomColorMode::CrystalOrientation, atoms);
    colourBy("Grain Boundaries", AtomColorMode::GrainBoundary, atoms);
    const auto style = [&](const char* label, AtomDisplayMode mode) {
        add("View > Atom Style", label, [this, mode, &structure, &updateBuffers] {
            if (atomDisplayMode == mode) return;
            atomDisplayMode = mode;
            updateBuffers(structure);
        }, atoms);
    };
    style("Balls", AtomDisplayMode::Balls);
    style("Ball and Stick", AtomDisplayMode::BallAndStick);
    style("Space Filling", AtomDisplayMode::SpaceFilling);
    style("Polyhedral", AtomDisplayMode::Polyhedral);
    add("View", "Lattice Planes Settings", [this] { showLatticePlanesDialog = true; }, structure.hasUnitCell);
    add("View", "Miller Directions Settings", [this] { showMillerDirectionsDialog = true; }, structure.hasUnitCell);
    add("View", "Polyhedral Viewer", [this] { showPolyhedralViewer = !showPolyhedralViewer; }, atoms);
    add("View", "Structure Info", [this] { requestStructureInfo = true; }, atoms);
    add("View", "Atom Info", [this] { requestAtomInfo = true; }, atoms);
    add("View", "Measure Distance", [this] { requestMeasureDistance = true; }, atoms);
    add("View", "Measure Angle", [this] { requestMeasureAngle = true; }, atoms);
    add("View > Theme", "Light theme", [] { setThemeMode(ThemeMode::Light); });
    add("View > Theme", "Dark theme", [] { setThemeMode(ThemeMode::Dark); });
    add("View > Theme", "Match the system theme", [] { setThemeMode(ThemeMode::System); });
    for (const auto& [label, scale] : std::vector<std::pair<const char*, float>>{
             {"Interface size: Automatic", 1.0f}, {"Interface size: Compact (85%)", .85f}, {"Interface size: Larger (125%)", 1.25f},
             {"Interface size: Large (150%)", 1.5f}, {"Interface size: Extra large (200%)", 2.0f}})
        add("Settings", label, [scale = scale] { setInterfaceScale(scale); });
    add("Settings", "Atomic Sizes", [&editMenuDialogs] { editMenuDialogs.openAtomicSizes(); });
    add("Settings", "Display Settings", [&editMenuDialogs] { editMenuDialogs.openDisplaySettings(); });
    add("Settings", "Polyhedral Settings", [this] { showPolyhedralSettingsDialog = true; });
    add("Help", "Manual", [this] { showManual = true; });
    add("Help", "Keyboard Shortcuts", [this] { showShortcuts = true; });
    add("Help", "About", [this] { showAbout = true; });

    // The query and the ranked matches.
    if (ImGui::IsWindowAppearing()) ImGui::SetKeyboardFocusHere();
    if (uiScript::take("palette.query")) std::snprintf(commandQuery, sizeof(commandQuery), "%s", "slab");
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::InputTextWithHint("##query", "Type a command, e.g. surface slab, rdf, dark theme", commandQuery, sizeof(commandQuery)))
        commandSelected = 0;
    std::vector<std::string> words;
    {
        std::istringstream stream(lower(commandQuery));
        for (std::string word; stream >> word;) words.push_back(word);
    }
    std::vector<std::pair<int, const Command*>> matches;
    for (const auto& command : commands)
        if (const int s = score(command, words); s >= 0) matches.push_back({s, &command});
    std::stable_sort(matches.begin(), matches.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

    const int count = static_cast<int>(matches.size());
    if (ImGui::IsKeyPressed(ImGuiKey_DownArrow)) commandSelected = std::min(commandSelected + 1, count - 1);
    if (ImGui::IsKeyPressed(ImGuiKey_UpArrow)) commandSelected = std::max(commandSelected - 1, 0);
    commandSelected = std::clamp(commandSelected, 0, std::max(0, count - 1));

    const Command* chosen = nullptr;
    if (count == 0) {
        ImGui::Spacing();
        ImGui::TextDisabled("No command matches \"%s\".", commandQuery);
    } else {
        const float rowHeight = ImGui::GetTextLineHeightWithSpacing();
        ImGui::BeginChild("##matches", ImVec2(0, std::min(count, 10) * rowHeight), ImGuiChildFlags_None);
        for (int i = 0; i < count; ++i) {
            const Command& command = *matches[static_cast<std::size_t>(i)].second;
            ImGui::PushID(i);
            ImGui::BeginDisabled(!command.enabled);
            if (ImGui::Selectable("##row", i == commandSelected, ImGuiSelectableFlags_AllowOverlap)) chosen = &command;
            if (i == commandSelected && (ImGui::IsKeyPressed(ImGuiKey_DownArrow) || ImGui::IsKeyPressed(ImGuiKey_UpArrow)))
                ImGui::SetScrollHereY();
            ImGui::SameLine(ImGui::GetStyle().ItemSpacing.x);
            ImGui::TextUnformatted(command.label.c_str());
            const float pathWidth = ImGui::CalcTextSize(command.path.c_str()).x;
            ImGui::SameLine(std::max(ImGui::GetContentRegionMax().x - pathWidth, ImGui::GetCursorPosX() + responsive::dp(16)));
            ImGui::TextDisabled("%s", command.path.c_str());
            ImGui::EndDisabled();
            ImGui::PopID();
        }
        ImGui::EndChild();
        if (ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter))
            chosen = matches[static_cast<std::size_t>(commandSelected)].second;
    }
    ImGui::Spacing();
    ImGui::TextDisabled("Up/Down to choose, Enter to run, Esc to close");
    if (chosen && chosen->enabled) {
        ImGui::CloseCurrentPopup();
        chosen->run();
    }
    ImGui::EndPopup();
}
