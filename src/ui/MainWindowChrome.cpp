// The main window around the 3D view: the toolbar under the menu bar, the
// status bar at the bottom and the welcome screen shown when nothing is open.
#include "FileBrowser.h"

#include "ImGuiSetup.h"
#include "ui/DialogLayout.h"
#include "ui/Fonts.h"
#include "ui/ResponsiveLayout.h"
#include "util/ApplicationPaths.h"

#include "imgui.h"
#include "imgui_internal.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

namespace
{
// A toolbar button: icon and/or text, with a tooltip.
bool toolButton(const char* icon, const char* text, const char* tooltip, bool enabled = true, bool iconOnly = false)
{
    char label[128];
    if (uiFonts::icons && icon) std::snprintf(label, sizeof(label), iconOnly ? "%s##%s" : "%s  %s", icon, text);
    else std::snprintf(label, sizeof(label), "%s", text);
    ImGui::BeginDisabled(!enabled);
    const bool pressed = ImGui::Button(label);
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled | ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", tooltip);
    return pressed;
}

// A toolbar button for a mode: highlighted while the mode is on.
bool toggleButton(const char* text, const char* tooltip, bool on, bool enabled = true)
{
    if (on) {
        ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_HeaderActive));
        ImGui::PushStyleColor(ImGuiCol_Border, dialogLayout::accent());
    }
    const bool pressed = toolButton(nullptr, text, tooltip, enabled);
    if (on) ImGui::PopStyleColor(2);
    return pressed;
}

// A muted caption before a group of toolbar controls.
void groupCaption(const char* text)
{
    ImGui::AlignTextToFramePadding();
    ImGui::TextDisabled("%s", text);
    ImGui::SameLine();
}

// Starts the next toolbar group on this row when it fits, else on a new row,
// with a thin vertical rule between groups on the same row.
void nextGroup(float width)
{
    ImGui::SameLine();
    if (ImGui::GetContentRegionAvail().x < width) { ImGui::NewLine(); return; }
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float h = ImGui::GetFrameHeight();
    ImGui::GetWindowDrawList()->AddLine(ImVec2(p.x + responsive::dp(2), p.y + h * 0.18f), ImVec2(p.x + responsive::dp(2), p.y + h * 0.82f),
                                        ImGui::GetColorU32(ImGuiCol_Separator), 1.0f);
    ImGui::Dummy(ImVec2(responsive::dp(5), h));
    ImGui::SameLine();
}

const char* colourModeName(AtomColorMode mode)
{
    switch (mode) {
    case AtomColorMode::CrystalOrientation: return "Crystal orientation";
    case AtomColorMode::GrainBoundary: return "Grain boundaries";
    case AtomColorMode::AtomProperty: return "Atom property";
    case AtomColorMode::ElementType: default: return "Element";
    }
}

std::string fileNameOf(const std::string& path)
{
    const auto slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}
}

void FileBrowser::drawToolbar(const Structure& structure, bool canUndo, bool canRedo)
{
    const float width = ImGui::GetIO().DisplaySize.x;
    ImGui::SetNextWindowPos(ImVec2(0.0f, ImGui::GetFrameHeight()), ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0), ImVec2(width, FLT_MAX));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, responsive::size(10, 6));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_MenuBarBg));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                                   ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoScrollbar |
                                   ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus;
    if (ImGui::Begin("##ViewToolbar", nullptr, flags)) {
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, responsive::size(5, 5));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, responsive::size(9, 4));
        const bool hasAtoms = !structure.atoms.empty();
        const bool hasCell = structure.hasUnitCell && hasAtoms;
        const float unit = ImGui::GetFrameHeight();

        if (toolButton(Icon::Open, "Open", "Open a structure file (Ctrl+O)")) openStructurePopup = true;
        ImGui::SameLine();
        if (toolButton(Icon::Save, "Save", "Save the structure (Ctrl+S)", hasAtoms)) requestSaveFile = true;

        nextGroup(unit * 4);
        const std::string undoTip = (undoDescription.empty() ? std::string("Undo") : "Undo " + undoDescription) + " (Ctrl+Z)";
        const std::string redoTip = (redoDescription.empty() ? std::string("Redo") : "Redo " + redoDescription) + " (Ctrl+Y)";
        if (toolButton(Icon::Undo, "Undo", undoTip.c_str(), canUndo, true)) requestUndo = true;
        ImGui::SameLine();
        if (toolButton(Icon::Redo, "Redo", redoTip.c_str(), canRedo, true)) requestRedo = true;

        nextGroup(unit * 9.5f);
        groupCaption("View");
        if (toolButton(nullptr, "X##axis", "Look along x (1)")) requestViewAxisX = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "Y##axis", "Look along y (2)")) requestViewAxisY = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "Z##axis", "Look along z (3)")) requestViewAxisZ = true;
        ImGui::SameLine(0, responsive::dp(9));
        if (toolButton(nullptr, "a##latview", "Look along the lattice vector a (Shift+1)", hasCell)) requestViewLatticeA = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "b##latview", "Look along the lattice vector b (Shift+2)", hasCell)) requestViewLatticeB = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "c##latview", "Look along the lattice vector c (Shift+3)", hasCell)) requestViewLatticeC = true;
        ImGui::SameLine(0, responsive::dp(9));
        if (toolButton(Icon::Home, "Reset", "Reset the view (R)", true, true)) requestResetDefaultView = true;

        nextGroup(unit * 8);
        groupCaption("Rotate");
        ImGui::SetNextItemWidth(responsive::dp(62));
        ImGui::DragFloat("##rotangle", &rotateCrystalAngle, 1.0f, -360.0f, 360.0f, "%.1f\xc2\xb0");
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Rotation angle (degrees); drag or double-click to type");
        ImGui::SameLine();
        if (toolButton(nullptr, "X##rot", "Rotate the crystal about x", hasAtoms)) requestRotateCrystalX = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "Y##rot", "Rotate the crystal about y", hasAtoms)) requestRotateCrystalY = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "Z##rot", "Rotate the crystal about z", hasAtoms)) requestRotateCrystalZ = true;

        nextGroup(unit * 8);
        groupCaption("Measure");
        if (toolButton(Icon::Ruler, "Distance##measure", "Measure the distance between two atoms (pick them in the view)", hasAtoms)) requestMeasureDistance = true;
        ImGui::SameLine();
        if (toolButton(nullptr, "Angle##measure", "Measure the angle between three atoms (pick them in the view)", hasAtoms)) requestMeasureAngle = true;

        nextGroup(unit * 7);
        groupCaption("Select");
        if (toggleButton("Box##select", "Box selection: drag a rectangle in the view to select atoms (Ctrl adds)", boxSelectMode, hasAtoms)) {
            boxSelectMode = !boxSelectMode;
            if (boxSelectMode) lassoSelectMode = false;
        }
        ImGui::SameLine();
        if (toggleButton("Lasso##select", "Lasso selection: draw around atoms in the view to select them (Ctrl adds)", lassoSelectMode, hasAtoms)) {
            lassoSelectMode = !lassoSelectMode;
            if (lassoSelectMode) boxSelectMode = false;
        }

        // Right-aligned: commands, the structure pipeline, the projection and the theme.
        const char* commandsLabel = "Commands";
        const char* pipelineLabel = "Pipeline";
        const char* projectionLabel = viewMode == ViewMode::Orthographic ? "Orthographic" : "Perspective";
        const float rightWidth = ImGui::CalcTextSize(commandsLabel).x + ImGui::CalcTextSize(pipelineLabel).x + ImGui::CalcTextSize("Orthographic").x
                                 + ImGui::CalcTextSize("Light theme").x + 8 * ImGui::GetStyle().FramePadding.x + 4 * ImGui::GetStyle().ItemSpacing.x
                                 + (uiFonts::icons ? unit * 2.4f : 0.0f);
        ImGui::SameLine();
        const float space = ImGui::GetContentRegionAvail().x - rightWidth;
        if (space > 0) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + space);
        else ImGui::NewLine();
        if (toolButton(Icon::Search, commandsLabel, "Find any command by name (Ctrl+K)"))
            showCommandPalette = true;
        ImGui::SameLine();
        if (toolButton(Icon::Pipeline, pipelineLabel, "Edit > Structure pipeline: non-destructive editing steps", hasAtoms || pipelineDialog.isActive()))
            pipelineDialog.show();
        ImGui::SameLine();
        if (toolButton(nullptr, projectionLabel, "Switch between perspective and orthographic projection")) {
            viewMode = viewMode == ViewMode::Orthographic ? ViewMode::Isometric : ViewMode::Orthographic;
            requestResetDefaultView = true;
        }
        ImGui::SameLine();
        const bool light = isLightThemeEnabled();
        if (toolButton(nullptr, light ? "Dark theme" : "Light theme", "Switch the colour theme (Settings > Theme can also follow the system)"))
            setThemeMode(light ? ThemeMode::Dark : ThemeMode::Light);
        ImGui::PopStyleVar(2);
    }
    ImGui::End();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar(3);
}

void FileBrowser::drawStatusBar(const Structure& structure)
{
    const ImGuiIO& io = ImGui::GetIO();
    const float height = ImGui::GetFrameHeight() + responsive::dp(4);
    ImGui::SetNextWindowPos(ImVec2(0.0f, io.DisplaySize.y - height), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, height), ImGuiCond_Always);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, responsive::size(12, 2));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, responsive::size(6, 0));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_MenuBarBg));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
                                   ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus |
                                   ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav;
    if (ImGui::Begin("##StatusBar", nullptr, flags)) {
        // A rule along the top edge.
        const ImVec2 p = ImGui::GetWindowPos();
        ImGui::GetWindowDrawList()->AddLine(p, ImVec2(p.x + ImGui::GetWindowWidth(), p.y), ImGui::GetColorU32(ImGuiCol_Separator), 1.0f);
        // Fields separated by bullets.
        bool first = true;
        const auto field = [&first](const char* text, ImVec4 colour) {
            if (!first) { ImGui::SameLine(); ImGui::TextDisabled("\xe2\x80\xa2"); ImGui::SameLine(); }
            first = false;
            ImGui::AlignTextToFramePadding();
            ImGui::TextColored(colour, "%s", text);
        };
        const ImVec4 normal = ImGui::GetStyleColorVec4(ImGuiCol_Text), muted = ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);
        char text[256];
        if (largeDatasetShown) field("Large dataset (view only)", normal);
        else if (structure.atoms.empty()) field("No structure open", muted);
        else {
            const std::string name = currentStructurePath.empty() ? std::string("Untitled") : fileNameOf(currentStructurePath);
            field(name.c_str(), normal);
            std::snprintf(text, sizeof(text), "%zu atoms", structure.atoms.size());
            field(text, muted);
            if (viewSelectionAtoms.size() == 1 && viewSelectionAtoms[0] >= 0 &&
                static_cast<std::size_t>(viewSelectionAtoms[0]) < structure.atoms.size()) {
                // One atom: which, and where.
                const AtomSite& atom = structure.atoms[static_cast<std::size_t>(viewSelectionAtoms[0])];
                std::snprintf(text, sizeof(text), "Selected %s #%d at (%.3f, %.3f, %.3f) \xc3\x85", atom.symbol.c_str(), viewSelectionAtoms[0],
                              atom.x, atom.y, atom.z);
                field(text, dialogLayout::accent());
            } else if (!viewSelectionAtoms.empty()) {
                // Several: how many of each element (the three most common).
                std::string composition;
                if (viewSelectionAtoms.size() <= 200000) {
                    std::vector<std::pair<std::string, int>> counts;
                    for (const int index : viewSelectionAtoms) {
                        if (index < 0 || static_cast<std::size_t>(index) >= structure.atoms.size()) continue;
                        const std::string& symbol = structure.atoms[static_cast<std::size_t>(index)].symbol;
                        auto it = std::find_if(counts.begin(), counts.end(), [&](const auto& c) { return c.first == symbol; });
                        if (it == counts.end()) counts.push_back({symbol, 1});
                        else ++it->second;
                    }
                    std::sort(counts.begin(), counts.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
                    for (std::size_t i = 0; i < counts.size() && i < 3; ++i)
                        composition += (i ? ", " : " (") + counts[i].first + " " + std::to_string(counts[i].second);
                    if (!composition.empty()) composition += counts.size() > 3 ? ", ...)" : ")";
                }
                std::snprintf(text, sizeof(text), "%zu selected%s", viewSelectionAtoms.size(), composition.c_str());
                field(text, dialogLayout::accent());
            }
            if (structure.hasUnitCell) {
                double length[3];
                for (int r = 0; r < 3; ++r)
                    length[r] = std::sqrt(structure.cellVectors[r][0] * structure.cellVectors[r][0] + structure.cellVectors[r][1] * structure.cellVectors[r][1] +
                                          structure.cellVectors[r][2] * structure.cellVectors[r][2]);
                std::snprintf(text, sizeof(text), "Cell %.3f \xc3\x97 %.3f \xc3\x97 %.3f \xc3\x85", length[0], length[1], length[2]);
                field(text, muted);
            } else field("No periodic cell", muted);
        }
        ImGui::SameLine();
        // Right side: projection and colouring, each clickable (projection
        // switches, colouring opens its choices).
        const char* projection = viewMode == ViewMode::Orthographic ? "Orthographic" : "Perspective";
        char colour[64];
        std::snprintf(colour, sizeof(colour), "Colour: %s", colourModeName(atomColorMode));
        const float gap = ImGui::CalcTextSize("   \xe2\x80\xa2   ").x;
        const float right = ImGui::CalcTextSize(projection).x + gap + ImGui::CalcTextSize(colour).x;
        ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() - right - responsive::dp(12)));
        const auto link = [](const char* label, const char* tooltip) {
            ImGui::AlignTextToFramePadding();
            ImGui::TextDisabled("%s", label);
            const bool hovered = ImGui::IsItemHovered();
            if (hovered) {
                const ImVec2 a = ImGui::GetItemRectMin(), b = ImGui::GetItemRectMax();
                ImGui::GetWindowDrawList()->AddLine(ImVec2(a.x, b.y - 1), ImVec2(b.x, b.y - 1), ImGui::GetColorU32(ImGuiCol_TextDisabled));
                ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
                ImGui::SetTooltip("%s", tooltip);
            }
            return hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left);
        };
        if (link(projection, "Switch between perspective and orthographic projection")) {
            viewMode = viewMode == ViewMode::Orthographic ? ViewMode::Isometric : ViewMode::Orthographic;
            requestResetDefaultView = true;
        }
        ImGui::SameLine(0, 0);
        ImGui::TextDisabled("   \xe2\x80\xa2   ");
        ImGui::SameLine(0, 0);
        if (link(colour, "Choose how atoms are coloured")) ImGui::OpenPopup("##statusColour");
        if (ImGui::BeginPopup("##statusColour")) {
            const bool atoms = !structure.atoms.empty();
            if (ImGui::MenuItem("Element", nullptr, atomColorMode == AtomColorMode::ElementType)) pendingColourMode = static_cast<int>(AtomColorMode::ElementType);
            if (ImGui::MenuItem("Crystal Orientation", nullptr, atomColorMode == AtomColorMode::CrystalOrientation, atoms))
                pendingColourMode = static_cast<int>(AtomColorMode::CrystalOrientation);
            if (ImGui::MenuItem("Grain Boundaries", nullptr, atomColorMode == AtomColorMode::GrainBoundary, atoms))
                pendingColourMode = static_cast<int>(AtomColorMode::GrainBoundary);
            ImGui::EndPopup();
        }
    }
    ImGui::End();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar(4);
}

void FileBrowser::drawWelcome(const Structure& structure)
{
    if (!structure.atoms.empty() || largeDatasetShown) return;
    const ImGuiIO& io = ImGui::GetIO();
    const float width = std::min(responsive::dp(560), io.DisplaySize.x - responsive::dp(40));
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.52f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0), ImVec2(width, io.DisplaySize.y * 0.8f));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, responsive::size(28, 24));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, responsive::dp(12));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_ChildBg));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
                                   ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus |
                                   ImGuiWindowFlags_NoFocusOnAppearing;
    if (ImGui::Begin("##Welcome", nullptr, flags)) {
        uiFonts::pushHeading(1.9f);
        ImGui::TextUnformatted("AtomForge");
        uiFonts::popHeading();
        ImGui::SameLine();
        ImGui::AlignTextToFramePadding();
#ifdef ATOMFORGE_VERSION
        ImGui::TextDisabled("Version %s", ATOMFORGE_VERSION);
#endif
        ImGui::TextDisabled("Build, edit, analyse and render atomic structures.");
        ImGui::Spacing();
        ImGui::Spacing();

        const float half = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
        const ImVec2 big(half, ImGui::GetFrameHeight() * 1.6f);
        char label[96];
        if (dialogLayout::primaryButton(uiFonts::withIcon(Icon::Open, "Open a structure", label, sizeof(label)), big)) openStructurePopup = true;
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("CIF, POSCAR, XYZ, PDB, LAMMPS and more (Ctrl+O)");
        ImGui::SameLine();
        if (dialogLayout::secondaryButton(uiFonts::withIcon(Icon::Cube, "Build a crystal", label, sizeof(label)), big)) bulkCrystalDialog.open();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Build > Bulk Crystal: a crystal from its space group and lattice");
        if (dialogLayout::secondaryButton(uiFonts::withIcon(Icon::Folder, "Open a project", label, sizeof(label)), big)) {
            projectAction = 1;
            projectPicker.open("Open AtomForge project", false, "project.afproject");
        }
        ImGui::SameLine();
        if (dialogLayout::secondaryButton(uiFonts::withIcon(Icon::Help, "Read the manual", label, sizeof(label)), big)) showManual = true;

        // Example structures, when installed (or in the source tree).
        static const std::filesystem::path examples = findExamplesDirectory(applicationDirectory());
        if (!examples.empty()) {
            struct Example { const char* file; const char* label; const char* about; };
            static const Example kExamples[] = {
                {"cu_fcc.cif", "Copper (fcc)", "The four-atom conventional cell of copper"},
                {"fe_bcc.cif", "Iron (bcc)", "The two-atom conventional cell of iron"},
                {"cu_ni_alloy.cif", "Cu-Ni alloy", "A 500-atom random Cu70Ni30 solid solution"},
                {"cu_sigma5.cif", "Grain boundary", "A copper Sigma-5 [001] tilt grain boundary"},
                {"cu_polycrystal.cif", "Polycrystal", "A four-grain copper polycrystal"},
                {"cu_dislocation.cif", "Dislocation", "A copper crystal containing an edge dislocation"},
                {"cu_sphere.xyz", "Nanoparticle", "A spherical copper nanoparticle"},
                {"sio2_pack.xyz", "Amorphous SiO2", "A 120-atom amorphous silica packing"}};
            dialogLayout::section("Examples");
            const float right = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;
            bool first = true;
            for (const Example& example : kExamples) {
                const auto path = examples / example.file;
                std::error_code error;
                if (!std::filesystem::is_regular_file(path, error)) continue;
                const float w = ImGui::CalcTextSize(example.label).x + 2 * ImGui::GetStyle().FramePadding.x;
                if (!first) {
                    ImGui::SameLine();
                    if (ImGui::GetCursorPosX() + w > right) ImGui::NewLine();
                }
                first = false;
                if (dialogLayout::secondaryButton(example.label)) {
                    pendingOpenPath = path.generic_u8string();
                    openStatusMsg[0] = '\0';
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s\n%s", example.about, example.file);
            }
        }

        const auto recent = recentFiles();
        if (!recent.empty()) {
            dialogLayout::section("Recent files");
            // Whether each file still exists, checked every few seconds (slow
            // network drives must not stall the frame).
            static std::vector<std::pair<std::string, bool>> present;
            static double checked = -10.0;
            if (ImGui::GetTime() - checked > 5.0) {
                checked = ImGui::GetTime();
                present.clear();
                for (std::size_t i = 0; i < recent.size() && i < 6; ++i) {
                    std::error_code error;
                    present.push_back({recent[i], std::filesystem::is_regular_file(std::filesystem::u8path(recent[i]), error)});
                }
            }
            int shown = 0;
            for (const std::string& path : recent) {
                if (++shown > 6) break;
                ImGui::PushID(path.c_str());
                const std::string name = fileNameOf(path);
                const auto found = std::find_if(present.begin(), present.end(), [&](const auto& p) { return p.first == path; });
                const bool missing = found != present.end() && !found->second;
                ImGui::BeginDisabled(missing);
                if (ImGui::Selectable(name.c_str(), false, ImGuiSelectableFlags_None, ImVec2(0, ImGui::GetFrameHeight()))) {
                    pendingOpenPath = path;
                    openStatusMsg[0] = '\0';
                }
                ImGui::EndDisabled();
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort | ImGuiHoveredFlags_AllowWhenDisabled))
                    ImGui::SetTooltip("%s%s", path.c_str(), missing ? "\n(the file is no longer there)" : "");
                // The folder, shortened from the left to fit.
                const float column = ImGui::GetContentRegionAvail().x * 0.42f;
                ImGui::SameLine(column);
                const float room = ImGui::GetContentRegionAvail().x;
                std::string folderText = path.substr(0, path.size() - name.size());
                if (ImGui::CalcTextSize(folderText.c_str()).x > room) {
                    while (folderText.size() > 4 && ImGui::CalcTextSize(("\xe2\x80\xa6" + folderText).c_str()).x > room) folderText.erase(0, 1);
                    folderText = "\xe2\x80\xa6" + folderText;
                }
                ImGui::AlignTextToFramePadding();
                ImGui::TextDisabled("%s", folderText.c_str());
                ImGui::PopID();
            }
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("Drop structure files here. Dumps of billions of atoms open with the large-data viewer.");
    }
    ImGui::End();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar(2);
}

// Help > Keyboard Shortcuts (F1): every key and mouse binding of the main view.
void FileBrowser::drawShortcuts()
{
    if (showShortcuts) {
        ImGui::OpenPopup("Keyboard Shortcuts");
        showShortcuts = false;
    }
    ImGui::SetNextWindowSizeConstraints(ImVec2(responsive::dp(560), 0), ImVec2(responsive::dp(560), FLT_MAX));
    bool open = true;
    if (!responsive::beginModal("Keyboard Shortcuts", &open, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize)) return;
    struct Row { const char* keys; const char* action; };
    const auto group = [](const char* title, std::initializer_list<Row> rows) {
        dialogLayout::section(title);
        if (!ImGui::BeginTable(title, 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) return;
        ImGui::TableSetupColumn("keys", ImGuiTableColumnFlags_WidthStretch, 0.38f);
        ImGui::TableSetupColumn("action", ImGuiTableColumnFlags_WidthStretch, 0.62f);
        for (const Row& row : rows) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(row.keys);
            ImGui::TableSetColumnIndex(1);
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
            ImGui::TextWrapped("%s", row.action);
            ImGui::PopStyleColor();
        }
        ImGui::EndTable();
    };
    group("Files", {{"Ctrl+O", "Open a structure"},
                    {"Ctrl+S", "Save"},
                    {"Ctrl+Shift+S", "Save as"},
                    {"Ctrl+Alt+S", "Export an image"},
                    {"Ctrl+W", "Close the structure"}});
    group("Editing", {{"Ctrl+Z", "Undo"},
                      {"Ctrl+Y, Ctrl+Shift+Z", "Redo"},
                      {"Delete", "Delete the selected atoms"},
                      {"S", "Substitute the selected atoms"},
                      {"I", "Insert an atom at the midpoint of the selection"},
                      {"G", "Move the selected atoms; X, Y or Z locks an axis, click or Enter confirms, Esc cancels"}});
    group("Selection", {{"Hover", "Name the atom under the mouse (element, index, position)"},
                        {"Click", "Select an atom"},
                        {"Ctrl+click", "Add or remove an atom"},
                        {"Ctrl+A", "Select all atoms"},
                        {"Esc, Ctrl+D", "Clear the selection"}});
    group("View", {{"Left drag", "Orbit"},
                   {"Right drag", "Pan (when no box or lasso selection is active)"},
                   {"Wheel", "Zoom"},
                   {"R", "Reset the view"},
                   {"1, 2, 3", "Look along x, y, z"},
                   {"Shift+1, 2, 3", "Look along the lattice vectors a, b, c"},
                   {"F1", "This list"},
                   {"Ctrl+K", "Find any command by name"}});
    group("Dialogs", {{"Esc", "Close the dialog (in a field: cancel the edit)"},
                      {"Tab, Shift+Tab", "Move between fields"}});
    ImGui::Spacing();
    if (dialogLayout::footer("Close", nullptr) == 1) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
}
