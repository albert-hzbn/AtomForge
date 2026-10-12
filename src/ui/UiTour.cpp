// ATOMFORGE_UI_TOUR=FOLDER: screenshots of the main window (welcome screen,
// a structure, every menu) and of the main dialogs in both themes, written to
// FOLDER as numbered PNGs, after which the application quits. Used to review
// the interface as a whole.
#include "FileBrowser.h"
#include "ImGuiSetup.h"
#include "ui/UiScript.h"
#include "ui/MenuParity.h"
#include "science/ScienceCatalog.h"

#include "third_party/stb_image_write.h"

#include "imgui.h"
#include "imgui_internal.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

struct FileBrowser::UiTour
{
    std::filesystem::path folder;
    std::size_t index = 0;
    int frame = 0;
    using Update = std::function<void(Structure&)>;
    // Steps get the current frame's structure, dialogs and update callback.
    struct Step { std::string name; std::function<void(Structure&, EditMenuDialogs&, const Update&)> setup; };
    std::vector<Step> steps;
};

namespace
{
void screenshotFrontBuffer(const std::filesystem::path& file)
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    const int w = viewport[2], h = viewport[3];
    if (w <= 0 || h <= 0) return;
    std::vector<unsigned char> pixels(static_cast<std::size_t>(w) * h * 3), flipped(pixels.size());
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    glReadBuffer(GL_BACK);
    for (int row = 0; row < h; ++row)
        std::copy_n(pixels.begin() + static_cast<std::ptrdiff_t>(row) * w * 3, w * 3, flipped.begin() + static_cast<std::ptrdiff_t>(h - 1 - row) * w * 3);
    stbi_write_png(file.u8string().c_str(), w, h, 3, flipped.data(), w * 3);
}

Structure tourCrystal()
{
    // 4 x 4 x 4 cells of an fcc Cu-Ni ordered alloy.
    Structure s;
    const double a = 3.6;
    const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                for (int b = 0; b < 4; ++b) {
                    AtomSite atom{};
                    atom.symbol = b == 0 ? "Ni" : "Cu";
                    atom.atomicNumber = b == 0 ? 28 : 29;
                    atom.x = (i + basis[b][0]) * a; atom.y = (j + basis[b][1]) * a; atom.z = (k + basis[b][2]) * a;
                    s.atoms.push_back(atom);
                }
    s.hasUnitCell = true;
    s.cellVectors = {{{4 * a, 0, 0}, {0, 4 * a, 0}, {0, 0, 4 * a}}};
    return s;
}
}

void FileBrowser::runUiTour(Structure& current, EditMenuDialogs& dialogs, const std::function<void(Structure&)>& update)
{
    if (!uiTour) {
        static bool started = false;
        const char* folder = std::getenv("ATOMFORGE_UI_TOUR");
        if (started || !folder || !*folder) return;
        started = true;
        uiTour = std::make_shared<UiTour>();
        auto& t = *uiTour;
        t.folder = std::filesystem::u8path(folder);
        std::filesystem::create_directories(t.folder);
        const auto load = [this](Structure& structure, EditMenuDialogs&, const UiTour::Update& updateBuffers) {
            if (structure.atoms.empty()) { structure = tourCrystal(); currentStructurePath = "CuNi-L12.cif"; updateBuffers(structure); requestResetDefaultView = true; }
        };
        // ATOMFORGE_UI_TOUR_DARK=1 takes the whole tour in the dark theme.
        const bool allDark = std::getenv("ATOMFORGE_UI_TOUR_DARK") && *std::getenv("ATOMFORGE_UI_TOUR_DARK") == '1';
        const auto theme = [this, allDark](bool light) { light = light && !allDark; if (light) applyLightTheme(); else applyDarkTheme(); };
        t.steps = {
            {"welcome", [theme](Structure&, EditMenuDialogs&, const UiTour::Update&) { theme(true); }},
            {"main-window", load},
        };
        for (const char* menu : {"File", "Edit", "Build", "View", "Analysis", "Settings", "Help"})
            t.steps.push_back({std::string("menu-") + menu, [this, menu](Structure&, EditMenuDialogs&, const UiTour::Update&) { tourMenu = menu; }});
        t.steps.push_back({"pipeline-panel", [this](Structure& structure, EditMenuDialogs&, const UiTour::Update&) {
            pipelineDialog.restore(R"({"format":"atomforge-pipeline","version":1,"active":true,"open":true,"modifiers":[
                {"type":"replicate","parameters":{"counts":[2,1,1]}},
                {"type":"select-expression","parameters":{"expression":"fz > 0.5"}},
                {"type":"assign-element","parameters":{"element":"Al","target":"selected"}},
                {"type":"wrap","enabled":false}]})", &structure);
            pipelineDialog.show();
        }});
        t.steps.push_back({"operation-replicate", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { operationDialog.open("replicate"); }});
        t.steps.push_back({"operation-select", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { operationDialog.open("select-expression"); }});
        t.steps.push_back({"builder-bulk", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { bulkCrystalDialog.open(); }});
        t.steps.push_back({"open-file", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { openStructurePopup = true; }});
        t.steps.push_back({"export-image", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            exportImagePopup = true;
            exportDir = openDir;
            exportDirHistory = dirHistory;
            exportHistoryIndex = historyIndex;
            if (exportFilename[0] == '\0') std::snprintf(exportFilename, sizeof(exportFilename), "%s", "structure.png");
        }});
        t.steps.push_back({"display-settings", [](Structure&, EditMenuDialogs& editMenuDialogs, const UiTour::Update&) { editMenuDialogs.openDisplaySettings(); }});
        t.steps.push_back({"atomic-sizes", [](Structure&, EditMenuDialogs& editMenuDialogs, const UiTour::Update&) { editMenuDialogs.openAtomicSizes(); }});
        t.steps.push_back({"edit-structure", [](Structure&, EditMenuDialogs& editMenuDialogs, const UiTour::Update&) { editMenuDialogs.openEditStructure(); }});
        t.steps.push_back({"structure-info", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { requestStructureInfo = true; }});
        t.steps.push_back({"science-atoms", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { scientificToolsDialog.open("structure-type"); }});
        t.steps.push_back({"science-atoms-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            scientificToolsDialog.open("structure-type");
            uiScript::request("science.run");
        }});
        t.steps.push_back({"science-plot-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            scientificToolsDialog.open("powder-xrd");
            uiScript::request("science.run");
        }});
        t.steps.push_back({"science-error-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            scientificToolsDialog.open("band-gap");
            uiScript::request("science.run");
        }});
        t.steps.push_back({"science-plot", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { scientificToolsDialog.open("powder-xrd"); }});
        t.steps.push_back({"science-simulation", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { scientificToolsDialog.open("relax"); }});
        t.steps.push_back({"analysis-cna-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { cnaDialog.open(); uiScript::request("cna.run"); }});
        t.steps.push_back({"analysis-cna", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { cnaDialog.open(); }});
        t.steps.push_back({"analysis-rdf", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { rdfDialog.open(); }});
        t.steps.push_back({"analysis-adf-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { angularDistributionDialog.open(); uiScript::request("adf.run"); }});
        t.steps.push_back({"analysis-adf", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { angularDistributionDialog.open(); }});
        t.steps.push_back({"analysis-sro-result", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { shortRangeOrderDialog.openRequested = true; uiScript::request("sro.run"); }});
        t.steps.push_back({"analysis-sro", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { shortRangeOrderDialog.openRequested = true; }});
        t.steps.push_back({"analysis-dislocations", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { dislocationAnalysisDialog.open(); }});
        t.steps.push_back({"trajectory-loaded", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            // A small generated trajectory: 4 atoms moving over 20 frames.
            const auto file = std::filesystem::temp_directory_path() / "atomforge-tour-trajectory.xyz";
            std::ofstream out(file);
            for (int f = 0; f < 20; ++f) {
                out << "4\nLattice=\"4 0 0 0 4 0 0 0 4\" Properties=species:S:1:pos:R:3\n";
                for (int a = 0; a < 4; ++a) out << (a % 2 ? "Ni " : "Cu ") << a + 0.05 * f << " " << (a % 2) * 2.0 << " 0\n";
            }
            out.close();
            trajectoryDialog.show();
            trajectoryDialog.load(file.u8string());
        }});
        t.steps.push_back({"trajectory", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { trajectoryDialog.show(); }});
        t.steps.push_back({"electronic", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { electronicDialog.open(); }});
        t.steps.push_back({"wannier", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { wannierDialog.open(); }});
        t.steps.push_back({"lobster", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { lobsterDialog.open(); }});
        t.steps.push_back({"lattice-planes", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showLatticePlanesDialog = true; }});
        t.steps.push_back({"miller-directions", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showMillerDirectionsDialog = true; }});
        t.steps.push_back({"polyhedral-settings", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showPolyhedralSettingsDialog = true; }});
        t.steps.push_back({"about", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showAbout = true; }});
        t.steps.push_back({"unsaved-close", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { closeStructure(); }});
        t.steps.push_back({"unsaved-quit", [](Structure&, EditMenuDialogs&, const UiTour::Update&) { uiScript::request("app.quit"); }});
        t.steps.push_back({"selection-one", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { requestedSelection = {100}; selectionRequested = true; }});
        t.steps.push_back({"selection-many", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            requestedSelection.clear();
            for (int i = 0; i <= 20; ++i) requestedSelection.push_back(i);
            selectionRequested = true;
        }});
        t.steps.push_back({"measure-distance", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            requestedSelection = {0, 1};
            selectionRequested = true;
            requestMeasureDistance = true;
        }});
        t.steps.push_back({"measure-angle", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            requestedSelection = {0, 1, 2};
            selectionRequested = true;
            requestMeasureAngle = true;
        }});
        t.steps.push_back({"menu-history", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { tourMenu = "Edit"; tourSubmenu = "History"; }});
        t.steps.push_back({"menu-open-recent", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { tourMenu = "File"; tourSubmenu = "Open Recent"; }});
        t.steps.push_back({"overlay-polyhedra", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showPolyhedralViewer = true; }});
        t.steps.push_back({"overlay-voronoi", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showPolyhedralViewer = false; showVoronoi = true; }});
        t.steps.push_back({"overlay-off", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showVoronoi = false; }});
        t.steps.push_back({"hover-atom", [](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            ImGui::GetIO().AddMousePosEvent(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f);
        }});
        t.steps.push_back({"structure-info-small", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { requestStructureInfo = true; }});
        t.steps.push_back({"periodic-table", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            requestedSelection = {100};
            selectionRequested = true;
            uiScript::request("atom.substitute");
            ImGui::GetIO().AddInputCharactersUTF8("cu");  // typed symbol picks Cu
        }});
        t.steps.push_back({"command-palette", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showCommandPalette = true; }});
        t.steps.push_back({"command-palette-search", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            showCommandPalette = true;
            uiScript::request("palette.query");
        }});
        t.steps.push_back({"shortcuts", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { showShortcuts = true; }});
        t.steps.push_back({"notifications", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            showNotification("Saved cu_host.cif (500 atoms).");
            showNotification("Error: failed to save (format may not support this structure).", true);
        }});
        // ATOMFORGE_UI_TOUR_ALL=1 adds every remaining dialog: the other
        // builders, every operation and every scientific tool.
        if (const char* all = std::getenv("ATOMFORGE_UI_TOUR_ALL"); all && *all == '1') {
            t.steps.push_back({"builder-interface", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { interfaceBuilderDialog.open(); }});
            t.steps.push_back({"builder-stacking-fault", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { stackingFaultDialog.open(); }});
            t.steps.push_back({"builder-merge", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { mergeStructuresDialog.open(); }});
            for (const MenuLocation& location : menuLocations())
                if (!location.dialog)
                    t.steps.push_back({std::string("operation-") + location.step, [this, step = std::string(location.step)](Structure&, EditMenuDialogs&, const UiTour::Update&) { operationDialog.open(step); }});
            for (const ScienceToolDef& tool : scienceToolCatalog())
                t.steps.push_back({std::string("science-") + tool.id, [this, id = std::string(tool.id)](Structure&, EditMenuDialogs&, const UiTour::Update&) { scientificToolsDialog.open(id); }});
        }
        t.steps.push_back({"dark-main-window", [theme](Structure&, EditMenuDialogs&, const UiTour::Update&) { theme(false); }});
        t.steps.push_back({"dark-menu-View", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { tourMenu = "View"; }});
        t.steps.push_back({"dark-pipeline-panel", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { pipelineDialog.show(); }});
        t.steps.push_back({"dark-builder-bulk", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { bulkCrystalDialog.open(); }});
        t.steps.push_back({"dark-operation-replicate", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { operationDialog.open("replicate"); }});
        t.steps.push_back({"dark-science-plot", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) { scientificToolsDialog.open("powder-xrd"); }});
        t.steps.push_back({"dark-notifications", [this](Structure&, EditMenuDialogs&, const UiTour::Update&) {
            showNotification("Saved cu_host.cif (500 atoms).");
            showNotification("Error: failed to save (format may not support this structure).", true);
        }});
        t.steps.push_back({"dark-welcome", [this](Structure& structure, EditMenuDialogs&, const UiTour::Update& updateBuffers) {
            pipelineDialog.restore(R"({"format":"atomforge-pipeline","version":1,"active":false,"open":false,"modifiers":[]})", nullptr);
            structure = Structure();
            currentStructurePath.clear();
            updateBuffers(structure);
        }});
    }
    auto& t = *uiTour;
    if (t.index >= t.steps.size()) {
        applyLightTheme();
        glfwSetWindowShouldClose(glfwGetCurrentContext(), true);
        uiTour.reset();
        return;
    }
    if (t.frame == 0) t.steps[t.index].setup(current, dialogs, update);
    const std::string& stepName = t.steps[t.index].name;
    if (++t.frame < (stepName == "hover-atom" || stepName.find("-result") != std::string::npos ? 60 : 12)) return;  // hover tooltips wait for the mouse to rest
    char name[256];
    std::snprintf(name, sizeof(name), "%02zu-%s.png", t.index + 1, t.steps[t.index].name.c_str());
    screenshotFrontBuffer(t.folder / name);
    // Close what the step opened before the next one.
    tourMenu = nullptr;
    if (!ImGui::GetCurrentContext()->OpenPopupStack.empty()) ImGui::ClosePopupToLevel(0, true);
    operationDialog.close();
    scientificToolsDialog.close();
    pipelineDialog.hide();
    showAbout = false;
    trajectoryDialog.hide();
    electronicDialog.close();
    wannierDialog.close();
    lobsterDialog.close();
    dislocationAnalysisDialog.close();
    ++t.index;
    t.frame = 0;
}

void FileBrowser::forceTourMenu(const char* menu)
{
    if (tourMenu && std::string(tourMenu) == menu && !ImGui::IsPopupOpen(menu)) ImGui::OpenPopup(menu);
}

void FileBrowser::forceTourSubmenu(const char* menu)
{
    if (tourSubmenu && std::string(tourSubmenu) == menu && !ImGui::IsPopupOpen(menu)) ImGui::OpenPopup(menu);
}
