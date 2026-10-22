// ATOMFORGE_TUTORIAL_CAPTURE=FOLDER: recaptures the screenshots of the
// manual's illustrated tutorials (docs/manual/tutorials.json) in the running
// application. Each capture sets up a state (structures from the manual's
// examples, an open dialog, scripted actions through ui/UiScript.h), waits for
// it to settle, saves a full-window JPEG and records the rectangle of the
// window it shows; FOLDER/bounds-*.json then hold every image's crop bounds
// (docs/manual/scripts/update_tutorial_captures.py applies them).
#include "FileBrowser.h"
#include "ImGuiSetup.h"
#include "ui/UiScript.h"
#include "ui/TrajectoryLinesOverlay.h"

#include "io/StructureLoader.h"
#include "graphics/FrameCapture.h"
#include "graphics/GraphicsBackend.h"
#include "third_party/stb_image_write.h"

#include "imgui.h"
#include "imgui_internal.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace
{
using Update = std::function<void(Structure&)>;
struct Context { Structure& structure; EditMenuDialogs& dialogs; const Update& update; };
using Action = std::function<void(Context&)>;

ImRect popupRect(bool topOnly)
{
    ImGuiContext& g = *ImGui::GetCurrentContext();
    ImRect r(FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX);
    for (int i = topOnly ? g.OpenPopupStack.Size - 1 : 0; i >= 0 && i < g.OpenPopupStack.Size; ++i)
        if (ImGuiWindow* w = g.OpenPopupStack[i].Window; w && w->WasActive) r.Add(w->Rect());
    return r;
}

ImRect windowRect(const char* name)
{
    if (ImGuiWindow* w = ImGui::FindWindowByName(name); w && w->WasActive) return w->Rect();
    return ImRect(FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX);
}

// The 3D view between the toolbar (and tab bar) and the status bar.
ImRect viewportRect()
{
    const ImGuiIO& io = ImGui::GetIO();
    float top = ImGui::GetFrameHeight();
    for (const char* name : {"##ViewToolbar", "##StructureTabBar"})
        if (ImGuiWindow* w = ImGui::FindWindowByName(name); w && w->WasActive) top = std::max(top, w->Pos.y + w->Size.y);
    return ImRect(0, top, io.DisplaySize.x, io.DisplaySize.y - ImGui::GetFrameHeight() - 4);
}

void saveJpeg(const std::filesystem::path& file)
{
    std::vector<unsigned char> pixels;
    int w = 0, h = 0;
    if (!graphics::lastFrame(pixels, w, h)) return;
    std::filesystem::create_directories(file.parent_path());
    if (file.extension() == ".png") stbi_write_png(file.u8string().c_str(), w, h, 3, pixels.data(), w * 3);
    else stbi_write_jpg(file.u8string().c_str(), w, h, 3, pixels.data(), 92);
}
}

struct FileBrowser::TutorialCapture
{
    struct Shot
    {
        std::string file;
        std::vector<std::pair<int, Action>> actions;  // frame -> action
        int settle = 14;                              // frames after the last action
        std::function<ImRect()> target;
    };
    std::filesystem::path folder, examples, chgcar;
    std::vector<Shot> shots;
    std::size_t index = 0;
    int frame = 0;
    std::string bounds = "{", first;
};

void FileBrowser::runTutorialCapture(Structure& current, EditMenuDialogs& dialogsRef, const std::function<void(Structure&)>& update)
{
    if (!tutorialCapture) {
        static bool started = false;
        const char* folder = std::getenv("ATOMFORGE_TUTORIAL_CAPTURE");
        if (started || !folder || !*folder) return;
        started = true;
        tutorialCapture = std::make_shared<TutorialCapture>();
        auto& c = *tutorialCapture;
        c.folder = std::filesystem::u8path(folder);
#ifdef ATOMFORGE_SOURCE_MANUAL
        c.examples = std::filesystem::u8path(ATOMFORGE_SOURCE_MANUAL).parent_path() / "examples";
        c.chgcar = std::filesystem::u8path(ATOMFORGE_SOURCE_MANUAL).parent_path().parent_path().parent_path() / "tmp" / "CHGCAR-electronic-demo";
#endif
        if (const char* e = std::getenv("ATOMFORGE_TUTORIAL_EXAMPLES")) c.examples = std::filesystem::u8path(e);
        if (const char* e = std::getenv("ATOMFORGE_TUTORIAL_CHGCAR")) c.chgcar = std::filesystem::u8path(e);
        const std::string host = (c.examples / "cu_host.cif").u8string(), alloy = (c.examples / "cu_ni_alloy.cif").u8string(),
                          cu4 = (c.examples / "cu_fcc.cif").u8string(), octahedron = (c.examples / "octahedron.obj").u8string();
        // The main structure for a capture.
        const auto main = [this](const std::string& path) -> Action {
            return [this, path](Context& x) {
                Structure s;
                std::string error;
                if (loadStructureFromFile(path, s, error)) {
                    x.structure = std::move(s);
                    currentStructurePath = path;
                    x.update(x.structure);
                    requestResetDefaultView = true;
                }
            };
        };
        const auto script = [](const char* action) -> Action { return [action](Context&) { uiScript::request(action); }; };
        // Empties the main view, so a builder's result replaces it instead of
        // opening a new tab (whose own FileBrowser would not run this capture).
        const Action clearMain = [this](Context& x) {
            x.structure = Structure();
            currentStructurePath.clear();
            x.update(x.structure);
        };
        // Selects the atom nearest the centre (it has no periodic images).
        const Action selectCentral = [this](Context& x) {
            const auto& atoms = x.structure.atoms;
            if (atoms.empty()) return;
            double cx = 0, cy = 0, cz = 0;
            for (const AtomSite& a : atoms) { cx += a.x; cy += a.y; cz += a.z; }
            cx /= atoms.size(); cy /= atoms.size(); cz /= atoms.size();
            int best = 0;
            double bestDistance = 1e300;
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const double d = (atoms[i].x - cx) * (atoms[i].x - cx) + (atoms[i].y - cy) * (atoms[i].y - cy) + (atoms[i].z - cz) * (atoms[i].z - cz);
                if (d < bestDistance) { bestDistance = d; best = static_cast<int>(i); }
            }
            requestedSelection = {best};
            selectionRequested = true;
        };
        const auto popups = [] { return popupRect(false); };
        // Menus with the menu bar they drop from.
        const auto menus = [] { ImRect r = popupRect(false); if (r.Min.y <= r.Max.y) { r.Min.y = 0; r.Min.x = std::max(0.0f, r.Min.x - 40.0f); } return r; };
        const auto topPopup = [] { return popupRect(true); };
        const auto menu = [this](const char* name, const char* sub) -> Action {
            return [this, name, sub](Context&) { tourMenu = name; tourSubmenu = sub; };
        };
        auto& s = c.shots;
        using Shot = TutorialCapture::Shot;
        // Menus.
        s.push_back(Shot{"tutorials/view-menu.jpg", {{0, main(host)}, {4, menu("View", nullptr)}}, 14, menus});
        s.push_back(Shot{"tutorials/atom-representation.jpg", {{0, menu("View", nullptr)}, {6, menu("View", "Atom Style")}}, 14, menus});
        s.push_back(Shot{"tutorials/structure-colouring.jpg", {{0, menu("View", nullptr)}, {6, menu("View", "Colour Atoms By")}}, 14, menus});
        s.push_back(Shot{"tutorials/theme-options.jpg", {{0, menu("View", nullptr)}, {6, menu("View", "Theme")}}, 14, menus});
        s.push_back(Shot{"tutorials/edit-menu.jpg", {{0, menu("Edit", nullptr)}}, 14, menus});
        s.push_back(Shot{"tutorials/build-menu.jpg", {{0, menu("Build", nullptr)}}, 14, menus});
        s.push_back(Shot{"tutorials/interface-size.jpg", {{0, menu("Settings", nullptr)}, {6, menu("Settings", "Interface Size")}}, 14, menus});
        // Builders.
        s.push_back(Shot{"tutorials/bulk-controls.jpg", {{0, [this](Context&) { bulkCrystalDialog.open(); }}}, 14, popups});
        s.push_back(Shot{"bulk-dialog.png", {{0, [this](Context&) { bulkCrystalDialog.open(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/bulk-systems.jpg", {{0, [this](Context&) { bulkCrystalDialog.open(); }}, {6, script("bulk.systems")}}, 14, popups});
        s.push_back(Shot{"tutorials/alloy-controls.jpg", {{0, [this](Context&) { substitutionalSolidSolutionDialog.open(); }},
                                                         {4, [this, host](Context&) { substitutionalSolidSolutionDialog.feedDroppedFile(host); }}}, 30, popups});
        s.push_back(Shot{"tutorials/element-picker.jpg", {{0, [this](Context&) { substitutionalSolidSolutionDialog.open(); }},
                                                         {10, script("sss.picker")}}, 20, topPopup});
        s.push_back(Shot{"tutorials/amorphous-controls.jpg", {{0, [this](Context&) { amorphousBuilderDialog.open(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/amorphous-advanced.jpg", {{0, [this](Context&) { amorphousBuilderDialog.open(); }}, {4, script("amorphous.advanced")},
                                                             {8, script("amorphous.scroll")}}, 14, popups});
        s.push_back(Shot{"tutorials/amorphous-manual.jpg", {{0, [this](Context&) { amorphousBuilderDialog.open(); }}, {4, script("amorphous.manual")}}, 14, popups});
        s.push_back(Shot{"tutorials/polycrystal-controls.jpg", {{0, [this](Context&) { polyCrystalDialog.open(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/polycrystal-loaded.jpg", {{0, [this](Context&) { polyCrystalDialog.open(); }},
                                                             {4, [this, cu4](Context&) { polyCrystalDialog.feedDroppedFile(cu4); }}}, 30, popups});
        s.push_back(Shot{"tutorials/polycrystal-modes.jpg", {{0, [this](Context&) { polyCrystalDialog.open(); }}, {6, script("poly.scroll")}, {10, script("poly.modes")}}, 14, popups});
        s.push_back(Shot{"tutorials/csl-loaded.jpg", {{0, [this](Context&) { cslDialog.open(); }}, {4, [this, cu4](Context&) { cslDialog.feedDroppedFile(cu4); }},
                                                     {10, script("csl.find")}}, 40, popups});
        s.push_back(Shot{"tutorials/csl-options.jpg", {{0, [this](Context&) { cslDialog.open(); }}, {6, script("csl.scroll")}}, 20, popups});
        s.push_back(Shot{"tutorials/custom-loaded.jpg", {{0, [this](Context&) { customStructureDialog.open(); }},
                                                        {4, [this, cu4, octahedron](Context&) { customStructureDialog.feedDroppedFile(cu4); customStructureDialog.feedDroppedFile(octahedron); }},
                                                        {8, script("custom.scale12")}, {14, script("custom.scroll")}}, 40, popups});
        s.push_back(Shot{"tutorials/custom-euler.jpg", {{0, [this](Context&) { customStructureDialog.open(); }}, {4, script("custom.type1")}, {8, script("custom.scroll")}}, 20, popups});
        s.push_back(Shot{"tutorials/custom-miller.jpg", {{0, [this](Context&) { customStructureDialog.open(); }}, {4, script("custom.type2")}, {8, script("custom.scroll")}}, 20, popups});
        s.push_back(Shot{"tutorials/custom-result.jpg", {{0, clearMain}, {1, [this](Context&) { customStructureDialog.open(); }}, {4, script("custom.plain")},
                                                        {8, script("custom.build")}, {90, [](Context&) { ImGui::ClosePopupToLevel(0, true); }}}, 60, viewportRect});
        s.push_back(Shot{"tutorials/dislocation-controls.jpg", {{0, main(host)}, {2, [this](Context&) { dislocationDialog.open(); }},
                                                               {6, [this, host](Context&) { dislocationDialog.feedDroppedFile(host); }}}, 40, popups});
        s.push_back(Shot{"tutorials/dislocation-preview.jpg", {{0, [this](Context&) { dislocationDialog.open(); }}, {6, script("dislocation.generate")}}, 240, popups});
        s.push_back(Shot{"tutorials/dislocation-summary.jpg", {{0, [this](Context&) { dislocationDialog.open(); }}, {6, script("dislocation.scroll")}}, 20, popups});
        s.push_back(Shot{"tutorials/interstitial-detection.jpg", {{0, [this](Context&) { interstitialAtomsDialog.open(); }},
                                                                 {4, [this, host](Context&) { interstitialAtomsDialog.feedDroppedFile(host); }}}, 40, popups});
        s.push_back(Shot{"tutorials/interstitial-voids.jpg", {{0, [this](Context&) { interstitialAtomsDialog.open(); }}, {6, script("interstitial.detect")}}, 180, popups});
        s.push_back(Shot{"tutorials/interstitial-placement.jpg", {{0, [this](Context&) { interstitialAtomsDialog.open(); }}, {6, script("interstitial.scroll")}}, 20, popups});
        s.push_back(Shot{"tutorials/interstitial-modes.jpg", {{0, [this](Context&) { interstitialAtomsDialog.open(); }}, {6, script("interstitial.scroll")},
                                                             {12, script("interstitial.modes")}}, 16, popups});
        s.push_back(Shot{"tutorials/interstitial-mesh.jpg", {{0, [this](Context&) { interstitialAtomsDialog.open(); }}, {6, script("interstitial.mesh")},
                                                            {10, script("interstitial.scroll")}}, 16, popups});
        s.push_back(Shot{"tutorials/sculptor-periodic.jpg", {{0, main(host)}, {2, [this](Context&) { cellSculptorDialog.open(); }}, {6, script("sculptor.useScene")},
                                                            {12, script("sculptor.addSlab")}}, 40, popups});
        s.push_back(Shot{"tutorials/sculptor-manual.jpg", {{0, [this](Context&) { cellSculptorDialog.open(); }}, {6, script("sculptor.manual")}}, 40, popups});
        s.push_back(Shot{"tutorials/nano-controls.jpg", {{0, [this](Context&) { nanoCrystalDialog.open(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/wulff-loaded.jpg", {{0, [this](Context&) { nanoCrystalDialog.open(); }}, {4, [this, cu4](Context&) { nanoCrystalDialog.feedDroppedFile(cu4); }},
                                                       {10, script("nano.wulff")}}, 60, popups});
        s.push_back(Shot{"tutorials/wulff-built.jpg", {{0, clearMain}, {1, [this](Context&) { nanoCrystalDialog.open(); }}, {6, script("nano.build")}}, 300, popups});
        // wulff-built already placed the particle in the main view.
        s.push_back(Shot{"tutorials/wulff-result.jpg", {{0, [](Context&) {}}}, 60, viewportRect});
        s.push_back(Shot{"tutorials/merge-controls.jpg", {{0, [this](Context&) { mergeStructuresDialog.open(); }}}, 20, popups});
        s.push_back(Shot{"tutorials/transform-controls.jpg", {{0, main(host)}, {2, [this](Context&) { transformDialog.open(); }}}, 14, popups});
        // Edit and display settings.
        s.push_back(Shot{"tutorials/edit-structure.jpg", {{0, main(host)}, {2, [](Context& x) { x.dialogs.openEditStructure(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/atomic-sizes.jpg", {{0, [](Context& x) { x.dialogs.openAtomicSizes(); }}}, 14, popups});
        s.push_back(Shot{"tutorials/display-style.jpg", {{0, main(host)}, {2, [](Context& x) { x.dialogs.openDisplaySettings(); }}, {6, script("display.style")}}, 14, popups});
        s.push_back(Shot{"tutorials/display-lighting.jpg", {{0, [](Context& x) { x.dialogs.openDisplaySettings(); }}, {4, script("display.lighting")}}, 14, popups});
        s.push_back(Shot{"tutorials/display-atoms.jpg", {{0, [](Context& x) { x.dialogs.openDisplaySettings(); }}, {4, script("display.atoms")}}, 14, popups});
        // The Publication look in the view, then back to the default look.
        s.push_back(Shot{"tutorials/display-publication.jpg", {{0, [](Context& x) {
            ImGui::ClosePopupToLevel(0, true);
            atomforge::display::applyPreset(x.dialogs.display, atomforge::display::presets().back());
        }}}, 20, viewportRect});
        s.push_back(Shot{"tutorials/display-colors.jpg", {{0, [](Context& x) { x.dialogs.display = atomforge::display::DisplaySettings(); x.dialogs.openDisplaySettings(); }},
                                                         {4, script("display.colors")}}, 14, popups});
        s.push_back(Shot{"tutorials/lattice-planes.jpg", {{0, [this](Context&) { showLatticePlanesDialog = true; }}, {6, script("lattice.add")}}, 14, popups});
        s.push_back(Shot{"tutorials/miller-directions.jpg", {{0, [this](Context&) { showMillerDirectionsDialog = true; }}, {6, script("miller.add")}}, 14, popups});
        s.push_back(Shot{"tutorials/polyhedral-settings.jpg", {{0, [this](Context&) { showPolyhedralSettingsDialog = true; }}}, 14, popups});
        s.push_back(Shot{"tutorials/export-image.jpg", {{0, [this, &c](Context&) {
            exportImagePopup = true;
            exportDir = c.examples.generic_u8string();
            exportDirHistory = {exportDir};
            exportHistoryIndex = 0;
            if (exportFilename[0] == '\0') std::snprintf(exportFilename, sizeof(exportFilename), "%s", "structure.png");
        }}}, 14, popups});
        s.push_back(Shot{"tutorials/file-browser.jpg", {{0, [this, &c](Context&) {
            openDir = c.examples.generic_u8string();
            std::snprintf(openFilename, sizeof(openFilename), "%s", "cu_ni_alloy.cif");
            openStructurePopup = true;
        }}}, 14, popups});
        // Information and analysis.
        s.push_back(Shot{"tutorials/structure-info.jpg", {{0, main(alloy)}, {4, [this](Context&) { requestStructureInfo = true; }}}, 30, popups});
        s.push_back(Shot{"tutorials/atom-info.jpg", {{0, main(host)}, {4, selectCentral},
                                                    {8, [this](Context&) { requestAtomInfo = true; }}}, 20, topPopup});
        s.push_back(Shot{"tutorials/atom-context.jpg", {{0, selectCentral}, {8, script("atom.context")}}, 14, topPopup});
        s.push_back(Shot{"tutorials/rdf-analysis.jpg", {{0, main(alloy)}, {2, [this](Context&) { rdfDialog.open(); }}, {6, script("rdf.run")}}, 180, popups});
        s.push_back(Shot{"tutorials/sro-warren.jpg", {{0, [this](Context&) { shortRangeOrderDialog.openRequested = true; }}, {6, script("sro.run")}}, 240, popups});
        s.push_back(Shot{"tutorials/sro-rao.jpg", {{0, [this](Context&) { shortRangeOrderDialog.openRequested = true; }}, {6, script("sro.rao")}}, 20, popups});
        // OVITO-style analyses: run in the scientific tool window, atoms coloured by the result.
        const std::string polycrystal = (c.examples / "cu_polycrystal.cif").u8string(), sphere = (c.examples / "cu_sphere.xyz").u8string();
        const auto tool = [this](const char* id) -> Action { return [this, id](Context&) { scientificToolsDialog.open(id); }; };
        s.push_back(Shot{"analysis/ptm.jpg", {{0, main(polycrystal)}, {4, tool("ptm")}, {10, script("science.run")}}, 240, {}});
        s.push_back(Shot{"analysis/grain-segmentation.jpg", {{0, tool("grain-segmentation")}, {6, script("science.run")}}, 300, {}});
        s.push_back(Shot{"analysis/voronoi.jpg", {{0, tool("voronoi")}, {6, script("science.run")}}, 300, {}});
        const std::string md = (c.examples / "cu_md.extxyz").u8string();
        s.push_back(Shot{"analysis/trajectory-lines.jpg", {{0, main(cu4)}, {2, [this, md](Context&) { trajectoryDialog.show(); trajectoryDialog.load(md); }},
                                                           {60, script("trajectory.lines")},
                                                           {80, [this](Context&) { trajectoryLinesOverlay().wholePath = true; trajectoryDialog.hide(); }}}, 30, {}});
        s.push_back(Shot{"analysis/surface-mesh.jpg", {{0, main(sphere)}, {4, tool("surface-mesh")}, {10, script("science.run")}, {200, script("science.surface")}}, 60, {}});
        // Electronic structure (needs the demonstration CHGCAR).
        if (std::filesystem::exists(c.chgcar)) {
            const std::string chgcar = c.chgcar.u8string();
            s.push_back(Shot{"electronic-volume.png", {{0, [this](Context&) { electronicDialog.open(); }}, {4, [this, chgcar](Context&) { electronicDialog.feedDroppedFile(chgcar); }}},
                             400, [] { return windowRect("###Electronic Post-processing"); }});
            s.push_back(Shot{"electronic-arbitrary.png", {{0, [this](Context&) { electronicDialog.open(); }}, {6, script("electronic.arbitrary")}},
                             120, [] { return windowRect("###Electronic Post-processing"); }});
        }
    }
    auto& c = *tutorialCapture;
    Context x{current, dialogsRef, update};
    // ATOMFORGE_TUTORIAL_FROM / _UNTIL (parts of image names) capture a part
    // of the list (in a fresh process each, for long sessions); each part
    // writes bounds-FIRST.json.
    if (c.first.empty()) {
        if (const char* from = std::getenv("ATOMFORGE_TUTORIAL_FROM"))
            while (c.index < c.shots.size() && c.shots[c.index].file.find(from) == std::string::npos) ++c.index;
        c.first = c.index < c.shots.size() ? std::filesystem::u8path(c.shots[c.index].file).stem().u8string() : "none";
    }
    const char* until = std::getenv("ATOMFORGE_TUTORIAL_UNTIL");
    if (c.index >= c.shots.size() || (c.frame == 0 && until && *until && c.shots[c.index].file.find(until) != std::string::npos)) {
        c.bounds += "\n}\n";
        std::ofstream(c.folder / ("bounds-" + c.first + ".json")) << c.bounds;
        glfwSetWindowShouldClose(graphics::mainWindow(), true);
        tutorialCapture.reset();
        return;
    }
    auto& shot = c.shots[c.index];
    int last = 0;
    if (c.frame % 60 == 0) std::fprintf(stderr, "[TutorialCapture] %s frame %d\n", shot.file.c_str(), c.frame);
    for (const auto& [at, action] : shot.actions) {
        if (at == c.frame) action(x);
        last = std::max(last, at);
    }
    if (++c.frame < last + shot.settle) return;
    // The window it shows, then the image.
    ImRect r = shot.target ? shot.target() : ImRect();
    const ImGuiIO& io = ImGui::GetIO();
    if (r.Min.x > r.Max.x) r = ImRect(0, 0, io.DisplaySize.x, io.DisplaySize.y);
    const float pad = 6.0f;
    r.Min.x = std::max(0.0f, r.Min.x - pad); r.Min.y = std::max(0.0f, r.Min.y - pad);
    r.Max.x = std::min(io.DisplaySize.x, r.Max.x + pad); r.Max.y = std::min(io.DisplaySize.y, r.Max.y + pad);
    saveJpeg(c.folder / shot.file);
    char entry[512];
    std::snprintf(entry, sizeof(entry), "%s\n  \"%s\": [%d, %d, %d, %d]", c.bounds.size() > 1 ? "," : "", shot.file.c_str(),
                  static_cast<int>(r.Min.x), static_cast<int>(r.Min.y), static_cast<int>(r.Max.x), static_cast<int>(r.Max.y));
    c.bounds += entry;
    // Reset before the next capture (dialog states are kept: later captures build on them).
    tourMenu = nullptr;
    tourSubmenu = nullptr;
    uiScript::clear();
    if (!ImGui::GetCurrentContext()->OpenPopupStack.empty()) ImGui::ClosePopupToLevel(0, true);
    operationDialog.close();
    scientificToolsDialog.close();
    electronicDialog.close();
    ++c.index;
    c.frame = 0;
}
