// ATOMFORGE_STEP_DIALOG_CHECK: end-to-end check of editing pipeline steps in
// their Build, Edit and Analysis dialogs, in the running application. Every
// step of ui/MenuParity is opened in its dialog (rendered in the real window
// for a few frames), its settings are read back from the dialog, and the step
// built from them must give the same structure as the original step.
#include "FileBrowser.h"
#include "ui/MenuParity.h"

#include "io/StructureLoader.h"
#include "pipeline/Pipeline.h"

#include "third_party/stb_image_write.h"

#include "imgui.h"
#include "imgui_internal.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

using atomforge::pipeline::Json;

struct FileBrowser::StepDialogCheck
{
    std::string report;               // output path
    std::filesystem::path folder;     // files the steps read
    Structure input;                  // the structure entering every step
    std::vector<atomforge::pipeline::Modifier> steps;
    std::size_t index = 0;
    int frame = 0;
    int failures = 0;
    Json results = Json::array();
};

namespace
{
// 3x3x3 conventional fcc copper cells (108 atoms).
Structure copperCells()
{
    Structure s;
    const double a = 3.615;
    const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (const auto& b : basis) {
                    AtomSite atom{};
                    atom.symbol = "Cu";
                    atom.atomicNumber = 29;
                    atom.x = (i + b[0]) * a; atom.y = (j + b[1]) * a; atom.z = (k + b[2]) * a;
                    s.atoms.push_back(atom);
                }
    s.hasUnitCell = true;
    s.cellVectors = {{{3 * a, 0, 0}, {0, 3 * a, 0}, {0, 0, 3 * a}}};
    return s;
}

std::string quoted(const std::filesystem::path& path) { return "\"" + path.generic_u8string() + "\""; }

// A representative setting of each step (others keep their defaults).
void sample(atomforge::pipeline::Modifier& m, const std::filesystem::path& folder)
{
    const auto options = [&](const std::string& text) { m.parameters["options"] = text; };
    const std::string& t = m.type;
    if (t == "build-bulk") options("--system hexagonal --spacegroup 194 --a 3.21 --c 5.21 --atom \"Mg 0.3333333333 0.6666666667 0.25\"");
    else if (t == "build-sss") options("--frac \"Cu=0.75,Ni=0.25\" --seed 3");
    else if (t == "build-gb") options("--axis \"0 0 1\" --sigma 5 --uca 2 --ucb 2 --vacuum 2");
    else if (t == "build-interface") options("--layerB " + quoted(folder / "layer-b.cif") + " --align --nmax 3 --mmax 4 --maxcells 9 --maxcellsB 16 --layersA 2 --layersB 2 --gap 2.5");
    else if (t == "build-nano") options("--shape sphere --radius 6 --vacuum 3");
    else if (t == "build-custom") options("--mesh " + quoted(folder / "cube.obj") + " --scale 1 --vacuum 3");
    else if (t == "build-poly") options("--sizex 20 --sizey 20 --sizez 20 --grains 3 --seed 5");
    else if (t == "build-amorphous") options("--element \"Si 24\" --density 2.33 --seed 4");
    else if (t == "build-dislocation") options("--character screw");
    else if (t == "build-sculpt") options("--slabs \"1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1\"" " --nx 2 --ny 2 --nz 1");
    else if (t == "build-vacancy") options("--count 3 --seed 2");
    else if (t == "add-atom") { m.parameters["element"] = std::string("Ni"); m.parameters["position"] = Json::array({1.0, 1.5, 2.0}); }
    else if (t == "set-cell") { m.parameters["a"] = Json::array({11.0, 0.0, 0.0}); m.parameters["b"] = Json::array({0.0, 11.0, 0.0}); m.parameters["c"] = Json::array({0.0, 0.0, 12.0}); }
    else if (t == "insert-interstitials") { m.parameters["kind"] = std::string("octahedral"); m.parameters["count"] = 4; m.parameters["element"] = std::string("C"); }
    else if (t == "supercell") m.parameters["matrix"] = std::string("1 1 0 -1 1 0 0 0 1");
    else if (t == "merge") { m.parameters["file"] = (folder / "layer-b.vasp").generic_u8string(); m.parameters["offset"] = Json::array({0.0, 0.0, 12.0}); }
    else if (t == "relax") { m.parameters["steps"] = 20; }
    else if (t == "nvt-dynamics" || t == "npt-dynamics") { m.parameters["steps"] = 10; m.parameters["seed"] = 7; }
    else if (t == "replicate") m.parameters["counts"] = Json::array({1.0, 1.0, 2.0});
    else if (t == "select-expression") m.parameters["expression"] = std::string("fz > 0.5");
}

// The last presented frame as a PNG (ATOMFORGE_STEP_DIALOG_SHOTS=FOLDER).
void screenshot(const std::filesystem::path& file)
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

// Whether two step results describe the same structure.
std::string difference(const atomforge::pipeline::PipelineData& a, const atomforge::pipeline::PipelineData& b)
{
    if (a.structure.atoms.size() != b.structure.atoms.size())
        return std::to_string(a.structure.atoms.size()) + " vs " + std::to_string(b.structure.atoms.size()) + " atoms";
    for (std::size_t i = 0; i < a.structure.atoms.size(); ++i) {
        const auto& p = a.structure.atoms[i];
        const auto& q = b.structure.atoms[i];
        if (p.symbol != q.symbol || std::abs(p.x - q.x) + std::abs(p.y - q.y) + std::abs(p.z - q.z) > 2e-3)
            return "atom " + std::to_string(i) + " differs";
    }
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k)
            if (std::abs(a.structure.cellVectors[r][k] - b.structure.cellVectors[r][k]) > 2e-3) return "cell differs";
    if (a.selected != b.selected) return "selection differs";
    return "";
}
}

Json FileBrowser::stepDialogSettings(const std::string& step, const Json& parameters, EditMenuDialogs& editMenuDialogs) const
{
    Json settings = parameters;
    const auto withOptions = [&](const std::string& options) { settings["options"] = options; return settings; };
    if (step == "build-bulk") return withOptions(bulkCrystalDialog.stepOptions());
    if (step == "build-gb") return withOptions(cslDialog.stepOptions());
    if (step == "build-nano") return withOptions(nanoCrystalDialog.stepOptions());
    if (step == "build-custom") return withOptions(customStructureDialog.stepOptions());
    if (step == "build-interface") return withOptions(interfaceBuilderDialog.stepOptions());
    if (step == "build-poly") return withOptions(polyCrystalDialog.stepOptions());
    if (step == "build-amorphous") return withOptions(amorphousBuilderDialog.stepOptions());
    if (step == "build-dislocation") return withOptions(dislocationDialog.stepOptions());
#if ATOMFORGE_ENABLE_SFE_BUILDER
    if (step == "build-stacking-fault") return withOptions(stackingFaultDialog.stepOptions());
#endif
#if ATOMFORGE_ENABLE_SSS_BUILDER
    if (step == "build-sss") return withOptions(substitutionalSolidSolutionDialog.stepOptions());
#endif
    if (step == "build-sculpt") return withOptions(cellSculptorDialog.stepOptions());
    if (step == "add-atom" || step == "set-cell") return editMenuDialogs.stepParameters();
    if (step == "insert-interstitials") return interstitialAtomsDialog.stepParameters();
    if (step == "supercell") return transformDialog.stepParameters();
    if (step == "merge") return mergeStructuresDialog.stepParameters();
    if (step == "relax" || step == "nvt-dynamics" || step == "npt-dynamics" || step == "standardize-cell") return scientificToolsDialog.stepParameters();
    return operationDialog.parameters();
}

void FileBrowser::runStepDialogCheck(EditMenuDialogs& editMenuDialogs)
{
    if (!stepDialogCheck) {
        static bool checked = false;
        const char* report = std::getenv("ATOMFORGE_STEP_DIALOG_CHECK");
        if (checked || !report || !*report) return;
        checked = true;
        stepDialogCheck = std::make_shared<StepDialogCheck>();
        auto& c = *stepDialogCheck;
        c.report = report;
        c.folder = std::filesystem::temp_directory_path() / "AtomForge-step-dialog-check";
        std::filesystem::create_directories(c.folder);
        c.input = copperCells();
        // Files read by the interface, custom-structure and merge steps.
        Structure layer;
        layer.hasUnitCell = true;
        layer.cellVectors = {{{3.52, 0, 0}, {0, 3.52, 0}, {0, 0, 3.52}}};
        for (const auto& b : {std::array<double, 3>{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}) {
            AtomSite atom{};
            atom.symbol = "Ni"; atom.atomicNumber = 28;
            atom.x = b[0] * 3.52; atom.y = b[1] * 3.52; atom.z = b[2] * 3.52;
            layer.atoms.push_back(atom);
        }
        saveStructure(layer, (c.folder / "layer-b.cif").u8string(), "cif");
        saveStructure(layer, (c.folder / "layer-b.vasp").u8string(), "vasp");
        std::ofstream cube(c.folder / "cube.obj");
        cube << "v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\nv 0 0 10\nv 10 0 10\nv 10 10 10\nv 0 10 10\n"
                "f 1 3 2\nf 1 4 3\nf 5 6 7\nf 5 7 8\nf 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\n";
        for (const auto& location : menuLocations()) {
            if (!atomforge::pipeline::findModifierType(location.step)) continue;
            auto m = atomforge::pipeline::makeModifier(location.step);
            sample(m, c.folder);
            c.steps.push_back(m);
        }
    }
    auto& c = *stepDialogCheck;
    if (c.index >= c.steps.size()) {
        Json document = Json::object();
        document["steps"] = c.results;
        document["failures"] = c.failures;
        std::ofstream(std::filesystem::u8path(c.report)) << document.dump(2) << "\n";
        glfwSetWindowShouldClose(glfwGetCurrentContext(), true);
        stepDialogCheck.reset();
        return;
    }
    const auto& step = c.steps[c.index];
    if (c.frame == 0) {
        // Open the step's dialog as "Edit in dialog" does.
        StepEdit request;
        request.step = step.type;
        request.parameters = step.parameters;
        request.input = c.input;
        request.commit = [](const Json&) {};
        operationDialog.close();
        scientificToolsDialog.close();
        openStepDialog(std::move(request), editMenuDialogs);
    } else if (c.frame == 6) {
        // The dialog has been drawn: what is on screen, and what it would write back.
        ImGuiContext& g = *ImGui::GetCurrentContext();
        std::string shown;
        if (!g.OpenPopupStack.empty() && g.OpenPopupStack.back().Window) shown = g.OpenPopupStack.back().Window->Name;
        else if (operationDialog.editingStep() && operationDialog.step() == step.type)
            shown = std::string("operation window: ") + atomforge::pipeline::findModifierType(step.type)->title;
        else if (scientificToolsDialog.isOpen()) shown = "scientific tool window";
        if (const char* shots = std::getenv("ATOMFORGE_STEP_DIALOG_SHOTS"); shots && *shots)
            screenshot(std::filesystem::u8path(shots) / (std::to_string(c.index + 1) + "-" + step.type + ".png"));
        Json entry = Json::object();
        entry["step"] = step.type;
        entry["menu"] = menuPathOf(step.type);
        entry["dialog"] = shown;
        std::string problem = shown.empty() ? "no dialog shown" : "";
        try {
            auto rebuilt = step;
            rebuilt.parameters = stepDialogSettings(step.type, step.parameters, editMenuDialogs);
            entry["read_back"] = rebuilt.parameters;
            atomforge::pipeline::Pipeline original, again;
            original.modifiers = {step};
            again.modifiers = {rebuilt};
            std::vector<atomforge::pipeline::StageResult> first, second;
            const auto a = original.evaluate(c.input, &first);
            const auto b = again.evaluate(c.input, &second);
            if (!first.front().error.empty()) entry["original_error"] = first.front().error;
            // A step that needs earlier steps (a selection, a property) fails alone in both forms.
            if (!second.front().error.empty() && second.front().error == first.front().error) entry["note"] = "fails without earlier steps, as the original";
            else if (!second.front().error.empty()) problem += (problem.empty() ? "" : "; ") + ("read-back step fails: " + second.front().error);
            else if (first.front().error.empty()) {
                const std::string diff = difference(a, b);
                if (!diff.empty()) problem += (problem.empty() ? "" : "; ") + ("different result: " + diff);
                entry["atoms"] = static_cast<int>(b.structure.atoms.size());
            }
        } catch (const std::exception& error) { problem += (problem.empty() ? "" : "; ") + std::string(error.what()); }
        entry["ok"] = problem.empty();
        if (!problem.empty()) { entry["problem"] = problem; ++c.failures; }
        c.results.push(entry);
        // Close the dialog (its step editing ends) before the next step.
        ImGui::ClosePopupToLevel(0, true);
    } else if (c.frame >= 8) {
        ++c.index;
        c.frame = 0;
        return;
    }
    ++c.frame;
}
