// Headless rendering of every scientific tool dialog. The dialogs run the
// native engine, so no Python interpreter or environment selection may appear.
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "ui/ScientificToolsDialog.h"
#include "ui/TrajectoryDialog.h"
#include "ui/PipelineDialog.h"
#include "ui/OperationDialog.h"
#include "ui/MenuParity.h"
#include "ui/SciencePlot.h"
#include <cmath>
#include "imgui.h"
#include "imgui_internal.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

int main()
{
    try {
        ImGui::CreateContext();
        auto& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(1280, 900);
        io.DeltaTime = 1.0f / 60;
        io.Fonts->AddFontDefault();
        io.Fonts->Build();
        Structure copper;
        copper.hasUnitCell = true;
        copper.cellVectors = {{{3.6, 0, 0}, {0, 3.6, 0}, {0, 0, 3.6}}};
        for (const auto& p : {std::array<double, 3>{0, 0, 0}, {0, 1.8, 1.8}, {1.8, 0, 1.8}, {1.8, 1.8, 0}})
            copper.atoms.push_back({"Cu", 29, p[0], p[1], p[2]});
        ScientificToolsDialog dialog;
        if (dialog.open("no-such-tool")) throw std::runtime_error("Unknown tool ids must be rejected");
        // Input-source rules: position inputs fill a cell that exists on the
        // same tool, and only trajectory tools accept the loaded trajectory.
        if (std::string(ScientificToolsDialog::companionCell("reference")) != "reference_cell" ||
            std::string(ScientificToolsDialog::companionCell("current")) != "current_cell" ||
            ScientificToolsDialog::companionCell("timestep_fs") != nullptr)
            throw std::runtime_error("Companion cell mapping");
        for (const auto& tool : scienceToolCatalog())
            for (const auto& parameter : tool.parameters) {
                const char* companion = ScientificToolsDialog::companionCell(parameter.name);
                // Only data inputs fill a companion cell; structure inputs carry their own.
                if (!companion || std::string(parameter.kind) != "data" || std::string(tool.id) == "structure-factor") continue;
                bool found = false;
                for (const auto& other : tool.parameters) found = found || std::string(other.name) == companion;
                if (!found) throw std::runtime_error(std::string(tool.id) + " lacks the companion " + companion);
            }
        if (!ScientificToolsDialog::acceptsTrajectory("msd", "positions") || !ScientificToolsDialog::acceptsTrajectory("vacf", "velocities") ||
            ScientificToolsDialog::acceptsTrajectory("centrosymmetry", "positions"))
            throw std::runtime_error("Trajectory input rules");
        dialog.setTrajectorySource("trajectory.xyz");
        int rendered = 0;
        int layoutCounts[6] = {0, 0, 0, 0, 0, 0};
        for (const auto& tool : scienceToolCatalog()) {
            if (!dialog.open(tool.id)) throw std::runtime_error(std::string("Cannot open ") + tool.id);
            for (int frame = 0; frame < 3; ++frame) {
                ImGui::NewFrame();
                dialog.draw(copper, [](Structure&) {}, [](const std::string&, const std::vector<double>&) {});
                ImGui::Render();
            }
            const auto layout = ScientificToolsDialog::layoutFor(tool.id);
            const std::string title = std::string(tool.title) + ScientificToolsDialog::windowId(layout);
            ImGuiWindow* window = ImGui::FindWindowByName(title.c_str());
            if (!window || !window->WasActive) throw std::runtime_error(std::string("Dialog not shown for ") + tool.id);
            ++layoutCounts[static_cast<int>(layout)];
            ++rendered;
        }
        // Tools are spread over distinct layouts rather than one shared form.
        using Layout = ScientificToolsDialog::Layout;
        for (int count : layoutCounts)
            if (count == 0) throw std::runtime_error("Every window layout must be used by at least one tool");
        if (*std::max_element(layoutCounts, layoutCounts + 6) * 3 > rendered) throw std::runtime_error("One layout serves too many tools");
        if (ScientificToolsDialog::layoutFor("centrosymmetry") != Layout::Atoms || ScientificToolsDialog::layoutFor("nvt") != Layout::Simulation ||
            ScientificToolsDialog::layoutFor("dft-inputs") != Layout::Generator || ScientificToolsDialog::layoutFor("msd") != Layout::Trajectory ||
            ScientificToolsDialog::layoutFor("elastic-tensor") != Layout::Properties || ScientificToolsDialog::layoutFor("powder-xrd") != Layout::Plot)
            throw std::runtime_error("Tool layout assignment");
        // The former shared Inputs/Method/Results tab bar is gone.
        for (auto& tabBar : ImGui::GetCurrentContext()->TabBars.Buf)
            for (auto& tab : tabBar.Tabs)
                if (ImGui::TabBarGetTabName(&tabBar, &tab) == std::string("Inputs") || ImGui::TabBarGetTabName(&tabBar, &tab) == std::string("Method"))
                    throw std::runtime_error("Generic Inputs/Method tabs are still drawn");
        // Plot widget edge cases: log axis with non-positive values, NaN gaps,
        // markers and reference lines, a single point, and no data.
        using atomforge::science::PlotSpec;
        PlotSpec logPlot;
        logPlot.title = "log"; logPlot.xLabel = "x"; logPlot.yLabel = "y"; logPlot.logY = true;
        logPlot.series.push_back({"a", {0, 1, 2, 3, 4}, {1, 1e-3, 0, -1, 1e-6}, false});
        logPlot.series.push_back({"b", {0, 1, 2, 3}, {2, std::nan(""), 3, 4}, true});
        logPlot.markers = {{1.5, "X"}};
        logPlot.horizontal = {1e-2};
        PlotSpec single = logPlot;
        single.logY = false;
        single.series = {{"one", {1}, {1}, true}};
        PlotSpec empty;
        for (int frame = 0; frame < 3; ++frame) {
            io.MousePos = ImVec2(400, 200);
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(ImVec2(800, 900));
            ImGui::Begin("Plots");
            uiPlot::drawSciencePlot("log", logPlot);
            uiPlot::drawSciencePlot("single", single);
            uiPlot::drawSciencePlot("empty", empty);
            ImGui::End();
            ImGui::Render();
        }
        // Project persistence: tool, inputs and the last result survive snapshot/restore.
        {
            using atomforge::science::Json;
            Json request = Json::parse(R"({"positions": [[0,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]],
                "cell": [[3.6,0,0],[0,3.6,0],[0,0,3.6]], "pbc": [true,true,true], "cutoff_A": 3.0})");
            const auto output = atomforge::science::runTool("centrosymmetry", request, ".");
            const Json document = atomforge::science::resultDocument("centrosymmetry", request, output);
            ScientificToolsDialog saved;
            if (!saved.open("centrosymmetry")) throw std::runtime_error("Cannot open centrosymmetry");
            Json state = Json::parse(saved.snapshot());
            if (state.at("tool").string() != "centrosymmetry" || !state.at("open").boolean() || state.contains("result"))
                throw std::runtime_error("Snapshot of a fresh dialog");
            state["fields"].items()[1]["value"] = "2.75";
            state["result"] = document;
            state["frames"] = "1\nframe\nCu 0 0 0\n";
            ScientificToolsDialog restored;
            restored.restore(state.dump());
            const Json again = Json::parse(restored.snapshot());
            if (again.at("tool").string() != "centrosymmetry" || !again.at("open").boolean())
                throw std::runtime_error("Restored tool or visibility");
            if (again.at("fields").items()[1].at("value").string() != "2.75") throw std::runtime_error("Restored input value");
            if (again.at("result").dump() != document.dump()) throw std::runtime_error("Restored result document");
            if (again.at("frames").string() != "1\nframe\nCu 0 0 0\n") throw std::runtime_error("Restored result structures");
            // The restored result is shown, including the per-atom colouring control.
            std::string coloured;
            for (int frame = 0; frame < 3; ++frame) {
                ImGui::NewFrame();
                restored.draw(copper, [](Structure&) {}, [&](const std::string& name, const std::vector<double>&) { coloured = name; });
                ImGui::Render();
            }
            // Kept runs overlay their plots on the current result's plots.
            const std::size_t baseSeries = restored.displayedPlot(0).series.size();
            restored.keepForComparison();
            restored.keepForComparison();
            if (restored.keptRuns() != 2) throw std::runtime_error("Kept runs");
            const auto overlaid = restored.displayedPlot(0);
            if (overlaid.series.size() != 3 * baseSeries || overlaid.series.front().name.rfind("Run 1: ", 0) != 0 ||
                overlaid.series.back().name.rfind("Current: ", 0) != 0)
                throw std::runtime_error("Overlay of kept runs");
            for (int frame = 0; frame < 2; ++frame) {
                ImGui::NewFrame();
                restored.draw(copper, [](Structure&) {}, [&](const std::string& name, const std::vector<double>&) { coloured = name; });
                ImGui::Render();
            }
            ImGuiWindow* window = ImGui::FindWindowByName("Centrosymmetry parameter###Atom analysis");
            if (!window || !window->WasActive) throw std::runtime_error("Restored dialog not shown");
            // Unknown tools (from a newer release) and malformed state leave the dialog closed.
            ScientificToolsDialog unknown;
            unknown.restore(R"({"open": true, "tool": "from-the-future"})");
            unknown.restore("not json");
            if (Json::parse(unknown.snapshot()).at("open").boolean()) throw std::runtime_error("Unknown tool restored as open");
        }
        // Every layout renders a real result: plots, tensors, generated files and frames.
        {
            using atomforge::science::Json;
            const std::string copperJson = R"({"symbols":["Cu","Cu","Cu","Cu"],"positions":[[0,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]],"cell":[[3.6,0,0],[0,3.6,0],[0,0,3.6]]})";
            const std::vector<std::pair<std::string, std::string>> cases = {
                {"msd", R"({"positions":[[[0,0,0],[1,1,1]],[[0.1,0,0],[1,1.1,1]],[[0.2,0,0],[1,1.2,1]],[[0.3,0,0],[1,1.3,1]]],"timestep_fs":1.0})"},
                {"elastic-tensor", R"({"strains":[[0.01,0,0,0,0,0],[0,0.01,0,0,0,0],[0,0,0.01,0,0,0],[0,0,0,0.01,0,0],[0,0,0,0,0.01,0],[0,0,0,0,0,0.01],[0,0,0,0,0,0]],
                    "stresses_GPa":[[1.7,1.2,1.2,0,0,0],[1.2,1.7,1.2,0,0,0],[1.2,1.2,1.7,0,0,0],[0,0,0,0.75,0,0],[0,0,0,0,0.75,0], [0,0,0,0,0,0.75],[0,0,0,0,0,0]]})"},
                {"powder-xrd", std::string(R"({"structure":)") + copperJson + "}"},
                {"relax", std::string(R"({"structure":)") + copperJson + R"(,"calculator":{"potential":"EMT"},"steps":5})"},
                {"lammps-export", std::string(R"({"structure":)") + copperJson + "}"},
            };
            for (const auto& [id, text] : cases) {
                const Json request = Json::parse(text);
                const auto output = atomforge::science::runTool(id, request, ".");
                ScientificToolsDialog view;
                if (!view.open(id)) throw std::runtime_error("Cannot open " + id);
                Json state = Json::parse(view.snapshot());
                state["result"] = atomforge::science::resultDocument(id, request, output);
                if (!output.frames.empty()) state["frames"] = "1\nframe\nCu 0 0 0\n";
                view.restore(state.dump());
                for (int frame = 0; frame < 3; ++frame) {
                    ImGui::NewFrame();
                    view.draw(copper, [](Structure&) {}, [](const std::string&, const std::vector<double>&) {});
                    ImGui::Render();
                }
                const auto& tool = *std::find_if(scienceToolCatalog().begin(), scienceToolCatalog().end(), [&](const auto& t) { return id == t.id; });
                const std::string title = std::string(tool.title) + ScientificToolsDialog::windowId(ScientificToolsDialog::layoutFor(id));
                ImGuiWindow* window = ImGui::FindWindowByName(title.c_str());
                if (!window || !window->WasActive) throw std::runtime_error("Result layout not shown for " + id);
            }
        }
        // Per-frame analyses of Trajectory playback colour each frame's atoms.
        {
            Structure supercell = copper;
            supercell.atoms.clear();
            supercell.cellVectors = {{{7.2, 0, 0}, {0, 7.2, 0}, {0, 0, 7.2}}};
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k)
                for (const auto& atom : copper.atoms) {
                    AtomSite site = atom;
                    site.x += 3.6 * i; site.y += 3.6 * j; site.z += 3.6 * k;
                    supercell.atoms.push_back(site);
                }
            const auto types = TrajectoryDialog::analyseFrame(supercell, "structure-type", 3.0);
            if (types.values.size() != 32) throw std::runtime_error("Per-frame structure types per atom");
            for (double t : types.values) if (t != 1) throw std::runtime_error("Perfect fcc frame is all fcc");
            const auto csp = TrajectoryDialog::analyseFrame(supercell, "centrosymmetry", 3.0);
            for (double v : csp.values) if (!(std::abs(v) < 1e-9)) throw std::runtime_error("Perfect fcc frame has zero centrosymmetry");
            const auto q6 = TrajectoryDialog::analyseFrame(supercell, "bond-order", 3.0);
            if (q6.values.size() != 32 || std::abs(q6.values[0] - 0.574524) > 1e-5) throw std::runtime_error("Per-frame q6 of fcc");
            if (TrajectoryDialog::frameAnalyses().size() != 3) throw std::runtime_error("Per-frame analysis list");
            for (const auto& [id, title] : TrajectoryDialog::frameAnalyses())
                if (!findScienceTool(id)) throw std::runtime_error("Per-frame analysis names an unknown tool " + id);
        }
        // Structure pipeline panel: an active pipeline shows its output and keeps its input.
        {
            Structure shown = copper;
            PipelineDialog panel;
            panel.restore(R"({"format":"atomforge-pipeline","version":1,"active":true,"open":true,"modifiers":[
                {"type":"replicate","parameters":{"counts":[2,2,1]}},
                {"type":"select-expression","parameters":{"expression":"fz > 0.4"}},
                {"type":"delete-selected"}]})", &copper);
            int updates = 0;
            for (int frame = 0; frame < 3; ++frame) {
                ImGui::NewFrame();
                panel.draw(shown, [&](Structure&) { ++updates; if (!panel.isUpdating()) throw std::runtime_error("Pipeline updates are flagged"); });
                ImGui::Render();
            }
            if (shown.atoms.size() != 8 || updates != 1) throw std::runtime_error("Pipeline output in the view (" + std::to_string(shown.atoms.size()) + " atoms)");
            if (panel.editor().input().atoms.size() != 4) throw std::runtime_error("Pipeline input kept");
            if (!ImGui::FindWindowByName("Structure pipeline") || !ImGui::FindWindowByName("Structure pipeline")->WasActive)
                throw std::runtime_error("Pipeline panel not shown");
            // Editing re-evaluates on the next frame; the project snapshot round-trips.
            panel.editor().setEnabled(2, false);
            ImGui::NewFrame();
            panel.draw(shown, [&](Structure&) { ++updates; });
            ImGui::Render();
            if (shown.atoms.size() != 16) throw std::runtime_error("Disabling a step restores its atoms");
            PipelineDialog reopened;
            reopened.restore(panel.snapshot(), &copper);
            if (reopened.editor().size() != 3 || reopened.editor().modifier(2).enabled) throw std::runtime_error("Pipeline snapshot round trip");
        }
        // Moving and removing steps with real mouse input on the drawn controls.
        {
            Structure shown = copper;
            PipelineDialog panel;
            panel.restore(R"({"format":"atomforge-pipeline","version":1,"active":true,"open":true,"modifiers":[
                {"type":"replicate","parameters":{"counts":[2,1,1]}},
                {"type":"wrap"},
                {"type":"select-element","parameters":{"elements":"Cu"}},
                {"type":"invert-selection"}]})", &copper);
            auto& io = ImGui::GetIO();
            const auto frame = [&] {
                ImGui::NewFrame();
                ImGui::SetNextWindowPos(ImVec2(0, 0));
                ImGui::SetNextWindowSize(ImVec2(700, 880));
                panel.draw(shown, [](Structure&) {});
                ImGui::Render();
            };
            const auto centre = [&](const std::string& name) {
                for (const auto& [control, rect] : panel.controls())
                    if (control == name) return ImVec2(0.5f * (rect[0] + rect[2]), 0.5f * (rect[1] + rect[3]));
                throw std::runtime_error("Control not drawn: " + name);
            };
            const auto order = [&] {
                std::string text;
                for (std::size_t i = 0; i < panel.editor().size(); ++i) text += (i ? "," : "") + panel.editor().modifier(i).type;
                return text;
            };
            const auto click = [&](const std::string& name) {
                frame();
                const ImVec2 at = centre(name);
                io.AddMousePosEvent(at.x, at.y); frame();
                io.AddMouseButtonEvent(0, true); frame();
                io.AddMouseButtonEvent(0, false); frame();
                frame();
            };
            for (int k = 0; k < 3; ++k) frame();
            click("down 0");
            if (order() != "wrap,replicate,select-element,invert-selection") throw std::runtime_error("Down arrow did not move the step: " + order());
            click("up 2");
            if (order() != "wrap,select-element,replicate,invert-selection") throw std::runtime_error("Up arrow did not move the step: " + order());
            click("remove 0");
            if (order() != "select-element,replicate,invert-selection") throw std::runtime_error("X did not remove the step: " + order());
            // The selected step's own buttons.
            click("step 2");
            click("move up");
            if (order() != "select-element,invert-selection,replicate") throw std::runtime_error("Move up did not move the selected step: " + order());
            click("remove");
            if (order() != "select-element,replicate") throw std::runtime_error("Remove did not remove the selected step: " + order());
            click("duplicate");
            if (order() != "select-element,replicate,replicate") throw std::runtime_error("Duplicate: " + order());
            // Drag the first step onto the last.
            frame();
            const ImVec2 from = centre("step 0"), to = centre("step 2");
            io.AddMousePosEvent(from.x, from.y); frame();
            io.AddMouseButtonEvent(0, true); frame();
            for (int k = 1; k <= 8; ++k) { io.AddMousePosEvent(from.x + (to.x - from.x) * k / 8, from.y + (to.y - from.y) * k / 8); frame(); }
            io.AddMouseButtonEvent(0, false); frame();
            frame();
            if (order() != "replicate,replicate,select-element") throw std::runtime_error("Dragging did not move the step: " + order());
            // The output follows the new order: two doublings of the 4-atom cell.
            if (shown.atoms.size() != 16) throw std::runtime_error("Output after reordering: " + std::to_string(shown.atoms.size()));
            io.AddMousePosEvent(-1, -1);
            frame();
        }
        // Steps added before a pipeline exists start it on the active structure and can be removed.
        {
            Structure shown = copper;
            PipelineDialog panel;
            panel.restore(R"({"format":"atomforge-pipeline","version":1,"active":false,"open":true,"modifiers":[{"type":"wrap"},{"type":"center"}]})", nullptr);
            ImGui::NewFrame();
            panel.draw(shown, [](Structure&) {});
            ImGui::Render();
            bool listed = false;
            for (const auto& control : panel.controls()) listed = listed || control.first == "remove 1";
            if (!listed) throw std::runtime_error("Steps are listed (and removable) before the pipeline starts");
            if (panel.isActive()) throw std::runtime_error("Not started yet");
        }
        // Menu parity: every pipeline step has a menu location, and menu operations work.
        {
            for (const auto& type : atomforge::pipeline::modifierTypes())
                if (menuPathOf(type.id).empty()) throw std::runtime_error(std::string("Pipeline step without a menu entry: ") + type.id);
            for (const auto& location : menuLocations())
                if (!location.dialog && std::string(location.step).rfind("build-", 0) != 0 && !atomforge::pipeline::findModifierType(location.step))
                    throw std::runtime_error(std::string("Menu entry for an unknown step: ") + location.step);
            Structure cell = copper;
            OperationDialog dialog;
            if (!dialog.open("delete-selected") || !OperationDialog::actsOnAtoms("delete-selected")) throw std::runtime_error("Delete operation");
            dialog.setTarget(OperationDialog::Condition, "fx > 0.25");
            dialog.apply(cell, {});
            if (cell.atoms.size() != 2) throw std::runtime_error("Delete by condition: " + std::to_string(cell.atoms.size()));
            cell = copper;
            dialog.open("assign-element");
            dialog.setParameter("element", atomforge::pipeline::Json(std::string("Ni")));
            dialog.setTarget(OperationDialog::ViewSelection);
            dialog.apply(cell, {1, 3});
            if (cell.atoms[1].symbol != "Ni" || cell.atoms[3].symbol != "Ni" || cell.atoms[0].symbol != "Cu") throw std::runtime_error("Assign element to the view selection");
            bool noMatch = false;
            dialog.setTarget(OperationDialog::Condition, "x > 100");
            try { dialog.apply(cell, {}); } catch (const std::exception&) { noMatch = true; }
            if (!noMatch) throw std::runtime_error("An empty condition must not change the structure");
            // Selection steps return the new view selection, combined with the current one.
            dialog.open("select-expression");
            dialog.setParameter("expression", atomforge::pipeline::Json(std::string("fz > 0.25")));
            dialog.setParameter("mode", atomforge::pipeline::Json(std::string("add")));
            const auto selected = dialog.apply(cell, {0});
            if (selected != std::vector<int>({0, 1, 2})) throw std::runtime_error("Selection step adds to the view selection");
            if (cell.atoms.size() != 4) throw std::runtime_error("Selection steps leave the structure unchanged");
            // Drawn with "Add to pipeline" wiring.
            std::vector<atomforge::pipeline::Modifier> added;
            OperationDialog::Callbacks callbacks;
            callbacks.update = [](Structure&) {};
            callbacks.select = [](const std::vector<int>&) {};
            callbacks.addToPipeline = [&](const std::vector<atomforge::pipeline::Modifier>& steps) { added = steps; };
            dialog.open("replicate");
            for (int frame = 0; frame < 2; ++frame) {
                ImGui::NewFrame();
                dialog.draw(cell, {}, callbacks);
                ImGui::Render();
            }
            if (!ImGui::FindWindowByName("Replicate###Structure operation") || !ImGui::FindWindowByName("Replicate###Structure operation")->WasActive)
                throw std::runtime_error("Operation dialog not shown");
            PipelineDialog panel;
            panel.addSteps({atomforge::pipeline::makeModifier("replicate")}, copper);
            if (!panel.isActive() || panel.editor().size() != 1) throw std::runtime_error("Add to pipeline starts the pipeline");
        }
        // Simulation steps open in their tool window, filled from the step.
        {
            ScientificToolsDialog tools;
            StepEdit edit;
            edit.step = "relax";
            edit.parameters = atomforge::pipeline::makeModifier("relax").parameters;
            edit.parameters["fmax"] = 0.02;
            edit.parameters["relax_cell"] = true;
            edit.parameters["calculator"] = std::string("{\"potential\": \"LennardJones\", \"epsilon\": 0.4, \"sigma\": 2.3, \"cutoff\": 6}");
            edit.input = copper;
            atomforge::pipeline::Json committed;
            edit.commit = [&](const atomforge::pipeline::Json& p) { committed = p; };
            if (!tools.editStep(std::move(edit))) throw std::runtime_error("Relax step opens in the tool window");
            for (int frame = 0; frame < 2; ++frame) {
                ImGui::NewFrame();
                tools.draw(copper, [](Structure&) {}, {});
                ImGui::Render();
            }
            const auto back = tools.stepParameters();
            if (std::abs(back.at("fmax").number() - 0.02) > 1e-12 || !back.at("relax_cell").boolean()) throw std::runtime_error("Relax settings round trip");
            const auto calculator = atomforge::pipeline::Json::parse(back.at("calculator").string());
            if (calculator.at("potential").string() != "LennardJones" || std::abs(calculator.at("sigma").number() - 2.3) > 1e-9)
                throw std::runtime_error("Potential round trip: " + back.at("calculator").string());
        }
        // A step opens in its menu dialog ("Edit in dialog"), and "Update step" writes it back.
        {
            Structure shown = copper;
            PipelineDialog panel;
            panel.restore(R"({"format":"atomforge-pipeline","version":1,"active":true,"open":true,"modifiers":[
                {"type":"replicate","parameters":{"counts":[2,1,1]}},
                {"type":"assign-element","parameters":{"element":"Ni","target":"all"}}]})", &copper);
            OperationDialog dialog;
            OperationDialog::Callbacks none;
            auto& io = ImGui::GetIO();
            const auto frame = [&] {
                ImGui::NewFrame();
                ImGui::SetNextWindowPos(ImVec2(0, 0));
                ImGui::SetNextWindowSize(ImVec2(700, 880));
                panel.draw(shown, [](Structure&) {});
                if (StepEdit request; panel.consumeDialogRequest(request)) dialog.editStep(std::move(request));
                ImGui::SetNextWindowPos(ImVec2(720, 0));
                dialog.draw(shown, {}, none);
                ImGui::Render();
            };
            const auto clickAt = [&](ImVec2 at) {
                io.AddMousePosEvent(at.x, at.y); frame();
                io.AddMouseButtonEvent(0, true); frame();
                io.AddMouseButtonEvent(0, false); frame();
                frame();
            };
            const auto centre = [&](const std::string& name) {
                for (const auto& [control, rect] : panel.controls())
                    if (control == name) return ImVec2(0.5f * (rect[0] + rect[2]), 0.5f * (rect[1] + rect[3]));
                throw std::runtime_error("Control not drawn: " + name);
            };
            for (int k = 0; k < 3; ++k) frame();
            // Double-clicking step 2 opens the operation window of Edit > Structure Operations > Assign Element.
            const ImVec2 row = centre("step 1");
            io.AddMousePosEvent(row.x, row.y); frame();
            for (int k = 0; k < 2; ++k) { io.AddMouseButtonEvent(0, true); frame(); io.AddMouseButtonEvent(0, false); frame(); }
            frame();
            if (!dialog.editingStep() || dialog.step() != "assign-element") throw std::runtime_error("Double-click opens the step's dialog");
            // Its settings come from the step; change one and press Update step.
            dialog.setParameter("element", atomforge::pipeline::Json(std::string("Zn")));
            frame();
            const auto& rect = dialog.updateButton();
            clickAt(ImVec2(0.5f * (rect[0] + rect[2]), 0.5f * (rect[1] + rect[3])));
            if (panel.editor().modifier(1).parameters.at("element").string() != "Zn") throw std::runtime_error("Update step writes the dialog's settings into the step");
            if (shown.atoms.size() != 8 || shown.atoms[0].symbol != "Zn") throw std::runtime_error("The pipeline output follows the updated step");
            io.AddMousePosEvent(-1, -1);
            frame();
        }
        const auto ticks = uiPlot::niceTicks(0.013, 0.98, 5);
        if (ticks.size() < 4 || ticks.front() < 0.013 || ticks.back() > 0.98) throw std::runtime_error("Tick generation");
        ImGui::DestroyContext();
        std::cout << "Rendered " << rendered << " native scientific tool dialogs and plot edge cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Scientific dialog test failed: " << error.what() << '\n';
        return 1;
    }
}
