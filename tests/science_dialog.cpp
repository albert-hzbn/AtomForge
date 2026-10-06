// Headless rendering of every scientific tool dialog. The dialogs run the
// native engine, so no Python interpreter or environment selection may appear.
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "ui/ScientificToolsDialog.h"
#include "ui/TrajectoryDialog.h"
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
