// Headless rendering of every scientific tool dialog. The dialogs run the
// native engine, so no Python interpreter or environment selection may appear.
#include "science/ScienceCatalog.h"
#include "ui/ScientificToolsDialog.h"
#include "imgui.h"
#include "imgui_internal.h"

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
        int rendered = 0;
        for (const auto& tool : scienceToolCatalog()) {
            if (!dialog.open(tool.id)) throw std::runtime_error(std::string("Cannot open ") + tool.id);
            for (int frame = 0; frame < 3; ++frame) {
                ImGui::NewFrame();
                dialog.draw(copper, [](Structure&) {});
                ImGui::Render();
            }
            const std::string title = std::string(tool.title) + "###Scientific analysis";
            ImGuiWindow* window = ImGui::FindWindowByName(title.c_str());
            if (!window || !window->WasActive) throw std::runtime_error(std::string("Dialog not shown for ") + tool.id);
            ++rendered;
        }
        ImGui::DestroyContext();
        std::cout << "Rendered " << rendered << " native scientific tool dialogs\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Scientific dialog test failed: " << error.what() << '\n';
        return 1;
    }
}
