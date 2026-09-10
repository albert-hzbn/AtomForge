#include "ui/ResponsiveLayout.h"
#include "ui/ImGuiSetup.h"
#include "imgui_internal.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}

void checkBounds()
{
    // ImGui measures a newly opened modal for one hidden frame before centering it.
    if (ImGui::GetCurrentWindow()->Hidden) return;
    const auto* viewport = ImGui::GetMainViewport();
    const auto position = ImGui::GetWindowPos();
    const auto extent = ImGui::GetWindowSize();
    if (position.x + extent.x > viewport->WorkPos.x + viewport->WorkSize.x + 1 ||
        position.y + extent.y > viewport->WorkPos.y + viewport->WorkSize.y + 1)
        std::cerr << ImGui::GetCurrentWindow()->Name << ": pos=" << position.x << ',' << position.y
                  << " size=" << extent.x << ',' << extent.y << " work=" << viewport->WorkSize.x
                  << ',' << viewport->WorkSize.y << '\n';
    require(position.x >= viewport->WorkPos.x - 1 && position.y >= viewport->WorkPos.y - 1,
            "Dialog starts outside the work area");
    require(position.x + extent.x <= viewport->WorkPos.x + viewport->WorkSize.x + 1 &&
            position.y + extent.y <= viewport->WorkPos.y + viewport->WorkSize.y + 1,
            "Dialog extends outside the work area");
}

void run(ImVec2 display, float dpi)
{
    ImGui::CreateContext();
    auto& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.DisplaySize = display;
    io.DeltaTime = 1.0f / 60;
    ImFontConfig config;
    config.SizePixels = 18;
    io.Fonts->AddFontDefault(&config);
    io.Fonts->Build();
    ImGui::GetStyle().ScaleAllSizes(dpi);
    ImGui::GetStyle().FontScaleMain = dpi;
    for (int frame = 0; frame < 4; ++frame)
    {
        if (frame == 2) io.DisplaySize = ImVec2(800, 600);
        ImGui::NewFrame();
        responsive::windowSize(ImVec2(1600, 960), ImGuiCond_FirstUseEver);
        responsive::windowConstraints(ImVec2(1100, 680), ImVec2(3200, 3200));
        responsive::begin("Large builder");
        checkBounds();
        // ImGui rounds fractional font sizes to device pixels.
        require(std::abs(responsive::dp(100) - 100 * dpi) < 3.0f, "Design size does not follow font scale");
        responsive::button("Generate and export structure", responsive::size(80, 12));
        const auto bounds = ImGui::GetItemRectSize();
        const auto text = ImGui::CalcTextSize("Generate and export structure");
        require(bounds.x >= text.x + 2 * ImGui::GetStyle().FramePadding.x - .1f &&
                bounds.y >= ImGui::GetFrameHeight() - .1f, "Button clips its label or padding");
        if (frame >= 2) require(responsive::stacked(), "Narrow layout should stack panels");
        // Simulate an old dialog position on a larger monitor.
        if (frame == 1) ImGui::SetWindowPos(ImVec2(1500, 900));
        ImGui::End();
        if (frame == 0) ImGui::OpenPopup("Large modal");
        responsive::windowSize(ImVec2(1200, 900), ImGuiCond_Appearing);
        if (responsive::beginModal("Large modal"))
        {
            checkBounds();
            ImGui::EndPopup();
        }
        ImGui::Render();
    }
    ImGui::DestroyContext();
}
}

int main()
{
    try
    {
        ImGui::CreateContext();
        auto& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(1280, 720);
        ImFontConfig baseFont;
        baseFont.SizePixels = 18;
        io.Fonts->AddFontDefault(&baseFont);
        applyLightTheme();
        for (float multiplier : {1.25f, 2.0f, .85f, 1.0f})
        {
            setInterfaceScale(multiplier);
            updateImGuiScale(nullptr);
            io.Fonts->Build();
            ImGui::NewFrame();
            if (std::abs(ImGui::GetFontSize() - 18 * multiplier) > 1.0f)
                std::cerr << "Multiplier=" << multiplier << " rendered=" << ImGui::GetFontSize()
                          << " base=" << ImGui::GetStyle().FontSizeBase << '\n';
            require(std::abs(ImGui::GetFontSize() - 18 * multiplier) <= 1.0f,
                    "Interface preference did not resize the rendered font");
            ImGui::Render();
            for (int i = 0; i < 3; ++i) { applyDarkTheme(); applyLightTheme(); }
            require(std::abs(ImGui::GetStyle().FramePadding.y - 7 * multiplier) <= 1.0f,
                    "Switching themes accumulated scaling");
        }
        ImGui::DestroyContext();
        for (const auto display : {ImVec2(1024,600), ImVec2(1280,720), ImVec2(1920,1080), ImVec2(3840,2160)})
            for (const float dpi : {1.0f, 1.25f, 1.5f, 2.0f}) run(display, dpi);
        std::cout << "Responsive layouts passed at 16 resolution/scaling combinations, including live resizing.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
