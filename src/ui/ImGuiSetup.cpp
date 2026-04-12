#include "ImGuiSetup.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>

#include <algorithm>

namespace
{
constexpr float kBaseFontSizePixels = 16.0f;

ImGuiStyle gBaseStyle;
float gUiScale = 1.0f;
bool gBaseStyleCaptured = false;
bool gImGuiBackendsReady = false;

float clampUiScale(float scale)
{
    return std::max(1.0f, std::min(scale, 4.0f));
}

float computeUiScale(GLFWwindow* window)
{
    if (!window)
        return 1.0f;

    float contentScaleX = 1.0f;
    float contentScaleY = 1.0f;
    glfwGetWindowContentScale(window, &contentScaleX, &contentScaleY);

    int windowWidth = 0;
    int windowHeight = 0;
    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);

    float framebufferScaleX = 1.0f;
    float framebufferScaleY = 1.0f;
    if (windowWidth > 0)
        framebufferScaleX = (float)framebufferWidth / (float)windowWidth;
    if (windowHeight > 0)
        framebufferScaleY = (float)framebufferHeight / (float)windowHeight;

    const float scale = std::max(
        std::max(contentScaleX, contentScaleY),
        std::max(framebufferScaleX, framebufferScaleY));
    return clampUiScale(scale);
}

void captureBaseStyle()
{
    gBaseStyle = ImGui::GetStyle();
    gBaseStyleCaptured = true;
}

void applyScaledStyle()
{
    ImGuiStyle scaled = gBaseStyle;
    scaled.ScaleAllSizes(gUiScale);
    scaled.DisplayWindowPadding = ImVec2(0.0f, 0.0f);
    scaled.DisplaySafeAreaPadding = ImVec2(0.0f, 0.0f);
    ImGui::GetStyle() = scaled;
}

void rebuildFonts()
{
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();

    ImFontConfig fontConfig;
    fontConfig.SizePixels = kBaseFontSizePixels * gUiScale;
    io.Fonts->AddFontDefault(&fontConfig);
    io.FontGlobalScale = 1.0f;

    if (gImGuiBackendsReady)
    {
        ImGui_ImplOpenGL3_DestroyDeviceObjects();
        io.Fonts->Build();
        ImGui_ImplOpenGL3_CreateDeviceObjects();
    }
}

void applyCommonStyle()
{
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowPadding = ImVec2(14.0f, 12.0f);
    style.FramePadding = ImVec2(10.0f, 7.0f);
    style.CellPadding = ImVec2(8.0f, 6.0f);
    style.ItemSpacing = ImVec2(10.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(8.0f, 6.0f);
    style.TouchExtraPadding = ImVec2(1.0f, 1.0f);
    style.IndentSpacing = 20.0f;
    style.ScrollbarSize = 14.0f;
    style.GrabMinSize = 10.0f;

    style.WindowRounding = 12.0f;
    style.ChildRounding = 10.0f;
    style.FrameRounding = 9.0f;
    style.PopupRounding = 12.0f;
    style.ScrollbarRounding = 12.0f;
    style.GrabRounding = 12.0f;
    style.TabRounding = 10.0f;

    style.WindowBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.PopupBorderSize = 1.0f;
    style.FrameBorderSize = 0.0f;
    style.TabBorderSize = 0.0f;
}
} // namespace

void applyDarkTheme()
{
    ImGuiStyle& style = ImGui::GetStyle();
    ImVec4* colors = style.Colors;

    ImGui::StyleColorsDark(&style);
    applyCommonStyle();

    colors[ImGuiCol_Text]                 = ImVec4(0.93f, 0.95f, 0.97f, 1.00f);
    colors[ImGuiCol_TextDisabled]         = ImVec4(0.53f, 0.58f, 0.63f, 1.00f);
    colors[ImGuiCol_WindowBg]             = ImVec4(0.08f, 0.10f, 0.13f, 1.00f);
    colors[ImGuiCol_ChildBg]              = ImVec4(0.10f, 0.12f, 0.16f, 1.00f);
    colors[ImGuiCol_PopupBg]              = ImVec4(0.09f, 0.11f, 0.15f, 1.00f);
    colors[ImGuiCol_Border]               = ImVec4(0.22f, 0.28f, 0.34f, 0.70f);
    colors[ImGuiCol_BorderShadow]         = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);

    colors[ImGuiCol_FrameBg]              = ImVec4(0.14f, 0.17f, 0.22f, 1.00f);
    colors[ImGuiCol_FrameBgHovered]       = ImVec4(0.18f, 0.23f, 0.29f, 1.00f);
    colors[ImGuiCol_FrameBgActive]        = ImVec4(0.21f, 0.28f, 0.35f, 1.00f);

    colors[ImGuiCol_TitleBg]              = ImVec4(0.07f, 0.09f, 0.12f, 1.00f);
    colors[ImGuiCol_TitleBgActive]        = ImVec4(0.10f, 0.13f, 0.17f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed]     = ImVec4(0.07f, 0.09f, 0.12f, 0.80f);
    colors[ImGuiCol_MenuBarBg]            = ImVec4(0.09f, 0.12f, 0.16f, 1.00f);

    colors[ImGuiCol_Button]               = ImVec4(0.15f, 0.39f, 0.46f, 1.00f);
    colors[ImGuiCol_ButtonHovered]        = ImVec4(0.20f, 0.52f, 0.60f, 1.00f);
    colors[ImGuiCol_ButtonActive]         = ImVec4(0.24f, 0.60f, 0.68f, 1.00f);

    colors[ImGuiCol_Header]               = ImVec4(0.14f, 0.32f, 0.38f, 0.85f);
    colors[ImGuiCol_HeaderHovered]        = ImVec4(0.20f, 0.46f, 0.53f, 0.92f);
    colors[ImGuiCol_HeaderActive]         = ImVec4(0.25f, 0.56f, 0.63f, 1.00f);

    colors[ImGuiCol_CheckMark]            = ImVec4(0.89f, 0.73f, 0.31f, 1.00f);
    colors[ImGuiCol_SliderGrab]           = ImVec4(0.31f, 0.74f, 0.77f, 0.95f);
    colors[ImGuiCol_SliderGrabActive]     = ImVec4(0.42f, 0.82f, 0.84f, 1.00f);

    captureBaseStyle();
    applyScaledStyle();
}

void applyLightTheme()
{
    ImGuiStyle& style = ImGui::GetStyle();
    ImVec4* colors = style.Colors;

    ImGui::StyleColorsLight(&style);
    applyCommonStyle();

    colors[ImGuiCol_Text]                 = ImVec4(0.10f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_TextDisabled]         = ImVec4(0.45f, 0.47f, 0.50f, 1.00f);
    colors[ImGuiCol_WindowBg]             = ImVec4(0.95f, 0.95f, 0.96f, 1.00f);
    colors[ImGuiCol_ChildBg]              = ImVec4(0.92f, 0.93f, 0.94f, 1.00f);
    colors[ImGuiCol_PopupBg]              = ImVec4(0.97f, 0.97f, 0.98f, 1.00f);
    colors[ImGuiCol_Border]               = ImVec4(0.72f, 0.74f, 0.78f, 0.70f);
    colors[ImGuiCol_BorderShadow]         = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);

    colors[ImGuiCol_FrameBg]              = ImVec4(0.88f, 0.89f, 0.91f, 1.00f);
    colors[ImGuiCol_FrameBgHovered]       = ImVec4(0.82f, 0.84f, 0.87f, 1.00f);
    colors[ImGuiCol_FrameBgActive]        = ImVec4(0.76f, 0.79f, 0.83f, 1.00f);

    colors[ImGuiCol_TitleBg]              = ImVec4(0.88f, 0.89f, 0.91f, 1.00f);
    colors[ImGuiCol_TitleBgActive]        = ImVec4(0.82f, 0.84f, 0.87f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed]     = ImVec4(0.88f, 0.89f, 0.91f, 0.80f);
    colors[ImGuiCol_MenuBarBg]            = ImVec4(0.90f, 0.91f, 0.93f, 1.00f);

    colors[ImGuiCol_Button]               = ImVec4(0.22f, 0.52f, 0.60f, 1.00f);
    colors[ImGuiCol_ButtonHovered]        = ImVec4(0.28f, 0.60f, 0.68f, 1.00f);
    colors[ImGuiCol_ButtonActive]         = ImVec4(0.18f, 0.45f, 0.52f, 1.00f);

    colors[ImGuiCol_Header]               = ImVec4(0.72f, 0.84f, 0.88f, 0.85f);
    colors[ImGuiCol_HeaderHovered]        = ImVec4(0.62f, 0.78f, 0.84f, 0.92f);
    colors[ImGuiCol_HeaderActive]         = ImVec4(0.52f, 0.72f, 0.78f, 1.00f);

    colors[ImGuiCol_CheckMark]            = ImVec4(0.18f, 0.45f, 0.52f, 1.00f);
    colors[ImGuiCol_SliderGrab]           = ImVec4(0.22f, 0.52f, 0.60f, 0.95f);
    colors[ImGuiCol_SliderGrabActive]     = ImVec4(0.18f, 0.45f, 0.52f, 1.00f);

    style.FrameBorderSize = 1.0f;

    captureBaseStyle();
    applyScaledStyle();
}

void initImGui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    applyDarkTheme();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
    gImGuiBackendsReady = true;

    updateImGuiScale(window);
}

void updateImGuiScale(GLFWwindow* window)
{
    if (!gBaseStyleCaptured)
        return;

    const float newScale = computeUiScale(window);
    if (newScale == gUiScale)
        return;

    gUiScale = newScale;
    applyScaledStyle();
    rebuildFonts();
}

void shutdownImGui()
{
    gImGuiBackendsReady = false;
    gBaseStyleCaptured = false;
    gUiScale = 1.0f;
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}
