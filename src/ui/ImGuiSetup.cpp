#include "ImGuiSetup.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <filesystem>
#include <cmath>

namespace
{
constexpr float kBaseFontSizePixels = 18.0f;

ImGuiStyle gBaseStyle;
float gUiScale = 1.0f;
float gInterfaceScale = 1.0f;
bool gBaseStyleCaptured = false;
bool gImGuiBackendsReady = false;

float clampUiScale(float scale)
{
    return std::clamp(scale, 0.85f, 4.0f);
}

float computeUiScale(GLFWwindow* window)
{
    if (!window)
        return gInterfaceScale;

    float contentScaleX = 1.0f;
    float contentScaleY = 1.0f;
    glfwGetWindowContentScale(window, &contentScaleX, &contentScaleY);

    int windowWidth = 0;
    int windowHeight = 0;
    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);

    // A minimized window can temporarily report a zero-sized framebuffer.
    if (windowWidth <= 0 || windowHeight <= 0 || framebufferWidth <= 0 || framebufferHeight <= 0)
        return gUiScale;

    float framebufferScaleX = 1.0f;
    float framebufferScaleY = 1.0f;
    if (windowWidth > 0)
        framebufferScaleX = (float)framebufferWidth / (float)windowWidth;
    if (windowHeight > 0)
        framebufferScaleY = (float)framebufferHeight / (float)windowHeight;

    // GLFW coordinates already include framebuffer scaling on Retina displays.
    // Apply only the remaining logical content scale to avoid doubling UI sizes.
    const float scale = std::max(
        contentScaleX / framebufferScaleX, contentScaleY / framebufferScaleY);
    return clampUiScale(scale * gInterfaceScale);
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
    // ImGui 1.92 bakes glyphs on demand at the requested size. Keep a stable
    // base size and use the style multiplier rather than recreating the atlas.
    scaled.FontSizeBase = kBaseFontSizePixels;
    scaled.FontScaleMain = gUiScale;
    scaled.DisplayWindowPadding = ImVec2(0.0f, 0.0f);
    scaled.DisplaySafeAreaPadding = ImVec2(0.0f, 0.0f);
    ImGui::GetStyle() = scaled;
}

void rebuildFonts()
{
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();

    ImFontConfig fontConfig;
    fontConfig.SizePixels = kBaseFontSizePixels;
    io.FontGlobalScale = 1.0f;

    // Prefer the platform UI font, with a built-in fallback on minimal systems.
#ifdef _WIN32
    static const char* kSystemFontCandidates[] = {
        "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/arial.ttf",
        nullptr
    };
#else
    static const char* kSystemFontCandidates[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        nullptr
    };
#endif
    ImFontConfig mergeConfig;
    mergeConfig.MergeMode   = false;
    mergeConfig.SizePixels  = kBaseFontSizePixels;
    mergeConfig.OversampleH = 2;
    mergeConfig.OversampleV = 1;
    bool loaded = false;
    static const ImWchar ranges[] = {0x20,0xFF,0x2190,0x2192,0};
    for (int i = 0; kSystemFontCandidates[i]; ++i)
    {
        std::error_code error;
        if (!std::filesystem::is_regular_file(kSystemFontCandidates[i],error)) continue;
        if (io.Fonts->AddFontFromFileTTF(
                kSystemFontCandidates[i],
                kBaseFontSizePixels,
                &mergeConfig,
                ranges))
        { loaded = true; break; }
    }
    if (!loaded) io.Fonts->AddFontDefault(&fontConfig);

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

    style.WindowPadding = ImVec2(18.0f, 16.0f);
    style.FramePadding = ImVec2(12.0f, 7.0f);
    style.CellPadding = ImVec2(8.0f, 6.0f);
    style.ItemSpacing = ImVec2(10.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(8.0f, 6.0f);
    style.TouchExtraPadding = ImVec2(1.0f, 1.0f);
    style.IndentSpacing = 20.0f;
    style.ScrollbarSize = 14.0f;
    style.GrabMinSize = 10.0f;

    style.WindowRounding = 8.0f;
    style.ChildRounding = 6.0f;
    style.FrameRounding = 5.0f;
    style.PopupRounding = 8.0f;
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
    style = ImGuiStyle();
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

    colors[ImGuiCol_Button]               = ImVec4(0.18f, 0.25f, 0.31f, 1.00f);
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
    style = ImGuiStyle();
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

    colors[ImGuiCol_Button]               = ImVec4(0.83f, 0.89f, 0.93f, 1.00f);
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
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr;

    applyLightTheme();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
    gImGuiBackendsReady = true;

    // Always build the font atlas here so font selection in
    // rebuildFonts() runs even on non-HiDPI displays where updateImGuiScale
    // would return early (newScale == gUiScale == 1.0).
    rebuildFonts();

    updateImGuiScale(window);
}

void updateImGuiScale(GLFWwindow* window)
{
    if (!gBaseStyleCaptured)
        return;

    const float newScale = computeUiScale(window);
    if (std::abs(newScale - gUiScale) < 0.01f)
        return;

    gUiScale = newScale;
    applyScaledStyle();
}

float interfaceScale() { return gInterfaceScale; }

void setInterfaceScale(float multiplier)
{
    if (std::isfinite(multiplier)) gInterfaceScale = std::clamp(multiplier, 0.85f, 2.0f);
}

void shutdownImGui()
{
    gImGuiBackendsReady = false;
    gBaseStyleCaptured = false;
    gUiScale = 1.0f;
    gInterfaceScale = 1.0f;
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}
