#include "ImGuiSetup.h"
#include "ui/Fonts.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "graphics/GraphicsBackend.h"
#ifdef ATOMFORGE_VULKAN
#include "imgui_impl_vulkan.h"
#include "graphics/vulkan/VulkanContext.h"
#endif

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace
{
constexpr float kBaseFontSizePixels = 18.0f;
#ifdef ATOMFORGE_VULKAN
VkInstance gVulkanInstance = VK_NULL_HANDLE;
#endif

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
    io.FontGlobalScale = 1.0f;
    uiFonts::heading = nullptr;
    uiFonts::icons = false;

    // The platform UI font (regular and semibold), with a built-in fallback.
#ifdef _WIN32
    static const char* kTextFonts[] = {"C:/Windows/Fonts/segoeui.ttf", "C:/Windows/Fonts/arial.ttf", nullptr};
    static const char* kHeadingFonts[] = {"C:/Windows/Fonts/seguisb.ttf", "C:/Windows/Fonts/segoeuib.ttf", "C:/Windows/Fonts/arialbd.ttf", nullptr};
    static const char* kIconFonts[] = {"C:/Windows/Fonts/SegoeIcons.ttf", "C:/Windows/Fonts/segmdl2.ttf", nullptr};
#else
    static const char* kTextFonts[] = {"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                                       "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
                                       "/usr/share/fonts/truetype/freefont/FreeSans.ttf", nullptr};
    static const char* kHeadingFonts[] = {"/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                                          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
                                          "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf", nullptr};
    static const char* kIconFonts[] = {nullptr};
#endif
    const auto firstExisting = [](const char* const* candidates) -> const char* {
        for (int i = 0; candidates[i]; ++i) {
            std::error_code error;
            if (std::filesystem::is_regular_file(candidates[i], error)) return candidates[i];
        }
        return nullptr;
    };
    static const ImWchar textRanges[] = {0x20, 0xFF, 0x2190, 0x2192, 0x2022, 0x2022, 0x2026, 0x2026, 0x00B0, 0x00B0, 0x212B, 0x212B, 0x0391, 0x03C9,
                                         0x2013, 0x2014, 0x2070, 0x209F, 0x1D62, 0x1D62, 0x2C7C, 0x2C7C, 0x2212, 0x2212, 0x2248, 0x2248, 0x2260, 0x2265, 0};
    static const ImWchar iconRanges[] = {0xE700, 0xF8FF, 0};
    ImFontConfig text;
    text.OversampleH = 2;
    text.OversampleV = 1;
    ImFont* body = nullptr;
    if (const char* path = firstExisting(kTextFonts)) body = io.Fonts->AddFontFromFileTTF(path, kBaseFontSizePixels, &text, textRanges);
    if (!body) {
        ImFontConfig fallback;
        fallback.SizePixels = kBaseFontSizePixels;
        io.Fonts->AddFontDefault(&fallback);
    } else if (const char* path = firstExisting(kIconFonts)) {
        // Icons share the text font, centred on its line.
        ImFontConfig icon;
        icon.MergeMode = true;
        icon.GlyphOffset = ImVec2(0.0f, 2.0f);
        icon.GlyphMinAdvanceX = kBaseFontSizePixels;
        uiFonts::icons = io.Fonts->AddFontFromFileTTF(path, kBaseFontSizePixels - 2.0f, &icon, iconRanges) != nullptr;
    }
    if (const char* path = firstExisting(kHeadingFonts)) {
        ImFontConfig headingConfig;
        headingConfig.OversampleH = 2;
        uiFonts::heading = io.Fonts->AddFontFromFileTTF(path, kBaseFontSizePixels, &headingConfig, textRanges);
    }

    // The Vulkan backend updates its font texture by itself (ImGui 1.92 textures).
    if (gImGuiBackendsReady && graphics::usingOpenGL())
    {
        ImGui_ImplOpenGL3_DestroyDeviceObjects();
        io.Fonts->Build();
        ImGui_ImplOpenGL3_CreateDeviceObjects();
    }
}

// Geometry shared by both themes: compact, evenly spaced, softly rounded.
void applyCommonStyle()
{
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowPadding = ImVec2(16.0f, 14.0f);
    style.FramePadding = ImVec2(10.0f, 6.0f);
    style.CellPadding = ImVec2(8.0f, 5.0f);
    style.ItemSpacing = ImVec2(8.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(7.0f, 5.0f);
    style.TouchExtraPadding = ImVec2(0.0f, 0.0f);
    style.IndentSpacing = 18.0f;
    style.ScrollbarSize = 12.0f;
    style.GrabMinSize = 8.0f;
    style.SeparatorTextBorderSize = 1.0f;
    style.SeparatorTextPadding = ImVec2(0.0f, 4.0f);

    style.WindowRounding = 10.0f;
    style.ChildRounding = 8.0f;
    style.FrameRounding = 6.0f;
    style.PopupRounding = 8.0f;
    style.ScrollbarRounding = 6.0f;
    style.GrabRounding = 6.0f;
    style.TabRounding = 6.0f;

    style.WindowBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.PopupBorderSize = 1.0f;
    style.FrameBorderSize = 1.0f;
    style.TabBorderSize = 0.0f;

    style.WindowTitleAlign = ImVec2(0.0f, 0.5f);
    style.WindowMenuButtonPosition = ImGuiDir_None;
    style.SelectableTextAlign = ImVec2(0.0f, 0.5f);
}

// Colours from a small palette: neutral surfaces, one accent, quiet borders.
struct Palette
{
    ImVec4 text, textMuted, window, surface, popup, border, frame, frameHover, frameActive,
           title, titleActive, menuBar, button, buttonHover, buttonActive, accent, accentStrong, separator, dim;
};

ImVec4 withAlpha(ImVec4 c, float a) { c.w = a; return c; }

void applyPalette(const Palette& p, bool light)
{
    ImVec4* c = ImGui::GetStyle().Colors;
    c[ImGuiCol_Text] = p.text;
    c[ImGuiCol_TextDisabled] = p.textMuted;
    c[ImGuiCol_WindowBg] = p.window;
    c[ImGuiCol_ChildBg] = p.surface;
    c[ImGuiCol_PopupBg] = p.popup;
    c[ImGuiCol_Border] = p.border;
    c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_FrameBg] = p.frame;
    c[ImGuiCol_FrameBgHovered] = p.frameHover;
    c[ImGuiCol_FrameBgActive] = p.frameActive;
    c[ImGuiCol_TitleBg] = p.title;
    c[ImGuiCol_TitleBgActive] = p.titleActive;
    c[ImGuiCol_TitleBgCollapsed] = p.title;
    c[ImGuiCol_MenuBarBg] = p.menuBar;
    c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ScrollbarGrab] = withAlpha(p.textMuted, 0.35f);
    c[ImGuiCol_ScrollbarGrabHovered] = withAlpha(p.textMuted, 0.55f);
    c[ImGuiCol_ScrollbarGrabActive] = withAlpha(p.textMuted, 0.75f);
    c[ImGuiCol_CheckMark] = p.accent;
    c[ImGuiCol_SliderGrab] = withAlpha(p.accent, 0.65f);  // the value text stays readable over it
    c[ImGuiCol_SliderGrabActive] = p.accentStrong;
    c[ImGuiCol_Button] = p.button;
    c[ImGuiCol_ButtonHovered] = p.buttonHover;
    c[ImGuiCol_ButtonActive] = p.buttonActive;
    c[ImGuiCol_Header] = withAlpha(p.accent, light ? 0.14f : 0.24f);
    c[ImGuiCol_HeaderHovered] = withAlpha(p.accent, light ? 0.22f : 0.34f);
    c[ImGuiCol_HeaderActive] = withAlpha(p.accent, light ? 0.32f : 0.46f);
    c[ImGuiCol_Separator] = p.separator;
    c[ImGuiCol_SeparatorHovered] = withAlpha(p.accent, 0.70f);
    c[ImGuiCol_SeparatorActive] = p.accent;
    c[ImGuiCol_ResizeGrip] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_ResizeGripHovered] = withAlpha(p.accent, 0.45f);
    c[ImGuiCol_ResizeGripActive] = withAlpha(p.accent, 0.80f);
    c[ImGuiCol_InputTextCursor] = p.text;
    c[ImGuiCol_Tab] = p.button;
    c[ImGuiCol_TabHovered] = p.buttonHover;
    c[ImGuiCol_TabSelected] = p.surface;
    c[ImGuiCol_TabSelectedOverline] = p.accent;
    c[ImGuiCol_TabDimmed] = p.button;
    c[ImGuiCol_TabDimmedSelected] = p.surface;
    c[ImGuiCol_TabDimmedSelectedOverline] = withAlpha(p.accent, 0.5f);
    c[ImGuiCol_PlotLines] = p.accent;
    c[ImGuiCol_PlotLinesHovered] = p.accentStrong;
    c[ImGuiCol_PlotHistogram] = p.accent;
    c[ImGuiCol_PlotHistogramHovered] = p.accentStrong;
    c[ImGuiCol_TableHeaderBg] = p.button;
    c[ImGuiCol_TableBorderStrong] = p.border;
    c[ImGuiCol_TableBorderLight] = p.separator;
    c[ImGuiCol_TableRowBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_TableRowBgAlt] = light ? ImVec4(0, 0, 0, 0.025f) : ImVec4(1, 1, 1, 0.03f);
    c[ImGuiCol_TextLink] = p.accent;
    c[ImGuiCol_TextSelectedBg] = withAlpha(p.accent, 0.28f);
    c[ImGuiCol_TreeLines] = p.separator;
    c[ImGuiCol_DragDropTarget] = p.accent;
    c[ImGuiCol_NavCursor] = p.accent;
    c[ImGuiCol_NavWindowingHighlight] = withAlpha(p.accent, 0.7f);
    c[ImGuiCol_NavWindowingDimBg] = p.dim;
    c[ImGuiCol_ModalWindowDimBg] = p.dim;
}
} // namespace

void applyDarkTheme()
{
    ImGuiStyle& style = ImGui::GetStyle();
    style = ImGuiStyle();
    ImGui::StyleColorsDark(&style);
    applyCommonStyle();
    Palette p;
    p.text = ImVec4(0.90f, 0.92f, 0.95f, 1.00f);
    p.textMuted = ImVec4(0.56f, 0.60f, 0.66f, 1.00f);
    p.window = ImVec4(0.105f, 0.115f, 0.135f, 1.00f);
    p.surface = ImVec4(0.125f, 0.137f, 0.160f, 1.00f);
    p.popup = ImVec4(0.135f, 0.147f, 0.172f, 1.00f);
    p.border = ImVec4(0.235f, 0.260f, 0.300f, 1.00f);
    p.frame = ImVec4(0.155f, 0.170f, 0.198f, 1.00f);
    p.frameHover = ImVec4(0.185f, 0.203f, 0.237f, 1.00f);
    p.frameActive = ImVec4(0.215f, 0.240f, 0.282f, 1.00f);
    p.title = ImVec4(0.115f, 0.126f, 0.148f, 1.00f);
    p.titleActive = ImVec4(0.140f, 0.153f, 0.180f, 1.00f);
    p.menuBar = ImVec4(0.120f, 0.131f, 0.154f, 1.00f);
    p.button = ImVec4(0.180f, 0.197f, 0.230f, 1.00f);
    p.buttonHover = ImVec4(0.225f, 0.248f, 0.290f, 1.00f);
    p.buttonActive = ImVec4(0.265f, 0.295f, 0.345f, 1.00f);
    p.accent = ImVec4(0.33f, 0.63f, 0.93f, 1.00f);
    p.accentStrong = ImVec4(0.45f, 0.72f, 0.98f, 1.00f);
    p.separator = ImVec4(0.215f, 0.238f, 0.275f, 1.00f);
    p.dim = ImVec4(0.02f, 0.03f, 0.05f, 0.55f);
    applyPalette(p, false);
    captureBaseStyle();
    applyScaledStyle();
}

void applyLightTheme()
{
    ImGuiStyle& style = ImGui::GetStyle();
    style = ImGuiStyle();
    ImGui::StyleColorsLight(&style);
    applyCommonStyle();
    Palette p;
    p.text = ImVec4(0.11f, 0.12f, 0.14f, 1.00f);
    p.textMuted = ImVec4(0.40f, 0.43f, 0.48f, 1.00f);  // at least 4.5:1 on the window (WCAG AA)
    p.window = ImVec4(0.965f, 0.969f, 0.976f, 1.00f);
    p.surface = ImVec4(1.000f, 1.000f, 1.000f, 1.00f);
    p.popup = ImVec4(1.000f, 1.000f, 1.000f, 1.00f);
    p.border = ImVec4(0.800f, 0.820f, 0.850f, 1.00f);
    p.frame = ImVec4(1.000f, 1.000f, 1.000f, 1.00f);
    p.frameHover = ImVec4(0.945f, 0.955f, 0.972f, 1.00f);
    p.frameActive = ImVec4(0.905f, 0.930f, 0.965f, 1.00f);
    p.title = ImVec4(0.930f, 0.937f, 0.950f, 1.00f);
    p.titleActive = ImVec4(0.900f, 0.912f, 0.932f, 1.00f);
    p.menuBar = ImVec4(0.975f, 0.978f, 0.984f, 1.00f);
    p.button = ImVec4(0.935f, 0.945f, 0.960f, 1.00f);
    p.buttonHover = ImVec4(0.885f, 0.905f, 0.935f, 1.00f);
    p.buttonActive = ImVec4(0.830f, 0.865f, 0.915f, 1.00f);
    p.accent = ImVec4(0.12f, 0.44f, 0.75f, 1.00f);
    p.accentStrong = ImVec4(0.08f, 0.35f, 0.64f, 1.00f);
    p.separator = ImVec4(0.865f, 0.880f, 0.900f, 1.00f);
    p.dim = ImVec4(0.10f, 0.12f, 0.16f, 0.30f);
    applyPalette(p, true);
    captureBaseStyle();
    applyScaledStyle();
}

// Where the window and docking layout is kept between sessions; null (none)
// in the automated modes, so that they always start from the default layout.
static const char* layoutFile()
{
    for (const char* mode : {"ATOMFORGE_UI_TOUR", "ATOMFORGE_TUTORIAL_CAPTURE", "ATOMFORGE_STEP_DIALOG_CHECK",
                             "ATOMFORGE_CLOUD_BENCHMARK", "ATOMFORGE_NO_LAYOUT"})
        if (const char* value = std::getenv(mode); value && *value) return nullptr;
    static std::string path;
    try {
        std::filesystem::path root;
#ifdef _WIN32
        if (const char* folder = std::getenv("LOCALAPPDATA")) root = std::filesystem::u8path(folder);
#else
        if (const char* folder = std::getenv("XDG_CONFIG_HOME")) root = std::filesystem::u8path(folder);
        else if (const char* home = std::getenv("HOME")) root = std::filesystem::u8path(home) / ".config";
#endif
        if (root.empty()) return nullptr;
        root /= "AtomForge";
        std::filesystem::create_directories(root);
        path = (root / "layout.ini").u8string();
    } catch (const std::exception&) {
        return nullptr;
    }
    return path.c_str();
}

void newImGuiRendererFrame()
{
#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        ImGui_ImplVulkan_SetMinImageCount(vkr::context().minImageCount());
        ImGui_ImplVulkan_NewFrame();
        return;
    }
#endif
    ImGui_ImplOpenGL3_NewFrame();
}

void renderImGuiDrawData()
{
#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        if (vkr::context().inFrame())
            ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), vkr::context().commandBuffer());
        return;
    }
#endif
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void initImGui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigWindowsMoveFromTitleBarOnly = true;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    // Tool windows (the pipeline, scientific tools, trajectory playback, ...)
    // dock to the sides of the view or into each other as tabs.
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigDockingWithShift = false;
    io.IniFilename = layoutFile();

    applyLightTheme();

#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        vkr::Context& ctx = vkr::context();
        ImGui_ImplGlfw_InitForVulkan(window, true);
        gVulkanInstance = ctx.instance();
        ImGui_ImplVulkan_LoadFunctions(VK_API_VERSION_1_1, [](const char* name, void* instance) {
            return vkload::getInstanceProcAddr()(*static_cast<VkInstance*>(instance), name);
        }, &gVulkanInstance);
        ImGui_ImplVulkan_InitInfo info{};
        info.ApiVersion = VK_API_VERSION_1_1;
        info.Instance = ctx.instance();
        info.PhysicalDevice = ctx.physicalDevice();
        info.Device = ctx.device();
        info.QueueFamily = ctx.queueFamily();
        info.Queue = ctx.queue();
        info.DescriptorPool = ctx.descriptorPool();
        info.MinImageCount = ctx.minImageCount();
        info.ImageCount = ctx.imageCount();
        info.PipelineInfoMain.RenderPass = ctx.mainRenderPass();
        info.PipelineInfoMain.MSAASamples = ctx.samples();
        ImGui_ImplVulkan_Init(&info);
    }
    else
#endif
    {
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 130");
    }
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

namespace
{
ThemeMode gThemeMode = ThemeMode::Light;
bool gPreferencesLoaded = false;

std::filesystem::path preferencesPath()
{
#ifdef _WIN32
    const char* home = std::getenv("USERPROFILE");
#else
    const char* home = std::getenv("HOME");
#endif
    return std::filesystem::u8path(home && *home ? home : ".") / ".atomforge_preferences.txt";
}

void savePreferences()
{
    // Only once the user's preferences are in effect, so that defaults
    // applied during start-up never overwrite them.
    if (!gPreferencesLoaded) return;
    std::ofstream out(preferencesPath(), std::ios::trunc);
    if (!out) return;
    out << "theme=" << (gThemeMode == ThemeMode::Dark ? "dark" : gThemeMode == ThemeMode::System ? "system" : "light") << '\n'
        << "interface_scale=" << gInterfaceScale << '\n';
}

// The operating system's app theme (light unless it asks for dark).
bool systemPrefersLight()
{
#ifdef _WIN32
    DWORD value = 1, size = sizeof(value);
    if (RegGetValueW(HKEY_CURRENT_USER, L"Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize",
                     L"AppsUseLightTheme", RRF_RT_REG_DWORD, nullptr, &value, &size) == ERROR_SUCCESS)
        return value != 0;
    return true;
#else
    const char* gtk = std::getenv("GTK_THEME");
    return !(gtk && std::string(gtk).find("dark") != std::string::npos);
#endif
}

void applyThemeMode()
{
    const bool light = gThemeMode == ThemeMode::System ? systemPrefersLight() : gThemeMode == ThemeMode::Light;
    if (light) applyLightTheme();
    else applyDarkTheme();
}
}

void setInterfaceScale(float multiplier)
{
    if (std::isfinite(multiplier)) gInterfaceScale = std::clamp(multiplier, 0.85f, 2.0f);
    savePreferences();
}

ThemeMode themeMode() { return gThemeMode; }

void setThemeMode(ThemeMode mode)
{
    gThemeMode = mode;
    applyThemeMode();
    savePreferences();
}

void loadAppearancePreferences()
{
    std::ifstream in(preferencesPath());
    for (std::string line; std::getline(in, line);) {
        const auto equals = line.find('=');
        if (equals == std::string::npos) continue;
        const std::string key = line.substr(0, equals), value = line.substr(equals + 1);
        if (key == "theme")
            gThemeMode = value == "dark" ? ThemeMode::Dark : value == "system" ? ThemeMode::System : ThemeMode::Light;
        else if (key == "interface_scale") {
            const float scale = std::strtof(value.c_str(), nullptr);
            if (std::isfinite(scale) && scale > 0) gInterfaceScale = std::clamp(scale, 0.85f, 2.0f);
        }
    }
    applyThemeMode();
    gPreferencesLoaded = true;
}

void updateSystemTheme()
{
    // A registry read every two seconds is enough to follow a change.
    if (gThemeMode != ThemeMode::System) return;
    static double lastCheck = -10.0;
    const double now = ImGui::GetTime();
    if (now - lastCheck < 2.0) return;
    lastCheck = now;
    const bool light = systemPrefersLight();
    const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    if (light != ((bg.x + bg.y + bg.z) / 3.0f > 0.5f)) applyThemeMode();
}

void shutdownImGui()
{
    gImGuiBackendsReady = false;
    gBaseStyleCaptured = false;
    gUiScale = 1.0f;
    gInterfaceScale = 1.0f;
#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        vkDeviceWaitIdle(vkr::context().device());
        ImGui_ImplVulkan_Shutdown();
    }
    else
#endif
        ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}
