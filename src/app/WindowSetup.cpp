#include "app/WindowSetup.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>

#include "Camera.h"

struct SplashScreen
{
    GLFWwindow* window = nullptr;
    int width = 560;
    int height = 320;
};

namespace {

void renderSplash(SplashScreen* splash, float progress);

const unsigned char kIconGrid[4] = { 46, 94, 119, 235 };
const unsigned char kIconNode[4] = { 63, 110, 134, 255 };
const unsigned char kIconCore[4] = { 18, 63, 87, 255 };
const unsigned char kIconPlus[4] = { 46, 94, 119, 255 };

std::array<unsigned char, 64 * 64 * 4> makeGridCoreIconPixels()
{
    std::array<unsigned char, 64 * 64 * 4> pixels{};

    auto putPixel = [&](int x, int y, const unsigned char color[4]) {
        if (x < 0 || x >= 64 || y < 0 || y >= 64)
            return;
        const size_t idx = static_cast<size_t>((y * 64 + x) * 4);
        pixels[idx + 0] = color[0];
        pixels[idx + 1] = color[1];
        pixels[idx + 2] = color[2];
        pixels[idx + 3] = color[3];
    };

    auto fillRect = [&](int x0, int y0, int x1, int y1, const unsigned char color[4]) {
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x)
                putPixel(x, y, color);
    };

    auto fillCircle = [&](int cx, int cy, int r, const unsigned char color[4]) {
        const int rr = r * r;
        for (int y = cy - r; y <= cy + r; ++y)
        {
            for (int x = cx - r; x <= cx + r; ++x)
            {
                const int dx = x - cx;
                const int dy = y - cy;
                if (dx * dx + dy * dy <= rr)
                    putPixel(x, y, color);
            }
        }
    };

    // Keep icon background fully transparent.

    // Thicker grid bonds so the icon stays legible at 16-32px.
    fillRect(8, 7, 56, 11, kIconGrid);
    fillRect(8, 30, 56, 34, kIconGrid);
    fillRect(8, 53, 56, 57, kIconGrid);

    fillRect(7, 8, 11, 56, kIconGrid);
    fillRect(30, 8, 34, 56, kIconGrid);
    fillRect(53, 8, 57, 56, kIconGrid);

    fillCircle(10, 10, 6, kIconNode);
    fillCircle(54, 10, 6, kIconNode);
    fillCircle(10, 54, 6, kIconNode);
    fillCircle(54, 54, 6, kIconNode);

    fillCircle(32, 32, 10, kIconCore);
    fillRect(16, 29, 48, 35, kIconPlus);
    fillRect(29, 16, 35, 48, kIconPlus);

    return pixels;
}

void applyWindowIcon(GLFWwindow* window)
{
    if (!window)
        return;

    static const std::array<unsigned char, 64 * 64 * 4> kIconPixels = makeGridCoreIconPixels();
    GLFWimage image{};
    image.width = 64;
    image.height = 64;
    image.pixels = const_cast<unsigned char*>(kIconPixels.data());
    glfwSetWindowIcon(window, 1, &image);
}

bool ensureGlfwInitialized()
{
    return glfwInit() == GLFW_TRUE;
}

void fillRect(GLFWwindow* window, int x, int y, int width, int height, float r, float g, float b)
{
    if (!window || width <= 0 || height <= 0)
        return;

    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
    if (framebufferWidth <= 0 || framebufferHeight <= 0)
        return;

    x = std::max(0, std::min(x, framebufferWidth));
    y = std::max(0, std::min(y, framebufferHeight));
    width = std::min(width, framebufferWidth - x);
    height = std::min(height, framebufferHeight - y);
    if (width <= 0 || height <= 0)
        return;

    glEnable(GL_SCISSOR_TEST);
    glScissor(x, y, width, height);
    glClearColor(r, g, b, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_SCISSOR_TEST);
}

void fillRectTopLeft(GLFWwindow* window, int x, int y, int width, int height, float r, float g, float b)
{
    if (!window)
        return;

    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
    fillRect(window, x, framebufferHeight - y - height, width, height, r, g, b);
}

void drawGlyph(GLFWwindow* window,
               int x,
               int y,
               int scale,
               const std::array<const char*, 7>& rows,
               float r,
               float g,
               float b)
{
    for (int row = 0; row < 7; ++row)
    {
        for (int col = 0; rows[row][col] != '\0'; ++col)
        {
            if (rows[row][col] == ' ')
                continue;

            fillRectTopLeft(window, x + col * scale, y + row * scale, scale - 1, scale - 1, r, g, b);
        }
    }
}

void drawAtomForgeLogo(GLFWwindow* window, int x, int y, int scale)
{
    static const std::array<const char*, 7> kA = {" ### ", "#   #", "#   #", "#####", "#   #", "#   #", "#   #"};
    static const std::array<const char*, 7> kT = {"#####", "  #  ", "  #  ", "  #  ", "  #  ", "  #  ", "  #  "};
    static const std::array<const char*, 7> kO = {" ### ", "#   #", "#   #", "#   #", "#   #", "#   #", " ### "};
    static const std::array<const char*, 7> kM = {"#   #", "## ##", "# # #", "#   #", "#   #", "#   #", "#   #"};
    static const std::array<const char*, 7> kF = {"#####", "#    ", "#    ", "#### ", "#    ", "#    ", "#    "};
    static const std::array<const char*, 7> kR = {"#### ", "#   #", "#   #", "#### ", "# #  ", "#  # ", "#   #"};
    static const std::array<const char*, 7> kG = {" ### ", "#   #", "#    ", "# ###", "#   #", "#   #", " ### "};
    static const std::array<const char*, 7> kE = {"#####", "#    ", "#    ", "#### ", "#    ", "#    ", "#####"};

    const std::array<const std::array<const char*, 7>*, 9> glyphs = {
        &kA, &kT, &kO, &kM, &kF, &kO, &kR, &kG, &kE
    };

    int cursorX = x;
    for (const auto* glyph : glyphs)
    {
        drawGlyph(window, cursorX, y, scale, *glyph, 0.95f, 0.96f, 0.98f);
        cursorX += scale * 6;
    }
}

void drawSplashBadge(GLFWwindow* window, int x, int y)
{
    fillRectTopLeft(window, x, y, 88, 88, 0.20f, 0.43f, 0.68f);
    fillRectTopLeft(window, x + 10, y + 10, 68, 68, 0.08f, 0.11f, 0.16f);
    fillRectTopLeft(window, x + 21, y + 16, 14, 56, 0.38f, 0.76f, 0.92f);
    fillRectTopLeft(window, x + 39, y + 16, 14, 56, 0.95f, 0.45f, 0.25f);
    fillRectTopLeft(window, x + 57, y + 16, 10, 56, 0.94f, 0.82f, 0.28f);
    fillRectTopLeft(window, x + 15, y + 42, 58, 8, 0.86f, 0.89f, 0.93f);
}

void renderSplash(SplashScreen* splash, float progress)
{
    if (!splash || !splash->window)
        return;

    const float clampedProgress = std::max(0.0f, std::min(progress, 1.0f));
    GLFWwindow* previousContext = glfwGetCurrentContext();
    glfwMakeContextCurrent(splash->window);

    int framebufferWidth = 0;
    int framebufferHeight = 0;
    glfwGetFramebufferSize(splash->window, &framebufferWidth, &framebufferHeight);
    glViewport(0, 0, framebufferWidth, framebufferHeight);
    glDisable(GL_DEPTH_TEST);
    glClearColor(0.07f, 0.09f, 0.13f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    fillRectTopLeft(splash->window, 0, 0, framebufferWidth, framebufferHeight, 0.07f, 0.09f, 0.13f);
    fillRectTopLeft(splash->window, 0, 0, framebufferWidth, 12, 0.19f, 0.45f, 0.76f);
    fillRectTopLeft(splash->window, 0, 12, framebufferWidth, 28, 0.11f, 0.14f, 0.20f);
    fillRectTopLeft(splash->window, 36, 86, framebufferWidth - 72, 2, 0.24f, 0.29f, 0.38f);
    fillRectTopLeft(splash->window, 36, framebufferHeight - 68, framebufferWidth - 72, 1, 0.18f, 0.21f, 0.29f);

    drawSplashBadge(splash->window, 38, 58);
    drawAtomForgeLogo(splash->window, 146, 74, 8);

    fillRectTopLeft(splash->window, 146, 150, 228, 6, 0.24f, 0.30f, 0.39f);
    fillRectTopLeft(splash->window, 146, 150, static_cast<int>(228.0f * clampedProgress), 6, 0.34f, 0.77f, 0.91f);
    fillRectTopLeft(splash->window, 38, 246, framebufferWidth - 76, 24, 0.10f, 0.13f, 0.19f);
    fillRectTopLeft(splash->window, 38, 246, static_cast<int>((framebufferWidth - 76) * clampedProgress), 24, 0.20f, 0.43f, 0.68f);

    glfwSwapBuffers(splash->window);
    glfwPollEvents();
    glfwMakeContextCurrent(previousContext);
}

GLFWwindow* createInitialWindow(int contextMajor, int contextMinor, int profile)
{
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, contextMajor);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, contextMinor);
    glfwWindowHint(GLFW_OPENGL_PROFILE, profile);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_SAMPLES, 4); // 4× MSAA

    glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);
    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    return glfwCreateWindow(1280, 800, "AtomForge", nullptr, nullptr);
}

} // namespace

SplashScreen* createSplashScreen()
{
    if (!ensureGlfwInitialized())
        return nullptr;

    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_FLOATING, GLFW_TRUE);

    SplashScreen* splash = new SplashScreen();
    splash->window = glfwCreateWindow(splash->width, splash->height, "AtomForge", nullptr, nullptr);
    if (!splash->window)
    {
        delete splash;
        return nullptr;
    }

    applyWindowIcon(splash->window);

    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
    if (primaryMonitor)
    {
        int workX = 0;
        int workY = 0;
        int workWidth = 0;
        int workHeight = 0;
        glfwGetMonitorWorkarea(primaryMonitor, &workX, &workY, &workWidth, &workHeight);
        if (workWidth > 0 && workHeight > 0)
        {
            const int posX = workX + (workWidth - splash->width) / 2;
            const int posY = workY + (workHeight - splash->height) / 2;
            glfwSetWindowPos(splash->window, posX, posY);
        }
    }

    updateSplashScreen(splash, 0.05f, "Starting");
    return splash;
}

void updateSplashScreen(SplashScreen* splash, float progress, const char* status)
{
    if (!splash || !splash->window)
        return;

    const std::string windowTitle = status && status[0]
        ? (std::string("AtomForge - ") + status)
        : std::string("AtomForge");
    glfwSetWindowTitle(splash->window, windowTitle.c_str());
    renderSplash(splash, progress);
}

void destroySplashScreen(SplashScreen* splash)
{
    if (!splash)
        return;

    if (splash->window)
        glfwDestroyWindow(splash->window);
    delete splash;
}

GLFWwindow* createMainWindow()
{
    if (!ensureGlfwInitialized())
        return nullptr;

    GLFWwindow* window = createInitialWindow(4, 6, GLFW_OPENGL_CORE_PROFILE);
    if (!window)
    {
        // Fallback for systems without OpenGL 4.6 support.
        window = createInitialWindow(4, 3, GLFW_OPENGL_CORE_PROFILE);
    }

    if (!window)
    {
        // Final fallback for older systems.
        window = createInitialWindow(3, 3, GLFW_OPENGL_COMPAT_PROFILE);
        if (!window)
            return nullptr;
    }

    applyWindowIcon(window);

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
    {
        glfwDestroyWindow(window);
        return nullptr;
    }

    if (GLEW_VERSION_4_6)
        std::cout << "[OpenGL] Running with OpenGL 4.6 (compute path available)." << std::endl;
    else if (GLEW_VERSION_4_3)
        std::cout << "[OpenGL] Running with OpenGL 4.3+ (compute path available)." << std::endl;
    else
        std::cout << "[OpenGL] OpenGL 4.3 not available, using legacy instancing path." << std::endl;

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_DEPTH_CLAMP);  // Prevents near-plane from slicing sphere geometry (flat-cutoff artifact)
    return window;
}

void showMainWindow(GLFWwindow* window)
{
    if (!window)
        return;

    glfwShowWindow(window);
    glfwFocusWindow(window);
}

void configureCameraCallbacks(GLFWwindow* window, Camera& camera)
{
    Camera::instance = &camera;
    glfwSetMouseButtonCallback(window, Camera::mouseButton);
    glfwSetCursorPosCallback(window, Camera::cursor);
    glfwSetScrollCallback(window, Camera::scroll);
}