#pragma once
// The graphics API the desktop draws with: Vulkan by default, OpenGL as the
// fallback (when Vulkan cannot start, or when chosen in Settings, with
// ATOMFORGE_RENDERER=opengl or with --renderer opengl).

#include <string>

struct GLFWwindow;

namespace graphics
{
enum class Backend { OpenGL, Vulkan };

// The backend in use (set once the main window is created).
Backend backend();
void setBackend(Backend backend);
inline bool usingVulkan() { return backend() == Backend::Vulkan; }
inline bool usingOpenGL() { return backend() == Backend::OpenGL; }
const char* backendName(Backend backend);

// True when this build has the Vulkan renderer.
bool vulkanBuilt();

// The backend asked for: the command line, then ATOMFORGE_RENDERER, then the
// saved preference, else Vulkan (when built).
Backend requestedBackend();
// The --renderer value of the command line ("" when absent).
void setCommandLineRenderer(const std::string& value);

// The saved choice, applied at the next start ("renderer=" in the preferences file).
Backend preferredBackend();
void savePreferredBackend(Backend backend);

// Why the requested backend was not used ("" when it was), for a notice.
const std::string& fallbackReason();
void setFallbackReason(const std::string& reason);

// The main window (with Vulkan it has no OpenGL context to find it by).
GLFWwindow* mainWindow();
void setMainWindow(GLFWwindow* window);
}
