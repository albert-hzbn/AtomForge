#include "graphics/GraphicsBackend.h"

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace graphics
{
namespace
{
Backend gBackend = Backend::OpenGL;
std::string gCommandLine;
std::string gFallbackReason;

std::filesystem::path preferencesFile()
{
#ifdef _WIN32
    const char* home = std::getenv("USERPROFILE");
#else
    const char* home = std::getenv("HOME");
#endif
    return std::filesystem::u8path(home && *home ? home : ".") / ".atomforge_preferences.txt";
}

bool parse(const std::string& text, Backend& out)
{
    std::string lower;
    for (char c : text) lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (lower == "vulkan" || lower == "vk") { out = Backend::Vulkan; return true; }
    if (lower == "opengl" || lower == "gl") { out = Backend::OpenGL; return true; }
    return false;
}
}

Backend backend() { return gBackend; }
void setBackend(Backend value) { gBackend = value; }
const char* backendName(Backend value) { return value == Backend::Vulkan ? "Vulkan" : "OpenGL"; }

bool vulkanBuilt()
{
#ifdef ATOMFORGE_VULKAN
    return true;
#else
    return false;
#endif
}

void setCommandLineRenderer(const std::string& value) { gCommandLine = value; }

Backend preferredBackend()
{
    Backend chosen = vulkanBuilt() ? Backend::Vulkan : Backend::OpenGL;
    std::ifstream in(preferencesFile());
    for (std::string line; std::getline(in, line);)
        if (line.rfind("renderer=", 0) == 0) parse(line.substr(9), chosen);
    return chosen;
}

void savePreferredBackend(Backend value)
{
    // Rewrites the preferences file keeping its other settings.
    std::vector<std::string> lines;
    {
        std::ifstream in(preferencesFile());
        for (std::string line; std::getline(in, line);)
            if (line.rfind("renderer=", 0) != 0) lines.push_back(line);
    }
    lines.push_back(std::string("renderer=") + (value == Backend::Vulkan ? "vulkan" : "opengl"));
    std::ofstream out(preferencesFile(), std::ios::trunc);
    for (const auto& line : lines) out << line << '\n';
}

Backend requestedBackend()
{
    Backend chosen = preferredBackend();
    if (const char* environment = std::getenv("ATOMFORGE_RENDERER"); environment && *environment) parse(environment, chosen);
    if (!gCommandLine.empty()) parse(gCommandLine, chosen);
    if (chosen == Backend::Vulkan && !vulkanBuilt()) chosen = Backend::OpenGL;
    return chosen;
}

const std::string& fallbackReason() { return gFallbackReason; }
void setFallbackReason(const std::string& reason) { gFallbackReason = reason; }

namespace
{
GLFWwindow* gMainWindow = nullptr;
}
GLFWwindow* mainWindow() { return gMainWindow; }
void setMainWindow(GLFWwindow* window) { gMainWindow = window; }
}
