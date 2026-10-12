#include "graphics/FrameCapture.h"

#include "graphics/GraphicsBackend.h"
#ifdef ATOMFORGE_VULKAN
#include "graphics/vulkan/VulkanContext.h"
#endif

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cstddef>

namespace graphics
{
void enableFrameCapture()
{
#ifdef ATOMFORGE_VULKAN
    if (usingVulkan()) vkr::context().setContinuousCapture(true);
#endif
}

void setVsync(bool on)
{
#ifdef ATOMFORGE_VULKAN
    if (usingVulkan())
    {
        vkr::context().setVsync(on);
        return;
    }
#endif
    glfwSwapInterval(on ? 1 : 0);
}

void waitIdle()
{
#ifdef ATOMFORGE_VULKAN
    if (usingVulkan())
    {
        if (vkr::context().ready()) vkDeviceWaitIdle(vkr::context().device());
        return;
    }
#endif
    glFinish();
}

bool lastFrame(std::vector<unsigned char>& rgb, int& width, int& height)
{
#ifdef ATOMFORGE_VULKAN
    if (usingVulkan())
    {
        std::vector<unsigned char> rgba;
        if (!vkr::context().lastFrame(rgba, width, height) || width <= 0 || height <= 0) return false;
        rgb.resize(static_cast<std::size_t>(width) * height * 3);
        for (std::size_t i = 0, n = static_cast<std::size_t>(width) * height; i < n; ++i)
            std::copy_n(rgba.begin() + static_cast<std::ptrdiff_t>(i * 4), 3, rgb.begin() + static_cast<std::ptrdiff_t>(i * 3));
        return true;
    }
#endif
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    width = viewport[2];
    height = viewport[3];
    if (width <= 0 || height <= 0) return false;
    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * height * 3);
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    glReadBuffer(GL_BACK);
    rgb.resize(pixels.size());
    const std::ptrdiff_t row = static_cast<std::ptrdiff_t>(width) * 3;
    for (int y = 0; y < height; ++y)
        std::copy_n(pixels.begin() + y * row, row, rgb.begin() + static_cast<std::ptrdiff_t>(height - 1 - y) * row);
    return true;
}
}
