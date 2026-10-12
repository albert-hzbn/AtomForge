#include "graphics/vulkan/VulkanLoader.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#define ATOMFORGE_VK_DEFINE(name) PFN_##name name = nullptr;
ATOMFORGE_VK_GLOBAL_FUNCTIONS(ATOMFORGE_VK_DEFINE)
ATOMFORGE_VK_INSTANCE_FUNCTIONS(ATOMFORGE_VK_DEFINE)
ATOMFORGE_VK_DEVICE_FUNCTIONS(ATOMFORGE_VK_DEFINE)
#undef ATOMFORGE_VK_DEFINE

namespace
{
PFN_vkGetInstanceProcAddr gGetInstanceProcAddr = nullptr;
}

namespace vkload
{
PFN_vkGetInstanceProcAddr getInstanceProcAddr() { return gGetInstanceProcAddr; }

bool loadGlobal()
{
    if (!glfwVulkanSupported()) return false;
    gGetInstanceProcAddr = reinterpret_cast<PFN_vkGetInstanceProcAddr>(glfwGetInstanceProcAddress(nullptr, "vkGetInstanceProcAddr"));
    if (!gGetInstanceProcAddr) return false;
#define ATOMFORGE_VK_LOAD(name) name = reinterpret_cast<PFN_##name>(gGetInstanceProcAddr(nullptr, #name));
    ATOMFORGE_VK_GLOBAL_FUNCTIONS(ATOMFORGE_VK_LOAD)
#undef ATOMFORGE_VK_LOAD
    return vkCreateInstance != nullptr;
}

void loadInstance(VkInstance instance)
{
#define ATOMFORGE_VK_LOAD(name) name = reinterpret_cast<PFN_##name>(gGetInstanceProcAddr(instance, #name));
    ATOMFORGE_VK_INSTANCE_FUNCTIONS(ATOMFORGE_VK_LOAD)
#undef ATOMFORGE_VK_LOAD
}

void loadDevice(VkDevice device)
{
#define ATOMFORGE_VK_LOAD(name) name = reinterpret_cast<PFN_##name>(vkGetDeviceProcAddr(device, #name));
    ATOMFORGE_VK_DEVICE_FUNCTIONS(ATOMFORGE_VK_LOAD)
#undef ATOMFORGE_VK_LOAD
}
}
