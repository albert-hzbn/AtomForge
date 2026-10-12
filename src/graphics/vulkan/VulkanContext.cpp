#include "graphics/vulkan/VulkanContext.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cstring>
#include <iostream>

namespace vkr
{
Context& context()
{
    static Context instance;
    return instance;
}

namespace
{
bool check(VkResult result, const char* what, std::string& error)
{
    if (result == VK_SUCCESS) return true;
    error = std::string(what) + " failed (VkResult " + std::to_string(static_cast<int>(result)) + ")";
    return false;
}

VkImageView makeView(VkDevice device, VkImage image, VkFormat format, VkImageAspectFlags aspect)
{
    VkImageViewCreateInfo info{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    info.image = image;
    info.viewType = VK_IMAGE_VIEW_TYPE_2D;
    info.format = format;
    info.subresourceRange = {aspect, 0, 1, 0, 1};
    VkImageView view = VK_NULL_HANDLE;
    vkCreateImageView(device, &info, nullptr, &view);
    return view;
}
}

bool Context::init(GLFWwindow* window, std::string& error)
{
    m_window = window;
    if (!vkload::loadGlobal())
    {
        error = "no Vulkan loader is installed";
        return false;
    }

    // Instance, with the extensions GLFW needs for a window surface.
    uint32_t glfwCount = 0;
    const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwCount);
    if (!glfwExtensions)
    {
        error = "GLFW found no Vulkan surface extensions";
        return false;
    }
    std::vector<const char*> extensions(glfwExtensions, glfwExtensions + glfwCount);
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "AtomForge";
    app.pEngineName = "AtomForge";
    app.apiVersion = VK_API_VERSION_1_1;
    VkInstanceCreateInfo instanceInfo{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    instanceInfo.pApplicationInfo = &app;
    instanceInfo.enabledExtensionCount = static_cast<uint32_t>(extensions.size());
    instanceInfo.ppEnabledExtensionNames = extensions.data();
    // ATOMFORGE_VULKAN_VALIDATION=1: the Khronos validation layer, when installed.
    const char* validation = "VK_LAYER_KHRONOS_validation";
    if (const char* wanted = std::getenv("ATOMFORGE_VULKAN_VALIDATION"); wanted && *wanted == '1')
    {
        uint32_t layerCount = 0;
        vkEnumerateInstanceLayerProperties(&layerCount, nullptr);
        std::vector<VkLayerProperties> layers(layerCount);
        vkEnumerateInstanceLayerProperties(&layerCount, layers.data());
        for (const auto& layer : layers)
            if (std::strcmp(layer.layerName, validation) == 0)
            {
                instanceInfo.enabledLayerCount = 1;
                instanceInfo.ppEnabledLayerNames = &validation;
            }
        std::cout << (instanceInfo.enabledLayerCount ? "[Vulkan] Validation layer enabled." : "[Vulkan] Validation layer not found.")
                  << std::endl;
    }
    if (!check(vkCreateInstance(&instanceInfo, nullptr, &m_instance), "vkCreateInstance", error)) return false;
    vkload::loadInstance(m_instance);
    if (!check(glfwCreateWindowSurface(m_instance, window, nullptr, &m_surface), "glfwCreateWindowSurface", error)) return false;

    // A device that can draw to the window: a discrete GPU first.
    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(m_instance, &deviceCount, nullptr);
    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(m_instance, &deviceCount, devices.data());
    int bestScore = -1;
    for (VkPhysicalDevice candidate : devices)
    {
        VkPhysicalDeviceProperties properties;
        vkGetPhysicalDeviceProperties(candidate, &properties);
        if (properties.apiVersion < VK_API_VERSION_1_1) continue;
        uint32_t extensionCount = 0;
        vkEnumerateDeviceExtensionProperties(candidate, nullptr, &extensionCount, nullptr);
        std::vector<VkExtensionProperties> available(extensionCount);
        vkEnumerateDeviceExtensionProperties(candidate, nullptr, &extensionCount, available.data());
        const bool swapchain = std::any_of(available.begin(), available.end(), [](const VkExtensionProperties& e) {
            return std::strcmp(e.extensionName, VK_KHR_SWAPCHAIN_EXTENSION_NAME) == 0;
        });
        if (!swapchain) continue;
        uint32_t familyCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(candidate, &familyCount, nullptr);
        std::vector<VkQueueFamilyProperties> families(familyCount);
        vkGetPhysicalDeviceQueueFamilyProperties(candidate, &familyCount, families.data());
        for (uint32_t family = 0; family < familyCount; ++family)
        {
            VkBool32 present = VK_FALSE;
            vkGetPhysicalDeviceSurfaceSupportKHR(candidate, family, m_surface, &present);
            if (!(families[family].queueFlags & VK_QUEUE_GRAPHICS_BIT) || !present) continue;
            const int score = properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU ? 3
                            : properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU ? 2 : 1;
            if (score > bestScore)
            {
                bestScore = score;
                m_physical = candidate;
                m_queueFamily = family;
                m_deviceName = properties.deviceName;
            }
            break;
        }
    }
    if (!m_physical)
    {
        error = "no Vulkan device can draw to the window";
        return false;
    }
    vkGetPhysicalDeviceMemoryProperties(m_physical, &m_memory);
    VkPhysicalDeviceFeatures supported{};
    vkGetPhysicalDeviceFeatures(m_physical, &supported);
    VkPhysicalDeviceProperties properties;
    vkGetPhysicalDeviceProperties(m_physical, &properties);
    VkPhysicalDeviceFeatures enabled{};
    enabled.wideLines = supported.wideLines;
    enabled.depthClamp = supported.depthClamp;
    // Atom clouds: point sprites larger than a pixel, occlusion queries that count samples.
    enabled.largePoints = supported.largePoints;
    enabled.occlusionQueryPrecise = supported.occlusionQueryPrecise;
    m_maxPointSize = supported.largePoints ? properties.limits.pointSizeRange[1] : 1.0f;
    m_preciseOcclusion = supported.occlusionQueryPrecise == VK_TRUE;
    m_timestampPeriod = properties.limits.timestampComputeAndGraphics ? properties.limits.timestampPeriod : 0.0f;
    m_wideLines = supported.wideLines == VK_TRUE;
    m_maxLineWidth = m_wideLines ? properties.limits.lineWidthRange[1] : 1.0f;
    m_depthClamp = supported.depthClamp == VK_TRUE;
    // 4x MSAA as the OpenGL window, when the device has it.
    const VkSampleCountFlags counts = properties.limits.framebufferColorSampleCounts & properties.limits.framebufferDepthSampleCounts;
    m_samples = (counts & VK_SAMPLE_COUNT_4_BIT) ? VK_SAMPLE_COUNT_4_BIT : VK_SAMPLE_COUNT_1_BIT;

    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queueInfo{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queueInfo.queueFamilyIndex = m_queueFamily;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;
    const char* deviceExtensions[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    VkDeviceCreateInfo deviceInfo{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
    deviceInfo.enabledExtensionCount = 1;
    deviceInfo.ppEnabledExtensionNames = deviceExtensions;
    deviceInfo.pEnabledFeatures = &enabled;
    if (!check(vkCreateDevice(m_physical, &deviceInfo, nullptr, &m_device), "vkCreateDevice", error)) return false;
    vkload::loadDevice(m_device);
    vkGetDeviceQueue(m_device, m_queueFamily, 0, &m_queue);

    // Descriptors for the interface's textures and the renderers' sets.
    const VkDescriptorPoolSize sizes[] = {
        {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 512},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 64},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 64},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 64},
        // The ImGui backend binds textures as a sampled image plus a sampler.
        {VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 512},
        {VK_DESCRIPTOR_TYPE_SAMPLER, 16},
    };
    VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    poolInfo.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    poolInfo.maxSets = 640;
    poolInfo.poolSizeCount = static_cast<uint32_t>(std::size(sizes));
    poolInfo.pPoolSizes = sizes;
    if (!check(vkCreateDescriptorPool(m_device, &poolInfo, nullptr, &m_descriptorPool), "vkCreateDescriptorPool", error)) return false;

    VkCommandPoolCreateInfo commandPoolInfo{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    commandPoolInfo.flags = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT | VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    commandPoolInfo.queueFamilyIndex = m_queueFamily;
    if (!check(vkCreateCommandPool(m_device, &commandPoolInfo, nullptr, &m_immediatePool), "vkCreateCommandPool", error)) return false;
    for (Frame& frame : m_frames)
    {
        if (!check(vkCreateCommandPool(m_device, &commandPoolInfo, nullptr, &frame.pool), "vkCreateCommandPool", error)) return false;
        VkCommandBufferAllocateInfo allocate{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        allocate.commandPool = frame.pool;
        allocate.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        allocate.commandBufferCount = 1;
        if (!check(vkAllocateCommandBuffers(m_device, &allocate, &frame.commandBuffer), "vkAllocateCommandBuffers", error)) return false;
        VkFenceCreateInfo fenceInfo{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
        fenceInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT;
        VkSemaphoreCreateInfo semaphoreInfo{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        if (!check(vkCreateFence(m_device, &fenceInfo, nullptr, &frame.fence), "vkCreateFence", error) ||
            !check(vkCreateSemaphore(m_device, &semaphoreInfo, nullptr, &frame.imageAvailable), "vkCreateSemaphore", error))
            return false;
    }

    // The window's colour format: plain 8-bit (as the OpenGL default framebuffer).
    uint32_t formatCount = 0;
    vkGetPhysicalDeviceSurfaceFormatsKHR(m_physical, m_surface, &formatCount, nullptr);
    std::vector<VkSurfaceFormatKHR> formats(formatCount);
    vkGetPhysicalDeviceSurfaceFormatsKHR(m_physical, m_surface, &formatCount, formats.data());
    if (formats.empty())
    {
        error = "the window offers no surface formats";
        return false;
    }
    m_colorFormat = formats[0].format == VK_FORMAT_UNDEFINED ? VK_FORMAT_B8G8R8A8_UNORM : formats[0].format;
    for (const auto& format : formats)
        if ((format.format == VK_FORMAT_B8G8R8A8_UNORM || format.format == VK_FORMAT_R8G8B8A8_UNORM) &&
            format.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR)
        {
            m_colorFormat = format.format;
            break;
        }

    // The main pass: MSAA colour and depth, resolved into the swapchain image.
    std::vector<VkAttachmentDescription> attachments;
    VkAttachmentDescription color{};
    color.format = m_colorFormat;
    color.samples = m_samples;
    color.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    color.storeOp = m_samples == VK_SAMPLE_COUNT_1_BIT ? VK_ATTACHMENT_STORE_OP_STORE : VK_ATTACHMENT_STORE_OP_DONT_CARE;
    color.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    color.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    color.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    color.finalLayout = m_samples == VK_SAMPLE_COUNT_1_BIT ? VK_IMAGE_LAYOUT_PRESENT_SRC_KHR : VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
    attachments.push_back(color);
    VkAttachmentDescription depth = color;
    depth.format = m_depthFormat;
    depth.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    depth.finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
    attachments.push_back(depth);
    if (m_samples != VK_SAMPLE_COUNT_1_BIT)
    {
        VkAttachmentDescription resolve = color;
        resolve.samples = VK_SAMPLE_COUNT_1_BIT;
        resolve.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        resolve.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        resolve.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
        attachments.push_back(resolve);
    }
    VkAttachmentReference colorRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};
    VkAttachmentReference resolveRef{2, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkSubpassDescription subpass{};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &colorRef;
    subpass.pDepthStencilAttachment = &depthRef;
    if (m_samples != VK_SAMPLE_COUNT_1_BIT) subpass.pResolveAttachments = &resolveRef;
    // The same dependencies as the canvases' pass (render passes are only
    // compatible, so share pipelines, when their dependencies are identical).
    const auto dependencies = passDependencies();
    VkRenderPassCreateInfo passInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    passInfo.attachmentCount = static_cast<uint32_t>(attachments.size());
    passInfo.pAttachments = attachments.data();
    passInfo.subpassCount = 1;
    passInfo.pSubpasses = &subpass;
    passInfo.dependencyCount = static_cast<uint32_t>(dependencies.size());
    passInfo.pDependencies = dependencies.data();
    if (!check(vkCreateRenderPass(m_device, &passInfo, nullptr, &m_mainPass), "vkCreateRenderPass", error)) return false;

    return createSwapchain(error);
}

bool Context::createSwapchain(std::string& error)
{
    VkSurfaceCapabilitiesKHR capabilities;
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(m_physical, m_surface, &capabilities);
    int width = 0, height = 0;
    glfwGetFramebufferSize(m_window, &width, &height);
    if (capabilities.currentExtent.width != UINT32_MAX)
        m_extent = capabilities.currentExtent;
    else
        m_extent = {static_cast<uint32_t>(std::max(width, 1)), static_cast<uint32_t>(std::max(height, 1))};
    m_extent.width = std::clamp(m_extent.width, capabilities.minImageExtent.width, std::max(capabilities.minImageExtent.width, capabilities.maxImageExtent.width));
    m_extent.height = std::clamp(m_extent.height, capabilities.minImageExtent.height, std::max(capabilities.minImageExtent.height, capabilities.maxImageExtent.height));
    if (m_extent.width == 0 || m_extent.height == 0) return true;  // minimised: made when it is shown again

    m_minImageCount = std::max(2u, capabilities.minImageCount);
    // Unthrottled presentation needs an image more, or acquiring waits for the display.
    uint32_t imageCount = capabilities.minImageCount + (m_vsync ? 1 : 2);
    if (capabilities.maxImageCount > 0) imageCount = std::min(imageCount, capabilities.maxImageCount);
    VkSwapchainCreateInfoKHR info{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    info.surface = m_surface;
    info.minImageCount = imageCount;
    info.imageFormat = m_colorFormat;
    info.imageColorSpace = VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
    info.imageExtent = m_extent;
    info.imageArrayLayers = 1;
    info.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                      (capabilities.supportedUsageFlags & VK_IMAGE_USAGE_TRANSFER_SRC_BIT);  // for screenshots
    info.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    info.preTransform = capabilities.currentTransform;
    info.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    if (!(capabilities.supportedCompositeAlpha & VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR))
        for (VkCompositeAlphaFlagBitsKHR mode : {VK_COMPOSITE_ALPHA_INHERIT_BIT_KHR, VK_COMPOSITE_ALPHA_PRE_MULTIPLIED_BIT_KHR, VK_COMPOSITE_ALPHA_POST_MULTIPLIED_BIT_KHR})
            if (capabilities.supportedCompositeAlpha & mode) { info.compositeAlpha = mode; break; }
    info.presentMode = VK_PRESENT_MODE_FIFO_KHR;  // vertical sync, always available
    if (!m_vsync)
    {
        // Unthrottled (benchmarks): immediate, else mailbox, when offered.
        uint32_t modeCount = 0;
        vkGetPhysicalDeviceSurfacePresentModesKHR(m_physical, m_surface, &modeCount, nullptr);
        std::vector<VkPresentModeKHR> modes(modeCount);
        vkGetPhysicalDeviceSurfacePresentModesKHR(m_physical, m_surface, &modeCount, modes.data());
        for (VkPresentModeKHR mode : {VK_PRESENT_MODE_IMMEDIATE_KHR, VK_PRESENT_MODE_MAILBOX_KHR})
            if (std::find(modes.begin(), modes.end(), mode) != modes.end()) { info.presentMode = mode; break; }
    }
    info.clipped = VK_TRUE;
    info.oldSwapchain = VK_NULL_HANDLE;
    if (!check(vkCreateSwapchainKHR(m_device, &info, nullptr, &m_swapchain), "vkCreateSwapchainKHR", error)) return false;

    uint32_t count = 0;
    vkGetSwapchainImagesKHR(m_device, m_swapchain, &count, nullptr);
    m_swapImages.resize(count);
    vkGetSwapchainImagesKHR(m_device, m_swapchain, &count, m_swapImages.data());
    m_swapViews.clear();
    for (VkImage image : m_swapImages) m_swapViews.push_back(makeView(m_device, image, m_colorFormat, VK_IMAGE_ASPECT_COLOR_BIT));

    if (m_samples != VK_SAMPLE_COUNT_1_BIT)
        m_msaaColor = createImage(m_extent.width, m_extent.height, m_colorFormat, VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT,
                                  m_samples, VK_IMAGE_ASPECT_COLOR_BIT);
    m_depth = createImage(m_extent.width, m_extent.height, m_depthFormat, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT,
                          m_samples, VK_IMAGE_ASPECT_DEPTH_BIT);
    m_framebuffers.clear();
    m_renderFinished.clear();
    for (std::size_t i = 0; i < m_swapImages.size(); ++i)
    {
        std::vector<VkImageView> views;
        if (m_samples == VK_SAMPLE_COUNT_1_BIT) views = {m_swapViews[i], m_depth.view};
        else views = {m_msaaColor.view, m_depth.view, m_swapViews[i]};
        VkFramebufferCreateInfo framebufferInfo{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        framebufferInfo.renderPass = m_mainPass;
        framebufferInfo.attachmentCount = static_cast<uint32_t>(views.size());
        framebufferInfo.pAttachments = views.data();
        framebufferInfo.width = m_extent.width;
        framebufferInfo.height = m_extent.height;
        framebufferInfo.layers = 1;
        VkFramebuffer framebuffer = VK_NULL_HANDLE;
        if (!check(vkCreateFramebuffer(m_device, &framebufferInfo, nullptr, &framebuffer), "vkCreateFramebuffer", error)) return false;
        m_framebuffers.push_back(framebuffer);
        VkSemaphoreCreateInfo semaphoreInfo{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        VkSemaphore semaphore = VK_NULL_HANDLE;
        vkCreateSemaphore(m_device, &semaphoreInfo, nullptr, &semaphore);
        m_renderFinished.push_back(semaphore);
    }
    return true;
}

void Context::destroySwapchain()
{
    for (VkFramebuffer framebuffer : m_framebuffers) vkDestroyFramebuffer(m_device, framebuffer, nullptr);
    for (VkSemaphore semaphore : m_renderFinished) vkDestroySemaphore(m_device, semaphore, nullptr);
    for (VkImageView view : m_swapViews) vkDestroyImageView(m_device, view, nullptr);
    m_framebuffers.clear();
    m_renderFinished.clear();
    m_swapViews.clear();
    m_swapImages.clear();
    // Not deferred: the device is idle here.
    for (Image* image : {&m_msaaColor, &m_depth})
        if (image->image)
        {
            vkDestroyImageView(m_device, image->view, nullptr);
            vkDestroyImage(m_device, image->image, nullptr);
            vkFreeMemory(m_device, image->memory, nullptr);
            *image = Image();
        }
    if (m_swapchain) vkDestroySwapchainKHR(m_device, m_swapchain, nullptr);
    m_swapchain = VK_NULL_HANDLE;
}

bool Context::recreateSwapchain()
{
    vkDeviceWaitIdle(m_device);
    destroySwapchain();
    std::string error;
    if (!createSwapchain(error))
    {
        std::cerr << "[Vulkan] " << error << '\n';
        return false;
    }
    m_resized = false;
    return m_swapchain != VK_NULL_HANDLE;
}

void Context::runDeletions(Frame& frame)
{
    for (auto& destroy : frame.deletions) destroy();
    frame.deletions.clear();
}

bool Context::beginFrame()
{
    if (!ready() || m_inFrame) return false;
    int width = 0, height = 0;
    glfwGetFramebufferSize(m_window, &width, &height);
    if (width == 0 || height == 0) return false;
    if (!m_swapchain || m_resized || static_cast<uint32_t>(width) != m_extent.width || static_cast<uint32_t>(height) != m_extent.height)
        if (!recreateSwapchain()) return false;

    Frame& frame = m_frames[m_frameIndex];
    vkWaitForFences(m_device, 1, &frame.fence, VK_TRUE, UINT64_MAX);
    runDeletions(frame);
    const VkResult acquired = vkAcquireNextImageKHR(m_device, m_swapchain, UINT64_MAX, frame.imageAvailable, VK_NULL_HANDLE, &m_imageIndex);
    if (acquired == VK_ERROR_OUT_OF_DATE_KHR)
    {
        m_resized = true;
        return false;
    }
    if (acquired != VK_SUCCESS && acquired != VK_SUBOPTIMAL_KHR) return false;
    vkResetFences(m_device, 1, &frame.fence);
    vkResetCommandPool(m_device, frame.pool, 0);
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(frame.commandBuffer, &begin);
    m_inFrame = true;
    m_captureRecorded = false;
    return true;
}

void Context::beginMainPass(const float clear[4])
{
    VkClearValue values[3];
    values[0].color = {{clear[0], clear[1], clear[2], clear[3]}};
    values[1].depthStencil = {1.0f, 0};
    values[2].color = values[0].color;
    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = m_mainPass;
    begin.framebuffer = m_framebuffers[m_imageIndex];
    begin.renderArea = {{0, 0}, m_extent};
    begin.clearValueCount = m_samples == VK_SAMPLE_COUNT_1_BIT ? 2u : 3u;
    begin.pClearValues = values;
    vkCmdBeginRenderPass(commandBuffer(), &begin, VK_SUBPASS_CONTENTS_INLINE);
}

void Context::endMainPass()
{
    VkCommandBuffer cmd = commandBuffer();
    vkCmdEndRenderPass(cmd);
    if (!m_captureRequested && !m_continuousCapture) return;
    // A copy of the presented image for a screenshot.
    const VkDeviceSize bytes = static_cast<VkDeviceSize>(m_extent.width) * m_extent.height * 4;
    if (m_captureBuffer.size < bytes)
    {
        destroyBuffer(m_captureBuffer);
        m_captureBuffer = createBuffer(bytes, VK_BUFFER_USAGE_TRANSFER_DST_BIT, true);
    }
    VkImageMemoryBarrier toTransfer{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    toTransfer.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
    toTransfer.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    toTransfer.oldLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    toTransfer.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    toTransfer.srcQueueFamilyIndex = toTransfer.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    toTransfer.image = m_swapImages[m_imageIndex];
    toTransfer.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 0, nullptr, 0, nullptr, 1, &toTransfer);
    VkBufferImageCopy region{};
    region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
    region.imageExtent = {m_extent.width, m_extent.height, 1};
    vkCmdCopyImageToBuffer(cmd, m_swapImages[m_imageIndex], VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, m_captureBuffer.buffer, 1, &region);
    VkImageMemoryBarrier toPresent = toTransfer;
    toPresent.srcAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    toPresent.dstAccessMask = 0;
    toPresent.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    toPresent.newLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, 0, 0, nullptr, 0, nullptr, 1, &toPresent);
    m_captureRecorded = true;
    m_captureRequested = false;
}

void Context::endFrame()
{
    if (!m_inFrame) return;
    Frame& frame = m_frames[m_frameIndex];
    vkEndCommandBuffer(frame.commandBuffer);
    const VkPipelineStageFlags wait = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.waitSemaphoreCount = 1;
    submit.pWaitSemaphores = &frame.imageAvailable;
    submit.pWaitDstStageMask = &wait;
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &frame.commandBuffer;
    submit.signalSemaphoreCount = 1;
    submit.pSignalSemaphores = &m_renderFinished[m_imageIndex];
    vkQueueSubmit(m_queue, 1, &submit, frame.fence);

    if (m_captureRecorded)
    {
        // The screenshot: wait for this frame, then read the copy.
        vkWaitForFences(m_device, 1, &frame.fence, VK_TRUE, UINT64_MAX);
        m_captureWidth = static_cast<int>(m_extent.width);
        m_captureHeight = static_cast<int>(m_extent.height);
        m_capturePixels.resize(static_cast<std::size_t>(m_captureWidth) * m_captureHeight * 4);
        std::memcpy(m_capturePixels.data(), m_captureBuffer.mapped, m_capturePixels.size());
        if (m_colorFormat == VK_FORMAT_B8G8R8A8_UNORM || m_colorFormat == VK_FORMAT_B8G8R8A8_SRGB)
            for (std::size_t i = 0; i < m_capturePixels.size(); i += 4) std::swap(m_capturePixels[i], m_capturePixels[i + 2]);
        for (std::size_t i = 3; i < m_capturePixels.size(); i += 4) m_capturePixels[i] = 255;
        m_captureReady = true;
        m_captureRecorded = false;
    }

    VkPresentInfoKHR present{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    present.waitSemaphoreCount = 1;
    present.pWaitSemaphores = &m_renderFinished[m_imageIndex];
    present.swapchainCount = 1;
    present.pSwapchains = &m_swapchain;
    present.pImageIndices = &m_imageIndex;
    const VkResult presented = vkQueuePresentKHR(m_queue, &present);
    if (presented == VK_ERROR_OUT_OF_DATE_KHR || presented == VK_SUBOPTIMAL_KHR) m_resized = true;
    m_inFrame = false;
    m_frameIndex = (m_frameIndex + 1) % kFramesInFlight;
}

bool Context::lastFrame(std::vector<unsigned char>& rgba, int& width, int& height) const
{
    if (!m_continuousCapture || m_capturePixels.empty()) return false;
    rgba = m_capturePixels;
    width = m_captureWidth;
    height = m_captureHeight;
    return true;
}

bool Context::takeCapture(std::vector<unsigned char>& rgba, int& width, int& height)
{
    if (!m_captureReady || m_continuousCapture) return false;
    rgba = std::move(m_capturePixels);
    width = m_captureWidth;
    height = m_captureHeight;
    m_captureReady = false;
    return true;
}

void Context::shutdown()
{
    if (!m_device)
    {
        if (m_surface) vkDestroySurfaceKHR(m_instance, m_surface, nullptr);
        if (m_instance) vkDestroyInstance(m_instance, nullptr);
        m_surface = VK_NULL_HANDLE;
        m_instance = VK_NULL_HANDLE;
        return;
    }
    vkDeviceWaitIdle(m_device);
    destroyBuffer(m_captureBuffer);
    for (Frame& frame : m_frames) runDeletions(frame);
    destroySwapchain();
    for (Frame& frame : m_frames)
    {
        vkDestroyFence(m_device, frame.fence, nullptr);
        vkDestroySemaphore(m_device, frame.imageAvailable, nullptr);
        vkDestroyCommandPool(m_device, frame.pool, nullptr);
        frame = Frame();
    }
    for (Frame& frame : m_frames) runDeletions(frame);
    vkDestroyCommandPool(m_device, m_immediatePool, nullptr);
    vkDestroyRenderPass(m_device, m_mainPass, nullptr);
    vkDestroyDescriptorPool(m_device, m_descriptorPool, nullptr);
    vkDestroyDevice(m_device, nullptr);
    vkDestroySurfaceKHR(m_instance, m_surface, nullptr);
    vkDestroyInstance(m_instance, nullptr);
    m_device = VK_NULL_HANDLE;
    m_surface = VK_NULL_HANDLE;
    m_instance = VK_NULL_HANDLE;
}

uint32_t Context::memoryType(uint32_t typeBits, VkMemoryPropertyFlags properties) const
{
    for (uint32_t i = 0; i < m_memory.memoryTypeCount; ++i)
        if ((typeBits & (1u << i)) && (m_memory.memoryTypes[i].propertyFlags & properties) == properties) return i;
    return UINT32_MAX;
}

Buffer Context::createBuffer(VkDeviceSize size, VkBufferUsageFlags usage, bool hostVisible)
{
    Buffer result;
    result.size = std::max<VkDeviceSize>(size, 4);
    VkBufferCreateInfo info{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    info.size = result.size;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    if (vkCreateBuffer(m_device, &info, nullptr, &result.buffer) != VK_SUCCESS) return Buffer();
    VkMemoryRequirements requirements;
    vkGetBufferMemoryRequirements(m_device, result.buffer, &requirements);
    VkMemoryAllocateInfo allocate{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocate.allocationSize = requirements.size;
    allocate.memoryTypeIndex = memoryType(requirements.memoryTypeBits,
        hostVisible ? (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT) : VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (allocate.memoryTypeIndex == UINT32_MAX || vkAllocateMemory(m_device, &allocate, nullptr, &result.memory) != VK_SUCCESS)
    {
        vkDestroyBuffer(m_device, result.buffer, nullptr);
        return Buffer();
    }
    vkBindBufferMemory(m_device, result.buffer, result.memory, 0);
    if (hostVisible) vkMapMemory(m_device, result.memory, 0, VK_WHOLE_SIZE, 0, &result.mapped);
    return result;
}

Buffer Context::createDeviceBuffer(const void* data, VkDeviceSize size, VkBufferUsageFlags usage)
{
    Buffer result = createBuffer(size, usage | VK_BUFFER_USAGE_TRANSFER_DST_BIT, false);
    if (!result.buffer || !data || size == 0) return result;
    Buffer staging = createBuffer(size, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, true);
    std::memcpy(staging.mapped, data, static_cast<std::size_t>(size));
    immediate([&](VkCommandBuffer cmd) {
        VkBufferCopy copy{0, 0, size};
        vkCmdCopyBuffer(cmd, staging.buffer, result.buffer, 1, &copy);
    });
    vkUnmapMemory(m_device, staging.memory);
    vkDestroyBuffer(m_device, staging.buffer, nullptr);
    vkFreeMemory(m_device, staging.memory, nullptr);
    return result;
}

void Context::destroyLater(std::function<void()> destroy)
{
    if (!m_device) return;
    m_frames[m_frameIndex].deletions.push_back(std::move(destroy));
}

void Context::destroyBuffer(Buffer& buffer)
{
    if (!buffer.buffer) return;
    const Buffer doomed = buffer;
    VkDevice device = m_device;
    destroyLater([device, doomed] {
        if (doomed.mapped) vkUnmapMemory(device, doomed.memory);
        vkDestroyBuffer(device, doomed.buffer, nullptr);
        vkFreeMemory(device, doomed.memory, nullptr);
    });
    buffer = Buffer();
}

std::array<VkSubpassDependency, 2> Context::passDependencies()
{
    std::array<VkSubpassDependency, 2> dependencies{};
    dependencies[0].srcSubpass = VK_SUBPASS_EXTERNAL;
    dependencies[0].dstSubpass = 0;
    dependencies[0].srcStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT |
                                   VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT;
    dependencies[0].dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    dependencies[0].srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
    dependencies[0].dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_READ_BIT |
                                    VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
    dependencies[1].srcSubpass = 0;
    dependencies[1].dstSubpass = VK_SUBPASS_EXTERNAL;
    dependencies[1].srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    dependencies[1].dstStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT;
    dependencies[1].srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
    dependencies[1].dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_TRANSFER_READ_BIT;
    return dependencies;
}

void Context::destroyBufferNow(Buffer& buffer)
{
    if (!buffer.buffer || !m_device) return;
    if (buffer.mapped) vkUnmapMemory(m_device, buffer.memory);
    vkDestroyBuffer(m_device, buffer.buffer, nullptr);
    vkFreeMemory(m_device, buffer.memory, nullptr);
    buffer = Buffer();
}

Image Context::createImage(uint32_t width, uint32_t height, VkFormat format, VkImageUsageFlags usage,
                           VkSampleCountFlagBits samples, VkImageAspectFlags aspect)
{
    Image result;
    result.format = format;
    result.width = width;
    result.height = height;
    VkImageCreateInfo info{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    info.imageType = VK_IMAGE_TYPE_2D;
    info.format = format;
    info.extent = {width, height, 1};
    info.mipLevels = 1;
    info.arrayLayers = 1;
    info.samples = samples;
    info.tiling = VK_IMAGE_TILING_OPTIMAL;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    if (vkCreateImage(m_device, &info, nullptr, &result.image) != VK_SUCCESS) return Image();
    VkMemoryRequirements requirements;
    vkGetImageMemoryRequirements(m_device, result.image, &requirements);
    VkMemoryAllocateInfo allocate{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocate.allocationSize = requirements.size;
    allocate.memoryTypeIndex = memoryType(requirements.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (allocate.memoryTypeIndex == UINT32_MAX || vkAllocateMemory(m_device, &allocate, nullptr, &result.memory) != VK_SUCCESS)
    {
        vkDestroyImage(m_device, result.image, nullptr);
        return Image();
    }
    vkBindImageMemory(m_device, result.image, result.memory, 0);
    result.view = makeView(m_device, result.image, format, aspect);
    return result;
}

void Context::destroyImage(Image& image)
{
    if (!image.image) return;
    const Image doomed = image;
    VkDevice device = m_device;
    destroyLater([device, doomed] {
        vkDestroyImageView(device, doomed.view, nullptr);
        vkDestroyImage(device, doomed.image, nullptr);
        vkFreeMemory(device, doomed.memory, nullptr);
    });
    image = Image();
}

void Context::immediate(const std::function<void(VkCommandBuffer)>& record)
{
    VkCommandBufferAllocateInfo allocate{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    allocate.commandPool = m_immediatePool;
    allocate.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocate.commandBufferCount = 1;
    VkCommandBuffer cmd = VK_NULL_HANDLE;
    vkAllocateCommandBuffers(m_device, &allocate, &cmd);
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(cmd, &begin);
    record(cmd);
    vkEndCommandBuffer(cmd);
    VkFenceCreateInfo fenceInfo{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    VkFence fence = VK_NULL_HANDLE;
    vkCreateFence(m_device, &fenceInfo, nullptr, &fence);
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &cmd;
    vkQueueSubmit(m_queue, 1, &submit, fence);
    vkWaitForFences(m_device, 1, &fence, VK_TRUE, UINT64_MAX);
    vkDestroyFence(m_device, fence, nullptr);
    vkFreeCommandBuffers(m_device, m_immediatePool, 1, &cmd);
}

VkShaderModule Context::shaderModule(const uint32_t* code, size_t bytes) const
{
    VkShaderModuleCreateInfo info{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    info.codeSize = bytes;
    info.pCode = code;
    VkShaderModule module = VK_NULL_HANDLE;
    vkCreateShaderModule(m_device, &info, nullptr, &module);
    return module;
}
}
