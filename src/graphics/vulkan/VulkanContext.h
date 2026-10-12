#pragma once
// The Vulkan device, window swapchain and frame loop of the desktop: one
// command buffer per frame records any off-screen passes (previews, shadow
// map) and then the main pass (the scene and the interface), which is
// resolved from 4x MSAA into the swapchain image and presented.

#include "graphics/vulkan/VulkanLoader.h"

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

struct GLFWwindow;

namespace vkr
{
struct Buffer
{
    VkBuffer buffer = VK_NULL_HANDLE;
    VkDeviceMemory memory = VK_NULL_HANDLE;
    VkDeviceSize size = 0;
    void* mapped = nullptr;  // host-visible buffers stay mapped
};

struct Image
{
    VkImage image = VK_NULL_HANDLE;
    VkDeviceMemory memory = VK_NULL_HANDLE;
    VkImageView view = VK_NULL_HANDLE;
    VkFormat format = VK_FORMAT_UNDEFINED;
    uint32_t width = 0, height = 0;
};

class Context
{
public:
    static constexpr uint32_t kFramesInFlight = 2;

    bool init(GLFWwindow* window, std::string& error);
    void shutdown();
    bool ready() const { return m_device != VK_NULL_HANDLE; }

    // A frame: beginFrame acquires the next swapchain image and starts the
    // frame's command buffer (false: nothing is drawn this frame, e.g. the
    // window is minimised); off-screen passes are recorded into
    // commandBuffer(); beginMainPass / endMainPass bracket the window pass;
    // endFrame submits and presents.
    bool beginFrame();
    bool inFrame() const { return m_inFrame; }
    VkCommandBuffer commandBuffer() const { return m_frames[m_frameIndex].commandBuffer; }
    // A pausable pass can be interrupted once (pauseMainPass: its depth is then
    // readable by compute shaders as depthImage()) and continued (resumeMainPass).
    void beginMainPass(const float clear[4], bool pausable = false);
    void pauseMainPass();
    void resumeMainPass();
    const Image& depthImage() const { return m_depth; }
    // Changes whenever the swapchain's images (and so depthImage()) are made again.
    uint64_t swapchainGeneration() const { return m_swapchainGeneration; }
    void endMainPass();
    void endFrame();

    // The pixels of the next presented frame (RGBA, top row first), once endFrame has run.
    void requestCapture() { m_captureRequested = true; }
    bool takeCapture(std::vector<unsigned char>& rgba, int& width, int& height);
    // Keeps a copy of every presented frame (screenshots of the review and
    // capture modes, which read the last frame as OpenGL's front buffer).
    void setContinuousCapture(bool on) { m_continuousCapture = on; }
    // Vertical sync (on by default); the swapchain is remade on the next frame.
    void setVsync(bool on) { if (on != m_vsync) { m_vsync = on; m_resized = true; } }
    bool lastFrame(std::vector<unsigned char>& rgba, int& width, int& height) const;

    // Resources. Host-visible buffers are mapped persistently; destroy* waits
    // until no frame in flight can still use the resource.
    Buffer createBuffer(VkDeviceSize size, VkBufferUsageFlags usage, bool hostVisible);
    Buffer createDeviceBuffer(const void* data, VkDeviceSize size, VkBufferUsageFlags usage);
    void destroyBuffer(Buffer& buffer);
    // At once: for buffers no recorded command uses any more (after immediate()).
    void destroyBufferNow(Buffer& buffer);
    Image createImage(uint32_t width, uint32_t height, VkFormat format, VkImageUsageFlags usage,
                      VkSampleCountFlagBits samples, VkImageAspectFlags aspect);
    void destroyImage(Image& image);
    void destroyLater(std::function<void()> destroy);
    // Records and runs commands at once, waiting for them (uploads, readbacks).
    void immediate(const std::function<void(VkCommandBuffer)>& record);
    VkShaderModule shaderModule(const uint32_t* code, size_t bytes) const;
    uint32_t memoryType(uint32_t typeBits, VkMemoryPropertyFlags properties) const;

    VkInstance instance() const { return m_instance; }
    VkPhysicalDevice physicalDevice() const { return m_physical; }
    VkDevice device() const { return m_device; }
    VkQueue queue() const { return m_queue; }
    uint32_t queueFamily() const { return m_queueFamily; }
    VkDescriptorPool descriptorPool() const { return m_descriptorPool; }
    VkRenderPass mainRenderPass() const { return m_mainPass; }
    VkSampleCountFlagBits samples() const { return m_samples; }
    VkFormat colorFormat() const { return m_colorFormat; }
    VkFormat depthFormat() const { return m_depthFormat; }
    VkExtent2D extent() const { return m_extent; }
    uint32_t imageCount() const { return static_cast<uint32_t>(m_swapImages.size()); }
    uint32_t minImageCount() const { return m_minImageCount; }
    uint32_t frameIndex() const { return m_frameIndex; }
    bool wideLines() const { return m_wideLines; }
    float maxLineWidth() const { return m_maxLineWidth; }
    bool depthClamp() const { return m_depthClamp; }
    float maxPointSize() const { return m_maxPointSize; }
    bool preciseOcclusion() const { return m_preciseOcclusion; }
    float timestampPeriod() const { return m_timestampPeriod; }  // ns per tick; 0: no timestamps
    // The subpass dependencies of every pass sharing the scene pipelines
    // (the window's and the canvases'): earlier reads or copies of the
    // attachments before drawing, and sampling or copying them afterwards.
    static std::array<VkSubpassDependency, 2> passDependencies();
    const std::string& deviceName() const { return m_deviceName; }

private:
    struct Frame
    {
        VkCommandPool pool = VK_NULL_HANDLE;
        VkCommandBuffer commandBuffer = VK_NULL_HANDLE;
        VkFence fence = VK_NULL_HANDLE;
        VkSemaphore imageAvailable = VK_NULL_HANDLE;
        std::vector<std::function<void()>> deletions;
    };

    bool createSwapchain(std::string& error);
    void destroySwapchain();
    bool recreateSwapchain();
    void runDeletions(Frame& frame);

    GLFWwindow* m_window = nullptr;
    VkInstance m_instance = VK_NULL_HANDLE;
    VkSurfaceKHR m_surface = VK_NULL_HANDLE;
    VkPhysicalDevice m_physical = VK_NULL_HANDLE;
    VkDevice m_device = VK_NULL_HANDLE;
    VkQueue m_queue = VK_NULL_HANDLE;
    uint32_t m_queueFamily = 0;
    VkDescriptorPool m_descriptorPool = VK_NULL_HANDLE;
    VkCommandPool m_immediatePool = VK_NULL_HANDLE;
    VkPhysicalDeviceMemoryProperties m_memory{};
    std::string m_deviceName;
    bool m_wideLines = false;
    float m_maxLineWidth = 1.0f;
    bool m_depthClamp = false;
    float m_maxPointSize = 1.0f;
    bool m_preciseOcclusion = false;
    float m_timestampPeriod = 0.0f;

    VkSwapchainKHR m_swapchain = VK_NULL_HANDLE;
    VkFormat m_colorFormat = VK_FORMAT_B8G8R8A8_UNORM;
    VkFormat m_depthFormat = VK_FORMAT_D32_SFLOAT;
    VkExtent2D m_extent{0, 0};
    uint32_t m_minImageCount = 2;
    VkSampleCountFlagBits m_samples = VK_SAMPLE_COUNT_4_BIT;
    VkRenderPass m_mainPass = VK_NULL_HANDLE;
    std::vector<VkImage> m_swapImages;
    std::vector<VkImageView> m_swapViews;
    std::vector<VkFramebuffer> m_framebuffers;
    std::vector<VkSemaphore> m_renderFinished;  // one per swapchain image
    Image m_msaaColor, m_depth;
    VkRenderPass m_firstHalfPass = VK_NULL_HANDLE, m_secondHalfPass = VK_NULL_HANDLE;
    bool m_pausable = false, m_paused = false;
    uint64_t m_swapchainGeneration = 0;

    std::array<Frame, kFramesInFlight> m_frames{};
    uint32_t m_frameIndex = 0;
    uint32_t m_imageIndex = 0;
    bool m_inFrame = false;
    bool m_resized = false;

    bool m_captureRequested = false;
    bool m_captureRecorded = false;
    Buffer m_captureBuffer;
    std::vector<unsigned char> m_capturePixels;
    int m_captureWidth = 0, m_captureHeight = 0;
    bool m_captureReady = false;
    bool m_continuousCapture = false;
    bool m_vsync = true;
};

Context& context();
}
