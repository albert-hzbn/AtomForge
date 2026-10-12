#pragma once
// An off-screen image the Vulkan scene renderer draws into: the dialogs'
// previews (shown with ImGui::Image) and image export (copied out). Drawing
// is recorded into the frame's command buffer, or into one given to begin()
// (image export runs it at once). Two images alternate so a frame still on
// the GPU never sees its preview overwritten.

#include "graphics/vulkan/VulkanSceneRenderer.h"

#include "imgui.h"

#include <glm/glm.hpp>

#include <vector>

namespace vkr
{
class Canvas
{
public:
    Canvas() = default;
    Canvas(const Canvas&) = delete;
    Canvas& operator=(const Canvas&) = delete;
    ~Canvas() { release(); }

    // Starts a width x height image cleared to `background`; false when
    // nothing can be drawn (no frame is being recorded).
    bool begin(int width, int height, const glm::vec3& background, VkCommandBuffer cmd = VK_NULL_HANDLE);
    // The camera and lights of the draws that follow.
    void setView(const SceneView& view, const SceneLighting& lighting);
    void drawScene(const SceneBuffers& buffers, const SceneFlags& flags, const glm::vec3& boxColor);
    void drawGhostAtoms(const SceneBuffers& buffers, float opacity);
    void drawLines(const std::vector<glm::vec3>& points, const glm::vec3& color, LineTopology topology, bool depthTest);
    void drawMesh(const std::vector<MeshVertex>& vertices, const glm::vec4& color, MeshStyle style, bool translucent,
                  const glm::vec3& lightDirection = glm::vec3(0.0f, 0.0f, 1.0f));
    void end();

    // After end(): records a copy of the image into `buffer` (host-visible,
    // width * height * 4 bytes, rows top first, in colorFormat() order).
    void copyTo(VkCommandBuffer cmd, const Buffer& buffer);
    VkFormat colorFormat() const;

    // The last image drawn (null before any), rows top first.
    ImTextureID texture();
    // ImGui::Image of the last image at `size` (nothing before any).
    void image(const ImVec2& size);
    bool drawn() const { return m_latest >= 0; }
    VkCommandBuffer commandBuffer() const { return m_cmd; }
    uint32_t frameOffset() const { return m_frameOffset; }
    VkExtent2D extent() const { return m_extent; }

    void release();

private:
    struct Target
    {
        Image color;    // multisampled (unused with one sample)
        Image depth;
        Image resolved; // sampled and copied
        VkFramebuffer framebuffer = VK_NULL_HANDLE;
        VkDescriptorSet texture = VK_NULL_HANDLE;
        int width = 0, height = 0;
    };
    bool prepareTarget(Target& target, int width, int height);
    void releaseTarget(Target& target);

    Target m_targets[2];
    int m_latest = -1, m_current = -1;
    VkCommandBuffer m_cmd = VK_NULL_HANDLE;
    uint32_t m_frameOffset = 0;
    VkExtent2D m_extent{0, 0};
};
}
