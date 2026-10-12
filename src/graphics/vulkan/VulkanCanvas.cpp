#include "graphics/vulkan/VulkanCanvas.h"

#include "imgui_impl_vulkan.h"

namespace vkr
{
bool Canvas::prepareTarget(Target& target, int width, int height)
{
    if (target.framebuffer && target.width == width && target.height == height) return true;
    releaseTarget(target);
    Context& ctx = context();
    const bool multisampled = ctx.samples() != VK_SAMPLE_COUNT_1_BIT;
    const VkImageUsageFlags finalUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT;
    const uint32_t w = static_cast<uint32_t>(width), h = static_cast<uint32_t>(height);
    if (multisampled)
        target.color = ctx.createImage(w, h, ctx.colorFormat(), VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT,
                                       ctx.samples(), VK_IMAGE_ASPECT_COLOR_BIT);
    target.depth = ctx.createImage(w, h, ctx.depthFormat(), VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT, ctx.samples(), VK_IMAGE_ASPECT_DEPTH_BIT);
    target.resolved = ctx.createImage(w, h, ctx.colorFormat(), finalUsage, VK_SAMPLE_COUNT_1_BIT, VK_IMAGE_ASPECT_COLOR_BIT);
    if (!target.depth.image || !target.resolved.image || (multisampled && !target.color.image))
    {
        releaseTarget(target);
        return false;
    }
    VkImageView views[3];
    uint32_t count = 0;
    if (multisampled)
    {
        views[count++] = target.color.view;
        views[count++] = target.depth.view;
        views[count++] = target.resolved.view;
    }
    else
    {
        views[count++] = target.resolved.view;
        views[count++] = target.depth.view;
    }
    VkFramebufferCreateInfo info{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
    info.renderPass = sceneRenderer().offscreenPass();
    info.attachmentCount = count;
    info.pAttachments = views;
    info.width = w;
    info.height = h;
    info.layers = 1;
    vkCreateFramebuffer(ctx.device(), &info, nullptr, &target.framebuffer);
    target.width = width;
    target.height = height;
    return target.framebuffer != VK_NULL_HANDLE;
}

void Canvas::releaseTarget(Target& target)
{
    Context& ctx = context();
    if (target.framebuffer || target.texture)
    {
        const VkFramebuffer framebuffer = target.framebuffer;
        const VkDescriptorSet texture = target.texture;
        const VkDevice device = ctx.device();
        ctx.destroyLater([device, framebuffer, texture] {
            if (framebuffer) vkDestroyFramebuffer(device, framebuffer, nullptr);
            if (texture && ImGui::GetCurrentContext()) ImGui_ImplVulkan_RemoveTexture(texture);
        });
    }
    ctx.destroyImage(target.color);
    ctx.destroyImage(target.depth);
    ctx.destroyImage(target.resolved);
    target = Target();
}

bool Canvas::begin(int width, int height, const glm::vec3& background, VkCommandBuffer cmd)
{
    Context& ctx = context();
    SceneRenderer& renderer = sceneRenderer();
    if (!ctx.ready() || !renderer.ready() || width <= 0 || height <= 0) return false;
    if (!cmd)
    {
        if (!ctx.inFrame()) return false;
        cmd = ctx.commandBuffer();
    }
    // The image not shown by the frame that may still be on the GPU.
    m_current = m_latest == 0 ? 1 : 0;
    Target& target = m_targets[m_current];
    if (!prepareTarget(target, width, height)) return false;
    m_cmd = cmd;
    m_extent = {static_cast<uint32_t>(width), static_cast<uint32_t>(height)};
    VkClearValue clears[3]{};
    clears[0].color = {{background.r, background.g, background.b, 1.0f}};
    clears[1].depthStencil = {1.0f, 0};
    clears[2].color = clears[0].color;
    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = renderer.offscreenPass();
    begin.framebuffer = target.framebuffer;
    begin.renderArea = {{0, 0}, m_extent};
    begin.clearValueCount = ctx.samples() != VK_SAMPLE_COUNT_1_BIT ? 3u : 2u;
    begin.pClearValues = clears;
    vkCmdBeginRenderPass(cmd, &begin, VK_SUBPASS_CONTENTS_INLINE);
    return true;
}

void Canvas::setView(const SceneView& view, const SceneLighting& lighting)
{
    m_frameOffset = sceneRenderer().prepare(view, lighting);
}

void Canvas::drawScene(const SceneBuffers& buffers, const SceneFlags& flags, const glm::vec3& boxColor)
{
    if (m_cmd) sceneRenderer().drawScene(m_cmd, buffers, m_frameOffset, flags, m_extent, boxColor);
}

void Canvas::drawGhostAtoms(const SceneBuffers& buffers, float opacity)
{
    if (m_cmd) sceneRenderer().drawGhostAtoms(m_cmd, buffers, m_frameOffset, opacity, m_extent);
}

void Canvas::drawLines(const std::vector<glm::vec3>& points, const glm::vec3& color, LineTopology topology, bool depthTest)
{
    if (m_cmd) sceneRenderer().drawLines(m_cmd, m_frameOffset, points.data(), points.size(), topology, color, depthTest, m_extent);
}

void Canvas::drawMesh(const std::vector<MeshVertex>& vertices, const glm::vec4& color, MeshStyle style, bool translucent,
                      const glm::vec3& lightDirection)
{
    if (m_cmd)
        sceneRenderer().drawMesh(m_cmd, m_frameOffset, vertices.data(), vertices.size(), color, style, translucent, lightDirection, m_extent);
}

void Canvas::end()
{
    if (!m_cmd) return;
    vkCmdEndRenderPass(m_cmd);
    m_latest = m_current;
    m_cmd = VK_NULL_HANDLE;
}

VkFormat Canvas::colorFormat() const
{
    return context().colorFormat();
}

void Canvas::copyTo(VkCommandBuffer cmd, const Buffer& buffer)
{
    if (m_latest < 0) return;
    const Target& target = m_targets[m_latest];
    VkImageMemoryBarrier barrier{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
    barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    barrier.oldLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
    barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    barrier.srcQueueFamilyIndex = barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.image = target.resolved.image;
    barrier.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
                         VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 0, nullptr, 0, nullptr, 1, &barrier);
    VkBufferImageCopy region{};
    region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
    region.imageExtent = {static_cast<uint32_t>(target.width), static_cast<uint32_t>(target.height), 1};
    vkCmdCopyImageToBuffer(cmd, target.resolved.image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, buffer.buffer, 1, &region);
    barrier.srcAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    barrier.newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 0, nullptr, 0, nullptr, 1, &barrier);
    VkBufferMemoryBarrier hostRead{VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER};
    hostRead.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    hostRead.dstAccessMask = VK_ACCESS_HOST_READ_BIT;
    hostRead.srcQueueFamilyIndex = hostRead.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    hostRead.buffer = buffer.buffer;
    hostRead.size = VK_WHOLE_SIZE;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_HOST_BIT, 0, 0, nullptr, 1, &hostRead, 0, nullptr);
}

ImTextureID Canvas::texture()
{
    if (m_latest < 0) return ImTextureID{};
    Target& target = m_targets[m_latest];
    if (!target.texture)
        target.texture = ImGui_ImplVulkan_AddTexture(sceneRenderer().linearSampler(), target.resolved.view,
                                                     VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL);
    return reinterpret_cast<ImTextureID>(target.texture);
}

void Canvas::image(const ImVec2& size)
{
    if (const ImTextureID id = texture()) ImGui::Image(id, size);
    else ImGui::Dummy(size);
}

void Canvas::release()
{
    if (!context().ready()) return;
    for (Target& target : m_targets) releaseTarget(target);
    m_latest = m_current = -1;
}
}
