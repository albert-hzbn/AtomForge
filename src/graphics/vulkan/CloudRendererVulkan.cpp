// The Vulkan side of CloudRenderer: chunk buffers streamed in with staging
// copies recorded before the render pass, point-sprite impostors, occlusion
// queries on the chunk boxes (read a few frames later, never stalling) and
// GPU timestamps; the same plan of chunks and atoms as with OpenGL.
#ifdef ATOMFORGE_VULKAN
#include "graphics/CloudRenderer.h"

#include "graphics/vulkan/VulkanCanvas.h"
#include "graphics/vulkan/VulkanSceneRenderer.h"

#include "vkshaders/cloud_box.vert.h"
#include "vkshaders/cloud_point.frag.h"
#include "vkshaders/cloud_point.vert.h"
#include "vkshaders/depth_only.frag.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <thread>

using atomforge::cloud::PackedAtom;

namespace
{
struct CloudUniforms
{
    glm::mat4 view, projection;
    glm::vec4 light, params, colouring;
};

struct ChunkPush
{
    glm::vec4 lower, scale;
};

constexpr VkShaderStageFlags kPushStages = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
}

struct CloudRenderer::VulkanState
{
    std::vector<vkr::Buffer> chunks;
    std::vector<std::uint64_t> issued;   // frame each pending occlusion query was recorded in
    std::vector<std::size_t> toQuery;    // chunks whose queries were reset for this frame
    std::vector<glm::vec3> outline;
    std::vector<vkr::Buffer> staging;    // uploads recorded, not yet freed
    vkr::Buffer species, box;
    VkDescriptorSetLayout setLayout = VK_NULL_HANDLE;
    VkPipelineLayout layout = VK_NULL_HANDLE;
    VkPipeline points = VK_NULL_HANDLE, boxes = VK_NULL_HANDLE;
    VkShaderModule modules[3]{};
    VkDescriptorSet sets[vkr::Context::kFramesInFlight]{};
    VkQueryPool occlusion = VK_NULL_HANDLE, timers = VK_NULL_HANDLE;
    std::uint64_t timerFrame[vkr::Context::kFramesInFlight]{};
    bool timerUsed[vkr::Context::kFramesInFlight]{};
    VkCommandBuffer cmd = VK_NULL_HANDLE;  // where uploads are recorded (in prepareVulkan)

    ~VulkanState()
    {
        vkr::Context& ctx = vkr::context();
        for (vkr::Buffer& buffer : chunks) ctx.destroyBuffer(buffer);
        for (vkr::Buffer& buffer : staging) ctx.destroyBuffer(buffer);
        ctx.destroyBuffer(species);
        ctx.destroyBuffer(box);
        const VkDevice device = ctx.device();
        const VkDescriptorPool pool = ctx.descriptorPool();
        std::vector<VkDescriptorSet> doomedSets;
        for (VkDescriptorSet set : sets)
            if (set) doomedSets.push_back(set);
        const VkDescriptorSetLayout doomedSetLayout = setLayout;
        const VkPipelineLayout doomedLayout = layout;
        const VkPipeline doomedPoints = points, doomedBoxes = boxes;
        const VkShaderModule m0 = modules[0], m1 = modules[1], m2 = modules[2];
        const VkQueryPool doomedOcclusion = occlusion, doomedTimers = timers;
        ctx.destroyLater([=] {
            if (!doomedSets.empty()) vkFreeDescriptorSets(device, pool, static_cast<uint32_t>(doomedSets.size()), doomedSets.data());
            for (VkPipeline pipeline : {doomedPoints, doomedBoxes})
                if (pipeline) vkDestroyPipeline(device, pipeline, nullptr);
            for (VkShaderModule module : {m0, m1, m2})
                if (module) vkDestroyShaderModule(device, module, nullptr);
            if (doomedLayout) vkDestroyPipelineLayout(device, doomedLayout, nullptr);
            if (doomedSetLayout) vkDestroyDescriptorSetLayout(device, doomedSetLayout, nullptr);
            if (doomedOcclusion) vkDestroyQueryPool(device, doomedOcclusion, nullptr);
            if (doomedTimers) vkDestroyQueryPool(device, doomedTimers, nullptr);
        });
    }
};

namespace
{
VkPipeline makeCloudPipeline(VkPipelineLayout layout, VkShaderModule vertex, VkShaderModule fragment, bool points)
{
    vkr::Context& ctx = vkr::context();
    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = vertex;
    stages[0].pName = "main";
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = fragment;
    stages[1].pName = "main";
    const VkVertexInputBindingDescription binding{0, points ? static_cast<uint32_t>(sizeof(PackedAtom)) : 12u, VK_VERTEX_INPUT_RATE_VERTEX};
    const VkVertexInputAttributeDescription attribute{0, 0, points ? VK_FORMAT_R16G16B16A16_UINT : VK_FORMAT_R32G32B32_SFLOAT, 0};
    VkPipelineVertexInputStateCreateInfo input{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    input.vertexBindingDescriptionCount = 1;
    input.pVertexBindingDescriptions = &binding;
    input.vertexAttributeDescriptionCount = 1;
    input.pVertexAttributeDescriptions = &attribute;
    VkPipelineInputAssemblyStateCreateInfo assembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology = points ? VK_PRIMITIVE_TOPOLOGY_POINT_LIST : VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    VkPipelineViewportStateCreateInfo viewport{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewport.viewportCount = 1;
    viewport.scissorCount = 1;
    VkPipelineRasterizationStateCreateInfo raster{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    raster.polygonMode = VK_POLYGON_MODE_FILL;
    raster.cullMode = VK_CULL_MODE_NONE;
    raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
    raster.lineWidth = 1.0f;
    VkPipelineMultisampleStateCreateInfo multisample{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples = ctx.samples();
    VkPipelineDepthStencilStateCreateInfo depth{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depth.depthTestEnable = VK_TRUE;
    depth.depthWriteEnable = points ? VK_TRUE : VK_FALSE;  // the query boxes only test
    depth.depthCompareOp = VK_COMPARE_OP_LESS;
    VkPipelineColorBlendAttachmentState colorBlend{};
    colorBlend.colorWriteMask = points ? (VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT) : 0u;
    VkPipelineColorBlendStateCreateInfo blend{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blend.attachmentCount = 1;
    blend.pAttachments = &colorBlend;
    const VkDynamicState dynamicStates[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    VkPipelineDynamicStateCreateInfo dynamic{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    dynamic.dynamicStateCount = 2;
    dynamic.pDynamicStates = dynamicStates;
    VkGraphicsPipelineCreateInfo info{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    info.stageCount = 2;
    info.pStages = stages;
    info.pVertexInputState = &input;
    info.pInputAssemblyState = &assembly;
    info.pViewportState = &viewport;
    info.pRasterizationState = &raster;
    info.pMultisampleState = &multisample;
    info.pDepthStencilState = &depth;
    info.pColorBlendState = &blend;
    info.pDynamicState = &dynamic;
    info.layout = layout;
    info.renderPass = ctx.mainRenderPass();
    VkPipeline pipeline = VK_NULL_HANDLE;
    vkCreateGraphicsPipelines(ctx.device(), VK_NULL_HANDLE, 1, &info, nullptr, &pipeline);
    return pipeline;
}

void transferBarrier(VkCommandBuffer cmd)
{
    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT | VK_ACCESS_TRANSFER_WRITE_BIT | VK_ACCESS_VERTEX_ATTRIBUTE_READ_BIT;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT | VK_PIPELINE_STAGE_VERTEX_INPUT_BIT,
                         0, 1, &barrier, 0, nullptr, 0, nullptr);
}
}

void CloudRenderer::openVulkan(const std::vector<glm::vec3>& outline)
{
    vkr::Context& ctx = vkr::context();
    vkr::SceneRenderer& scene = vkr::sceneRenderer();
    if (!ctx.ready() || !scene.ready()) throw std::runtime_error("The Vulkan renderer is not ready for the atom cloud");
    const auto& info = m_file->info();
    auto state = std::make_shared<VulkanState>();
    state->chunks.assign(info.chunks.size(), vkr::Buffer());
    state->issued.assign(info.chunks.size(), 0);
    state->outline = outline;
    const VkDevice device = ctx.device();

    // Species colours and radii.
    std::vector<glm::vec4> species(std::max<std::size_t>(1, info.species.size()), glm::vec4(0.7f));
    for (std::size_t i = 0; i < info.species.size(); ++i)
        species[i] = glm::vec4(info.species[i].color[0], info.species[i].color[1], info.species[i].color[2], info.species[i].radius);
    state->species = ctx.createDeviceBuffer(species.data(), species.size() * sizeof(glm::vec4), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    // Unit cube (triangles) for the chunk occlusion queries.
    {
        static const float corners[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
        static const int faces[12][3] = {{0, 2, 1}, {1, 2, 3}, {4, 5, 6}, {5, 7, 6}, {0, 1, 4}, {1, 5, 4},
                                         {2, 6, 3}, {3, 6, 7}, {0, 4, 2}, {2, 4, 6}, {1, 3, 5}, {3, 7, 5}};
        std::vector<float> triangles;
        for (const auto& f : faces)
            for (int v : f) triangles.insert(triangles.end(), corners[v], corners[v] + 3);
        state->box = ctx.createDeviceBuffer(triangles.data(), triangles.size() * sizeof(float), VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
    }

    VkDescriptorSetLayoutBinding bindings[2]{};
    bindings[0].binding = 0;
    bindings[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC;
    bindings[0].descriptorCount = 1;
    bindings[0].stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
    bindings[1].binding = 1;
    bindings[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings[1].descriptorCount = 1;
    bindings[1].stageFlags = VK_SHADER_STAGE_VERTEX_BIT;
    VkDescriptorSetLayoutCreateInfo setInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    setInfo.bindingCount = 2;
    setInfo.pBindings = bindings;
    vkCreateDescriptorSetLayout(device, &setInfo, nullptr, &state->setLayout);
    const VkPushConstantRange push{kPushStages, 0, sizeof(glm::mat4)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.setLayoutCount = 1;
    layoutInfo.pSetLayouts = &state->setLayout;
    layoutInfo.pushConstantRangeCount = 1;
    layoutInfo.pPushConstantRanges = &push;
    vkCreatePipelineLayout(device, &layoutInfo, nullptr, &state->layout);
    state->modules[0] = ctx.shaderModule(kSpv_cloud_point_vert, sizeof(kSpv_cloud_point_vert));
    state->modules[1] = ctx.shaderModule(kSpv_cloud_point_frag, sizeof(kSpv_cloud_point_frag));
    state->modules[2] = ctx.shaderModule(kSpv_cloud_box_vert, sizeof(kSpv_cloud_box_vert));
    VkShaderModule depthOnly = ctx.shaderModule(kSpv_depth_only_frag, sizeof(kSpv_depth_only_frag));
    state->points = makeCloudPipeline(state->layout, state->modules[0], state->modules[1], true);
    state->boxes = makeCloudPipeline(state->layout, state->modules[2], depthOnly, false);
    vkDestroyShaderModule(device, depthOnly, nullptr);

    for (uint32_t i = 0; i < vkr::Context::kFramesInFlight; ++i)
    {
        VkDescriptorSetAllocateInfo allocate{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        allocate.descriptorPool = ctx.descriptorPool();
        allocate.descriptorSetCount = 1;
        allocate.pSetLayouts = &state->setLayout;
        vkAllocateDescriptorSets(device, &allocate, &state->sets[i]);
        const VkDescriptorBufferInfo uniforms{scene.uniformRing(i), 0, sizeof(CloudUniforms)};
        const VkDescriptorBufferInfo speciesInfo{state->species.buffer, 0, VK_WHOLE_SIZE};
        VkWriteDescriptorSet writes[2]{};
        writes[0].sType = writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].dstSet = writes[1].dstSet = state->sets[i];
        writes[0].dstBinding = 0;
        writes[0].descriptorCount = 1;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC;
        writes[0].pBufferInfo = &uniforms;
        writes[1].dstBinding = 1;
        writes[1].descriptorCount = 1;
        writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[1].pBufferInfo = &speciesInfo;
        vkUpdateDescriptorSets(device, 2, writes, 0, nullptr);
    }

    VkQueryPoolCreateInfo occlusion{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    occlusion.queryType = VK_QUERY_TYPE_OCCLUSION;
    occlusion.queryCount = static_cast<uint32_t>(std::max<std::size_t>(1, info.chunks.size()));
    vkCreateQueryPool(device, &occlusion, nullptr, &state->occlusion);
    if (ctx.timestampPeriod() > 0.0f)
    {
        VkQueryPoolCreateInfo timers{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
        timers.queryType = VK_QUERY_TYPE_TIMESTAMP;
        timers.queryCount = 2 * vkr::Context::kFramesInFlight;
        vkCreateQueryPool(device, &timers, nullptr, &state->timers);
    }
    // Queries start reset (Vulkan requires a reset before first use).
    ctx.immediate([&](VkCommandBuffer cmd) {
        vkCmdResetQueryPool(cmd, state->occlusion, 0, occlusion.queryCount);
        if (state->timers) vkCmdResetQueryPool(cmd, state->timers, 0, 2 * vkr::Context::kFramesInFlight);
    });
    if (!state->points || !state->boxes || !state->occlusion) throw std::runtime_error("Cannot create the atom cloud pipelines");
    m_vk = state;
}

void CloudRenderer::closeVulkan()
{
    m_vk.reset();
}

void CloudRenderer::releaseVulkanChunk(std::size_t chunk)
{
    if (m_vk && chunk < m_vk->chunks.size()) vkr::context().destroyBuffer(m_vk->chunks[chunk]);
}

bool CloudRenderer::uploadVulkan(Load& load, std::uint64_t capacity)
{
    vkr::Context& ctx = vkr::context();
    ChunkState& c = m_chunks[load.chunk];
    vkr::Buffer& buffer = m_vk->chunks[load.chunk];
    VkCommandBuffer cmd = m_vk->cmd;
    if (capacity != c.capacity || !buffer.buffer)
    {
        vkr::Buffer grown = ctx.createBuffer(capacity * sizeof(PackedAtom),
                                             VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT, false);
        if (!grown.buffer) return false;  // out of GPU memory: the chunk keeps what it has
        if (buffer.buffer && c.resident)
        {
            const VkBufferCopy region{0, 0, c.resident * sizeof(PackedAtom)};
            vkCmdCopyBuffer(cmd, buffer.buffer, grown.buffer, 1, &region);
            transferBarrier(cmd);
        }
        ctx.destroyBuffer(buffer);
        buffer = grown;
    }
    const VkDeviceSize bytes = load.count * sizeof(PackedAtom);
    vkr::Buffer staging = ctx.createBuffer(bytes, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, true);
    if (!staging.mapped) return false;
    std::memcpy(staging.mapped, load.atoms.data(), static_cast<std::size_t>(bytes));
    const VkBufferCopy region{0, load.first * sizeof(PackedAtom), bytes};
    vkCmdCopyBuffer(cmd, staging.buffer, buffer.buffer, 1, &region);
    transferBarrier(cmd);
    m_vk->staging.push_back(staging);
    return true;
}

void CloudRenderer::readVulkanQueries()
{
    // Occlusion results of queries recorded at least two frames ago (their
    // frames have finished), as soon as they are available: no stall.
    const bool precise = vkr::context().preciseOcclusion();
    const VkDevice device = vkr::context().device();
    for (std::size_t i = 0; i < m_queryPending.size(); ++i) {
        if (!m_queryPending[i] || m_frame < m_vk->issued[i] + vkr::Context::kFramesInFlight) continue;
        std::uint64_t result[2] = {0, 0};
        const VkResult status = vkGetQueryPoolResults(device, m_vk->occlusion, static_cast<uint32_t>(i), 1, sizeof(result), result, sizeof(result),
                                                      VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WITH_AVAILABILITY_BIT);
        if ((status != VK_SUCCESS && status != VK_NOT_READY) || !result[1]) continue;
        m_occluded[i] = precise ? result[0] < 16 : result[0] == 0;
        m_queryPending[i] = 0;
    }
}

void CloudRenderer::prepareVulkan(VkCommandBuffer cmd, const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic,
                                  bool immediate)
{
    if (!m_file || !m_vk || width <= 0 || height <= 0) return;
    ++m_frame;
    m_stats.uploadedMBThisFrame = 0;
    m_vk->cmd = cmd;
    consumeLoads();
    m_vk->cmd = VK_NULL_HANDLE;
    if (!immediate) {
        // Freed once the frame has run.
        for (vkr::Buffer& buffer : m_vk->staging) vkr::context().destroyBuffer(buffer);
        m_vk->staging.clear();
    }
    readVulkanQueries();
    plan(projection, view, width, height, orthographic);

    // GPU time of the slot's last draw (its frame has finished), then this frame's.
    vkr::Context& ctx = vkr::context();
    const uint32_t slot = ctx.frameIndex();
    if (m_vk->timers)
    {
        if (m_vk->timerUsed[slot])
        {
            std::uint64_t stamps[4] = {};
            const VkResult status = vkGetQueryPoolResults(ctx.device(), m_vk->timers, slot * 2, 2, sizeof(stamps), stamps, 2 * sizeof(std::uint64_t),
                                                          VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WITH_AVAILABILITY_BIT);
            if (status == VK_SUCCESS && stamps[1] && stamps[3] && stamps[2] >= stamps[0])
                m_stats.gpuMilliseconds = static_cast<double>(stamps[2] - stamps[0]) * ctx.timestampPeriod() * 1e-6;
        }
        vkCmdResetQueryPool(cmd, m_vk->timers, slot * 2, 2);
    }
    // The occlusion queries recorded this frame: reset before the render pass.
    m_vk->toQuery.clear();
    if (settings.occlusionCulling)
        for (std::size_t i : m_tested)
            if (!m_queryPending[i]) {
                vkCmdResetQueryPool(cmd, m_vk->occlusion, static_cast<uint32_t>(i), 1);
                m_vk->toQuery.push_back(i);
            }
}

void CloudRenderer::drawVulkan(VkCommandBuffer cmd, int width, int height, bool lightTheme)
{
    if (!m_file || !m_vk || width <= 0 || height <= 0) return;
    vkr::Context& ctx = vkr::context();
    vkr::SceneRenderer& scene = vkr::sceneRenderer();
    const auto& info = m_file->info();
    const uint32_t slot = ctx.frameIndex();
    const VkExtent2D extent{static_cast<uint32_t>(width), static_cast<uint32_t>(height)};
    const VkViewport viewport{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f, 1.0f};
    const VkRect2D scissor{{0, 0}, extent};
    vkCmdSetViewport(cmd, 0, 1, &viewport);
    vkCmdSetScissor(cmd, 0, 1, &scissor);
    if (m_vk->timers)
    {
        vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, m_vk->timers, slot * 2);
        m_vk->timerUsed[slot] = true;
    }

    const glm::mat4 projection = vkr::clipCorrection(true) * m_projection;
    CloudUniforms uniforms{};
    uniforms.view = m_view;
    uniforms.projection = projection;
    uniforms.light = glm::vec4(glm::normalize(glm::vec3(0.45f, 0.6f, 0.65f)), 0.38f);
    uniforms.params = glm::vec4(m_pixelsPerUnit, m_orthographic ? 1.0f : 0.0f, settings.radiusScale, std::max(1.0f, ctx.maxPointSize()));
    uniforms.colouring = glm::vec4(static_cast<float>(colourAxis), 0.0f, 1.0f, 0.0f);
    if (colourAxis > 0) {
        const std::size_t k = static_cast<std::size_t>(colourAxis - 1);
        uniforms.colouring.y = static_cast<float>(info.lower[k]);
        uniforms.colouring.z = static_cast<float>(info.upper[k]);
    }
    const uint32_t offset = scene.allocateUniform(&uniforms, sizeof(uniforms));
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_vk->points);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_vk->layout, 0, 1, &m_vk->sets[slot], 1, &offset);
    std::uint64_t drawn = 0;
    for (const auto& v : m_visible) {
        const ChunkState& c = m_chunks[v.chunk];
        const std::uint64_t n = std::min<std::uint64_t>(c.resident, static_cast<std::uint64_t>(std::ceil(v.want)));
        if (!n || !m_vk->chunks[v.chunk].buffer) continue;
        const auto& chunk = info.chunks[v.chunk];
        // A subsample is drawn with larger points so that the chunk still looks solid.
        const double fraction = static_cast<double>(n) / static_cast<double>(chunk.count);
        ChunkPush push;
        push.lower = glm::vec4(static_cast<float>(chunk.lower[0]), static_cast<float>(chunk.lower[1]), static_cast<float>(chunk.lower[2]),
                               static_cast<float>(std::clamp(std::cbrt(1.0 / fraction), 1.0, 2.5)));
        push.scale = glm::vec4(static_cast<float>((chunk.upper[0] - chunk.lower[0]) / 65535.0), static_cast<float>((chunk.upper[1] - chunk.lower[1]) / 65535.0),
                               static_cast<float>((chunk.upper[2] - chunk.lower[2]) / 65535.0), 0.0f);
        vkCmdPushConstants(cmd, m_vk->layout, kPushStages, 0, sizeof(push), &push);
        const VkDeviceSize zero = 0;
        vkCmdBindVertexBuffers(cmd, 0, 1, &m_vk->chunks[v.chunk].buffer, &zero);
        vkCmdDraw(cmd, static_cast<uint32_t>(n), 1, 0, 0);
        drawn += n;
    }

    // Occlusion queries: each chunk's box against the atoms drawn now.
    std::size_t occluded = 0;
    if (settings.occlusionCulling) {
        const glm::mat4 viewProjection = projection * m_view;
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_vk->boxes);
        const VkDeviceSize zero = 0;
        vkCmdBindVertexBuffers(cmd, 0, 1, &m_vk->box.buffer, &zero);
        for (std::size_t i : m_tested)
            if (m_occluded[i]) ++occluded;
        const VkQueryControlFlags control = ctx.preciseOcclusion() ? VK_QUERY_CONTROL_PRECISE_BIT : 0;
        for (std::size_t i : m_vk->toQuery) {
            const auto& chunk = info.chunks[i];
            const glm::vec3 a = glm::vec3(chunk.lower[0], chunk.lower[1], chunk.lower[2]) - glm::vec3(m_maxRadius);
            const glm::vec3 b = glm::vec3(chunk.upper[0], chunk.upper[1], chunk.upper[2]) + glm::vec3(m_maxRadius);
            glm::mat4 model(1.0f);
            model[0][0] = b.x - a.x; model[1][1] = b.y - a.y; model[2][2] = b.z - a.z;
            model[3] = glm::vec4(a, 1.0f);
            const glm::mat4 mvp = viewProjection * model;
            vkCmdPushConstants(cmd, m_vk->layout, kPushStages, 0, sizeof(mvp), &mvp);
            vkCmdBeginQuery(cmd, m_vk->occlusion, static_cast<uint32_t>(i), control);
            vkCmdDraw(cmd, 36, 1, 0, 0);
            vkCmdEndQuery(cmd, m_vk->occlusion, static_cast<uint32_t>(i));
            m_queryPending[i] = 1;
            m_vk->issued[i] = m_frame;
        }
        m_vk->toQuery.clear();
    }

    if (settings.showBox) {
        vkr::SceneView view;
        view.projection = m_projection;
        view.view = m_view;
        const uint32_t frame = scene.prepare(view, vkr::SceneLighting());
        scene.drawLines(cmd, frame, m_vk->outline.data(), m_vk->outline.size(), vkr::LineTopology::List,
                        lightTheme ? glm::vec3(0.25f) : glm::vec3(0.85f), true, extent, 1.0f);
    }
    if (m_vk->timers) vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, m_vk->timers, slot * 2 + 1);
    finishStats(drawn, occluded);
}

bool CloudRenderer::captureVulkan(const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic, bool lightTheme,
                                  const glm::vec3& background, std::vector<unsigned char>& rgba, std::string& error)
{
    vkr::Context& ctx = vkr::context();
    if (!m_file || !m_vk) { error = "No atom cloud is open."; return false; }
    const VkDeviceSize bytes = static_cast<VkDeviceSize>(width) * static_cast<VkDeviceSize>(height) * 4;
    vkr::Buffer readback = ctx.createBuffer(bytes, VK_BUFFER_USAGE_TRANSFER_DST_BIT, true);
    if (!readback.mapped) { error = "Not enough memory to export an image of this size."; return false; }
    vkr::Canvas canvas;
    bool drawn = false;
    const double start = glfwGetTime();
    for (int pass = 0;; ++pass) {
        // The cloud at the export resolution: redraw while its atoms stream in (at most 30 s).
        drawn = false;
        ctx.immediate([&](VkCommandBuffer cmd) {
            prepareVulkan(cmd, projection, view, width, height, orthographic, true);
            if (!canvas.begin(width, height, background, cmd)) return;
            drawVulkan(cmd, width, height, lightTheme);
            canvas.end();
            canvas.copyTo(cmd, readback);
            drawn = true;
        });
        for (vkr::Buffer& buffer : m_vk->staging) ctx.destroyBufferNow(buffer);
        m_vk->staging.clear();
        if (!drawn || (pass > 2 && m_stats.complete) || glfwGetTime() - start > 30.0) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    if (drawn) {
        rgba.resize(static_cast<std::size_t>(bytes));
        std::memcpy(rgba.data(), readback.mapped, rgba.size());
        if (canvas.colorFormat() == VK_FORMAT_B8G8R8A8_UNORM || canvas.colorFormat() == VK_FORMAT_B8G8R8A8_SRGB)
            for (std::size_t i = 0; i < rgba.size(); i += 4) std::swap(rgba[i], rgba[i + 2]);
        for (std::size_t i = 3; i < rgba.size(); i += 4) rgba[i] = 255;
    } else
        error = "Failed to create the off-screen image for export (too large?).";
    ctx.destroyBuffer(readback);
    return drawn;
}
#endif
