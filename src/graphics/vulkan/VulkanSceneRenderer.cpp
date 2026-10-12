#include "graphics/vulkan/VulkanSceneRenderer.h"

#include "graphics/BillboardMesh.h"
#include "graphics/CylinderMesh.h"
#include "graphics/LowPolyMesh.h"
#include "graphics/SceneBuffers.h"
#include "graphics/SphereMesh.h"

#include "vkshaders/atom.frag.h"
#include "vkshaders/atom.vert.h"
#include "vkshaders/billboard.frag.h"
#include "vkshaders/billboard.vert.h"
#include "vkshaders/bond.frag.h"
#include "vkshaders/bond.vert.h"
#include "vkshaders/depth_only.frag.h"
#include "vkshaders/halo.frag.h"
#include "vkshaders/halo.vert.h"
#include "vkshaders/line.frag.h"
#include "vkshaders/line.vert.h"
#include "vkshaders/mesh.frag.h"
#include "vkshaders/mesh.vert.h"
#include "vkshaders/shadow_atom.vert.h"
#include "vkshaders/shadow_billboard.frag.h"
#include "vkshaders/shadow_billboard.vert.h"
#include "vkshaders/shadow_bond.vert.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

namespace vkr
{
SceneRenderer& sceneRenderer()
{
    static SceneRenderer instance;
    return instance;
}

glm::mat4 clipCorrection(bool flipY)
{
    glm::mat4 c(1.0f);
    if (flipY) c[1][1] = -1.0f;
    c[2][2] = 0.5f;
    c[3][2] = 0.5f;
    return c;
}

SceneResources::~SceneResources()
{
    Context& ctx = context();
    for (Buffer* buffer : {&positions, &colors, &scales, &shininess, &bondStart, &bondEnd, &bondColorA, &bondColorB,
                           &bondRadius, &bondShininessA, &bondShininessB, &boxLines, &loop})
        ctx.destroyBuffer(*buffer);
}

namespace
{
// The frame uniforms (std140, as common.glsl).
struct FrameUniforms
{
    glm::mat4 projection, view, lightMVP, lightClip;
    glm::vec4 frustumPlanes[6];
    glm::vec4 lightPos, viewPos, params0, params1;
};

struct Push
{
    glm::vec4 color;
    glm::vec4 extra;
};

constexpr VkShaderStageFlags kPushStages = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
constexpr VkDeviceSize kStreamMinimum = 1024 * 1024;

std::array<glm::vec4, 6> frustumPlanes(const glm::mat4& vp)
{
    const glm::vec4 row0(vp[0][0], vp[1][0], vp[2][0], vp[3][0]);
    const glm::vec4 row1(vp[0][1], vp[1][1], vp[2][1], vp[3][1]);
    const glm::vec4 row2(vp[0][2], vp[1][2], vp[2][2], vp[3][2]);
    const glm::vec4 row3(vp[0][3], vp[1][3], vp[2][3], vp[3][3]);
    std::array<glm::vec4, 6> planes = {row3 + row0, row3 - row0, row3 + row1, row3 - row1, row3 + row2, row3 - row2};
    for (glm::vec4& plane : planes)
    {
        const float length = std::sqrt(plane.x * plane.x + plane.y * plane.y + plane.z * plane.z);
        if (length > 1e-6f) plane /= length;
    }
    return planes;
}

VkVertexInputBindingDescription binding(uint32_t index, uint32_t stride, bool instance)
{
    return {index, stride, instance ? VK_VERTEX_INPUT_RATE_INSTANCE : VK_VERTEX_INPUT_RATE_VERTEX};
}

VkVertexInputAttributeDescription attribute(uint32_t location, uint32_t bindingIndex, VkFormat format, uint32_t offset = 0)
{
    return {location, bindingIndex, format, offset};
}

// Per-atom streams: vertex positions, then instance position, colour, scale and shininess.
std::vector<VkVertexInputBindingDescription> atomBindings()
{
    return {binding(0, 12, false), binding(1, 16, true), binding(2, 16, true), binding(3, 4, true), binding(4, 4, true)};
}

std::vector<VkVertexInputBindingDescription> bondBindings()
{
    return {binding(0, 12, false), binding(1, 12, true), binding(2, 12, true), binding(3, 12, true), binding(4, 12, true),
            binding(5, 4, true), binding(6, 4, true), binding(7, 4, true)};
}

Buffer hostBuffer(const void* data, std::size_t bytes, VkBufferUsageFlags usage)
{
    Buffer buffer = context().createBuffer(std::max<std::size_t>(bytes, 16), usage, true);
    if (buffer.mapped && data && bytes) std::memcpy(buffer.mapped, data, bytes);
    return buffer;
}

template <typename T>
Buffer deviceBuffer(const std::vector<T>& values, VkBufferUsageFlags usage)
{
    if (values.empty()) return Buffer();
    return context().createDeviceBuffer(values.data(), values.size() * sizeof(T), usage);
}
}

enum class Blend { None, Alpha, ConstantAlpha };

struct SceneRenderer::PipelineSpec
{
    const uint32_t* vertexCode = nullptr;
    size_t vertexBytes = 0;
    const uint32_t* fragmentCode = nullptr;
    size_t fragmentBytes = 0;
    std::vector<VkVertexInputBindingDescription> bindings;
    std::vector<VkVertexInputAttributeDescription> attributes;
    VkPrimitiveTopology topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    bool shadowPass = false;  // the shadow map's pass (else the window's and the canvases')
    bool depthTest = true;
    bool depthWrite = true;
    VkCompareOp depthOp = VK_COMPARE_OP_LESS;
    Blend blend = Blend::None;
    bool lines = false;       // line width set per draw
    bool polygonOffset = false;
};

VkPipeline SceneRenderer::makePipeline(const PipelineSpec& spec)
{
    Context& ctx = context();
    VkShaderModule vertex = ctx.shaderModule(spec.vertexCode, spec.vertexBytes);
    VkShaderModule fragment = ctx.shaderModule(spec.fragmentCode, spec.fragmentBytes);
    m_modules.push_back(vertex);
    m_modules.push_back(fragment);
    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = vertex;
    stages[0].pName = "main";
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = fragment;
    stages[1].pName = "main";

    VkPipelineVertexInputStateCreateInfo input{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    input.vertexBindingDescriptionCount = static_cast<uint32_t>(spec.bindings.size());
    input.pVertexBindingDescriptions = spec.bindings.data();
    input.vertexAttributeDescriptionCount = static_cast<uint32_t>(spec.attributes.size());
    input.pVertexAttributeDescriptions = spec.attributes.data();
    VkPipelineInputAssemblyStateCreateInfo assembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology = spec.topology;
    VkPipelineViewportStateCreateInfo viewport{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewport.viewportCount = 1;
    viewport.scissorCount = 1;
    VkPipelineRasterizationStateCreateInfo raster{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    raster.polygonMode = VK_POLYGON_MODE_FILL;
    raster.lineWidth = 1.0f;
    if (spec.shadowPass)
    {
        // As the OpenGL shadow pass: back faces only, with a polygon offset.
        // Its projection is not flipped in y, which reverses the winding
        // Vulkan sees compared with OpenGL's: OpenGL's front faces are clockwise here.
        raster.frontFace = VK_FRONT_FACE_CLOCKWISE;
        raster.cullMode = VK_CULL_MODE_FRONT_BIT;
        raster.depthBiasEnable = VK_TRUE;
        raster.depthBiasConstantFactor = 6.0f;
        raster.depthBiasSlopeFactor = 3.0f;
    }
    else
    {
        // The projection's y flip keeps OpenGL's winding (gl_FrontFacing as in OpenGL).
        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
        raster.cullMode = VK_CULL_MODE_NONE;
        // Keeps the near plane from slicing spheres (GL_DEPTH_CLAMP in OpenGL).
        raster.depthClampEnable = ctx.depthClamp() ? VK_TRUE : VK_FALSE;
        if (spec.polygonOffset)
        {
            raster.depthBiasEnable = VK_TRUE;
            raster.depthBiasConstantFactor = 1.0f;
            raster.depthBiasSlopeFactor = 1.0f;
        }
    }
    VkPipelineMultisampleStateCreateInfo multisample{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples = spec.shadowPass ? VK_SAMPLE_COUNT_1_BIT : ctx.samples();
    VkPipelineDepthStencilStateCreateInfo depth{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depth.depthTestEnable = spec.depthTest ? VK_TRUE : VK_FALSE;
    depth.depthWriteEnable = spec.depthWrite ? VK_TRUE : VK_FALSE;
    depth.depthCompareOp = spec.depthTest ? spec.depthOp : VK_COMPARE_OP_ALWAYS;
    VkPipelineColorBlendAttachmentState colorBlend{};
    colorBlend.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    if (spec.blend != Blend::None)
    {
        const bool constant = spec.blend == Blend::ConstantAlpha;
        colorBlend.blendEnable = VK_TRUE;
        colorBlend.srcColorBlendFactor = constant ? VK_BLEND_FACTOR_CONSTANT_ALPHA : VK_BLEND_FACTOR_SRC_ALPHA;
        colorBlend.dstColorBlendFactor = constant ? VK_BLEND_FACTOR_ONE_MINUS_CONSTANT_ALPHA : VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlend.colorBlendOp = VK_BLEND_OP_ADD;
        colorBlend.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        colorBlend.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        colorBlend.alphaBlendOp = VK_BLEND_OP_ADD;
    }
    VkPipelineColorBlendStateCreateInfo blendState{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blendState.attachmentCount = spec.shadowPass ? 0u : 1u;
    blendState.pAttachments = &colorBlend;
    std::vector<VkDynamicState> dynamicStates = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    if (spec.lines && ctx.wideLines()) dynamicStates.push_back(VK_DYNAMIC_STATE_LINE_WIDTH);
    if (spec.blend == Blend::ConstantAlpha) dynamicStates.push_back(VK_DYNAMIC_STATE_BLEND_CONSTANTS);
    VkPipelineDynamicStateCreateInfo dynamic{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    dynamic.dynamicStateCount = static_cast<uint32_t>(dynamicStates.size());
    dynamic.pDynamicStates = dynamicStates.data();

    VkGraphicsPipelineCreateInfo info{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    info.stageCount = 2;
    info.pStages = stages;
    info.pVertexInputState = &input;
    info.pInputAssemblyState = &assembly;
    info.pViewportState = &viewport;
    info.pRasterizationState = &raster;
    info.pMultisampleState = &multisample;
    info.pDepthStencilState = &depth;
    info.pColorBlendState = &blendState;
    info.pDynamicState = &dynamic;
    info.layout = m_layout;
    info.renderPass = spec.shadowPass ? m_shadowPass : ctx.mainRenderPass();
    VkPipeline pipeline = VK_NULL_HANDLE;
    vkCreateGraphicsPipelines(ctx.device(), VK_NULL_HANDLE, 1, &info, nullptr, &pipeline);
    if (pipeline) m_pipelines.push_back(pipeline);
    return pipeline;
}

bool SceneRenderer::init(const SphereMesh& sphere, const LowPolyMesh& lowPoly, const BillboardMesh& billboard, const CylinderMesh& cylinder)
{
    Context& ctx = context();
    if (!ctx.ready()) return false;
    VkDevice device = ctx.device();
    VkPhysicalDeviceProperties properties;
    vkGetPhysicalDeviceProperties(ctx.physicalDevice(), &properties);
    m_uniformAlign = static_cast<uint32_t>(std::max<VkDeviceSize>(properties.limits.minUniformBufferOffsetAlignment, 16));

    // Geometry shared by every scene.
    const VkBufferUsageFlags vertexUsage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, indexUsage = VK_BUFFER_USAGE_INDEX_BUFFER_BIT;
    m_sphereVertices = deviceBuffer(sphere.vertices, vertexUsage);
    m_sphereIndices = deviceBuffer(sphere.indices, indexUsage);
    m_sphereIndexCount = static_cast<uint32_t>(sphere.indices.size());
    m_lowPolyVertices = deviceBuffer(lowPoly.vertices, vertexUsage);
    m_lowPolyIndices = deviceBuffer(lowPoly.indices, indexUsage);
    m_lowPolyIndexCount = static_cast<uint32_t>(lowPoly.indices.size());
    m_billboardVertices = deviceBuffer(billboard.vertices, vertexUsage);
    m_billboardIndices = deviceBuffer(billboard.indices, indexUsage);
    m_billboardIndexCount = static_cast<uint32_t>(billboard.indices.size());
    m_cylinderVertices = deviceBuffer(cylinder.vertices, vertexUsage);
    m_cylinderVertexCount = static_cast<uint32_t>(cylinder.vertices.size() / 3);
    const std::vector<float> corners = {-1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1};
    m_haloCorners = deviceBuffer(corners, vertexUsage);

    // Frame uniforms (a ring per frame in flight) and the shadow map sampler.
    VkDescriptorSetLayoutBinding layoutBindings[2]{};
    layoutBindings[0].binding = 0;
    layoutBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC;
    layoutBindings[0].descriptorCount = 1;
    layoutBindings[0].stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
    layoutBindings[1].binding = 1;
    layoutBindings[1].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    layoutBindings[1].descriptorCount = 1;
    layoutBindings[1].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
    VkDescriptorSetLayoutCreateInfo setInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    setInfo.bindingCount = 2;
    setInfo.pBindings = layoutBindings;
    vkCreateDescriptorSetLayout(device, &setInfo, nullptr, &m_setLayout);
    VkPushConstantRange pushRange{kPushStages, 0, sizeof(Push)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.setLayoutCount = 1;
    layoutInfo.pSetLayouts = &m_setLayout;
    layoutInfo.pushConstantRangeCount = 1;
    layoutInfo.pPushConstantRanges = &pushRange;
    vkCreatePipelineLayout(device, &layoutInfo, nullptr, &m_layout);

    // Shadow map: depth only, sampled with comparison (sampler2DShadow), white outside.
    m_shadowMap = ctx.createImage(kShadowSize, kShadowSize, VK_FORMAT_D32_SFLOAT,
                                  VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT,
                                  VK_SAMPLE_COUNT_1_BIT, VK_IMAGE_ASPECT_DEPTH_BIT);
    VkAttachmentDescription shadowAttachment{};
    shadowAttachment.format = VK_FORMAT_D32_SFLOAT;
    shadowAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
    shadowAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    shadowAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    shadowAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    shadowAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    shadowAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    shadowAttachment.finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_READ_ONLY_OPTIMAL;
    VkAttachmentReference shadowRef{0, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};
    VkSubpassDescription shadowSubpass{};
    shadowSubpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    shadowSubpass.pDepthStencilAttachment = &shadowRef;
    VkSubpassDependency shadowDependencies[2]{};
    shadowDependencies[0].srcSubpass = VK_SUBPASS_EXTERNAL;
    shadowDependencies[0].dstSubpass = 0;
    shadowDependencies[0].srcStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    shadowDependencies[0].dstStageMask = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT;
    shadowDependencies[0].srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
    shadowDependencies[0].dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_READ_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
    shadowDependencies[1].srcSubpass = 0;
    shadowDependencies[1].dstSubpass = VK_SUBPASS_EXTERNAL;
    shadowDependencies[1].srcStageMask = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT;
    shadowDependencies[1].dstStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    shadowDependencies[1].srcAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
    shadowDependencies[1].dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    VkRenderPassCreateInfo shadowPassInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    shadowPassInfo.attachmentCount = 1;
    shadowPassInfo.pAttachments = &shadowAttachment;
    shadowPassInfo.subpassCount = 1;
    shadowPassInfo.pSubpasses = &shadowSubpass;
    shadowPassInfo.dependencyCount = 2;
    shadowPassInfo.pDependencies = shadowDependencies;
    vkCreateRenderPass(device, &shadowPassInfo, nullptr, &m_shadowPass);
    VkFramebufferCreateInfo shadowFramebuffer{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
    shadowFramebuffer.renderPass = m_shadowPass;
    shadowFramebuffer.attachmentCount = 1;
    shadowFramebuffer.pAttachments = &m_shadowMap.view;
    shadowFramebuffer.width = shadowFramebuffer.height = kShadowSize;
    shadowFramebuffer.layers = 1;
    vkCreateFramebuffer(device, &shadowFramebuffer, nullptr, &m_shadowFramebuffer);
    // Readable before any shadow pass has run, and cleared to "no shadow".
    ctx.immediate([&](VkCommandBuffer cmd) {
        VkClearValue clear{};
        clear.depthStencil = {1.0f, 0};
        VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
        begin.renderPass = m_shadowPass;
        begin.framebuffer = m_shadowFramebuffer;
        begin.renderArea = {{0, 0}, {kShadowSize, kShadowSize}};
        begin.clearValueCount = 1;
        begin.pClearValues = &clear;
        vkCmdBeginRenderPass(cmd, &begin, VK_SUBPASS_CONTENTS_INLINE);
        vkCmdEndRenderPass(cmd);
    });
    VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    samplerInfo.magFilter = samplerInfo.minFilter = VK_FILTER_LINEAR;
    samplerInfo.addressModeU = samplerInfo.addressModeV = samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_BORDER;
    samplerInfo.borderColor = VK_BORDER_COLOR_FLOAT_OPAQUE_WHITE;
    samplerInfo.compareEnable = VK_TRUE;
    samplerInfo.compareOp = VK_COMPARE_OP_LESS_OR_EQUAL;
    samplerInfo.maxLod = 1.0f;
    vkCreateSampler(device, &samplerInfo, nullptr, &m_shadowSampler);
    VkSamplerCreateInfo linearInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    linearInfo.magFilter = linearInfo.minFilter = VK_FILTER_LINEAR;
    linearInfo.addressModeU = linearInfo.addressModeV = linearInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    linearInfo.maxLod = 1.0f;
    vkCreateSampler(device, &linearInfo, nullptr, &m_linearSampler);

    for (uint32_t i = 0; i < Context::kFramesInFlight; ++i)
    {
        m_uniforms[i] = ctx.createBuffer(kUniformRing, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT, true);
        m_streams[i] = ctx.createBuffer(kStreamMinimum, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, true);
        VkDescriptorSetAllocateInfo allocate{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        allocate.descriptorPool = ctx.descriptorPool();
        allocate.descriptorSetCount = 1;
        allocate.pSetLayouts = &m_setLayout;
        vkAllocateDescriptorSets(device, &allocate, &m_sets[i]);
        VkDescriptorBufferInfo bufferInfo{m_uniforms[i].buffer, 0, sizeof(FrameUniforms)};
        VkDescriptorImageInfo imageInfo{m_shadowSampler, m_shadowMap.view, VK_IMAGE_LAYOUT_DEPTH_STENCIL_READ_ONLY_OPTIMAL};
        VkWriteDescriptorSet writes[2]{};
        writes[0].sType = writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].dstSet = writes[1].dstSet = m_sets[i];
        writes[0].dstBinding = 0;
        writes[0].descriptorCount = 1;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC;
        writes[0].pBufferInfo = &bufferInfo;
        writes[1].dstBinding = 1;
        writes[1].descriptorCount = 1;
        writes[1].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        writes[1].pImageInfo = &imageInfo;
        vkUpdateDescriptorSets(device, 2, writes, 0, nullptr);
    }

    // An off-screen pass compatible with the window's (same formats and
    // samples), ending in a sampled image that can also be copied out.
    {
        const bool multisampled = ctx.samples() != VK_SAMPLE_COUNT_1_BIT;
        VkAttachmentDescription attachments[3]{};
        attachments[0].format = ctx.colorFormat();
        attachments[0].samples = ctx.samples();
        attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        attachments[0].storeOp = multisampled ? VK_ATTACHMENT_STORE_OP_DONT_CARE : VK_ATTACHMENT_STORE_OP_STORE;
        attachments[0].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        attachments[0].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        attachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        attachments[0].finalLayout = multisampled ? VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL : VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        attachments[1] = attachments[0];
        attachments[1].format = ctx.depthFormat();
        attachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        attachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
        attachments[2] = attachments[0];
        attachments[2].samples = VK_SAMPLE_COUNT_1_BIT;
        attachments[2].loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        attachments[2].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        attachments[2].finalLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        VkAttachmentReference colorRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
        VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};
        VkAttachmentReference resolveRef{2, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
        VkSubpassDescription subpass{};
        subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
        subpass.colorAttachmentCount = 1;
        subpass.pColorAttachments = &colorRef;
        subpass.pDepthStencilAttachment = &depthRef;
        if (multisampled) subpass.pResolveAttachments = &resolveRef;
        VkSubpassDependency dependencies[2]{};
        dependencies[0].srcSubpass = VK_SUBPASS_EXTERNAL;
        dependencies[0].dstSubpass = 0;
        dependencies[0].srcStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT |
                                       VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT;
        dependencies[0].dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
        dependencies[0].srcAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_TRANSFER_READ_BIT;
        dependencies[0].dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        dependencies[1].srcSubpass = 0;
        dependencies[1].dstSubpass = VK_SUBPASS_EXTERNAL;
        dependencies[1].srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependencies[1].dstStageMask = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT;
        dependencies[1].srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        dependencies[1].dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_TRANSFER_READ_BIT;
        VkRenderPassCreateInfo passInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
        passInfo.attachmentCount = multisampled ? 3u : 2u;
        passInfo.pAttachments = attachments;
        passInfo.subpassCount = 1;
        passInfo.pSubpasses = &subpass;
        passInfo.dependencyCount = 2;
        passInfo.pDependencies = dependencies;
        vkCreateRenderPass(device, &passInfo, nullptr, &m_offscreenPass);
    }

    // Pipelines.
    const auto atomAttributes = std::vector<VkVertexInputAttributeDescription>{
        attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT), attribute(1, 1, VK_FORMAT_R32G32B32_SFLOAT),
        attribute(2, 2, VK_FORMAT_R32G32B32_SFLOAT), attribute(3, 3, VK_FORMAT_R32_SFLOAT), attribute(4, 4, VK_FORMAT_R32_SFLOAT)};
    const auto bondAttributes = std::vector<VkVertexInputAttributeDescription>{
        attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT), attribute(1, 1, VK_FORMAT_R32G32B32_SFLOAT),
        attribute(2, 2, VK_FORMAT_R32G32B32_SFLOAT), attribute(3, 3, VK_FORMAT_R32G32B32_SFLOAT),
        attribute(4, 4, VK_FORMAT_R32G32B32_SFLOAT), attribute(5, 5, VK_FORMAT_R32_SFLOAT),
        attribute(6, 6, VK_FORMAT_R32_SFLOAT), attribute(7, 7, VK_FORMAT_R32_SFLOAT)};
    const auto shadowAtomAttributes = std::vector<VkVertexInputAttributeDescription>{
        attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT), attribute(1, 1, VK_FORMAT_R32G32B32_SFLOAT), attribute(3, 3, VK_FORMAT_R32_SFLOAT)};
    const auto shadowBondAttributes = std::vector<VkVertexInputAttributeDescription>{
        attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT), attribute(1, 1, VK_FORMAT_R32G32B32_SFLOAT),
        attribute(2, 2, VK_FORMAT_R32G32B32_SFLOAT), attribute(5, 5, VK_FORMAT_R32_SFLOAT)};
#define SPV(vs, fs) spec.vertexCode = vs, spec.vertexBytes = sizeof(vs), spec.fragmentCode = fs, spec.fragmentBytes = sizeof(fs)
    {
        PipelineSpec spec;
        SPV(kSpv_atom_vert, kSpv_atom_frag);
        spec.bindings = atomBindings();
        spec.attributes = atomAttributes;
        m_atom = makePipeline(spec);
        // Ghost atoms: blended at a constant opacity over what is drawn, keeping its depth.
        spec.depthWrite = false;
        spec.blend = Blend::ConstantAlpha;
        m_ghostAtom = makePipeline(spec);
    }
    {
        PipelineSpec spec;
        SPV(kSpv_billboard_vert, kSpv_billboard_frag);
        spec.bindings = atomBindings();
        spec.attributes = atomAttributes;
        m_billboard = makePipeline(spec);
    }
    {
        PipelineSpec spec;
        SPV(kSpv_bond_vert, kSpv_bond_frag);
        spec.bindings = bondBindings();
        spec.attributes = bondAttributes;
        m_bond = makePipeline(spec);
    }
    {
        PipelineSpec spec;
        SPV(kSpv_line_vert, kSpv_line_frag);
        spec.bindings = {binding(0, 12, false)};
        spec.attributes = {attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT)};
        spec.lines = true;
        spec.topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST;
        m_lineList = makePipeline(spec);
        spec.topology = VK_PRIMITIVE_TOPOLOGY_LINE_STRIP;
        m_lineStrip = makePipeline(spec);
        // Overlays seen through the atoms (the dislocation loop, cutting-plane edges).
        spec.depthTest = spec.depthWrite = false;
        m_lineStripOverlay = makePipeline(spec);
        spec.topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST;
        m_lineListOverlay = makePipeline(spec);
    }
    {
        PipelineSpec spec;
        SPV(kSpv_mesh_vert, kSpv_mesh_frag);
        spec.bindings = {binding(0, sizeof(MeshVertex), false)};
        spec.attributes = {attribute(0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0),
                           attribute(1, 0, VK_FORMAT_R32G32B32_SFLOAT, static_cast<uint32_t>(sizeof(glm::vec3)))};
        m_mesh = makePipeline(spec);
        spec.depthWrite = false;
        spec.blend = Blend::Alpha;
        spec.polygonOffset = true;
        m_meshTranslucent = makePipeline(spec);
    }
    {
        // Halos: faint where hidden behind other atoms, full where visible; no depth writes.
        PipelineSpec spec;
        SPV(kSpv_halo_vert, kSpv_halo_frag);
        spec.bindings = {binding(0, 8, false), binding(1, 16, true)};
        spec.attributes = {attribute(0, 0, VK_FORMAT_R32G32_SFLOAT), attribute(1, 1, VK_FORMAT_R32G32B32A32_SFLOAT)};
        spec.depthWrite = false;
        spec.blend = Blend::Alpha;
        spec.depthOp = VK_COMPARE_OP_GREATER;
        m_haloHidden = makePipeline(spec);
        spec.depthOp = VK_COMPARE_OP_LESS_OR_EQUAL;
        m_haloVisible = makePipeline(spec);
    }
    {
        PipelineSpec spec;
        spec.shadowPass = true;
        SPV(kSpv_shadow_atom_vert, kSpv_depth_only_frag);
        spec.bindings = atomBindings();
        spec.attributes = shadowAtomAttributes;
        m_shadowAtom = makePipeline(spec);
        SPV(kSpv_shadow_billboard_vert, kSpv_shadow_billboard_frag);
        m_shadowBillboard = makePipeline(spec);
        SPV(kSpv_shadow_bond_vert, kSpv_depth_only_frag);
        spec.bindings = bondBindings();
        spec.attributes = shadowBondAttributes;
        m_shadowBond = makePipeline(spec);
    }
#undef SPV
    m_ready = m_pipelines.size() == 15;
    return m_ready;
}

void SceneRenderer::shutdown()
{
    Context& ctx = context();
    if (!ctx.ready()) return;
    VkDevice device = ctx.device();
    vkDeviceWaitIdle(device);
    for (VkPipeline pipeline : m_pipelines) vkDestroyPipeline(device, pipeline, nullptr);
    for (VkShaderModule module : m_modules) vkDestroyShaderModule(device, module, nullptr);
    for (Buffer* buffer : {&m_sphereVertices, &m_sphereIndices, &m_lowPolyVertices, &m_lowPolyIndices, &m_billboardVertices,
                           &m_billboardIndices, &m_cylinderVertices, &m_haloCorners})
        ctx.destroyBuffer(*buffer);
    for (uint32_t i = 0; i < Context::kFramesInFlight; ++i)
    {
        ctx.destroyBuffer(m_uniforms[i]);
        ctx.destroyBuffer(m_streams[i]);
    }
    if (m_shadowFramebuffer) vkDestroyFramebuffer(device, m_shadowFramebuffer, nullptr);
    if (m_shadowPass) vkDestroyRenderPass(device, m_shadowPass, nullptr);
    if (m_offscreenPass) vkDestroyRenderPass(device, m_offscreenPass, nullptr);
    if (m_shadowSampler) vkDestroySampler(device, m_shadowSampler, nullptr);
    if (m_linearSampler) vkDestroySampler(device, m_linearSampler, nullptr);
    ctx.destroyImage(m_shadowMap);
    if (m_layout) vkDestroyPipelineLayout(device, m_layout, nullptr);
    if (m_setLayout) vkDestroyDescriptorSetLayout(device, m_setLayout, nullptr);
    *this = SceneRenderer();
}

void SceneRenderer::upload(SceneBuffers& buffers, const std::vector<glm::vec4>& positions, const std::vector<glm::vec4>& colors,
                           const std::vector<float>& scales, const std::vector<float>& shininess,
                           const std::vector<glm::vec3>& bondStart, const std::vector<glm::vec3>& bondEnd,
                           const std::vector<glm::vec3>& bondColorA, const std::vector<glm::vec3>& bondColorB,
                           const std::vector<float>& bondRadius, const std::vector<float>& bondShininessA,
                           const std::vector<float>& bondShininessB, const std::vector<glm::vec3>& boxLines)
{
    // A new set of buffers; the old ones go when no frame in flight uses them.
    auto scene = std::make_shared<SceneResources>();
    const VkBufferUsageFlags vertex = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT;
    scene->atomCount = static_cast<uint32_t>(std::min({positions.size(), colors.size(), scales.size(), shininess.size()}));
    if (scene->atomCount)
    {
        scene->positions = hostBuffer(positions.data(), positions.size() * sizeof(glm::vec4), vertex);
        scene->colors = hostBuffer(colors.data(), colors.size() * sizeof(glm::vec4), vertex);
        scene->scales = hostBuffer(scales.data(), scales.size() * sizeof(float), vertex);
        scene->shininess = hostBuffer(shininess.data(), shininess.size() * sizeof(float), vertex);
    }
    scene->bondCount = static_cast<uint32_t>(bondStart.size());
    if (scene->bondCount)
    {
        scene->bondStart = deviceBuffer(bondStart, vertex);
        scene->bondEnd = deviceBuffer(bondEnd, vertex);
        scene->bondColorA = deviceBuffer(bondColorA, vertex);
        scene->bondColorB = deviceBuffer(bondColorB, vertex);
        scene->bondRadius = deviceBuffer(bondRadius, vertex);
        scene->bondShininessA = deviceBuffer(bondShininessA, vertex);
        scene->bondShininessB = deviceBuffer(bondShininessB, vertex);
    }
    scene->boxVertexCount = static_cast<uint32_t>(boxLines.size());
    if (scene->boxVertexCount) scene->boxLines = hostBuffer(boxLines.data(), boxLines.size() * sizeof(glm::vec3), vertex);
    if (buffers.vulkan && buffers.vulkan->loopVertexCount)
    {
        // The dislocation loop is uploaded on its own; keep it.
        scene->loop = buffers.vulkan->loop;
        scene->loopVertexCount = buffers.vulkan->loopVertexCount;
        buffers.vulkan->loop = Buffer();
    }
    buffers.vulkan = scene;
}

void SceneRenderer::setInstanceColor(SceneBuffers& buffers, int index, const glm::vec3& color)
{
    SceneResources* scene = buffers.vulkan.get();
    if (!scene || index < 0 || static_cast<uint32_t>(index) >= scene->atomCount || !scene->colors.mapped) return;
    const glm::vec4 packed(color, 1.0f);
    std::memcpy(static_cast<char*>(scene->colors.mapped) + static_cast<std::size_t>(index) * sizeof(glm::vec4), &packed, sizeof(packed));
}

void SceneRenderer::setInstancePosition(SceneBuffers& buffers, int index, const glm::vec3& position)
{
    SceneResources* scene = buffers.vulkan.get();
    if (!scene || index < 0 || static_cast<uint32_t>(index) >= scene->atomCount || !scene->positions.mapped) return;
    const glm::vec4 packed(position, 1.0f);
    std::memcpy(static_cast<char*>(scene->positions.mapped) + static_cast<std::size_t>(index) * sizeof(glm::vec4), &packed, sizeof(packed));
}

void SceneRenderer::uploadLoop(SceneBuffers& buffers, const std::vector<glm::vec3>& points)
{
    if (!buffers.vulkan) buffers.vulkan = std::make_shared<SceneResources>();
    SceneResources& scene = *buffers.vulkan;
    context().destroyBuffer(scene.loop);
    scene.loopVertexCount = 0;
    if (points.size() < 2) return;
    // A closed strip: the first point again at the end.
    std::vector<glm::vec3> strip = points;
    strip.push_back(points.front());
    scene.loop = hostBuffer(strip.data(), strip.size() * sizeof(glm::vec3), VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
    scene.loopVertexCount = static_cast<uint32_t>(strip.size());
}

void SceneRenderer::newFrame()
{
    m_uniformCursor = 0;
    m_streamCursor = 0;
}

uint32_t SceneRenderer::allocateUniform(const void* data, uint32_t bytes)
{
    Buffer& ring = m_uniforms[context().frameIndex()];
    const uint32_t size = (bytes + m_uniformAlign - 1) / m_uniformAlign * m_uniformAlign;
    if (m_uniformCursor + size > kUniformRing) m_uniformCursor = 0;  // more views than expected: reuse (no frame uses that many)
    const uint32_t offset = m_uniformCursor;
    m_uniformCursor += size;
    std::memcpy(static_cast<char*>(ring.mapped) + offset, data, bytes);
    return offset;
}

uint32_t SceneRenderer::prepare(const SceneView& view, const SceneLighting& lighting)
{
    FrameUniforms uniforms{};
    uniforms.projection = clipCorrection(true) * view.projection;
    uniforms.view = view.view;
    uniforms.lightMVP = view.lightMVP;
    uniforms.lightClip = clipCorrection(false) * view.lightMVP;
    const auto planes = frustumPlanes(view.projection * view.view);
    for (int i = 0; i < 6; ++i) uniforms.frustumPlanes[i] = planes[static_cast<std::size_t>(i)];
    uniforms.lightPos = glm::vec4(view.lightPosition, 1.0f);
    uniforms.viewPos = glm::vec4(view.cameraPosition, 1.0f);
    uniforms.params0 = glm::vec4(lighting.ambient, lighting.saturation, lighting.contrast, lighting.specularIntensity);
    uniforms.params1 = glm::vec4(lighting.shadowStrength, lighting.shininessScale, lighting.shininessFloor, 1.0f / static_cast<float>(kShadowSize));
    return allocateUniform(&uniforms, sizeof(uniforms));
}

VkDeviceSize SceneRenderer::stream(const void* data, VkDeviceSize bytes, VkBuffer& buffer)
{
    Context& ctx = context();
    Buffer& ring = m_streams[ctx.frameIndex()];
    const VkDeviceSize offset = (m_streamCursor + 15) / 16 * 16;
    if (offset + bytes > ring.size)
    {
        // Grow: the old stream stays until this frame is done with it.
        const VkDeviceSize size = std::max(ring.size * 2, (offset + bytes) * 2);
        Buffer grown = ctx.createBuffer(size, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, true);
        ctx.destroyBuffer(ring);
        ring = grown;
        m_streamCursor = 0;
        return stream(data, bytes, buffer);
    }
    std::memcpy(static_cast<char*>(ring.mapped) + offset, data, static_cast<std::size_t>(bytes));
    m_streamCursor = offset + bytes;
    buffer = ring.buffer;
    return offset;
}

void SceneRenderer::bindFrame(VkCommandBuffer cmd, uint32_t frameOffset)
{
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_layout, 0, 1, &m_sets[context().frameIndex()], 1, &frameOffset);
}

void SceneRenderer::setViewport(VkCommandBuffer cmd, VkExtent2D extent)
{
    const VkViewport viewport{0.0f, 0.0f, static_cast<float>(extent.width), static_cast<float>(extent.height), 0.0f, 1.0f};
    const VkRect2D scissor{{0, 0}, extent};
    vkCmdSetViewport(cmd, 0, 1, &viewport);
    vkCmdSetScissor(cmd, 0, 1, &scissor);
}

void SceneRenderer::setLineWidth(VkCommandBuffer cmd, float width)
{
    // 3.5 px as in OpenGL by default, where the device draws wide lines.
    if (context().wideLines()) vkCmdSetLineWidth(cmd, std::min(width, context().maxLineWidth()));
}

uint32_t SceneRenderer::bindAtoms(VkCommandBuffer cmd, const SceneResources& scene, bool billboards, bool lowPoly)
{
    const Buffer& vertices = billboards ? m_billboardVertices : lowPoly ? m_lowPolyVertices : m_sphereVertices;
    const Buffer& indices = billboards ? m_billboardIndices : lowPoly ? m_lowPolyIndices : m_sphereIndices;
    const VkBuffer streams[] = {vertices.buffer, scene.positions.buffer, scene.colors.buffer, scene.scales.buffer, scene.shininess.buffer};
    const VkDeviceSize offsets[5] = {};
    vkCmdBindVertexBuffers(cmd, 0, 5, streams, offsets);
    vkCmdBindIndexBuffer(cmd, indices.buffer, 0, VK_INDEX_TYPE_UINT32);
    return billboards ? m_billboardIndexCount : lowPoly ? m_lowPolyIndexCount : m_sphereIndexCount;
}

void SceneRenderer::bindBonds(VkCommandBuffer cmd, const SceneResources& scene)
{
    const VkBuffer streams[] = {m_cylinderVertices.buffer, scene.bondStart.buffer, scene.bondEnd.buffer, scene.bondColorA.buffer,
                                scene.bondColorB.buffer, scene.bondRadius.buffer, scene.bondShininessA.buffer, scene.bondShininessB.buffer};
    const VkDeviceSize offsets[8] = {};
    vkCmdBindVertexBuffers(cmd, 0, 8, streams, offsets);
}

void SceneRenderer::renderShadowMap(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags)
{
    const SceneResources* scene = buffers.vulkan.get();
    VkClearValue clear{};
    clear.depthStencil = {1.0f, 0};
    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = m_shadowPass;
    begin.framebuffer = m_shadowFramebuffer;
    begin.renderArea = {{0, 0}, {kShadowSize, kShadowSize}};
    begin.clearValueCount = 1;
    begin.pClearValues = &clear;
    vkCmdBeginRenderPass(cmd, &begin, VK_SUBPASS_CONTENTS_INLINE);
    if (scene && flags.atoms)
    {
        setViewport(cmd, {kShadowSize, kShadowSize});
        bindFrame(cmd, frameOffset);
        if (scene->atomCount)
        {
            const bool billboards = !flags.spheresOnly && buffers.renderMode == RenderingMode::BillboardImposters;
            const bool lowPoly = !flags.spheresOnly && buffers.renderMode == RenderingMode::LowPolyInstancing;
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, billboards ? m_shadowBillboard : m_shadowAtom);
            vkCmdDrawIndexed(cmd, bindAtoms(cmd, *scene, billboards, lowPoly), scene->atomCount, 0, 0, 0);
        }
        if (flags.bonds && scene->bondCount)
        {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_shadowBond);
            bindBonds(cmd, *scene);
            vkCmdDraw(cmd, m_cylinderVertexCount, scene->bondCount, 0, 0);
        }
    }
    vkCmdEndRenderPass(cmd);
}

void SceneRenderer::drawScene(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags,
                              VkExtent2D target, const glm::vec3& boxColor)
{
    const SceneResources* scene = buffers.vulkan.get();
    if (!scene) return;
    setViewport(cmd, target);
    bindFrame(cmd, frameOffset);
    if (flags.atoms && flags.bonds && scene->bondCount)
    {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_bond);
        bindBonds(cmd, *scene);
        vkCmdDraw(cmd, m_cylinderVertexCount, scene->bondCount, 0, 0);
    }
    if (flags.atoms && scene->atomCount)
    {
        const bool billboards = !flags.spheresOnly && buffers.renderMode == RenderingMode::BillboardImposters;
        const bool lowPoly = !flags.spheresOnly && buffers.renderMode == RenderingMode::LowPolyInstancing;
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, billboards ? m_billboard : m_atom);
        vkCmdDrawIndexed(cmd, bindAtoms(cmd, *scene, billboards, lowPoly), scene->atomCount, 0, 0, 0);
    }
    if (flags.box && scene->boxVertexCount)
    {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_lineList);
        setLineWidth(cmd);
        const Push push{glm::vec4(boxColor, 1.0f), glm::vec4(0.0f)};
        vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(cmd, 0, 1, &scene->boxLines.buffer, &offset);
        vkCmdDraw(cmd, scene->boxVertexCount, 1, 0, 0);
    }
    if (flags.dislocation && scene->loopVertexCount >= 3)
    {
        // The dislocation loop shows through the atoms, as in OpenGL.
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_lineStripOverlay);
        setLineWidth(cmd);
        const Push push{glm::vec4(1.0f, 0.25f, 0.1f, 1.0f), glm::vec4(0.0f)};
        vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(cmd, 0, 1, &scene->loop.buffer, &offset);
        vkCmdDraw(cmd, scene->loopVertexCount, 1, 0, 0);
    }
}

void SceneRenderer::drawGhostAtoms(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, float opacity, VkExtent2D target)
{
    const SceneResources* scene = buffers.vulkan.get();
    if (!scene || !scene->atomCount) return;
    setViewport(cmd, target);
    bindFrame(cmd, frameOffset);
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_ghostAtom);
    const float constants[4] = {0.0f, 0.0f, 0.0f, opacity};
    vkCmdSetBlendConstants(cmd, constants);
    vkCmdDrawIndexed(cmd, bindAtoms(cmd, *scene, false, false), scene->atomCount, 0, 0, 0);
}

void SceneRenderer::drawLines(VkCommandBuffer cmd, uint32_t frameOffset, const glm::vec3* points, std::size_t count,
                              LineTopology topology, const glm::vec3& color, bool depthTest, VkExtent2D target, float width)
{
    if (count < 2) return;
    VkBuffer buffer = VK_NULL_HANDLE;
    const VkDeviceSize offset = stream(points, count * sizeof(glm::vec3), buffer);
    setViewport(cmd, target);
    bindFrame(cmd, frameOffset);
    const bool list = topology == LineTopology::List;
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS,
                      depthTest ? (list ? m_lineList : m_lineStrip) : (list ? m_lineListOverlay : m_lineStripOverlay));
    setLineWidth(cmd, width);
    const Push push{glm::vec4(color, 1.0f), glm::vec4(0.0f)};
    vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
    vkCmdBindVertexBuffers(cmd, 0, 1, &buffer, &offset);
    vkCmdDraw(cmd, static_cast<uint32_t>(count), 1, 0, 0);
}

void SceneRenderer::drawMesh(VkCommandBuffer cmd, uint32_t frameOffset, const MeshVertex* vertices, std::size_t count,
                             const glm::vec4& color, MeshStyle style, bool translucent, const glm::vec3& lightDirection,
                             VkExtent2D target)
{
    if (count < 3) return;
    VkBuffer buffer = VK_NULL_HANDLE;
    const VkDeviceSize offset = stream(vertices, count * sizeof(MeshVertex), buffer);
    setViewport(cmd, target);
    bindFrame(cmd, frameOffset);
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, translucent ? m_meshTranslucent : m_mesh);
    const Push push{color, glm::vec4(static_cast<float>(static_cast<int>(style)), lightDirection)};
    vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
    vkCmdBindVertexBuffers(cmd, 0, 1, &buffer, &offset);
    vkCmdDraw(cmd, static_cast<uint32_t>(count / 3 * 3), 1, 0, 0);
}

void SceneRenderer::drawHalos(VkCommandBuffer cmd, uint32_t frameOffset, const std::vector<glm::vec3>& positions,
                              const std::vector<float>& radii, VkExtent2D target)
{
    if (positions.empty()) return;
    std::vector<glm::vec4> instances;
    instances.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) instances.emplace_back(positions[i], i < radii.size() ? radii[i] : 1.0f);
    VkBuffer buffer = VK_NULL_HANDLE;
    const VkDeviceSize offset = stream(instances.data(), instances.size() * sizeof(glm::vec4), buffer);
    setViewport(cmd, target);
    bindFrame(cmd, frameOffset);
    const VkBuffer streams[] = {m_haloCorners.buffer, buffer};
    const VkDeviceSize offsets[2] = {0, offset};
    vkCmdBindVertexBuffers(cmd, 0, 2, streams, offsets);
    constexpr float kExtent = 1.6f;  // quad half-size in silhouette radii
    const glm::vec4 color(1.0f, 0.66f, 0.12f, 1.0f);
    // Hidden parts faintly first (selected atoms behind others stay findable), then visible ones.
    Push push{color, glm::vec4(kExtent, 0.3f, 0.0f, 0.0f)};
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_haloHidden);
    vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
    vkCmdDraw(cmd, 6, static_cast<uint32_t>(instances.size()), 0, 0);
    push.extra.y = 1.0f;
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_haloVisible);
    vkCmdPushConstants(cmd, m_layout, kPushStages, 0, sizeof(push), &push);
    vkCmdDraw(cmd, 6, static_cast<uint32_t>(instances.size()), 0, 0);
}
}
