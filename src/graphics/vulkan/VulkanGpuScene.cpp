// The GPU-driven main view of the Vulkan scene renderer. Each frame a compute
// shader (cull.comp) tests every atom and bond against the view frustum,
// gives it a level of detail from its radius on screen and appends it to a
// compact list drawn with indirect draws. Occlusion culling is two-phase
// (Haar & Aaltonen, "GPU-Driven Rendering Pipelines", SIGGRAPH 2015): atoms
// drawn last frame are drawn first, a hierarchical depth buffer (Hi-Z) is
// built from them, and the other atoms and the bonds are drawn only if not
// hidden behind it, so nothing pops in while the camera moves. Atoms go from
// a 3200-triangle sphere close up through icospheres of 1280, 320, 80 and 20
// triangles to point sprites a few pixels wide; bonds from 32 to 5 sides,
// and below a quarter of a pixel they are not drawn. The shadow pass is
// culled against the light's frustum with coarse meshes.
#include "graphics/vulkan/VulkanSceneRenderer.h"
#include "graphics/vulkan/VulkanPipelineSpec.h"

#include "graphics/CylinderMesh.h"
#include "graphics/SceneBuffers.h"
#include "graphics/SphereMesh.h"

#include "vkshaders/atom.frag.h"
#include "vkshaders/atom_gpu.frag.h"
#include "vkshaders/atom_gpu.vert.h"
#include "vkshaders/atom_point.frag.h"
#include "vkshaders/atom_point.vert.h"
#include "vkshaders/billboard.frag.h"
#include "vkshaders/billboard_depth.frag.h"
#include "vkshaders/billboard_gpu.frag.h"
#include "vkshaders/billboard_gpu.vert.h"
#include "vkshaders/bond.frag.h"
#include "vkshaders/bond_gpu.vert.h"
#include "vkshaders/cull.comp.h"
#include "vkshaders/cull_basic.comp.h"
#include "vkshaders/depth_only.frag.h"
#include "vkshaders/hiz_down.comp.h"
#include "vkshaders/point_depth.frag.h"
#include "vkshaders/hiz_first.comp.h"
#include "vkshaders/shadow_atom_gpu.vert.h"
#include "vkshaders/shadow_bond_gpu.vert.h"
#include "vkshaders/shadow_point.frag.h"
#include "vkshaders/shadow_point.vert.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>

namespace vkr
{
namespace
{
// The indirect draws of a frame.
constexpr uint32_t kShadowAtoms = 0;   // 3 levels
constexpr uint32_t kFirstAtoms = 3;    // 6 levels: atoms drawn last frame
constexpr uint32_t kSecondAtoms = 9;   // 6 levels: the others, after occlusion culling
constexpr uint32_t kBonds = 15;        // 3 levels
constexpr uint32_t kShadowBonds = 18;  // 1 level
constexpr uint32_t kCommandCount = 19;
constexpr uint32_t kAtomLevels = 6;
// GPU timestamps of a frame: start, shadow culled, shadow drawn, first phase
// culled, first phase drawn, Hi-Z built, second phase culled, end.
constexpr uint32_t kStamps = 8;

// Radius on screen (pixels) from which each atom level is used, finest first:
// the 3200-triangle sphere, icospheres of 1280 and 320 triangles, ray-cast
// sphere impostors (2 triangles; two levels, as the lists come in pairs),
// and below the last a point sprite (1 vertex).
constexpr float kAtomThresholds[kAtomLevels - 1] = {80.0f, 32.0f, 12.0f, 5.0f, 2.2f};
constexpr uint32_t kFirstImpostorLevel = 3;
// Bond radius on screen: 32, 12 and 5 sides; thinner than this, none.
constexpr float kBondThresholds[2] = {10.0f, 2.5f};
constexpr float kBondCullRadius = 0.25f;
// Atom radius in shadow-map texels: icospheres of 320 and 80 triangles, then
// a disc (a point sprite) at the sphere's side facing the light.
constexpr float kShadowThresholds[2] = {24.0f, 10.0f};
// Atoms from this level on (a few pixels wide) take one shadow-map tap.
constexpr uint32_t kFastShadowLevel = 3;
constexpr float kShadowBondCullRadius = 0.5f;

struct DrawCommand
{
    uint32_t indexCount, instanceCount, firstIndex;
    int32_t vertexOffset;
    uint32_t firstInstance;
};

struct CullUniforms
{
    glm::mat4 viewProjection;
    glm::vec4 planes[6];
    glm::vec4 screen;
    glm::vec4 scale;
    glm::vec4 thresholds[2];
    glm::vec4 inflation[2];
    glm::vec4 eye;
    glm::vec4 forward;
    glm::uvec4 slots[8];
    glm::uvec4 counts;
};

struct GpuPush
{
    glm::uvec4 list;
    glm::vec4 extra;
};

constexpr VkShaderStageFlags kPushStages = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;

std::array<glm::vec4, 6> planesOf(const glm::mat4& vp)
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

// Pixels per world unit at clip w = 1 for a viewport `height` pixels tall
// (the projection's y scale: rotations in the view keep its row's length).
float pixelScaleOf(const glm::mat4& vp, uint32_t height)
{
    return glm::length(glm::vec3(vp[0][1], vp[1][1], vp[2][1])) * static_cast<float>(height) * 0.5f;
}

// An icosphere subdivided `level` times (20 * 4^level triangles), unit radius.
void icosphere(int level, std::vector<float>& vertices, std::vector<uint32_t>& indices)
{
    const float t = 1.6180339887f;
    std::vector<glm::vec3> points = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                                     {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    for (glm::vec3& p : points) p = glm::normalize(p);
    std::vector<uint32_t> faces = {0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
                                   3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1};
    for (int pass = 0; pass < level; ++pass)
    {
        std::map<std::pair<uint32_t, uint32_t>, uint32_t> middles;
        const auto middle = [&](uint32_t a, uint32_t b) {
            const auto key = std::minmax(a, b);
            const auto found = middles.find(key);
            if (found != middles.end()) return found->second;
            points.push_back(glm::normalize(points[a] + points[b]));
            return middles[key] = static_cast<uint32_t>(points.size() - 1);
        };
        std::vector<uint32_t> next;
        next.reserve(faces.size() * 4);
        for (std::size_t f = 0; f < faces.size(); f += 3)
        {
            const uint32_t a = faces[f], b = faces[f + 1], c = faces[f + 2];
            const uint32_t ab = middle(a, b), bc = middle(b, c), ca = middle(c, a);
            next.insert(next.end(), {a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca});
        }
        faces.swap(next);
    }
    for (const glm::vec3& p : points) vertices.insert(vertices.end(), {p.x, p.y, p.z});
    indices = faces;
}

VkDescriptorSetLayout setLayout(VkDevice device, const std::vector<std::pair<VkDescriptorType, VkShaderStageFlags>>& bindings)
{
    std::vector<VkDescriptorSetLayoutBinding> list;
    for (std::size_t i = 0; i < bindings.size(); ++i)
    {
        VkDescriptorSetLayoutBinding b{};
        b.binding = static_cast<uint32_t>(i);
        b.descriptorType = bindings[i].first;
        b.descriptorCount = 1;
        b.stageFlags = bindings[i].second;
        list.push_back(b);
    }
    VkDescriptorSetLayoutCreateInfo info{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    info.bindingCount = static_cast<uint32_t>(list.size());
    info.pBindings = list.data();
    VkDescriptorSetLayout layout = VK_NULL_HANDLE;
    vkCreateDescriptorSetLayout(device, &info, nullptr, &layout);
    return layout;
}

VkPipelineLayout pipelineLayout(VkDevice device, const std::vector<VkDescriptorSetLayout>& sets, VkShaderStageFlags pushStages, uint32_t pushBytes)
{
    const VkPushConstantRange push{pushStages, 0, pushBytes};
    VkPipelineLayoutCreateInfo info{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    info.setLayoutCount = static_cast<uint32_t>(sets.size());
    info.pSetLayouts = sets.data();
    info.pushConstantRangeCount = 1;
    info.pPushConstantRanges = &push;
    VkPipelineLayout layout = VK_NULL_HANDLE;
    vkCreatePipelineLayout(device, &info, nullptr, &layout);
    return layout;
}

VkPipeline computePipeline(VkDevice device, VkPipelineLayout layout, const uint32_t* code, size_t bytes)
{
    VkShaderModule module = context().shaderModule(code, bytes);
    VkComputePipelineCreateInfo info{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    info.stage.module = module;
    info.stage.pName = "main";
    info.layout = layout;
    VkPipeline pipeline = VK_NULL_HANDLE;
    vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline);
    vkDestroyShaderModule(device, module, nullptr);
    return pipeline;
}

VkDescriptorSet allocateSet(VkDescriptorSetLayout layout)
{
    VkDescriptorSetAllocateInfo info{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    info.descriptorPool = context().descriptorPool();
    info.descriptorSetCount = 1;
    info.pSetLayouts = &layout;
    VkDescriptorSet set = VK_NULL_HANDLE;
    vkAllocateDescriptorSets(context().device(), &info, &set);
    return set;
}

void writeBuffers(VkDescriptorSet set, const std::vector<VkBuffer>& buffers)
{
    std::vector<VkDescriptorBufferInfo> infos;
    for (VkBuffer buffer : buffers) infos.push_back({buffer, 0, VK_WHOLE_SIZE});
    std::vector<VkWriteDescriptorSet> writes(infos.size());
    for (std::size_t i = 0; i < infos.size(); ++i)
    {
        writes[i] = VkWriteDescriptorSet{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        writes[i].dstSet = set;
        writes[i].dstBinding = static_cast<uint32_t>(i);
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &infos[i];
    }
    vkUpdateDescriptorSets(context().device(), static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);
}

void barrier(VkCommandBuffer cmd, VkPipelineStageFlags src, VkAccessFlags srcAccess, VkPipelineStageFlags dst, VkAccessFlags dstAccess)
{
    VkMemoryBarrier memory{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    memory.srcAccessMask = srcAccess;
    memory.dstAccessMask = dstAccess;
    vkCmdPipelineBarrier(cmd, src, dst, 0, 1, &memory, 0, nullptr, 0, nullptr);
}

constexpr VkPipelineStageFlags kDrawReads = VK_PIPELINE_STAGE_DRAW_INDIRECT_BIT | VK_PIPELINE_STAGE_VERTEX_SHADER_BIT;
constexpr VkAccessFlags kDrawReadAccess = VK_ACCESS_INDIRECT_COMMAND_READ_BIT | VK_ACCESS_SHADER_READ_BIT;
constexpr VkAccessFlags kComputeAccess = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
}

void SceneRenderer::initGpu()
{
    Context& ctx = context();
    const VkDevice device = ctx.device();
    VkPhysicalDeviceProperties properties;
    vkGetPhysicalDeviceProperties(ctx.physicalDevice(), &properties);
    if (properties.limits.maxPerStageDescriptorStorageBuffers < 8 || properties.limits.maxBoundDescriptorSets < 2)
        return;  // the vertex-input path draws everything instead
    // Occlusion culling reads the window's multisampled depth in a compute shader.
    VkFormatProperties depthFormat;
    vkGetPhysicalDeviceFormatProperties(ctx.physicalDevice(), ctx.depthFormat(), &depthFormat);
    m_occlusion = (depthFormat.optimalTilingFeatures & VK_FORMAT_FEATURE_SAMPLED_IMAGE_BIT) &&
                  (properties.limits.sampledImageDepthSampleCounts & ctx.samples());

    const VkShaderStageFlags vertex = VK_SHADER_STAGE_VERTEX_BIT, compute = VK_SHADER_STAGE_COMPUTE_BIT;
    const VkDescriptorType storage = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    m_drawSetLayout = setLayout(device, {{storage, vertex}, {storage, vertex}, {storage, vertex}, {storage, vertex},
                                         {storage, vertex}, {storage, vertex}, {storage, vertex}});
    m_gpuLayout = pipelineLayout(device, {m_setLayout, m_drawSetLayout}, kPushStages, sizeof(GpuPush));
    m_cullFrameLayout = setLayout(device, {{VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, compute}, {storage, compute},
                                           {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, compute}});
    m_cullSceneLayout = setLayout(device, {{storage, compute}, {storage, compute}, {storage, compute}, {storage, compute}, {storage, compute}});
    m_cullLayout = pipelineLayout(device, {m_cullFrameLayout, m_cullSceneLayout}, compute, sizeof(uint32_t));
    m_hizSetLayout = setLayout(device, {{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, compute}, {VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, compute},
                                        {VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, compute}});
    m_hizLayout = pipelineLayout(device, {m_hizSetLayout}, compute, sizeof(glm::ivec4));
    // One atomic per subgroup and level where subgroup ballots are supported in compute shaders.
    VkPhysicalDeviceSubgroupProperties subgroup{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceProperties2 properties2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    properties2.pNext = &subgroup;
    if (vkGetPhysicalDeviceProperties2) vkGetPhysicalDeviceProperties2(ctx.physicalDevice(), &properties2);
    const bool ballot = (subgroup.supportedStages & VK_SHADER_STAGE_COMPUTE_BIT) &&
                        (subgroup.supportedOperations & VK_SUBGROUP_FEATURE_BALLOT_BIT) &&
                        (subgroup.supportedOperations & VK_SUBGROUP_FEATURE_BASIC_BIT);
    m_cullPipeline = ballot ? computePipeline(device, m_cullLayout, kSpv_cull_comp, sizeof(kSpv_cull_comp))
                            : computePipeline(device, m_cullLayout, kSpv_cull_basic_comp, sizeof(kSpv_cull_basic_comp));
    m_hizFirstPipeline = computePipeline(device, m_hizLayout, kSpv_hiz_first_comp, sizeof(kSpv_hiz_first_comp));
    m_hizDownPipeline = computePipeline(device, m_hizLayout, kSpv_hiz_down_comp, sizeof(kSpv_hiz_down_comp));

    {
        PipelineSpec spec;
        spec.gpu = true;
        spec.bindings = {{0, 12, VK_VERTEX_INPUT_RATE_VERTEX}};
        spec.attributes = {{0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0}};
        spec.vertexCode = kSpv_atom_gpu_vert, spec.vertexBytes = sizeof(kSpv_atom_gpu_vert);
        spec.fragmentCode = kSpv_atom_gpu_frag, spec.fragmentBytes = sizeof(kSpv_atom_gpu_frag);
        m_gpuAtom = makePipeline(spec);
        spec.vertexCode = kSpv_billboard_gpu_vert, spec.vertexBytes = sizeof(kSpv_billboard_gpu_vert);
        spec.fragmentCode = kSpv_billboard_gpu_frag, spec.fragmentBytes = sizeof(kSpv_billboard_gpu_frag);
        m_gpuBillboard = makePipeline(spec);
        // The depth pre-pass of the atoms drawn last frame, then their shading
        // where they are in front (each pixel shaded about once).
        PipelineSpec prepass = spec;
        prepass.colorWrite = false;
        prepass.vertexCode = kSpv_atom_gpu_vert, prepass.vertexBytes = sizeof(kSpv_atom_gpu_vert);
        prepass.fragmentCode = kSpv_depth_only_frag, prepass.fragmentBytes = sizeof(kSpv_depth_only_frag);
        m_gpuAtomDepth = makePipeline(prepass);
        prepass.vertexCode = kSpv_billboard_gpu_vert, prepass.vertexBytes = sizeof(kSpv_billboard_gpu_vert);
        prepass.fragmentCode = kSpv_billboard_depth_frag, prepass.fragmentBytes = sizeof(kSpv_billboard_depth_frag);
        m_gpuBillboardDepth = makePipeline(prepass);
        PipelineSpec shading = spec;
        shading.depthWrite = false;
        shading.depthOp = VK_COMPARE_OP_LESS_OR_EQUAL;
        shading.vertexCode = kSpv_atom_gpu_vert, shading.vertexBytes = sizeof(kSpv_atom_gpu_vert);
        shading.fragmentCode = kSpv_atom_gpu_frag, shading.fragmentBytes = sizeof(kSpv_atom_gpu_frag);
        m_gpuAtomShade = makePipeline(shading);
        shading.vertexCode = kSpv_billboard_gpu_vert, shading.vertexBytes = sizeof(kSpv_billboard_gpu_vert);
        shading.fragmentCode = kSpv_billboard_gpu_frag, shading.fragmentBytes = sizeof(kSpv_billboard_gpu_frag);
        m_gpuBillboardShade = makePipeline(shading);
        spec.vertexCode = kSpv_bond_gpu_vert, spec.vertexBytes = sizeof(kSpv_bond_gpu_vert);
        spec.fragmentCode = kSpv_bond_frag, spec.fragmentBytes = sizeof(kSpv_bond_frag);
        m_gpuBond = makePipeline(spec);
        spec.shadowPass = true;
        spec.vertexCode = kSpv_shadow_atom_gpu_vert, spec.vertexBytes = sizeof(kSpv_shadow_atom_gpu_vert);
        spec.fragmentCode = kSpv_depth_only_frag, spec.fragmentBytes = sizeof(kSpv_depth_only_frag);
        m_gpuShadowAtom = makePipeline(spec);
        spec.vertexCode = kSpv_shadow_bond_gpu_vert, spec.vertexBytes = sizeof(kSpv_shadow_bond_gpu_vert);
        m_gpuShadowBond = makePipeline(spec);
        PipelineSpec points;
        points.gpu = true;
        points.topology = VK_PRIMITIVE_TOPOLOGY_POINT_LIST;
        points.vertexCode = kSpv_atom_point_vert, points.vertexBytes = sizeof(kSpv_atom_point_vert);
        points.fragmentCode = kSpv_atom_point_frag, points.fragmentBytes = sizeof(kSpv_atom_point_frag);
        m_gpuPoint = makePipeline(points);
        PipelineSpec pointDepth = points;
        pointDepth.colorWrite = false;
        pointDepth.fragmentCode = kSpv_point_depth_frag, pointDepth.fragmentBytes = sizeof(kSpv_point_depth_frag);
        m_gpuPointDepth = makePipeline(pointDepth);
        PipelineSpec pointShade = points;
        pointShade.depthWrite = false;
        pointShade.depthOp = VK_COMPARE_OP_LESS_OR_EQUAL;
        m_gpuPointShade = makePipeline(pointShade);
        points.shadowPass = true;
        points.vertexCode = kSpv_shadow_point_vert, points.vertexBytes = sizeof(kSpv_shadow_point_vert);
        points.fragmentCode = kSpv_shadow_point_frag, points.fragmentBytes = sizeof(kSpv_shadow_point_frag);
        m_gpuShadowPoint = makePipeline(points);
    }

    // The level-of-detail meshes, in one vertex and one index buffer each for
    // atoms and bonds (an indirect draw picks its range).
    {
        std::vector<float> vertices;
        std::vector<uint32_t> indices;
        const auto add = [&](const std::vector<float>& v, const std::vector<uint32_t>& i) {
            MeshRange range;
            range.firstIndex = static_cast<uint32_t>(indices.size());
            range.indexCount = static_cast<uint32_t>(i.size());
            range.vertexOffset = static_cast<int32_t>(vertices.size() / 3);
            range.triangles = range.indexCount / 3;
            vertices.insert(vertices.end(), v.begin(), v.end());
            indices.insert(indices.end(), i.begin(), i.end());
            return range;
        };
        const SphereMesh sphere(40, 40);
        m_atomMeshes[0] = add(sphere.vertices, std::vector<uint32_t>(sphere.indices.begin(), sphere.indices.end()));
        for (int level = 3, k = 1; level >= 0; --level, ++k)
        {
            std::vector<float> v;
            std::vector<uint32_t> i;
            icosphere(level, v, i);
            // As large as the sphere on average: the faces just inside it, the
            // corners just outside (an inscribed coarse mesh looks smaller, and
            // one around it larger, than the atom).
            float inner = 1.0f;
            for (std::size_t f = 0; f < i.size(); f += 3)
            {
                const glm::vec3 a(v[3 * i[f]], v[3 * i[f] + 1], v[3 * i[f] + 2]);
                const glm::vec3 b(v[3 * i[f + 1]], v[3 * i[f + 1] + 1], v[3 * i[f + 1] + 2]);
                const glm::vec3 c(v[3 * i[f + 2]], v[3 * i[f + 2] + 1], v[3 * i[f + 2] + 2]);
                inner = std::min(inner, std::abs(glm::dot(glm::normalize(glm::cross(b - a, c - a)), a)));
            }
            const float scale = 1.0f / std::sqrt(inner);
            for (float& x : v) x *= scale;
            m_atomMeshes[k] = add(v, i);
            m_atomInflation[k] = scale;  // the corners' distance
        }
        m_quadMesh = add({-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0}, {0, 1, 2, 0, 2, 3});
        m_atomMeshes[5] = m_quadMesh;  // a point: one index, no vertex input
        m_atomMeshes[5].indexCount = 1;
        m_atomMeshes[5].triangles = 0;
        m_lodVertices = ctx.createDeviceBuffer(vertices.data(), vertices.size() * sizeof(float), VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
        m_lodIndices = ctx.createDeviceBuffer(indices.data(), indices.size() * sizeof(uint32_t), VK_BUFFER_USAGE_INDEX_BUFFER_BIT);
    }
    {
        std::vector<float> vertices;
        std::vector<uint32_t> indices;
        int k = 0;
        for (int sides : {32, 12, 5, 6})
        {
            const CylinderMesh cylinder(sides);
            MeshRange range;
            range.firstIndex = static_cast<uint32_t>(indices.size());
            range.vertexOffset = static_cast<int32_t>(vertices.size() / 3);
            range.indexCount = static_cast<uint32_t>(cylinder.vertices.size() / 3);
            range.triangles = range.indexCount / 3;
            for (uint32_t i = 0; i < range.indexCount; ++i) indices.push_back(i);
            vertices.insert(vertices.end(), cylinder.vertices.begin(), cylinder.vertices.end());
            m_bondMeshes[k++] = range;
        }
        m_cylinderLodVertices = ctx.createDeviceBuffer(vertices.data(), vertices.size() * sizeof(float), VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
        m_cylinderLodIndices = ctx.createDeviceBuffer(indices.data(), indices.size() * sizeof(uint32_t), VK_BUFFER_USAGE_INDEX_BUFFER_BIT);
    }
    const float zeros[4] = {};
    m_dummy = ctx.createDeviceBuffer(zeros, sizeof(zeros), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);

    VkSamplerCreateInfo nearest{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    nearest.magFilter = nearest.minFilter = VK_FILTER_NEAREST;
    nearest.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
    nearest.addressModeU = nearest.addressModeV = nearest.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    nearest.maxLod = VK_LOD_CLAMP_NONE;
    vkCreateSampler(device, &nearest, nullptr, &m_nearestSampler);

    for (uint32_t i = 0; i < Context::kFramesInFlight; ++i)
    {
        m_commands[i] = ctx.createBuffer(kCommandCount * sizeof(DrawCommand),
                                         VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT | VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                                             VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT, false);
        m_statsReadback[i] = ctx.createBuffer(kCommandCount * sizeof(DrawCommand), VK_BUFFER_USAGE_TRANSFER_DST_BIT, true);
        m_cullFrameSets[i] = allocateSet(m_cullFrameLayout);
        const VkDescriptorBufferInfo uniforms{uniformRing(i), 0, sizeof(CullUniforms)};
        const VkDescriptorBufferInfo commands{m_commands[i].buffer, 0, VK_WHOLE_SIZE};
        VkWriteDescriptorSet writes[2]{};
        writes[0].sType = writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].dstSet = writes[1].dstSet = m_cullFrameSets[i];
        writes[0].dstBinding = 0;
        writes[0].descriptorCount = 1;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC;
        writes[0].pBufferInfo = &uniforms;
        writes[1].dstBinding = 1;
        writes[1].descriptorCount = 1;
        writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[1].pBufferInfo = &commands;
        vkUpdateDescriptorSets(device, 2, writes, 0, nullptr);
    }
    if (ctx.timestampPeriod() > 0.0f)
    {
        VkQueryPoolCreateInfo timers{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
        timers.queryType = VK_QUERY_TYPE_TIMESTAMP;
        timers.queryCount = kStamps * Context::kFramesInFlight;
        vkCreateQueryPool(device, &timers, nullptr, &m_timers);
        ctx.immediate([&](VkCommandBuffer cmd) { vkCmdResetQueryPool(cmd, m_timers, 0, timers.queryCount); });
    }
    m_gpuReady = m_cullPipeline && m_hizFirstPipeline && m_hizDownPipeline && m_gpuAtom && m_gpuBillboard && m_gpuBond &&
                 m_gpuShadowAtom && m_gpuShadowBond && m_gpuPoint && m_gpuShadowPoint && m_gpuAtomDepth && m_gpuBillboardDepth &&
                 m_gpuPointDepth && m_gpuAtomShade && m_gpuBillboardShade && m_gpuPointShade && m_lodVertices.buffer && m_cylinderLodVertices.buffer &&
                 m_commands[0].buffer && m_commands[1].buffer && m_statsReadback[0].mapped && m_statsReadback[1].mapped;
}

void SceneRenderer::shutdownGpu()
{
    Context& ctx = context();
    const VkDevice device = ctx.device();
    releaseHiz();
    for (VkPipeline pipeline : {m_cullPipeline, m_hizFirstPipeline, m_hizDownPipeline})
        if (pipeline) vkDestroyPipeline(device, pipeline, nullptr);
    for (VkPipelineLayout layout : {m_gpuLayout, m_cullLayout, m_hizLayout})
        if (layout) vkDestroyPipelineLayout(device, layout, nullptr);
    for (VkDescriptorSetLayout layout : {m_drawSetLayout, m_cullFrameLayout, m_cullSceneLayout, m_hizSetLayout})
        if (layout) vkDestroyDescriptorSetLayout(device, layout, nullptr);
    for (Buffer* buffer : {&m_lodVertices, &m_lodIndices, &m_cylinderLodVertices, &m_cylinderLodIndices, &m_dummy})
        ctx.destroyBuffer(*buffer);
    for (uint32_t i = 0; i < Context::kFramesInFlight; ++i)
    {
        ctx.destroyBuffer(m_commands[i]);
        ctx.destroyBuffer(m_statsReadback[i]);
    }
    if (m_nearestSampler) vkDestroySampler(device, m_nearestSampler, nullptr);
    if (m_timers) vkDestroyQueryPool(device, m_timers, nullptr);
}

void SceneRenderer::createGpuResources(SceneResources& scene)
{
    if (!m_gpuReady) return;
    Context& ctx = context();
    const VkDeviceSize n = scene.atomCount, m = scene.bondCount;
    // Visible lists: three regions of n atoms (two levels each, one growing
    // from each end) and two of m bonds.
    scene.lists = ctx.createBuffer(std::max<VkDeviceSize>((3 * n + 2 * m) * sizeof(uint32_t), 16), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, false);
    scene.visibility = ctx.createBuffer(std::max<VkDeviceSize>(n * sizeof(uint32_t), 16),
                                        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT, false);
    if (!scene.lists.buffer || !scene.visibility.buffer) return;
    const VkBuffer visibility = scene.visibility.buffer;
    ctx.immediate([visibility](VkCommandBuffer cmd) { vkCmdFillBuffer(cmd, visibility, 0, VK_WHOLE_SIZE, 0); });
    const auto or_ = [this](const Buffer& buffer) { return buffer.buffer ? buffer.buffer : m_dummy.buffer; };
    scene.drawSet = allocateSet(m_drawSetLayout);
    scene.cullSet = allocateSet(m_cullSceneLayout);
    if (!scene.drawSet || !scene.cullSet) return;
    writeBuffers(scene.drawSet, {or_(scene.positions), or_(scene.colors), scene.lists.buffer, or_(scene.bondStart), or_(scene.bondEnd),
                                 or_(scene.bondColorA), or_(scene.bondColorB)});
    writeBuffers(scene.cullSet, {or_(scene.positions), or_(scene.bondStart), or_(scene.bondEnd), scene.lists.buffer, visibility});
}

void SceneRenderer::releaseHiz()
{
    Context& ctx = context();
    const VkDevice device = ctx.device();
    if (!m_hiz.image) return;
    vkDeviceWaitIdle(device);
    if (!m_hiz.sets.empty())
        vkFreeDescriptorSets(device, ctx.descriptorPool(), static_cast<uint32_t>(m_hiz.sets.size()), m_hiz.sets.data());
    for (VkImageView view : m_hiz.levels) vkDestroyImageView(device, view, nullptr);
    vkDestroyImageView(device, m_hiz.view, nullptr);
    vkDestroyImage(device, m_hiz.image, nullptr);
    vkFreeMemory(device, m_hiz.memory, nullptr);
    m_hiz = Hiz();
}

void SceneRenderer::ensureHiz(VkExtent2D extent)
{
    Context& ctx = context();
    const uint32_t width = std::max(1u, extent.width), height = std::max(1u, extent.height);  // level 0: one texel per pixel
    if (m_hiz.image && m_hiz.width == width && m_hiz.height == height && m_hiz.generation == ctx.swapchainGeneration()) return;
    releaseHiz();  // waits for the GPU: the sets are rewritten below (rare: on resizes)
    const VkDevice device = ctx.device();
    uint32_t levels = 1;
    while ((std::max(width, height) >> levels) > 0) ++levels;
    VkImageCreateInfo info{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    info.imageType = VK_IMAGE_TYPE_2D;
    info.format = VK_FORMAT_R32_SFLOAT;
    info.extent = {width, height, 1};
    info.mipLevels = levels;
    info.arrayLayers = 1;
    info.samples = VK_SAMPLE_COUNT_1_BIT;
    info.tiling = VK_IMAGE_TILING_OPTIMAL;
    info.usage = VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
    info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    if (vkCreateImage(device, &info, nullptr, &m_hiz.image) != VK_SUCCESS) { m_hiz = Hiz(); return; }
    VkMemoryRequirements requirements;
    vkGetImageMemoryRequirements(device, m_hiz.image, &requirements);
    VkMemoryAllocateInfo allocate{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocate.allocationSize = requirements.size;
    allocate.memoryTypeIndex = ctx.memoryType(requirements.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    vkAllocateMemory(device, &allocate, nullptr, &m_hiz.memory);
    vkBindImageMemory(device, m_hiz.image, m_hiz.memory, 0);
    VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    viewInfo.image = m_hiz.image;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = VK_FORMAT_R32_SFLOAT;
    viewInfo.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, levels, 0, 1};
    vkCreateImageView(device, &viewInfo, nullptr, &m_hiz.view);
    for (uint32_t level = 0; level < levels; ++level)
    {
        viewInfo.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, level, 1, 0, 1};
        VkImageView view = VK_NULL_HANDLE;
        vkCreateImageView(device, &viewInfo, nullptr, &view);
        m_hiz.levels.push_back(view);
    }
    m_hiz.width = width;
    m_hiz.height = height;
    m_hiz.generation = ctx.swapchainGeneration();
    m_hiz.fresh = true;

    // Building: level 0 from the depth, each level from the one above.
    const VkDescriptorImageInfo depth{m_nearestSampler, ctx.depthImage().view, VK_IMAGE_LAYOUT_DEPTH_STENCIL_READ_ONLY_OPTIMAL};
    for (uint32_t level = 0; level < levels; ++level)
    {
        VkDescriptorSet set = allocateSet(m_hizSetLayout);
        m_hiz.sets.push_back(set);
        const VkDescriptorImageInfo src{VK_NULL_HANDLE, m_hiz.levels[level == 0 ? 0 : level - 1], VK_IMAGE_LAYOUT_GENERAL};
        const VkDescriptorImageInfo dst{VK_NULL_HANDLE, m_hiz.levels[level], VK_IMAGE_LAYOUT_GENERAL};
        VkWriteDescriptorSet writes[3]{};
        for (VkWriteDescriptorSet& write : writes)
        {
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = set;
            write.descriptorCount = 1;
        }
        writes[0].dstBinding = 0;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        writes[0].pImageInfo = &depth;
        writes[1].dstBinding = 1;
        writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        writes[1].pImageInfo = &src;
        writes[2].dstBinding = 2;
        writes[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        writes[2].pImageInfo = &dst;
        vkUpdateDescriptorSets(device, 3, writes, 0, nullptr);
    }
    // Culling samples every level.
    const VkDescriptorImageInfo pyramid{m_nearestSampler, m_hiz.view, VK_IMAGE_LAYOUT_GENERAL};
    for (VkDescriptorSet set : m_cullFrameSets)
    {
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = set;
        write.dstBinding = 2;
        write.descriptorCount = 1;
        write.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        write.pImageInfo = &pyramid;
        vkUpdateDescriptorSets(device, 1, &write, 0, nullptr);
    }
    // In the general layout from the start (it is sampled before it is first built).
    ctx.immediate([this, levels](VkCommandBuffer cmd) {
        VkImageMemoryBarrier toGeneral{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
        toGeneral.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        toGeneral.newLayout = VK_IMAGE_LAYOUT_GENERAL;
        toGeneral.srcQueueFamilyIndex = toGeneral.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        toGeneral.image = m_hiz.image;
        toGeneral.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, levels, 0, 1};
        toGeneral.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, nullptr, 0, nullptr, 1, &toGeneral);
    });
    m_hiz.fresh = false;
}

void SceneRenderer::dispatchCull(VkCommandBuffer cmd, const SceneResources& scene, const CullSpec& spec)
{
    const uint32_t count = spec.kind == 0 ? scene.atomCount : scene.bondCount;
    if (!count) return;
    CullUniforms uniforms{};
    uniforms.viewProjection = clipCorrection(spec.flipY) * spec.viewProjection;
    const auto planes = planesOf(spec.viewProjection);
    for (int i = 0; i < 6; ++i) uniforms.planes[i] = planes[static_cast<std::size_t>(i)];
    uniforms.screen = glm::vec4(static_cast<float>(spec.target.width), static_cast<float>(spec.target.height),
                                static_cast<float>(std::max(1u, m_hiz.width)), static_cast<float>(std::max(1u, m_hiz.height)));
    uniforms.scale = glm::vec4(pixelScaleOf(spec.viewProjection, spec.target.height), static_cast<float>(std::max<std::size_t>(1, m_hiz.levels.size())),
                               spec.cullRadius, 0.0f);
    for (uint32_t i = 0; i < 8; ++i)
    {
        uniforms.thresholds[i >> 2][i & 3] = spec.thresholds[i];
        uniforms.inflation[i >> 2][i & 3] = spec.inflation[i];
        uniforms.slots[i] = spec.slots[i];
    }
    uniforms.counts = glm::uvec4(count, spec.levels, spec.phase, spec.occlusion && m_hiz.image ? 1u : 0u);
    const glm::mat4 inverseView = glm::inverse(spec.view);
    uniforms.eye = glm::vec4(glm::vec3(inverseView[3]), spec.orthographic ? 1.0f : 0.0f);
    uniforms.forward = glm::vec4(-glm::normalize(glm::vec3(inverseView[2])), 0.0f);
    const uint32_t offset = allocateUniform(&uniforms, sizeof(uniforms));
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_cullPipeline);
    const VkDescriptorSet sets[2] = {m_cullFrameSets[context().frameIndex()], scene.cullSet};
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_cullLayout, 0, 2, sets, 1, &offset);
    vkCmdPushConstants(cmd, m_cullLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &spec.kind);
    vkCmdDispatch(cmd, (count + 255) / 256, 1, 1);
}

void SceneRenderer::stamp(VkCommandBuffer cmd, uint32_t index)
{
    const uint32_t slot = context().frameIndex();
    if (m_timers && m_timerWritten[slot]) vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, m_timers, slot * kStamps + index);
}

void SceneRenderer::readStats()
{
    // The commands of this frame slot's last use (its frame has finished).
    Context& ctx = context();
    const uint32_t slot = ctx.frameIndex();
    if (!m_statsPending[slot]) return;
    m_statsPending[slot] = false;
    DrawCommand commands[kCommandCount];
    std::memcpy(commands, m_statsReadback[slot].mapped, sizeof(commands));
    SceneStats stats;
    for (uint32_t level = 0; level < kAtomLevels; ++level)
    {
        stats.atomLevels[level] = commands[kFirstAtoms + level].instanceCount + commands[kSecondAtoms + level].instanceCount;
        const bool impostor = (m_billboardFrame[slot] || level >= kFirstImpostorLevel) && level + 1 < kAtomLevels;
        stats.triangles += static_cast<uint64_t>(stats.atomLevels[level]) * (impostor ? 2u : m_atomMeshes[level].triangles);
    }
    for (uint32_t level = 0; level < 3; ++level)
    {
        stats.bondLevels[level] = commands[kBonds + level].instanceCount;
        stats.triangles += static_cast<uint64_t>(stats.bondLevels[level]) * m_bondMeshes[level].triangles;
    }
    stats.shadowAtoms = commands[kShadowAtoms].instanceCount + commands[kShadowAtoms + 1].instanceCount + commands[kShadowAtoms + 2].instanceCount;
    stats.shadowBonds = commands[kShadowBonds].instanceCount;
    for (uint32_t level = 0; level < kAtomLevels; ++level) stats.firstPhaseAtoms += commands[kFirstAtoms + level].instanceCount;
    stats.occlusion = m_occlusion;
    stats.atoms = m_statsAtoms[slot];
    stats.bonds = m_statsBonds[slot];
    stats.gpuMilliseconds = m_stats.gpuMilliseconds;
    if (m_timers && m_timerWritten[slot])
    {
        // Stages not run this frame have no timestamp: they take no time.
        uint64_t stamps[2 * kStamps] = {};
        const VkResult result = vkGetQueryPoolResults(ctx.device(), m_timers, slot * kStamps, kStamps, sizeof(stamps), stamps,
                                                      2 * sizeof(uint64_t), VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WITH_AVAILABILITY_BIT);
        if ((result == VK_SUCCESS || result == VK_NOT_READY) && stamps[1] && stamps[2 * (kStamps - 1) + 1])
        {
            const double period = ctx.timestampPeriod() * 1e-6;
            uint64_t previous = stamps[0];
            for (uint32_t k = 1; k < kStamps; ++k)
            {
                if (!stamps[2 * k + 1]) continue;
                const uint64_t now = std::max(stamps[2 * k], previous);
                stats.stageMilliseconds[k - 1] = static_cast<double>(now - previous) * period;
                previous = now;
            }
            stats.gpuMilliseconds = static_cast<double>(previous - stamps[0]) * period;
        }
    }
    stats.valid = true;
    m_stats = stats;
}

void SceneRenderer::beginGpuScene(VkCommandBuffer cmd, const SceneBuffers& buffers, const SceneView& view, uint32_t frameOffset,
                                  const SceneFlags& flags, VkExtent2D target)
{
    Context& ctx = context();
    const uint32_t slot = ctx.frameIndex();
    readStats();
    ensureHiz(target);
    flushPatches(cmd);
    const SceneResources* scene = buffers.vulkan.get();
    const bool billboards = buffers.renderMode == RenderingMode::BillboardImposters;
    m_billboardFrame[slot] = billboards;
    m_statsAtoms[slot] = scene ? scene->atomCount : 0;
    m_statsBonds[slot] = scene && flags.bonds ? scene->bondCount : 0;

    // This frame's indirect draws, with no instances yet.
    DrawCommand commands[kCommandCount]{};
    const auto set = [&](uint32_t index, const MeshRange& mesh) {
        commands[index] = {mesh.indexCount, 0, mesh.firstIndex, mesh.vertexOffset, 0};
    };
    set(kShadowAtoms, m_atomMeshes[2]);
    set(kShadowAtoms + 1, m_atomMeshes[3]);
    set(kShadowAtoms + 2, m_atomMeshes[5]);  // a point
    for (uint32_t level = 0; level < kAtomLevels; ++level)
    {
        const bool impostor = level >= kFirstImpostorLevel || billboards;
        const MeshRange& mesh = impostor && level + 1 < kAtomLevels ? m_quadMesh : m_atomMeshes[level];
        set(kFirstAtoms + level, mesh);
        set(kSecondAtoms + level, mesh);
    }
    for (uint32_t level = 0; level < 3; ++level) set(kBonds + level, m_bondMeshes[level]);
    set(kShadowBonds, m_bondMeshes[3]);
    // After every earlier use of the commands, lists and visibility (earlier frames' draws and culling).
    barrier(cmd, kDrawReads | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT,
            kDrawReadAccess | kComputeAccess | VK_ACCESS_TRANSFER_READ_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT | kComputeAccess);
    vkCmdUpdateBuffer(cmd, m_commands[slot].buffer, 0, sizeof(commands), commands);
    barrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_DRAW_INDIRECT_BIT, kComputeAccess | VK_ACCESS_INDIRECT_COMMAND_READ_BIT);
    if (m_timers)
    {
        vkCmdResetQueryPool(cmd, m_timers, slot * kStamps, kStamps);
        vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, m_timers, slot * kStamps);
        m_timerWritten[slot] = true;
    }

    const uint32_t n = scene ? scene->atomCount : 0, m = scene ? scene->bondCount : 0;
    const bool atoms = scene && flags.atoms && n && scene->cullSet;
    const bool bonds = atoms && flags.bonds && m;

    // The shadow map: culled against the light, coarse meshes. The light is
    // fixed to the structure, so the map stays valid until the light, the
    // atoms' positions or sizes, or what is shown change (cached shadow map).
    const ShadowKey key{view.lightMVP, scene ? scene->id : 0, scene ? scene->revision : 0, atoms, bonds};
    if (!(m_shadowValid && key == m_shadowKey))
    {
        if (atoms)
        {
            CullSpec shadow;
            shadow.viewProjection = view.lightMVP;
            shadow.flipY = false;
            shadow.target = {kShadowSize, kShadowSize};
            shadow.levels = 3;
            shadow.thresholds[0] = kShadowThresholds[0];
            shadow.thresholds[1] = kShadowThresholds[1];
            shadow.inflation[0] = m_atomInflation[2];
            shadow.inflation[1] = m_atomInflation[3];
            shadow.inflation[2] = 1.0f;
            shadow.slots[0] = {kShadowAtoms, 0, 1, 0};
            shadow.slots[1] = {kShadowAtoms + 1, n - 1, 0, 0};
            shadow.slots[2] = {kShadowAtoms + 2, n, 1, 0};
            dispatchCull(cmd, *scene, shadow);
            if (bonds)
            {
                CullSpec shadowBonds = shadow;
                shadowBonds.kind = 1;
                shadowBonds.levels = 1;
                shadowBonds.cullRadius = kShadowBondCullRadius;
                shadowBonds.slots[0] = {kShadowBonds, 3 * n, 1, 0};
                for (float& f : shadowBonds.inflation) f = 1.0f;
                dispatchCull(cmd, *scene, shadowBonds);
            }
            stamp(cmd, 1);
            barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT, kDrawReads, kDrawReadAccess);
        }
        VkClearValue clear{};
        clear.depthStencil = {1.0f, 0};
        VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
        begin.renderPass = m_shadowPass;
        begin.framebuffer = m_shadowFramebuffer;
        begin.renderArea = {{0, 0}, {kShadowSize, kShadowSize}};
        begin.clearValueCount = 1;
        begin.pClearValues = &clear;
        vkCmdBeginRenderPass(cmd, &begin, VK_SUBPASS_CONTENTS_INLINE);
        if (atoms)
        {
            setViewport(cmd, {kShadowSize, kShadowSize});
            const VkDescriptorSet sets[2] = {m_sets[slot], scene->drawSet};
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuLayout, 0, 2, sets, 1, &frameOffset);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuShadowAtom);
            const VkDeviceSize zero = 0;
            vkCmdBindVertexBuffers(cmd, 0, 1, &m_lodVertices.buffer, &zero);
            vkCmdBindIndexBuffer(cmd, m_lodIndices.buffer, 0, VK_INDEX_TYPE_UINT32);
            const glm::uvec4 lists[3] = {{0, 1, 0, 0}, {n - 1, 0, 0, 0}, {n, 1, 0, 0}};
            const float shadowScale = pixelScaleOf(view.lightMVP, kShadowSize);
            for (uint32_t level = 0; level < 3; ++level)
            {
                if (level == 2) vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuShadowPoint);
                const GpuPush push{lists[level], glm::vec4(shadowScale, 0.0f, 0.0f, 0.0f)};
                vkCmdPushConstants(cmd, m_gpuLayout, kPushStages, 0, sizeof(push), &push);
                vkCmdDrawIndexedIndirect(cmd, m_commands[slot].buffer, (kShadowAtoms + level) * sizeof(DrawCommand), 1, sizeof(DrawCommand));
            }
            if (bonds)
            {
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuShadowBond);
                vkCmdBindVertexBuffers(cmd, 0, 1, &m_cylinderLodVertices.buffer, &zero);
                vkCmdBindIndexBuffer(cmd, m_cylinderLodIndices.buffer, 0, VK_INDEX_TYPE_UINT32);
                const GpuPush push{glm::uvec4(3 * n, 1, 0, 0), glm::vec4(0.0f)};
                vkCmdPushConstants(cmd, m_gpuLayout, kPushStages, 0, sizeof(push), &push);
                vkCmdDrawIndexedIndirect(cmd, m_commands[slot].buffer, kShadowBonds * sizeof(DrawCommand), 1, sizeof(DrawCommand));
            }
        }
        vkCmdEndRenderPass(cmd);
        m_shadowKey = key;
        m_shadowValid = true;
    }
    stamp(cmd, 2);
    if (!atoms) return;

    // The atoms drawn last frame (all visible ones without occlusion culling);
    // the shadow pass's lists are reused.
    barrier(cmd, kDrawReads, kDrawReadAccess, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, kComputeAccess);
    CullSpec main;
    main.viewProjection = view.projection * view.view;
    main.view = view.view;
    main.orthographic = view.projection[3][3] == 1.0f;
    main.target = target;
    main.phase = m_occlusion ? 0u : 2u;
    main.levels = kAtomLevels;
    const float bias = buffers.renderMode == RenderingMode::LowPolyInstancing ? 2.5f : 1.0f;
    for (uint32_t level = 0; level + 1 < kAtomLevels; ++level) main.thresholds[level] = kAtomThresholds[level] * bias;
    for (uint32_t level = 0; level < kAtomLevels; ++level)
        main.inflation[level] = level >= kFirstImpostorLevel ? 1.0f : m_atomInflation[level];  // impostors are exact
    // Without large points, the farthest level is one pixel: only for atoms under it.
    if (context().maxPointSize() <= 1.0f) main.thresholds[kAtomLevels - 2] = 0.5f;
    for (uint32_t level = 0; level < kAtomLevels; ++level)
    {
        const uint32_t region = level / 2;
        const bool forward = level % 2 == 0;
        main.slots[level] = {kFirstAtoms + level, forward ? region * n : (region + 1) * n - 1, forward ? 1u : 0u, 0u};
    }
    dispatchCull(cmd, *scene, main);
    if (!m_occlusion && bonds)
    {
        CullSpec bondSpec = main;
        bondSpec.kind = 1;
        bondSpec.phase = 2;
        bondSpec.levels = 3;
        bondSpec.cullRadius = kBondCullRadius;
        for (float& f : bondSpec.inflation) f = 1.0f;
        bondSpec.thresholds[0] = kBondThresholds[0];
        bondSpec.thresholds[1] = kBondThresholds[1];
        bondSpec.slots[0] = {kBonds, 3 * n, 1, 0};
        bondSpec.slots[1] = {kBonds + 1, 3 * n + m - 1, 0, 0};
        bondSpec.slots[2] = {kBonds + 2, 3 * n + m, 1, 0};
        dispatchCull(cmd, *scene, bondSpec);
    }
    if (!m_occlusion)
    {
        // The statistics: what was counted (no second phase).
        barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_READ_BIT);
        const VkBufferCopy copy{0, 0, sizeof(commands)};
        vkCmdCopyBuffer(cmd, m_commands[slot].buffer, m_statsReadback[slot].buffer, 1, &copy);
        m_statsPending[slot] = true;
    }
    stamp(cmd, 3);
    barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT, kDrawReads, kDrawReadAccess);
    m_mainCull = main;
}

void SceneRenderer::drawAtomLevels(VkCommandBuffer cmd, const SceneResources& scene, uint32_t firstCommand, bool billboards, float pixelScale,
                                   AtomPass pass)
{
    const VkPipeline mesh = pass == AtomPass::Depth ? m_gpuAtomDepth : pass == AtomPass::Shade ? m_gpuAtomShade : m_gpuAtom;
    const VkPipeline impostor = pass == AtomPass::Depth ? m_gpuBillboardDepth : pass == AtomPass::Shade ? m_gpuBillboardShade : m_gpuBillboard;
    const VkPipeline point = pass == AtomPass::Depth ? m_gpuPointDepth : pass == AtomPass::Shade ? m_gpuPointShade : m_gpuPoint;
    const uint32_t slot = context().frameIndex();
    const uint32_t n = scene.atomCount;
    const VkDeviceSize zero = 0;
    vkCmdBindVertexBuffers(cmd, 0, 1, &m_lodVertices.buffer, &zero);
    vkCmdBindIndexBuffer(cmd, m_lodIndices.buffer, 0, VK_INDEX_TYPE_UINT32);
    for (uint32_t level = 0; level < kAtomLevels; ++level)
    {
        if (level == 0) vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, billboards ? impostor : mesh);
        if (level == kFirstImpostorLevel && !billboards) vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, impostor);
        if (level + 1 == kAtomLevels) vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, point);
        const uint32_t region = level / 2;
        const bool forward = level % 2 == 0;
        const GpuPush push{glm::uvec4(forward ? region * n : (region + 1) * n - 1, forward ? 1u : 0u, 0u, 0u),
                           glm::vec4(pixelScale, level >= kFastShadowLevel ? 1.0f : 0.0f, 0, 0)};
        vkCmdPushConstants(cmd, m_gpuLayout, kPushStages, 0, sizeof(push), &push);
        vkCmdDrawIndexedIndirect(cmd, m_commands[slot].buffer, (firstCommand + level) * sizeof(DrawCommand), 1, sizeof(DrawCommand));
    }
}

void SceneRenderer::drawGpuFirst(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags, VkExtent2D target)
{
    const SceneResources* scene = buffers.vulkan.get();
    if (!scene || !flags.atoms || !scene->atomCount || !scene->drawSet) return;
    setViewport(cmd, target);
    const VkDescriptorSet sets[2] = {m_sets[context().frameIndex()], scene->drawSet};
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuLayout, 0, 2, sets, 1, &frameOffset);
    const bool billboards = buffers.renderMode == RenderingMode::BillboardImposters;
    const float pixelScale = pixelScaleOf(m_mainCull.viewProjection, target.height);
    drawAtomLevels(cmd, *scene, kFirstAtoms, billboards, pixelScale, AtomPass::Depth);
    drawAtomLevels(cmd, *scene, kFirstAtoms, billboards, pixelScale, AtomPass::Shade);
    stamp(cmd, 4);
}

void SceneRenderer::cullGpuSecond(VkCommandBuffer cmd, const SceneBuffers& buffers, const SceneFlags& flags, VkExtent2D target)
{
    const SceneResources* scene = buffers.vulkan.get();
    if (!m_occlusion || !scene || !flags.atoms || !scene->atomCount || !scene->cullSet || !m_hiz.image) return;
    const uint32_t slot = context().frameIndex();
    const uint32_t n = scene->atomCount, m = scene->bondCount;

    // The Hi-Z pyramid of the atoms drawn so far: the farthest depth per texel.
    VkImageMemoryBarrier pyramid{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    pyramid.srcAccessMask = VK_ACCESS_SHADER_READ_BIT;
    pyramid.dstAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    pyramid.oldLayout = pyramid.newLayout = VK_IMAGE_LAYOUT_GENERAL;
    pyramid.srcQueueFamilyIndex = pyramid.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    pyramid.image = m_hiz.image;
    pyramid.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, static_cast<uint32_t>(m_hiz.levels.size()), 0, 1};
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, nullptr, 0, nullptr, 1, &pyramid);
    const VkExtent2D depthSize = context().extent();
    for (uint32_t level = 0; level < m_hiz.levels.size(); ++level)
    {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, level == 0 ? m_hizFirstPipeline : m_hizDownPipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_hizLayout, 0, 1, &m_hiz.sets[level], 0, nullptr);
        const glm::ivec4 push(static_cast<int>(depthSize.width), static_cast<int>(depthSize.height), static_cast<int>(context().samples()), 0);
        vkCmdPushConstants(cmd, m_hizLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push), &push);
        const uint32_t w = std::max(1u, m_hiz.width >> level), h = std::max(1u, m_hiz.height >> level);
        vkCmdDispatch(cmd, (w + 7) / 8, (h + 7) / 8, 1);
        barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_READ_BIT);
    }
    stamp(cmd, 5);
    // The first phase's draws read the lists that are now rewritten.
    barrier(cmd, kDrawReads, kDrawReadAccess, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, kComputeAccess);

    // The other atoms in view, unless hidden; then the bonds.
    CullSpec second = m_mainCull;
    second.phase = 1;
    second.occlusion = true;
    for (uint32_t level = 0; level < kAtomLevels; ++level) second.slots[level].x = kSecondAtoms + level;
    second.target = target;
    dispatchCull(cmd, *scene, second);
    if (flags.bonds && m)
    {
        CullSpec bondSpec = second;
        bondSpec.kind = 1;
        bondSpec.phase = 2;
        bondSpec.levels = 3;
        bondSpec.cullRadius = kBondCullRadius;
        for (float& f : bondSpec.inflation) f = 1.0f;
        bondSpec.thresholds[0] = kBondThresholds[0];
        bondSpec.thresholds[1] = kBondThresholds[1];
        bondSpec.slots[0] = {kBonds, 3 * n, 1, 0};
        bondSpec.slots[1] = {kBonds + 1, 3 * n + m - 1, 0, 0};
        bondSpec.slots[2] = {kBonds + 2, 3 * n + m, 1, 0};
        dispatchCull(cmd, *scene, bondSpec);
    }
    stamp(cmd, 6);
    barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT, kDrawReads | VK_PIPELINE_STAGE_TRANSFER_BIT,
            kDrawReadAccess | VK_ACCESS_TRANSFER_READ_BIT);
    const VkBufferCopy copy{0, 0, kCommandCount * sizeof(DrawCommand)};
    vkCmdCopyBuffer(cmd, m_commands[slot].buffer, m_statsReadback[slot].buffer, 1, &copy);
    m_statsPending[slot] = true;
}

void SceneRenderer::drawGpuSecond(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags,
                                  VkExtent2D target, const glm::vec3& boxColor)
{
    const SceneResources* scene = buffers.vulkan.get();
    if (!scene) return;
    const uint32_t slot = context().frameIndex();
    if (flags.atoms && scene->atomCount && scene->drawSet)
    {
        setViewport(cmd, target);
        const VkDescriptorSet sets[2] = {m_sets[slot], scene->drawSet};
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuLayout, 0, 2, sets, 1, &frameOffset);
        if (m_occlusion)
            drawAtomLevels(cmd, *scene, kSecondAtoms, buffers.renderMode == RenderingMode::BillboardImposters,
                           pixelScaleOf(m_mainCull.viewProjection, target.height), AtomPass::Full);
        if (flags.bonds && scene->bondCount)
        {
            const uint32_t n = scene->atomCount, m = scene->bondCount;
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, m_gpuBond);
            const VkDeviceSize zero = 0;
            vkCmdBindVertexBuffers(cmd, 0, 1, &m_cylinderLodVertices.buffer, &zero);
            vkCmdBindIndexBuffer(cmd, m_cylinderLodIndices.buffer, 0, VK_INDEX_TYPE_UINT32);
            const glm::uvec4 lists[3] = {{3 * n, 1, 0, 0}, {3 * n + m - 1, 0, 0, 0}, {3 * n + m, 1, 0, 0}};
            for (uint32_t level = 0; level < 3; ++level)
            {
                const GpuPush push{lists[level], glm::vec4(0.0f)};
                vkCmdPushConstants(cmd, m_gpuLayout, kPushStages, 0, sizeof(push), &push);
                vkCmdDrawIndexedIndirect(cmd, m_commands[slot].buffer, (kBonds + level) * sizeof(DrawCommand), 1, sizeof(DrawCommand));
            }
        }
    }
    drawOverlays(cmd, buffers, frameOffset, flags, target, boxColor);
    stamp(cmd, 7);
}
}
