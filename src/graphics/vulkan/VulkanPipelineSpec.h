#pragma once
// How the scene renderer's graphics pipelines are described (internal to
// VulkanSceneRenderer.cpp and VulkanGpuScene.cpp).

#include "graphics/vulkan/VulkanSceneRenderer.h"

#include <vector>

namespace vkr
{
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
    bool gpu = false;         // GPU-culled: instances from storage buffers (set 1)
    bool colorWrite = true;   // false: a depth pre-pass
};
}
