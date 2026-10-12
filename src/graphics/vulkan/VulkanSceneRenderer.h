#pragma once
// The atomic scene in Vulkan, as the OpenGL Renderer draws it: instanced
// atoms (spheres, low-poly spheres or billboards), bonds, the cell box, the
// dislocation loop and selection halos, with the key light's shadow map; and
// the lines and lit triangles the dialog previews add to it.

#include "graphics/vulkan/VulkanContext.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

struct SceneBuffers;
class SphereMesh;
class LowPolyMesh;
class BillboardMesh;
class CylinderMesh;

namespace vkr
{
// The GPU copies of one SceneBuffers (held by it while Vulkan draws).
struct SceneResources
{
    Buffer positions, colors, scales, shininess;  // per atom (host-visible: edits patch them)
    Buffer bondStart, bondEnd, bondColorA, bondColorB, bondRadius, bondShininessA, bondShininessB;
    Buffer boxLines, loop;
    uint32_t atomCount = 0, bondCount = 0, boxVertexCount = 0, loopVertexCount = 0;
    ~SceneResources();
};

struct SceneLighting
{
    float ambient = 0.40f, saturation = 1.55f, contrast = 1.25f, shadowStrength = 0.20f;
    float specularIntensity = 0.65f, shininessScale = 1.5f, shininessFloor = 32.0f;
};

struct SceneView
{
    glm::mat4 projection{1.0f};  // OpenGL convention (the renderer corrects it for Vulkan)
    glm::mat4 view{1.0f};
    glm::mat4 lightMVP{1.0f};
    glm::vec3 lightPosition{0.0f};
    glm::vec3 cameraPosition{0.0f};
};

struct SceneFlags
{
    bool atoms = true, bonds = true, box = true, dislocation = false, lightTheme = true;
    bool spheresOnly = false;  // previews: full spheres whatever the buffers' rendering mode
};

struct MeshVertex
{
    glm::vec3 position;
    glm::vec3 normal;
};

enum class LineTopology { List, Strip };

// The shading of lit triangles (as the OpenGL dialog shaders).
enum class MeshStyle
{
    Lit = 0,           // two-sided, lit from a direction (custom shapes)
    WulffPlane = 1,    // facets of the nanocrystal builder
    CuttingPlane = 2,  // the cell sculptor's planes
};

class SceneRenderer
{
public:
    bool init(const SphereMesh& sphere, const LowPolyMesh& lowPoly, const BillboardMesh& billboard, const CylinderMesh& cylinder);
    void shutdown();
    bool ready() const { return m_ready; }

    // Scene data (SceneBuffers calls these when Vulkan draws).
    void upload(SceneBuffers& buffers, const std::vector<glm::vec4>& positions, const std::vector<glm::vec4>& colors,
                const std::vector<float>& scales, const std::vector<float>& shininess,
                const std::vector<glm::vec3>& bondStart, const std::vector<glm::vec3>& bondEnd,
                const std::vector<glm::vec3>& bondColorA, const std::vector<glm::vec3>& bondColorB,
                const std::vector<float>& bondRadius, const std::vector<float>& bondShininessA,
                const std::vector<float>& bondShininessB, const std::vector<glm::vec3>& boxLines);
    void setInstanceColor(SceneBuffers& buffers, int index, const glm::vec3& color);
    void setInstancePosition(SceneBuffers& buffers, int index, const glm::vec3& position);
    void uploadLoop(SceneBuffers& buffers, const std::vector<glm::vec3>& points);

    // A frame. prepare() returns the uniform offset the draws of one view share;
    // the shadow map is rendered outside a render pass, the rest inside the
    // current pass (the window's, or a canvas's).
    void newFrame();
    uint32_t prepare(const SceneView& view, const SceneLighting& lighting);
    void renderShadowMap(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags);
    void drawScene(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags,
                   VkExtent2D target, const glm::vec3& boxColor);
    void drawHalos(VkCommandBuffer cmd, uint32_t frameOffset, const std::vector<glm::vec3>& positions,
                   const std::vector<float>& radii, VkExtent2D target);
    // Atoms blended over what is drawn at `opacity`, without writing depth.
    void drawGhostAtoms(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, float opacity, VkExtent2D target);
    void drawLines(VkCommandBuffer cmd, uint32_t frameOffset, const glm::vec3* points, std::size_t count, LineTopology topology,
                   const glm::vec3& color, bool depthTest, VkExtent2D target, float width = 3.5f);
    // Triangles (three vertices each); translucent ones blend by the colour's
    // alpha, keep the depth buffer and sit slightly behind coincident surfaces.
    void drawMesh(VkCommandBuffer cmd, uint32_t frameOffset, const MeshVertex* vertices, std::size_t count,
                  const glm::vec4& color, MeshStyle style, bool translucent, const glm::vec3& lightDirection, VkExtent2D target);

    // A render pass compatible with the scene pipelines that renders into a
    // sampled image (canvases: previews, image export).
    VkRenderPass offscreenPass() const { return m_offscreenPass; }
    VkSampler linearSampler() const { return m_linearSampler; }
    // Room in this frame's uniform ring for renderers with descriptor sets of
    // their own (bound as a dynamic uniform buffer on uniformRing(frame)).
    uint32_t allocateUniform(const void* data, uint32_t bytes);
    VkBuffer uniformRing(uint32_t frame) const { return m_uniforms[frame].buffer; }

private:
    struct PipelineSpec;
    VkPipeline makePipeline(const PipelineSpec& spec);
    void bindFrame(VkCommandBuffer cmd, uint32_t frameOffset);
    void setViewport(VkCommandBuffer cmd, VkExtent2D extent);
    void setLineWidth(VkCommandBuffer cmd, float width = 3.5f);
    // Binds the atom geometry and instance streams; returns the index count.
    uint32_t bindAtoms(VkCommandBuffer cmd, const SceneResources& scene, bool billboards, bool lowPoly);
    void bindBonds(VkCommandBuffer cmd, const SceneResources& scene);
    // Copies transient vertices into this frame's stream (host-visible);
    // returns their offset in `buffer`.
    VkDeviceSize stream(const void* data, VkDeviceSize bytes, VkBuffer& buffer);

    bool m_ready = false;
    VkDescriptorSetLayout m_setLayout = VK_NULL_HANDLE;
    VkPipelineLayout m_layout = VK_NULL_HANDLE;
    VkDescriptorSet m_sets[Context::kFramesInFlight]{};
    Buffer m_uniforms[Context::kFramesInFlight];
    uint32_t m_uniformCursor = 0;
    uint32_t m_uniformAlign = 256;
    static constexpr uint32_t kUniformRing = 2 * 1024 * 1024;
    Buffer m_streams[Context::kFramesInFlight];
    VkDeviceSize m_streamCursor = 0;

    Buffer m_sphereVertices, m_sphereIndices, m_lowPolyVertices, m_lowPolyIndices;
    Buffer m_billboardVertices, m_billboardIndices, m_cylinderVertices, m_haloCorners;
    uint32_t m_sphereIndexCount = 0, m_lowPolyIndexCount = 0, m_billboardIndexCount = 0, m_cylinderVertexCount = 0;

    Image m_shadowMap;
    VkRenderPass m_shadowPass = VK_NULL_HANDLE;
    VkFramebuffer m_shadowFramebuffer = VK_NULL_HANDLE;
    VkSampler m_shadowSampler = VK_NULL_HANDLE;
    static constexpr uint32_t kShadowSize = 2048;

    VkRenderPass m_offscreenPass = VK_NULL_HANDLE;
    VkSampler m_linearSampler = VK_NULL_HANDLE;

    VkPipeline m_atom = VK_NULL_HANDLE, m_billboard = VK_NULL_HANDLE, m_bond = VK_NULL_HANDLE, m_ghostAtom = VK_NULL_HANDLE;
    VkPipeline m_lineList = VK_NULL_HANDLE, m_lineStrip = VK_NULL_HANDLE;
    VkPipeline m_lineListOverlay = VK_NULL_HANDLE, m_lineStripOverlay = VK_NULL_HANDLE;
    VkPipeline m_mesh = VK_NULL_HANDLE, m_meshTranslucent = VK_NULL_HANDLE;
    VkPipeline m_haloHidden = VK_NULL_HANDLE, m_haloVisible = VK_NULL_HANDLE;
    VkPipeline m_shadowAtom = VK_NULL_HANDLE, m_shadowBillboard = VK_NULL_HANDLE, m_shadowBond = VK_NULL_HANDLE;
    std::vector<VkPipeline> m_pipelines;  // all of the above, for shutdown
    std::vector<VkShaderModule> m_modules;
};

SceneRenderer& sceneRenderer();

// The Vulkan clip-space correction: y down, depth in [0, 1].
glm::mat4 clipCorrection(bool flipY);
}
