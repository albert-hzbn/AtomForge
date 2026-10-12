#pragma once
// The atomic scene in Vulkan, as the OpenGL Renderer draws it: instanced
// atoms (spheres, low-poly spheres or billboards), bonds, the cell box, the
// dislocation loop and selection halos, with the key light's shadow map; and
// the lines and lit triangles the dialog previews add to it.

#include "graphics/vulkan/VulkanContext.h"

#include <glm/glm.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <cstdint>
#include <memory>
#include <string>
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
    Buffer positions, colors;  // per atom: xyz and radius, rgb and shininess
    Buffer bondStart, bondEnd, bondColorA, bondColorB;  // xyz and radius, xyz and shininess A, rgb and shininess B, rgb
    Buffer boxLines, loop;
    uint32_t atomCount = 0, bondCount = 0, boxVertexCount = 0, loopVertexCount = 0;
    // GPU culling of the main view: the visible lists, whether each atom was
    // drawn last frame (two-phase occlusion culling) and the descriptor sets.
    Buffer lists, visibility;
    VkDescriptorSet drawSet = VK_NULL_HANDLE, cullSet = VK_NULL_HANDLE;
    // Edits waiting for the next frame (copied in a batch outside render passes).
    struct Patch { bool color; VkDeviceSize offset; uint32_t count; float values[3]; };
    std::vector<Patch> patches;
    uint64_t revision = 0;  // counts edits of positions and radii (the shadow map's cache)
    const uint64_t id = nextId();  // unique (an address can be reused)
    static uint64_t nextId() { static uint64_t next = 0; return ++next; }
    ~SceneResources();
};

// The main view's GPU culling of the last finished frame.
struct SceneStats
{
    uint32_t atomLevels[6] = {};      // atoms drawn per level of detail (finest first)
    uint32_t bondLevels[3] = {};
    uint32_t shadowAtoms = 0, shadowBonds = 0;
    uint32_t firstPhaseAtoms = 0;     // of the atoms drawn, those drawn before occlusion culling
    bool occlusion = false;
    uint32_t atoms = 0, bonds = 0;    // in the scene
    uint64_t triangles = 0;           // triangles drawn by atoms and bonds (points as none)
    double gpuMilliseconds = 0.0;     // scene work on the GPU (culling, shadow, atoms, bonds)
    // Its stages: shadow culling, shadow pass, first-phase culling, first-phase
    // drawing, Hi-Z, second-phase culling (with bonds), second-phase drawing.
    double stageMilliseconds[7] = {};
    bool valid = false;
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
    void setInstanceRadius(SceneBuffers& buffers, int index, float radius);
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
    // The cell box and the dislocation loop (drawScene draws them too).
    void drawOverlays(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags,
                      VkExtent2D target, const glm::vec3& boxColor);

    // The main view, GPU-driven (VulkanGpuScene.cpp): per-atom and per-bond
    // frustum culling, two-phase Hi-Z occlusion culling, levels of detail from
    // the size on screen (meshes from 3200 to 20 triangles, then point
    // sprites) and small-feature culling; the shadow pass culled likewise.
    // Order in a frame: beginGpuScene (before the window's pass, which is then
    // begun pausable), drawGpuFirst, pauseMainPass, cullGpuSecond,
    // resumeMainPass, drawGpuSecond.
    bool gpuDriven() const { return m_gpuReady; }
    void beginGpuScene(VkCommandBuffer cmd, const SceneBuffers& buffers, const SceneView& view, uint32_t frameOffset,
                       const SceneFlags& flags, VkExtent2D target);
    void drawGpuFirst(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags, VkExtent2D target);
    void cullGpuSecond(VkCommandBuffer cmd, const SceneBuffers& buffers, const SceneFlags& flags, VkExtent2D target);
    void drawGpuSecond(VkCommandBuffer cmd, const SceneBuffers& buffers, uint32_t frameOffset, const SceneFlags& flags,
                       VkExtent2D target, const glm::vec3& boxColor);
    bool occlusionCulling() const { return m_occlusion; }
    const SceneStats& stats() const { return m_stats; }
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
    void queuePatch(SceneBuffers& buffers, int index, bool color, uint32_t first, uint32_t count, const float* values);
    void flushPatches(VkCommandBuffer cmd);
    std::vector<std::weak_ptr<SceneResources>> m_patched;

    // GPU-driven main view (VulkanGpuScene.cpp).
    void initGpu();
    void shutdownGpu();
    void createGpuResources(SceneResources& scene);
    void ensureHiz(VkExtent2D extent);
    void releaseHiz();
    // One culling dispatch: what is culled, for which view, into which lists.
    struct CullSpec
    {
        uint32_t kind = 0;               // 0 atoms, 1 bonds
        glm::mat4 viewProjection{1.0f};  // OpenGL convention
        glm::mat4 view{1.0f};            // its view part (the camera, for occlusion culling)
        bool orthographic = false;
        bool flipY = true;               // the window's clip correction (not the shadow map's)
        VkExtent2D target{1, 1};
        uint32_t phase = 2;              // 0 first, 1 second, 2 single
        bool occlusion = false;
        float cullRadius = 0.0f;         // pixels; smaller instances are not drawn
        uint32_t levels = 1;
        float thresholds[8] = {};        // radius on screen from which level i is used
        float inflation[8] = {1, 1, 1, 1, 1, 1, 1, 1};  // level i's mesh reaches out to this many radii
        glm::uvec4 slots[8] = {};        // per level: command, list start, forward
    };
    CullSpec m_mainCull;
    void dispatchCull(VkCommandBuffer cmd, const SceneResources& scene, const CullSpec& spec);
    enum class AtomPass { Full, Depth, Shade };  // depth and shading (two passes: Depth, then Shade)
    void drawAtomLevels(VkCommandBuffer cmd, const SceneResources& scene, uint32_t firstCommand, bool billboards, float pixelScale,
                        AtomPass pass);
    void readStats();
    void stamp(VkCommandBuffer cmd, uint32_t index);
    bool m_gpuReady = false, m_occlusion = false;
    VkDescriptorSetLayout m_drawSetLayout = VK_NULL_HANDLE, m_cullFrameLayout = VK_NULL_HANDLE, m_cullSceneLayout = VK_NULL_HANDLE,
                          m_hizSetLayout = VK_NULL_HANDLE;
    VkPipelineLayout m_gpuLayout = VK_NULL_HANDLE, m_cullLayout = VK_NULL_HANDLE, m_hizLayout = VK_NULL_HANDLE;
    VkPipeline m_cullPipeline = VK_NULL_HANDLE, m_hizFirstPipeline = VK_NULL_HANDLE, m_hizDownPipeline = VK_NULL_HANDLE;
    VkPipeline m_gpuAtom = VK_NULL_HANDLE, m_gpuBillboard = VK_NULL_HANDLE, m_gpuPoint = VK_NULL_HANDLE, m_gpuBond = VK_NULL_HANDLE;
    VkPipeline m_gpuShadowAtom = VK_NULL_HANDLE, m_gpuShadowBond = VK_NULL_HANDLE, m_gpuShadowPoint = VK_NULL_HANDLE;
    VkPipeline m_gpuAtomDepth = VK_NULL_HANDLE, m_gpuBillboardDepth = VK_NULL_HANDLE, m_gpuPointDepth = VK_NULL_HANDLE;
    VkPipeline m_gpuAtomShade = VK_NULL_HANDLE, m_gpuBillboardShade = VK_NULL_HANDLE, m_gpuPointShade = VK_NULL_HANDLE;
    VkDescriptorSet m_cullFrameSets[Context::kFramesInFlight]{};
    Buffer m_commands[Context::kFramesInFlight];       // indirect draws of each frame in flight
    Buffer m_statsReadback[Context::kFramesInFlight];  // their copy for the statistics (host-visible)
    Buffer m_lodVertices, m_lodIndices, m_cylinderLodVertices, m_cylinderLodIndices, m_dummy;
    struct MeshRange { uint32_t firstIndex = 0, indexCount = 0; int32_t vertexOffset = 0; uint32_t triangles = 0; };
    MeshRange m_atomMeshes[6], m_quadMesh, m_bondMeshes[4];
    float m_atomInflation[6] = {1, 1, 1, 1, 1, 1};
    VkSampler m_nearestSampler = VK_NULL_HANDLE;
    struct Hiz
    {
        VkImage image = VK_NULL_HANDLE;
        VkDeviceMemory memory = VK_NULL_HANDLE;
        VkImageView view = VK_NULL_HANDLE;   // every level (culling)
        std::vector<VkImageView> levels;     // one level each (building)
        std::vector<VkDescriptorSet> sets;   // per level: source and destination
        uint32_t width = 0, height = 0;
        uint64_t generation = 0;
        bool fresh = true;                   // not yet in the general layout
    } m_hiz;
    VkQueryPool m_timers = VK_NULL_HANDLE;
    bool m_timerWritten[Context::kFramesInFlight]{};
    bool m_statsPending[Context::kFramesInFlight]{};
    SceneStats m_stats;
    uint32_t m_statsAtoms[Context::kFramesInFlight]{}, m_statsBonds[Context::kFramesInFlight]{};
    bool m_billboardFrame[Context::kFramesInFlight]{};
    // The cached shadow map: what it was drawn for.
    struct ShadowKey
    {
        glm::mat4 lightMVP{0.0f};
        uint64_t scene = 0;  // SceneResources::id
        uint64_t revision = 0;
        bool atoms = false, bonds = false;
        bool operator==(const ShadowKey& o) const
        {
            return lightMVP == o.lightMVP && scene == o.scene && revision == o.revision && atoms == o.atoms && bonds == o.bonds;
        }
    };
    ShadowKey m_shadowKey;
    bool m_shadowValid = false;
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
