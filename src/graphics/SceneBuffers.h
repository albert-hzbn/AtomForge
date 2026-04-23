#pragma once

#include "StructureInstanceBuilder.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <array>
#include <vector>

// Rendering mode selection based on atom count
enum class RenderingMode
{
    StandardInstancing,     // N > 100,000 - indexed sphere with glDrawElementsInstanced
    LowPolyInstancing,      // 100k < N < 10M - 12-facet icosahedron
    BillboardImposters      // N >= 10M - quad billboards
};

// Owns the GPU buffers for instanced atom rendering and the bounding-box lines.
struct SceneBuffers
{
    GLuint instanceVBO = 0;
    GLuint colorVBO    = 0;
    GLuint scaleVBO    = 0;
    GLuint shininessVBO = 0;
    GLuint bondStartVBO = 0;
    GLuint bondEndVBO = 0;
    GLuint bondColorAVBO = 0;
    GLuint bondColorBVBO = 0;
    GLuint bondRadiusVBO = 0;
    GLuint bondShininessAVBO = 0;
    GLuint bondShininessBVBO = 0;
    GLuint lineVAO     = 0;
    GLuint lineVBO     = 0;

    size_t         atomCount   = 0;
    size_t         bondCount   = 0;
    RenderingMode  renderMode  = RenderingMode::StandardInstancing;  // Selected rendering mode
    glm::vec3      orbitCenter = glm::vec3(0.0f);
    std::vector<glm::vec3> boxLines;

    // CPU-side copies used for ray picking.
    // DISABLED for large structures (>100k atoms) to save memory.
    std::vector<glm::vec3> atomPositions;
    std::vector<glm::vec3> atomColors;   // base colours (no highlight)
    std::vector<float>     atomRadii;
    std::vector<float>     atomShininess;
    std::vector<int>       atomIndices;

    // CPU-side bond caches used for overlay tooling and SVG export.
    std::vector<glm::vec3> bondStarts;
    std::vector<glm::vec3> bondEnds;
    std::vector<glm::vec3> bondColorsA;
    std::vector<glm::vec3> bondColorsB;
    std::vector<float>     bondRadiiCpu;
    
    // Flag: if true, CPU caches are disabled (large structure)
    bool cpuCachesDisabled = false;

    // Per-tab geometry VAOs (owned by this SceneBuffers, NOT the shared mesh VAOs).
    // Each tab builds its own VAOs by pairing shared mesh geometry VBOs with its own
    // instance VBOs so that tabs never clobber each other's VAO attribute state.
    GLuint tabSphereVAO    = 0;
    GLuint tabLowPolyVAO   = 0;
    GLuint tabBillboardVAO = 0;
    GLuint tabCylinderVAO  = 0;

    // Geometry draw counts (cached from the mesh objects at init time).
    int tabSphereIndexCount    = 0;
    int tabLowPolyIndexCount   = 0;
    int tabBillboardIndexCount = 0;
    int tabCylinderVertexCount = 0;

    // Full init: creates per-tab VAOs by combining shared mesh geometry VBOs with
    // this tab's own instance VBOs.  Call once after a valid GL context is ready.
    void init(GLuint sphereVbo,    GLuint sphereEbo,    int sphereIndexCount,
              GLuint lowPolyVbo,   GLuint lowPolyEbo,   int lowPolyIndexCount,
              GLuint billboardVbo, GLuint billboardEbo, int billboardIndexCount,
              GLuint cylinderVbo,  int cylinderVertexCount);

    // Simplified init for preview renderers (standard instancing + bonds only).
    void init(GLuint sphereVbo, GLuint sphereEbo, int sphereIndexCount,
              GLuint cylinderVbo, int cylinderVertexCount);

    // Release all GPU resources owned by this SceneBuffers.
    void destroy();

    // Upload a new StructureInstanceData set to the GPU and cache derived data.
    void upload(const StructureInstanceData& data,
                bool bondElementFilterEnabled,
                const std::array<bool, 119>& bondElementFilterMask);

    // Patch the colour of one instance (e.g. to highlight a selected atom).
    void highlightAtom(int instanceIdx, glm::vec3 color);

    // Restore original (pre-highlight) colour for one instance.
    void restoreAtomColor(int instanceIdx);

    // Update the world-space position of one instance on the GPU.
    void updateAtomPosition(int instanceIdx, const glm::vec3& position);
};
