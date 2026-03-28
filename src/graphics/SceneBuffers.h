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

    // Allocate GPU objects and wire instance attributes into mesh VAOs.
    void init(GLuint sphereVAO, GLuint lowPolyVAO, GLuint billboardVAO, GLuint cylinderVAO);

    // Simplified init for preview renderers (standard instancing only).
    void init(GLuint sphereVAO, GLuint cylinderVAO);

    // Upload a new StructureInstanceData set to the GPU and cache derived data.
    void upload(const StructureInstanceData& data,
                bool bondElementFilterEnabled,
                const std::array<bool, 119>& bondElementFilterMask);

    // Patch the colour of one instance (e.g. to highlight a selected atom).
    void highlightAtom(int instanceIdx, glm::vec3 color);

    // Restore original (pre-highlight) colour for one instance.
    void restoreAtomColor(int instanceIdx);
};
