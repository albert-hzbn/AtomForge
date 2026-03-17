#pragma once

#include "StructureInstanceBuilder.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <vector>

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
    GLuint lineVAO     = 0;
    GLuint lineVBO     = 0;

    size_t         atomCount   = 0;
    size_t         bondCount   = 0;
    glm::vec3      orbitCenter = glm::vec3(0.0f);
    std::vector<glm::vec3> boxLines;

    // CPU-side copies used for ray picking.
    std::vector<glm::vec3> atomPositions;
    std::vector<glm::vec3> atomColors;   // base colours (no highlight)
    std::vector<float>     atomRadii;
    std::vector<float>     atomShininess;
    std::vector<int>       atomIndices;

    // Allocate GPU objects and wire instance attributes into sphereVAO.
    void init(GLuint sphereVAO, GLuint cylinderVAO);

    // Upload a new StructureInstanceData set to the GPU and cache derived data.
    void upload(const StructureInstanceData& data);

    // Patch the colour of one instance (e.g. to highlight a selected atom).
    void highlightAtom(int instanceIdx, glm::vec3 color);

    // Restore original (pre-highlight) colour for one instance.
    void restoreAtomColor(int instanceIdx);
};
