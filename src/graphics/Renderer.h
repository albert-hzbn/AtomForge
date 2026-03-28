#pragma once

#include "ShadowMap.h"
#include "CylinderMesh.h"
#include "SphereMesh.h"
#include "LowPolyMesh.h"
#include "BillboardMesh.h"
#include "SceneBuffers.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <cstddef>

// Owns the GLSL programs and exposes per-pass draw methods.
struct Renderer
{
    GLuint atomProgram   = 0;
    GLuint atomLowPolyProgram = 0;
    GLuint atomBillboardProgram = 0;
    GLuint bondProgram   = 0;
    GLuint shadowProgram = 0;
    GLuint shadowLowPolyProgram = 0;
    GLuint shadowBillboardProgram = 0;
    GLuint bondShadowProgram = 0;
    GLuint lineProgram   = 0;

    // Compile and link all shader programs.
    void init();

    // Select the appropriate rendering mode based on atom count
    RenderingMode selectRenderingMode(size_t atomCount) const;

    // Render all atoms into the shadow map (depth-only pass).
    void drawShadowPass(const ShadowMap& shadow,
                        const SphereMesh& sphere,
                        const glm::mat4& lightMVP,
                        size_t atomCount);

    void drawShadowPassLowPoly(const ShadowMap& shadow,
                               const LowPolyMesh& mesh,
                               const glm::mat4& lightMVP,
                               size_t atomCount);

    void drawShadowPassBillboard(const ShadowMap& shadow,
                                 const BillboardMesh& mesh,
                                 const glm::mat4& lightMVP,
                                 const glm::mat4& view,
                                 size_t atomCount);

    void drawBondShadowPass(const ShadowMap& shadow,
                            const CylinderMesh& cylinder,
                            const glm::mat4& lightMVP,
                            size_t bondCount);

    // Render all atoms into the colour buffer with shadow sampling.
    void drawAtoms(const glm::mat4& projection,
                   const glm::mat4& view,
                   const glm::mat4& lightMVP,
                   const glm::vec3& lightPos,
                   const glm::vec3& viewPos,
                   const ShadowMap& shadow,
                   const SphereMesh& sphere,
                   size_t atomCount);

    void drawAtomsLowPoly(const glm::mat4& projection,
                          const glm::mat4& view,
                          const glm::mat4& lightMVP,
                          const glm::vec3& lightPos,
                          const glm::vec3& viewPos,
                          const ShadowMap& shadow,
                          const LowPolyMesh& mesh,
                          size_t atomCount);

    void drawAtomsBillboard(const glm::mat4& projection,
                            const glm::mat4& view,
                            const glm::mat4& lightMVP,
                            const glm::vec3& lightPos,
                            const glm::vec3& viewPos,
                            const ShadowMap& shadow,
                            const BillboardMesh& mesh,
                            size_t atomCount);

    void drawBonds(const glm::mat4& projection,
                   const glm::mat4& view,
                   const glm::vec3& lightPos,
                   const glm::vec3& viewPos,
                   const CylinderMesh& cylinder,
                   size_t bondCount);

    // Render the bounding-box / lattice wireframe.
    void drawBoxLines(const glm::mat4& projection,
                      const glm::mat4& view,
                      GLuint lineVAO,
                      size_t lineVertexCount);
};
