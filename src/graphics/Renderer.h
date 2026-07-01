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
#include <vector>

// Owns the GLSL programs and exposes per-pass draw methods.
struct Renderer
{
    // Lighting parameters — set each frame before draw calls.
    float lightAmbient           = 0.40f;
    float lightSaturation        = 1.55f;
    float lightContrast          = 1.25f;
    float lightShadowStrength    = 0.20f;

    // Material parameters — set each frame before draw calls.
    float materialSpecularIntensity = 0.65f;
    float materialShininessScale    = 1.5f;
    float materialShininessFloor    = 32.0f;

    GLuint atomProgram   = 0;
    GLuint atomLowPolyProgram = 0;
    GLuint atomBillboardProgram = 0;
    GLuint bondProgram   = 0;
    GLuint shadowProgram = 0;
    GLuint shadowLowPolyProgram = 0;
    GLuint shadowBillboardProgram = 0;
    GLuint bondShadowProgram = 0;
    GLuint lineProgram      = 0;
    GLuint selWireProgram   = 0;
    GLuint selWireVAO       = 0;
    GLuint selWireVBO       = 0;
    int    selWireLineVtxCount = 0;

    // GPU-driven indirect rendering (GL 4.3+).
    // One compute cull pass per frame populates SceneBuffers::visibleIndexSSBO and
    // SceneBuffers::drawIndirectBuffer; the three SSBO programs then read per-instance
    // data from the existing VBOs (reused as SSBOs) via gl_InstanceID indirection.
    GLuint cullProgram              = 0;
    GLuint atomSSBOProgram          = 0;
    GLuint atomLowPolySSBOProgram   = 0;
    GLuint atomBillboardSSBOProgram = 0;

    // Depth prepass programs: write only depth (no color output).
    // The legacy program reads per-instance attributes from the VAO;
    // the SSBO program reads from visibleIndexSSBO + instance SSBOs (GL 4.3+).
    GLuint depthPrepassProgram     = 0;
    GLuint depthPrepassSSBOProgram = 0;

    // Shadow SSBO programs: shadow-pass variants that read instance data from SSBOs
    // instead of vertex attributes, enabling indirect shadow-pass draws (GL 4.3+).
    GLuint shadowSSBOProgram          = 0;
    GLuint shadowBillboardSSBOProgram = 0;

    // Compile and link all shader programs.
    void init();

    // Select the appropriate rendering mode based on atom count
    RenderingMode selectRenderingMode(size_t atomCount) const;

    // Render all atoms into the shadow map (depth-only pass).
    void drawShadowPass(const ShadowMap& shadow,
                        GLuint sphereVAO, int sphereIndexCount,
                        const glm::mat4& lightMVP,
                        size_t atomCount);

    void drawShadowPassLowPoly(const ShadowMap& shadow,
                               GLuint lowPolyVAO, int lowPolyIndexCount,
                               const glm::mat4& lightMVP,
                               size_t atomCount);

    void drawShadowPassBillboard(const ShadowMap& shadow,
                                 GLuint billboardVAO, int billboardIndexCount,
                                 const glm::mat4& lightMVP,
                                 const glm::mat4& view,
                                 size_t atomCount);

    void drawBondShadowPass(const ShadowMap& shadow,
                            GLuint cylinderVAO, int cylinderVertexCount,
                            const glm::mat4& lightMVP,
                            size_t bondCount);

    // Render all atoms into the colour buffer with shadow sampling.
    void drawAtoms(const glm::mat4& projection,
                   const glm::mat4& view,
                   const glm::mat4& lightMVP,
                   const glm::vec3& lightPos,
                   const glm::vec3& viewPos,
                   const ShadowMap& shadow,
                   GLuint sphereVAO, int sphereIndexCount,
                   size_t atomCount);

    void drawAtomsLowPoly(const glm::mat4& projection,
                          const glm::mat4& view,
                          const glm::mat4& lightMVP,
                          const glm::vec3& lightPos,
                          const glm::vec3& viewPos,
                          const ShadowMap& shadow,
                          GLuint lowPolyVAO, int lowPolyIndexCount,
                          size_t atomCount);

    void drawAtomsBillboard(const glm::mat4& projection,
                            const glm::mat4& view,
                            const glm::mat4& lightMVP,
                            const glm::vec3& lightPos,
                            const glm::vec3& viewPos,
                            const ShadowMap& shadow,
                            GLuint billboardVAO, int billboardIndexCount,
                            size_t atomCount);

    void drawBonds(const glm::mat4& projection,
                   const glm::mat4& view,
                   const glm::vec3& lightPos,
                   const glm::vec3& viewPos,
                   GLuint cylinderVAO, int cylinderVertexCount,
                   size_t bondCount);

    // Render the bounding-box / lattice wireframe.
    void drawBoxLines(const glm::mat4& projection,
                      const glm::mat4& view,
                      GLuint lineVAO,
                      size_t lineVertexCount,
                      const glm::vec3& color = glm::vec3(0.85f));

    // Render a closed polygon outline with GL_LINE_LOOP (e.g. dislocation loop).
    void drawLineLoop(const glm::mat4& projection,
                      const glm::mat4& view,
                      GLuint lineVAO,
                      size_t vertexCount,
                      const glm::vec3& color);

    // Draw a yellow low-poly wireframe sphere around each selected atom.
    void drawSelectionWireframes(const glm::mat4& projection,
                                 const glm::mat4& view,
                                 const std::vector<glm::vec3>& positions,
                                 const std::vector<float>& radii);

    // GPU-driven frustum cull: populate buf.visibleIndexSSBO with the indices of
    // atoms that pass the view frustum test, and reset buf.drawIndirectBuffer so
    // that instanceCount starts at 0 and is incremented atomically by the compute.
    // indexCount is the element count for the current render-mode mesh.
    void cullAtoms(const SceneBuffers& buf,
                   const glm::mat4& projection,
                   const glm::mat4& view,
                   int indexCount);

    // Light-frustum cull: same compute, but fills buf.shadowVisibleIndexSSBO and
    // buf.shadowDrawIndirectBuffer, enabling indirect shadow-pass draws.
    void cullAtomsForShadow(const SceneBuffers& buf,
                            const glm::mat4& lightMVP,
                            int indexCount);

    // Depth prepass: writes only the depth buffer (color writes disabled by caller).
    // Separate projection + view (not a pre-multiplied MVP) must be passed so the
    // GPU executes the identical fp operations as the colour-pass shaders.
    void drawDepthPrepass(const glm::mat4& projection, const glm::mat4& view,
                          GLuint vao, int indexCount, size_t atomCount);

    // Depth prepass SSBO path: reads from visibleIndexSSBO + per-instance SSBOs.
    void drawDepthPrepassIndirect(const glm::mat4& projection, const glm::mat4& view,
                                  const SceneBuffers& buf, GLuint vao);

    // Indirect shadow draw (Standard + LowPoly share kShadowSSBOVS).
    void drawShadowPassIndirect(const ShadowMap& shadow,
                                const glm::mat4& lightMVP,
                                const SceneBuffers& buf,
                                GLuint vao);

    // Indirect shadow draw for billboard mode.
    void drawShadowPassBillboardIndirect(const ShadowMap& shadow,
                                         const glm::mat4& lightMVP,
                                         const glm::vec3& viewPos,
                                         const SceneBuffers& buf,
                                         GLuint vao);

    // Indirect-draw variants: read per-instance data from SSBOs via visibleIndices.
    void drawAtomsIndirect(const glm::mat4& projection,
                           const glm::mat4& view,
                           const glm::mat4& lightMVP,
                           const glm::vec3& lightPos,
                           const glm::vec3& viewPos,
                           const ShadowMap& shadow,
                           const SceneBuffers& buf);

    void drawAtomsLowPolyIndirect(const glm::mat4& projection,
                                  const glm::mat4& view,
                                  const glm::mat4& lightMVP,
                                  const glm::vec3& lightPos,
                                  const glm::vec3& viewPos,
                                  const ShadowMap& shadow,
                                  const SceneBuffers& buf);

    void drawAtomsBillboardIndirect(const glm::mat4& projection,
                                    const glm::mat4& view,
                                    const glm::mat4& lightMVP,
                                    const glm::vec3& lightPos,
                                    const glm::vec3& viewPos,
                                    const ShadowMap& shadow,
                                    const SceneBuffers& buf);
};
