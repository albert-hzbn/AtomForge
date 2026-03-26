#pragma once

#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <functional>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct CSLGrainBoundaryDialog
{
    CSLGrainBoundaryDialog();
    ~CSLGrainBoundaryDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);
    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }
    void feedDroppedFile(const std::string& path);

private:
    bool m_openRequested = false;
    bool m_isOpen        = false;
    std::string m_pendingDropPath;

    // 3-D preview GL resources
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;
    CylinderMesh* m_previewCylinder = nullptr;
    SceneBuffers  m_previewBuffers;
    ShadowMap     m_previewShadow   = {};

    GLuint m_previewFBO      = 0;
    GLuint m_previewColorTex = 0;
    GLuint m_previewDepthRbo = 0;
    int    m_previewW        = 0;
    int    m_previewH        = 0;

    bool  m_glReady         = false;
    bool  m_previewBufDirty = true;

    float m_camYaw      = 45.0f;
    float m_camPitch    = 35.0f;
    float m_camDistance  = 10.0f;

    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers(const Structure& s,
                               const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void renderPreviewToFBO(int w, int h);
    void autoFitPreviewCamera();
};
