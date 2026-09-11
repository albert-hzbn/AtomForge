#pragma once

#include "ui/PathPicker.h"

#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct PolyCrystalBuilderDialog
{
    PolyCrystalBuilderDialog();
    ~PolyCrystalBuilderDialog();

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
    PathPicker m_sourcePicker;
    bool m_openRequested = false;
    bool m_isOpen        = false;

    // Reference structure for tiling
    Structure   m_reference;
    std::string m_referenceFilename; // just the basename

    // Path queued by feedDroppedFile(); consumed at start of drawDialog().
    std::string m_pendingDropPath;

    char m_statusMsg[256];

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

    bool tryLoadFile(const std::string& path,
                     const std::vector<float>& radii,
                     const std::vector<float>& shininess);
    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers(const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void renderPreviewToFBO(int w, int h);
    void autoFitPreviewCamera();
};
