#pragma once

#include "algorithms/NanoCrystalBuilder.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct WulffPreviewBatch
{
    glm::vec3 color = glm::vec3(1.0f);
    glm::vec3 normal = glm::vec3(0.0f, 0.0f, 1.0f);
    std::vector<glm::vec3> triangleVertices;
    std::vector<glm::vec3> vertices;
};

struct NanoCrystalBuilderDialog
{
    NanoCrystalBuilderDialog();
    ~NanoCrystalBuilderDialog();

    // Call once after the GL context and renderer are initialised.
    // Creates the preview sphere/cylinder meshes, SceneBuffers, shadow map,
    // and FBO owned by this dialog.
    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);

    // elementRadii / elementShininess are used for the 3-D preview and for
    // assigning colours when the nanocrystal is built.
    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }

    // Called by the drop handler to redirect a file drop into this dialog's
    // reference loader when the dialog is open.
    void feedDroppedFile(const std::string& path);

private:
    bool m_openRequested = false;
    bool m_isOpen        = false;

    // Reference structure – the source for tiling and shape carving.
    Structure   m_reference;

    // Path queued by feedDroppedFile(); consumed at the start of drawDialog().
    std::string m_pendingDropPath;

    // -----------------------------------------------------------------------
    // Embedded file-browser state
    // -----------------------------------------------------------------------
    std::string m_browsDir;
    std::vector<std::pair<std::string, bool>> m_browsEntries; // (name, isDir)
    bool  m_browsEntryDirty  = true;
    bool  m_browsShowPanel   = false;
    char  m_browsFilename[1024];
    char  m_browsStatusMsg[256];

    // -----------------------------------------------------------------------
    // 3-D preview GL resources  (created in initRenderResources)
    // -----------------------------------------------------------------------
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;  // owned
    CylinderMesh* m_previewCylinder = nullptr;  // owned
    SceneBuffers  m_previewBuffers;
    ShadowMap     m_previewShadow   = {};

    GLuint m_previewFBO      = 0;
    GLuint m_previewColorTex = 0;
    GLuint m_previewDepthRbo = 0;
    int    m_previewW        = 0;
    int    m_previewH        = 0;

    GLuint m_wulffPreviewFBO      = 0;
    GLuint m_wulffPreviewColorTex = 0;
    GLuint m_wulffPreviewDepthRbo = 0;
    int    m_wulffPreviewW        = 0;
    int    m_wulffPreviewH        = 0;

    GLuint m_wulffLineVAO = 0;
    GLuint m_wulffLineVBO = 0;
    GLuint m_wulffPlaneProgram = 0;
    GLuint m_wulffPlaneVAO = 0;
    GLuint m_wulffPlaneVBO = 0;

    bool m_glReady           = false;
    bool m_previewBufDirty   = true; // rebuild SceneBuffers from m_reference
    bool m_wulffPreviewDirty = true;

    WulffPreview m_wulffPreviewData;
    std::vector<WulffPreviewBatch> m_wulffPreviewBatches;
    std::vector<glm::vec3> m_wulffFamilyColors;
    std::size_t m_wulffPreviewSignature = 0;

    // Preview orbit camera
    float m_camYaw      = 45.0f;
    float m_camPitch    = 35.0f;
    float m_camDistance = 10.0f;

    float m_wulffCamYaw      = 35.0f;
    float m_wulffCamPitch    = 30.0f;
    float m_wulffCamDistance = 18.0f;
    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------
    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers(const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void renderPreviewToFBO(int w, int h);
    void ensureWulffPreviewFBO(int w, int h);
    void ensureWulffFamilyColors(std::size_t familyCount);
    void rebuildWulffPreviewGeometry(const WulffPreview& preview);
    void renderWulffPreviewToFBO(int w, int h);
    void autoFitPreviewCamera();
    void autoFitWulffPreviewCamera();
    void refreshBrowserEntries();
    bool tryLoadFile(const std::string& path,
                     const std::vector<float>& radii,
                     const std::vector<float>& shininess);
};
