#pragma once

#include "algorithms/CellSculptorAlgo.h"
#include "app/SceneView.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"

#include <array>
#include <functional>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct ImDrawList;
struct Renderer;
struct SphereMesh;
struct CylinderMesh;

// ---------------------------------------------------------------------------
// ScFaceBatch – one rendered face polygon (lower or upper face of a slab).
// Mirrors WulffPreviewBatch in NanoCrystalBuilderDialog.
// ---------------------------------------------------------------------------
struct ScFaceBatch
{
    std::vector<glm::vec3> triangleVerts; // fan-triangulated polygon
    std::vector<glm::vec3> edgeVerts;     // line-segment pairs for wireframe
    glm::vec3              normal;
    glm::vec3              color;
};

// ---------------------------------------------------------------------------
// CellSculptorDialog
// ---------------------------------------------------------------------------
// Two-panel dialog:
//   Left  — source structure with lit plane-face overlays (drag & drop)
//   Right — real-time result preview (atoms after slicing)
// Controls live in the right column below the result preview.
// ---------------------------------------------------------------------------
struct CellSculptorDialog
{
    CellSculptorDialog();
    ~CellSculptorDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);

    void drawDialog(Structure& structure,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen()    const { return m_isOpen; }
    void feedDroppedFile(const std::string& path);

private:
    // -----------------------------------------------------------------------
    // Source / supercell / slab state
    // -----------------------------------------------------------------------
    Structure   m_source;
    std::string m_sourceName;
    Structure   m_supercell;     // tiled source — kept in sync with m_sourceBuffers
    bool        m_hasSource = false;

    int m_nx = 1, m_ny = 1, m_nz = 1;

    std::vector<CellSlabPlane> m_slabs;

    // -----------------------------------------------------------------------
    // Dialog UI state
    // -----------------------------------------------------------------------
    bool        m_openRequested  = false;
    bool        m_isOpen         = false;
    std::string m_status;
    std::vector<std::string> m_pendingDropPaths;

    // New-slab input fields
    int   m_newH  = 1;
    int   m_newK  = 0;
    int   m_newL  = 0;
    float m_newD1 = -5.0f;
    float m_newD2 =  5.0f;

    // -----------------------------------------------------------------------
    // Plane face geometry (rebuilt whenever slabs change)
    // -----------------------------------------------------------------------
    std::vector<ScFaceBatch> m_faceBatches;
    bool                     m_facesDirty  = true;

    // -----------------------------------------------------------------------
    // Result
    // -----------------------------------------------------------------------
    bool      m_resultDirty  = true;
    Structure m_previewResult;

    // Per-atom inside/outside mask for applySlabDimming.
    // m_insideMask[i] = true iff m_supercell.atoms[i] passes all slab constraints.
    // Computed in rebuildResult() using the EXACT same loop as cscApplySlabs so
    // the dimming always matches the actual cut — no possibility of divergence.
    std::vector<bool> m_insideMask;

    // Fit-camera flags
    bool m_srcCameraFit = true;
    bool m_resCameraFit = true;

    // -----------------------------------------------------------------------
    // GL resources – source preview (left panel: atoms + plane faces)
    // -----------------------------------------------------------------------
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;
    CylinderMesh* m_previewCylinder = nullptr;
    ShadowMap     m_previewShadow   = {};
    bool          m_glReady         = false;

    // Source atoms FBO
    SceneBuffers m_sourceBuffers;
    SceneBuffers m_ghostBuffers;  // outside atoms only (scale=0 for inside) for transparency pass
    float        m_ghostAlpha     = 0.28f; // 0=invisible, 1=fully opaque (outside atoms)
    bool         m_sourceBufDirty = true;
    GLuint m_sourceFBO      = 0;
    GLuint m_sourceColorTex = 0;
    GLuint m_sourceDepthRbo = 0;
    int    m_sourceW        = 0;
    int    m_sourceH        = 0;

    // Result atoms FBO
    SceneBuffers m_resultBuffers;
    GLuint m_resultFBO      = 0;
    GLuint m_resultColorTex = 0;
    GLuint m_resultDepthRbo = 0;
    int    m_resultW        = 0;
    int    m_resultH        = 0;

    // Plane face GL objects (lit polygon rendering, like Wulff planes in NanoCrystal)
    GLuint m_planeProgram = 0;
    GLuint m_planeVAO     = 0;
    GLuint m_planeVBO     = 0;
    GLuint m_edgeVAO      = 0;
    GLuint m_edgeVBO      = 0;

    // -----------------------------------------------------------------------
    // Preview cameras
    // -----------------------------------------------------------------------
    float     m_srcYaw   = 45.0f;
    float     m_srcPitch = 35.0f;
    float     m_srcDist  = 12.0f;
    FrameView m_srcLastFrame;

    float     m_resYaw   = 45.0f;
    float     m_resPitch = 35.0f;
    float     m_resDist  = 12.0f;
    FrameView m_resLastFrame;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------
    bool loadSourceFromPath(const std::string& path);
    void loadSourceFromScene(const Structure& s, const std::string& name = "scene");

    void rebuildSourceBuffers();
    void rebuildResult();
    // Dims atoms in the source preview that lie outside all current slabs.
    void applySlabDimming();
    void rebuildFaces();

    // Upload a structure's atoms to a SceneBuffers for preview rendering.
    void uploadToPreview(const Structure& s, SceneBuffers& buffers);

    // Allocate (or reallocate) a colour+depth FBO to the given size.
    void ensureFBO(GLuint& fbo, GLuint& colorTex, GLuint& depthRbo,
                   int& storedW, int& storedH, int w, int h);

    // Fit camera distance/angles to the contents of the given buffers.
    void autoFitCamera(float& yaw, float& pitch, float& dist,
                       const SceneBuffers& buffers);

    // Render atoms (and optionally cutting-plane faces) into an FBO.
    void renderToFBO(int w, int h, GLuint fbo,
                     SceneBuffers& buffers,
                     float yaw, float pitch, float dist,
                     FrameView& lastFrame, bool drawPlanes);

    void renderSourceToFBO(int w, int h);
    void renderResultToFBO(int w, int h);
};
