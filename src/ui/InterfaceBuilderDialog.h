#pragma once

#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <functional>
#include <string>
#include <vector>
#include <array>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

// One match candidate produced by the strain-matching search.
struct InterfaceCandidate
{
    int idx;
    int nMatrix[2][2];
    int mMatrix[2][2];
    int detN;
    int detM;
    float thetaDeg;
    float orTargetDeg;
    float orMisfitDeg;
    float exx;
    float eyy;
    float exy;
    float meanAbsStrain;
    float elasticDensity;
    int   interfaceAtoms;
};

struct InterfaceBuilderDialog
{
    InterfaceBuilderDialog();
    ~InterfaceBuilderDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);

    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }

    // Called by drop handler to forward file drops when dialog is open.
    // slotHint: 0 = auto (fills A then B), 1 = force A, 2 = force B.
    void feedDroppedFile(const std::string& path, int slotHint = 0);

private:
    bool m_openRequested = false;
    bool m_isOpen        = false;

    // -- Two input structures -------------------------------------------------
    Structure m_structureA;
    Structure m_structureB;
    std::string m_pendingDropPathA;
    std::string m_pendingDropPathB;
    std::string m_statusA;
    std::string m_statusB;

    // Embedded file-browser state (shared for both slots)
    std::string m_browsDir;
    std::vector<std::pair<std::string, bool>> m_browsEntries;
    bool m_browsEntryDirty = true;

    // -- Interface Builder parameters (mirrors Python script defaults) --------
    int   m_nmax           = 8;
    int   m_mmax           = 8;
    int   m_maxCellsA      = 20;
    int   m_maxCellsB      = 20;
    int   m_layersA        = 12;
    int   m_layersB        = 12;
    float m_maxMeanStrain  = 0.1f;
    float m_maxRotationDeg = 10.0f;
    float m_orTolDeg       = 2.0f;
    float m_zGap           = 2.5f;
    float m_vacuum         = 15.0f;
    int   m_repeatX        = 2;
    int   m_repeatY        = 2;

    // Orientation relationship (optional)
    bool  m_useOrAngle     = false;
    float m_orAngleDeg     = 0.0f;

    // OR via plane + direction (Miller indices)
    bool  m_useOrPlaneDir  = true;
    float m_planeA[3]      = {1.0f, 1.0f, 0.0f};  // (h k l) for A
    float m_dirA[3]        = {1.0f, 1.0f, 1.0f};  // [u v w] for A
    float m_planeB[3]      = {0.0f, 0.0f, 1.0f};  // (h k l) for B
    float m_dirB[3]        = {1.0f, 1.0f, 0.0f};  // [u v w] for B

    // Stiffness matrix: user picks which structure (A or B) it belongs to
    int   m_stiffnessTarget = 1; // 0 = A, 1 = B
    float m_stiffness[6][6];     // 6x6 Voigt elastic stiffness (GPa)

    // -- Candidate search results ---------------------------------------------
    std::vector<InterfaceCandidate> m_candidates;
    std::vector<Structure>          m_interfaceStructures;
    bool m_searchDone  = false;
    int  m_selectedIdx = -1; // index into m_candidates of the picked candidate
    int  m_lastRenderedResultIdx = -2; // tracks when to rebuild result preview

    // -- 2D scatter plot state ------------------------------------------------
    float m_plotScrollX = 0.0f;
    float m_plotScrollY = 0.0f;
    float m_plotZoom    = 1.0f;
    int   m_hoveredCandidate = -1;

    // -- 3D preview for the bottom region (interface result) ------------------
    // Each preview needs its own sphere/cylinder mesh so that SceneBuffers::init
    // wires instance VBOs into a separate VAO (sharing a VAO breaks instancing).
    Renderer*     m_renderer          = nullptr;
    SphereMesh*   m_previewSphereA    = nullptr;
    CylinderMesh* m_previewCylinderA  = nullptr;
    SphereMesh*   m_previewSphereB    = nullptr;
    CylinderMesh* m_previewCylinderB  = nullptr;
    SphereMesh*   m_previewSphereR    = nullptr;
    CylinderMesh* m_previewCylinderR  = nullptr;
    SceneBuffers  m_previewBufA;
    SceneBuffers  m_previewBufB;
    SceneBuffers  m_previewBufResult;
    ShadowMap     m_previewShadow     = {};

    GLuint m_fboA         = 0, m_texA         = 0, m_rboA         = 0;
    GLuint m_fboB         = 0, m_texB         = 0, m_rboB         = 0;
    GLuint m_fboResult    = 0, m_texResult    = 0, m_rboResult    = 0;
    int    m_fboAW = 0, m_fboAH = 0;
    int    m_fboBW = 0, m_fboBH = 0;
    int    m_fboRW = 0, m_fboRH = 0;
    bool   m_glReady = false;

    // Preview cameras
    float m_camAYaw = 45.f, m_camAPitch = 35.f, m_camADist = 10.f;
    float m_camBYaw = 45.f, m_camBPitch = 35.f, m_camBDist = 10.f;
    float m_camRYaw = 45.f, m_camRPitch = 35.f, m_camRDist = 10.f;

    // -- Private helpers ------------------------------------------------------
    void ensureFBO(GLuint& fbo, GLuint& tex, GLuint& rbo, int& curW, int& curH, int w, int h);
    void rebuildPreviewBuffers(SceneBuffers& buf, const Structure& src,
                               const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void renderToFBO(GLuint fbo, int w, int h,
                     SceneBuffers& buf,
                     float yaw, float pitch, float dist);
    void autoFitCamera(const SceneBuffers& buf, float& dist);
    bool tryLoadFile(const std::string& path, Structure& dest, std::string& status,
                     SceneBuffers& buf,
                     const std::vector<float>& radii,
                     const std::vector<float>& shininess,
                     float& camDist);
    void runSearch();
    void buildInterfaceStructure(int candidateIdx);
    void drawDropZone(const char* label, Structure& struc, const std::string& status,
                      GLuint tex, SceneBuffers& buf,
                      float& yaw, float& pitch, float& dist,
                      const std::vector<float>& radii,
                      const std::vector<float>& shininess,
                      float width, float height, int slot);
    void draw2DPlot(float width, float height);
    void initDemoStiffness();
    void refreshBrowserEntries();
};
