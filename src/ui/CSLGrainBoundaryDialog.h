#pragma once

#include "ui/PathPicker.h"
#include "ui/StepEditing.h"

#include "algorithms/CSLComputation.h"

#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"

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

    // Opens the dialog on a pipeline step: filled from the step's settings, with
    // "Update step" writing them back instead of building.
    bool editStep(StepEdit edit);
    // The dialog's current settings as the step's options (command-line flags of
    // AtomForge --build gb).
    std::string stepOptions() const;

private:
    PathPicker m_sourcePicker;
    bool m_openRequested = false;
    bool m_isOpen        = false;
    std::string m_pendingDropPath;
    StepEdit m_step;

    // Dialog state
    Structure m_inputStructure;
    std::string m_referencePath;              // file the reference came from ("" = none or the step input)
    char m_statusMsg[256] = "(no structure loaded)";
    char m_loadedFileName[256] = "(none)";
    int m_axis[3] = {0, 0, 1};
    int m_sigmaMax = 200;
    std::vector<SigmaCandidate> m_sigmaCandidates;
    int m_sigmaSelection = 0;
    int m_planeSelection = 0;
    int m_lastAxisForSigma[3] = {0, 0, 0};
    int m_lastSigmaMaxForSigma = 0;
    int m_ucA = 1;
    int m_ucB = 1;
    float m_vacuumPadding = 0.0f;
    float m_gapDist = 0.0f;
    float m_overlapDist = 0.0f;
    bool m_conventionalCell = false;
    CSLBuildResult m_lastResult{};

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
    bool  m_fitPending      = false;   // fit the camera once the preview is rebuilt

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
