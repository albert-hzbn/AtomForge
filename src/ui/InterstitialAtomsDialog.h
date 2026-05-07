#pragma once

#include "app/SceneView.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <functional>
#include <string>
#include <vector>
#include <array>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct ImDrawList;
struct ImVec2;
struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct InterstitialAtomsDialog
{
    InterstitialAtomsDialog();
    ~InterstitialAtomsDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);

    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }
    void feedDroppedFile(const std::string& path);

    void drawVoidOverlay(ImDrawList* drawList,
                         const glm::mat4& projection,
                         const glm::mat4& view,
                         int framebufferWidth,
                         int framebufferHeight,
                         const Structure& activeStructure) const;

private:
    enum class VoidType { Tetrahedral = 0, Octahedral, Irregular };

    struct VoidCandidate
    {
        glm::vec3 position   = glm::vec3(0.0f);
        float clearance      = 0.0f;
        int coordination     = 0;
        float volume         = 0.0f;
        int volumeBin        = 0;
        VoidType type        = VoidType::Irregular;
        std::vector<int> neighborAtomIndices;
        std::vector<std::array<int, 3>> faces;
        std::vector<std::array<glm::vec3, 3>> faceVertices;
    };

    struct VoidTypeOption
    {
        std::string label;
        int key   = 0;
        int count = 0;
    };

    enum class PlacementMode { Random = 0, MeshRegion, ManualSelection };
    enum class TargetMode    { Percent = 0, Count };

    // ------------------------------------------------------------------
    // Mesh gizmo interaction
    // ------------------------------------------------------------------
    enum class MeshGizmoMode { Translate, Rotate, Scale };
    enum class MeshGizmoAxis { None = -1, X = 0, Y = 1, Z = 2 };

    struct MeshGizmoDragState
    {
        bool          active           = false;
        MeshGizmoMode mode             = MeshGizmoMode::Translate;
        MeshGizmoAxis axis             = MeshGizmoAxis::None;
        glm::vec3     startTranslation = glm::vec3(0.0f);
        glm::vec3     startRotationDeg = glm::vec3(0.0f);
        float         startScale       = 1.0f;
        float         startMouseX      = 0.0f;
        float         startMouseY      = 0.0f;
    };

    bool detectVoids();
    int  applyPlacement(Structure& structure, const std::vector<glm::vec3>& elementColors);
    void rebuildVoidTypeOptions();
    bool matchesSelectedVoidFilter(const VoidCandidate& v) const;

    bool loadStructureFromPath(const std::string& path,
                               const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    bool loadMeshFromPath(const std::string& path);

    bool pointInsideTransformedMesh(const glm::vec3& point) const;
    bool acceptsStructurePath(const std::string& path) const;
    bool acceptsMeshPath(const std::string& path) const;

    // GL preview
    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers(const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void rebuildMeshSurface();
    void renderMeshSurfaceToPreview();
    void renderPreviewToFBO(int w, int h);
    void autoFitPreviewCamera();
    void drawPreviewOverlays(ImDrawList* dl,
                             const ImVec2& canvasMin,
                             const ImVec2& canvasSize);

    // Mesh gizmo helpers
    glm::vec3 getMeshCenter() const;
    glm::vec2 projectToCanvas(const glm::vec3& wp,
                              const ImVec2& canvasMin,
                              const ImVec2& canvasSize) const;
    void drawMeshGizmo(ImDrawList* dl,
                       const ImVec2& canvasMin,
                       const ImVec2& canvasSize,
                       float gizmoScale);
    std::pair<MeshGizmoMode, MeshGizmoAxis> hitTestMeshGizmo(float mx, float my,
                                                              const ImVec2& canvasMin,
                                                              const ImVec2& canvasSize,
                                                              float gizmoScale) const;

    bool m_openRequested = false;
    bool m_isOpen        = false;

    std::vector<std::string> m_pendingDropPaths;

    Structure   m_source;
    bool        m_sourceLoaded = false;
    std::string m_sourceLabel;

    std::vector<VoidCandidate>  m_voids;
    std::vector<VoidTypeOption> m_voidTypeOptions;
    std::vector<char>           m_voidTypeSelections;
    std::vector<int>            m_selectedVoidIndices;

    int   m_selectedElement   = 6;
    int   m_selectedVoidIndex = -1;
    int   m_gridResolution    = 14;
    int   m_maxVoids          = 400;
    float m_minClearance      = 0.6f;
    float m_minSeparation     = 0.8f;
    bool  m_showVoidOverlay   = true;

    PlacementMode m_placementMode = PlacementMode::Random;
    TargetMode    m_targetMode    = TargetMode::Percent;
    float m_targetPercent     = 25.0f;
    int   m_targetCount       = 20;
    int   m_seed              = 42;
    float m_minInsertDistance = 0.6f;

    std::vector<glm::vec3>    m_meshVertices;
    std::vector<unsigned int> m_meshIndices;
    std::string m_meshLabel;
    bool      m_meshLoaded        = false;
    float     m_meshScale         = 1.0f;
    glm::vec3 m_meshTranslation   = glm::vec3(0.0f);
    glm::vec3 m_meshRotationDeg   = glm::vec3(0.0f);
    GLuint    m_meshProgram       = 0;
    GLuint    m_meshVAO           = 0;
    GLuint    m_meshVBO           = 0;
    int       m_meshTriCount      = 0;

    MeshGizmoDragState m_meshGizmoDrag;
    MeshGizmoMode      m_meshGizmoMode = MeshGizmoMode::Translate;

    std::string m_statusMsg;
    bool        m_statusIsError = false;

    // ---- GL preview resources ----
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;
    CylinderMesh* m_previewCylinder = nullptr;
    SceneBuffers  m_previewBuffers;
    ShadowMap     m_previewShadow   = {};

    GLuint m_previewFBO      = 0;
    GLuint m_previewColorTex = 0;
    GLuint m_previewDepthRbo = 0;
    int    m_previewFboW     = 0;
    int    m_previewFboH     = 0;

    bool m_glReady         = false;
    bool m_previewBufDirty = false;

    // ---- Preview camera ----
    float m_camYaw      = 45.0f;
    float m_camPitch    = 35.0f;
    float m_camDistance = 10.0f;

    // Last rendered frame matrices (used by drawPreviewOverlays)
    FrameView m_lastPreviewFrame;
};
