#pragma once

#include "app/SceneView.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"

#include <functional>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct ImDrawList;
struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct MergeStructuresDialog
{
    MergeStructuresDialog();
    ~MergeStructuresDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);

    void drawDialog(Structure& structure,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }
    void feedDroppedFile(const std::string& path);

private:
    // ------------------------------------------------------------------
    // Per-entry data
    // ------------------------------------------------------------------
    struct MergeEntry
    {
        Structure structure;
        std::string name;
        glm::vec3 translation = glm::vec3(0.0f);
        glm::vec3 rotationDeg = glm::vec3(0.0f);
        glm::vec3 pivot       = glm::vec3(0.0f); // centroid in original space
        bool      enabled     = true;
    };

    // ------------------------------------------------------------------
    // Gizmo interaction
    // ------------------------------------------------------------------
    enum class GizmoMode { Translate, Rotate };
    enum class GizmoAxis { None = -1, X = 0, Y = 1, Z = 2 };

    struct GizmoDragState
    {
        bool      active   = false;
        GizmoMode mode     = GizmoMode::Translate;
        GizmoAxis axis     = GizmoAxis::None;
        // Reference values at drag start
        glm::vec3 startTranslation = glm::vec3(0.0f);
        glm::vec3 startRotationDeg = glm::vec3(0.0f);
        // Mouse position at drag start (viewport-local)
        float     startMouseX = 0.0f;
        float     startMouseY = 0.0f;
    };

    // ------------------------------------------------------------------
    // Dialog state
    // ------------------------------------------------------------------
    bool m_openRequested = false;
    bool m_isOpen        = false;
    std::vector<std::string> m_pendingDropPaths;
    std::vector<MergeEntry>  m_entries;
    std::string              m_status;
    int m_selectedIndex = -1;

    // ------------------------------------------------------------------
    // GL resources
    // ------------------------------------------------------------------
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;
    CylinderMesh* m_previewCylinder = nullptr;
    SceneBuffers  m_previewBuffers;
    ShadowMap     m_previewShadow   = {};
    bool m_glReady      = false;
    bool m_previewDirty = true;
    bool m_autoFitOnRebuild = true;
    bool m_showPreviewBoundingBox = false;

    GLuint m_previewFBO      = 0;
    GLuint m_previewColorTex = 0;
    GLuint m_previewDepthRbo = 0;
    int m_previewW = 0;
    int m_previewH = 0;

    // ------------------------------------------------------------------
    // Preview camera
    // ------------------------------------------------------------------
    float m_camYaw      = 45.0f;
    float m_camPitch    = 35.0f;
    float m_camDistance = 12.0f;

    // Last rendered frame matrices (used for gizmo projection & ray picking)
    FrameView m_lastFrame;

    // ------------------------------------------------------------------
    // Gizmo drag state
    // ------------------------------------------------------------------
    GizmoDragState m_drag;

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------
    bool addStructureFromPath(const std::string& path);
    void addStructureFromScene(const Structure& structure);
    static glm::vec3 computeCentroid(const Structure& structure);
    static std::string baseName(const std::string& path);

    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers();
    void autoFitPreviewCamera();
    void renderPreviewToFBO(int w, int h);
    Structure buildCombinedPreviewStructure() const;

    // Gizmo helpers
    glm::vec3 getTransformedCentroid(int idx) const;
    glm::vec2 projectToViewport(const glm::vec3& worldPos, float vpX, float vpY, float vpW, float vpH) const;
    void drawGizmoOverlay(ImDrawList* dl, float vpX, float vpY, float vpW, float vpH,
                          float gizmoScale, GizmoMode mode);
    // Returns the handle hit by a viewport-local point, or {None,-1}
    std::pair<GizmoMode, GizmoAxis> hitTestGizmo(float lx, float ly,
                                                  float vpX, float vpY, float vpW, float vpH,
                                                  float gizmoScale) const;
    // Click-to-select: pick the closest entry atom under cursor
    int pickEntryAtCursor(float lx, float ly, float vpW, float vpH) const;
};
