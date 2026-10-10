#pragma once

#include "app/SceneView.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"
#include "ui/StepEditing.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <functional>
#include <string>
#include <vector>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

// ---------------------------------------------------------------------------
// SubstitutionalSolidSolutionDialog
//
// Allows the user to:
//   1. Supply a host structure (copied from the current scene or drag-and-dropped).
//      The loaded structure is shown in a 3-D preview on the left.
//   2. Define a target composition (only visible once a structure is loaded)
//      by specifying per-element integer atom counts that sum to the total.
//   3. Build the solid solution, which randomises element assignments while
//      preserving all atom positions and the unit cell.
// ---------------------------------------------------------------------------
struct SubstitutionalSolidSolutionDialog
{
    SubstitutionalSolidSolutionDialog();
    ~SubstitutionalSolidSolutionDialog();

    // Call once after the GL context and renderer are initialised.
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
    // AtomForge --build sss).
    std::string stepOptions() const;

private:
    // -----------------------------------------------------------------------
    // One row in the composition table
    // -----------------------------------------------------------------------
    struct CompositionEntry
    {
        int   atomicNumber = 26; // Fe by default
        float percent      = 0.0f; // user-specified atomic percent (any scale; auto-normalised)
        int   count        = 0;    // derived via largest-remainder from percent
    };

    // -----------------------------------------------------------------------
    // Dialog state
    // -----------------------------------------------------------------------
    bool m_openRequested = false;
    bool m_isOpen        = false;
    StepEdit m_step;

    // The host structure used as lattice template.
    Structure   m_source;
    bool        m_sourceLoaded = false;
    std::string m_sourceLabel;
    // File the host was loaded from (empty for the scene or the step's input).
    std::string m_sourcePath;
    // The preview camera is fitted on the next draw (host set by editStep).
    bool        m_refitCamera = false;

    // Path queued by feedDroppedFile(); consumed at the start of drawDialog().
    std::string m_pendingDropPath;

    // Composition rows.
    std::vector<CompositionEntry> m_entries;

    // Index of the row currently opening a periodic-table picker (-1 = none).
    int m_pickerTarget = -1;

    // RNG seed.
    int m_seed = 12345;

    // Status / result message shown at the bottom.
    std::string m_statusMsg;
    bool        m_statusIsError = false;

    // -----------------------------------------------------------------------
    // GL preview resources
    // -----------------------------------------------------------------------
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

    bool m_glReady      = false;
    bool m_previewDirty = true;

    // Preview orbit camera state.
    float m_camYaw      = 45.0f;
    float m_camPitch    = 35.0f;
    float m_camDistance = 10.0f;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------
    bool loadFromPath(const std::string& path,
                      const std::vector<float>& radii,
                      const std::vector<float>& shininess);
    void loadFromScene(const Structure& scene,
                       const std::vector<float>& radii,
                       const std::vector<float>& shininess);
    void onStructureLoaded(const std::vector<float>& radii,
                           const std::vector<float>& shininess);

    // Initialise composition entries with one row per unique element,
    // counts proportional to current occurrence in the host structure.
    void initCompositionFromSource();

    // Recompute every entry's count from its percent using largest-remainder,
    // so all counts always sum to the total atom count regardless of whether
    // the percents add to 100.
    void recomputeCountsFromPercents();

    // When the user changes entry changedIdx's percent to newPct,
    // scale all other entries proportionally so they still sum to 100,
    // then recompute integer counts.
    void redistributePercent(int changedIdx, float newPct);

    // When the user changes entry changedIdx's percent to newPct (where
    // changedIdx is NOT the last entry), keep all other non-last entries
    // untouched and only adjust the last entry so the total stays at 100.
    void redistributePercentLastFree(int changedIdx, float newPct);

    // When the user changes count of entry changedIdx to newCount,
    // adjust other entries proportionally so the total stays constant.
    void redistributeCount(int changedIdx, int newCount);

    void addEntry();
    void removeEntry(int idx);

    // GL preview helpers.
    void ensurePreviewFBO(int w, int h);
    void rebuildPreviewBuffers(const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void renderPreviewToFBO(int w, int h);
    void autoFitPreviewCamera();
};
