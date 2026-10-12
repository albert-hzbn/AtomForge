#pragma once

#include "ui/PathPicker.h"
#include "ui/StepEditing.h"

#include "algorithms/InterfaceMatching.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "model/Structure.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;
namespace vkr { class Canvas; }

// Build > Interface Builder: two crystals, each cut along a plane, matched
// by their surface lattices (algorithms/InterfaceMatching.h); the matches are
// listed and plotted (strain against atoms), and the selected one is shown
// and built.
struct InterfaceBuilderDialog
{
    InterfaceBuilderDialog() = default;
    ~InterfaceBuilderDialog();

    void initRenderResources(Renderer& renderer);

    void drawMenuItem(bool enabled);
    // Opens the dialog, as its menu item does.
    void open() { m_step.finish(); m_openRequested = true; }

    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }

    // Open, or opening this frame (it may build and replace the structure).

    bool active() const { return m_isOpen || m_openRequested; }

    // A dropped file: slotHint 0 fills A, then B; 1 forces A, 2 forces B.
    void feedDroppedFile(const std::string& path, int slotHint = 0);

    // Opens the dialog on a pipeline step: filled from the step's settings, with
    // "Update Step" writing them back instead of building.
    bool editStep(StepEdit edit);
    // The dialog's settings as the step's options (flags of AtomForge --build interface).
    std::string stepOptions() const;
    // Opens on an example (UI tour): fcc Cu (1 1 1) under fcc Ni (1 1 1);
    // the tour then presses Find Matches and Build ("interface.find",
    // "interface.build").
    // With cubeOnHexagonal, Cu (1 0 0) under Ni (1 1 1) instead: a larger,
    // rotated match.
    void openExample(bool cubeOnHexagonal = false);

private:
    struct Input
    {
        Structure crystal;
        std::string path;    // the file it came from ("" for the active structure or the step input)
        std::string status;  // what it is, for the panel
        bool cleave = true;
        int hkl[3] = {1, 1, 1};
        int layers = 4;
    };

    bool m_openRequested = false;
    bool m_isOpen = false;
    StepEdit m_step;
    bool m_stepSearch = false;  // search once the step's crystals are loaded
    int m_stepPick = -1;        // the step's --pick

    Input m_a, m_b;
    std::string m_pendingDropA, m_pendingDropB;
    PathPicker m_picker;
    int m_pickerSlot = 0;

    atomforge::interfaces::SearchParams m_search;
    atomforge::interfaces::Assembly m_assembly;
    int m_repeatX = 1, m_repeatY = 1;

    // Results of the last search, and the settings it was run with.
    Structure m_slabA, m_slabB;
    std::vector<atomforge::interfaces::Match> m_matches;
    bool m_searched = false;
    bool m_stale = false;  // settings changed since the search
    std::string m_searchMessage;
    int m_selected = -1;
    int m_suggested = -1;
    bool m_scrollToSelected = false;

    // The selected match as built, and its preview.
    Structure m_preview;
    bool m_previewDirty = true;
    Renderer* m_renderer = nullptr;
    SphereMesh* m_sphere = nullptr;
    CylinderMesh* m_cylinder = nullptr;
    SceneBuffers m_previewBuffers;
    ShadowMap m_shadow = {};
    GLuint m_fbo = 0, m_texture = 0, m_depth = 0;
    int m_fboWidth = 0, m_fboHeight = 0;
    bool m_glReady = false;
    std::shared_ptr<vkr::Canvas> m_canvas;  // the preview with Vulkan
    float m_yaw = 30.0f, m_pitch = 15.0f, m_distance = 30.0f;

    bool loadFile(const std::string& path, Input& input);
    void useStructure(const Structure& structure, Input& input, const char* what);
    void runSearch();
    void select(int index);
    void rebuildPreview(const std::vector<float>& radii, const std::vector<float>& shininess);
    void drawInput(const char* title, Input& input, int slot, const Structure& active);
    void drawSearchSettings();
    void drawResultsTable(float height);
    void drawPlot(float width, float height);
    void drawPreview(float width, float height);
    void drawAssembly();
};
