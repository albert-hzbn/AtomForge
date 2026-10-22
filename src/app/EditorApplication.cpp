#include "app/EditorApplication.h"

#include "app/FileDropHandler.h"
#include "app/EditorOps.h"
#include "app/SceneView.h"
#include "app/EditorState.h"
#include "app/ImageExport.h"
#include "app/InteractionHandlers.h"
#include "app/LargeData.h"
#include "app/SceneBenchmark.h"
#include "app/StructureFileService.h"
#include "app/WindowSetup.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"
#include "imgui_internal.h"

#include "Camera.h"
#include "CylinderMesh.h"
#include "Renderer.h"
#include "ShadowMap.h"
#include "SphereMesh.h"
#include "LowPolyMesh.h"
#include "BillboardMesh.h"
#include "graphics/FrameCapture.h"
#include "graphics/GraphicsBackend.h"
#ifdef ATOMFORGE_VULKAN
#include "graphics/vulkan/VulkanSceneRenderer.h"
#endif
#include "ui/DialogLayout.h"
#include "ui/ImGuiSetup.h"
#include "ui/UiScript.h"
#include "ui/LatticePlaneOverlay.h"
#include "ui/MillerDirectionOverlay.h"
#include "ui/PolyhedralOverlay.h"
#include "ui/SurfaceMeshOverlay.h"
#include "ui/TrajectoryLinesOverlay.h"
#include "ui/VoronoiOverlay.h"

#include <glm/glm.hpp>

#include <cmath>
#include <filesystem>
#include <cstdlib>
#include "model/Workspace.h"
#include <iostream>
#include <string>
#include <vector>

namespace
{
// "Opened cu_host.cif (500 atoms)." and, when they were saved with it, that
// the crystal orientations came back too.
std::string loadedMessage(const std::string& path, const Structure& structure)
{
    const std::size_t slash = path.find_last_of("/\\");
    std::string message = "Opened " + (slash == std::string::npos ? path : path.substr(slash + 1)) + " (" +
                          std::to_string(structure.atoms.size()) + (structure.atoms.size() == 1 ? " atom)." : " atoms).");
    if (structure.ipfLoadStatus.rfind("IPF loaded", 0) == 0)
        message += " Crystal orientations restored from the saved metadata.";
    if (!structure.loadNote.empty())
        message += " " + structure.loadNote;
    return message;
}

const glm::vec4 kDarkBackground(0.09f, 0.11f, 0.15f, 1.0f);
const glm::vec4 kLightBackground(0.92f, 0.93f, 0.94f, 1.0f);

// Dynamically selects rendering mode based on measured frame performance.
// Starts from an initial estimate (set by SceneBuffers::upload based on atom
// count) and adjusts up or down if the GPU can handle a better mode or is
// struggling with the current one.
// Seconds the main loop has slept waiting for events since the adaptive
// renderer last measured (idle time is not drawing time).
double g_idleSeconds = 0.0;

// The Vulkan main view culls and picks a level of detail per atom on the GPU
// (unless ATOMFORGE_GPU_CULLING=0), so the whole-scene modes below are not used.
bool gpuDrivenScene()
{
#ifdef ATOMFORGE_VULKAN
    static const bool enabled = [] {
        const char* value = std::getenv("ATOMFORGE_GPU_CULLING");
        return !(value && value[0] == '0');
    }();
    return enabled && graphics::usingVulkan() && vkr::sceneRenderer().gpuDriven();
#else
    return false;
#endif
}

struct AdaptiveRenderState
{
    static constexpr int    kWindowSize          = 30;     // rolling window (frames)
    static constexpr int    kCooldownFrames      = 60;     // pause after a mode switch
    static constexpr double kDowngradeThresholdMs = 33.0;  // > 33 ms ≈ < 30 fps
    static constexpr double kUpgradeThresholdMs   = 20.0;  // < 20 ms ≈ > 50 fps

    double frameTimes[30] = {};
    int    frameIndex   = 0;
    int    sampleCount  = 0;
    int    cooldown     = 0;
    double lastTime     = 0.0;
    size_t lastAtomCount = 0;
    bool   initialized  = false;

    void update(SceneBuffers& buffers)
    {
        double now = glfwGetTime();
        if (!initialized)
        {
            lastTime     = now;
            lastAtomCount = buffers.atomCount;
            initialized  = true;
            return;
        }

        double dtMs = (now - lastTime - g_idleSeconds) * 1000.0;
        g_idleSeconds = 0.0;
        lastTime = now;

        // Reset measurements when the structure changes.
        if (buffers.atomCount != lastAtomCount)
        {
            lastAtomCount = buffers.atomCount;
            sampleCount   = 0;
            frameIndex    = 0;
            cooldown      = kCooldownFrames;
            return;
        }

        frameTimes[frameIndex] = dtMs;
        frameIndex = (frameIndex + 1) % kWindowSize;
        if (sampleCount < kWindowSize)
            sampleCount++;

        if (cooldown > 0) { cooldown--; return; }
        if (sampleCount < kWindowSize) return;

        // Compute rolling average frame time.
        double avg = 0.0;
        for (int i = 0; i < kWindowSize; i++)
            avg += frameTimes[i];
        avg /= kWindowSize;

        RenderingMode current = buffers.renderMode;

        if (avg > kDowngradeThresholdMs)
        {
            // Frame rate too low – switch to a cheaper mode.
            if (current == RenderingMode::StandardInstancing)
            {
                buffers.renderMode = RenderingMode::LowPolyInstancing;
                resetSamples();
                std::cout << "[Adaptive] Downgraded to LowPoly (avg " << avg << " ms)\n";
            }
            else if (current == RenderingMode::LowPolyInstancing)
            {
                buffers.renderMode = RenderingMode::BillboardImposters;
                resetSamples();
                std::cout << "[Adaptive] Downgraded to Billboard (avg " << avg << " ms)\n";
            }
        }
        else if (avg < kUpgradeThresholdMs)
        {
            // Plenty of headroom – try a higher-quality mode.
            if (current == RenderingMode::BillboardImposters)
            {
                buffers.renderMode = RenderingMode::LowPolyInstancing;
                resetSamples();
                std::cout << "[Adaptive] Upgraded to LowPoly (avg " << avg << " ms)\n";
            }
            else if (current == RenderingMode::LowPolyInstancing)
            {
                buffers.renderMode = RenderingMode::StandardInstancing;
                resetSamples();
                std::cout << "[Adaptive] Upgraded to Standard (avg " << avg << " ms)\n";
            }
        }
    }

    void resetSamples()
    {
        cooldown    = kCooldownFrames;
        sampleCount = 0;
        frameIndex  = 0;
    }
};

// ---- Per-tab state -------------------------------------------------------

struct TabCameraState
{
    float     yaw      = 45.0f;
    float     pitch    = 35.2643897f;
    float     roll     = 0.0f;
    float     distance = 10.0f;
    glm::vec3 panOffset = glm::vec3(0.0f);
};

struct StructureTab
{
    EditorState         state;
    TabCameraState      cameraState;
    AdaptiveRenderState adaptiveRender;
    char                title[256] = "Untitled";
    bool                pendingClose = false;
    int                 closeAfterSaveFrames = 0;  // > 0: close once saved (Save in the close prompt)
};

// Closing the window: the close button only asks; the main loop confirms
// unsaved changes first (programmatic closes, e.g. the UI tour, still quit).
bool g_userCloseRequested = false;

void saveCameraToTab(const Camera& cam, StructureTab& tab)
{
    tab.cameraState.yaw       = cam.yaw;
    tab.cameraState.pitch     = cam.pitch;
    tab.cameraState.roll      = cam.roll;
    tab.cameraState.distance  = cam.distance;
    tab.cameraState.panOffset = cam.panOffset;
}

void restoreCameraFromTab(Camera& cam, const StructureTab& tab)
{
    cam.yaw       = tab.cameraState.yaw;
    cam.pitch     = tab.cameraState.pitch;
    cam.roll      = tab.cameraState.roll;
    cam.distance  = tab.cameraState.distance;
    cam.panOffset = tab.cameraState.panOffset;
}

// Called after every load and save: the tab now shows that file, unmodified.
// Three right-aligned buttons for the unsaved-changes prompts.
int unsavedPromptButtons(const char* danger, const char* primary)
{
    const float minimum = ImGui::GetFontSize() * 5.5f;
    dialogLayout::beginFooter(primary ? dialogLayout::buttonRowWidth({danger, "Cancel", primary}, minimum)
                                      : dialogLayout::buttonRowWidth({danger, "Cancel"}, minimum));
    const auto width = [minimum](const char* label) {
        return ImVec2(std::max(minimum, ImGui::CalcTextSize(label).x + 2 * ImGui::GetStyle().FramePadding.x), 0);
    };
    int result = 0;
    if (dialogLayout::dangerButton(danger, width(danger))) result = 1;
    ImGui::SameLine();
    if (dialogLayout::secondaryButton("Cancel", width("Cancel")) || ImGui::IsKeyPressed(ImGuiKey_Escape)) result = 2;
    if (primary) {
        ImGui::SameLine();
        if (dialogLayout::primaryButton(primary, width(primary))) result = 3;
    }
    return result;
}

void drawUnsavedChangesPrompts(GLFWwindow* window, std::vector<std::unique_ptr<StructureTab>>& tabs, int activeTabIdx,
                               int& pendingTabSwitch, int& closePromptTab, bool& quitPrompt)
{
    // "file.cif \u2022 \u2014 AtomForge" (bullet: unsaved changes).
    {
        static std::string shownTitle;
        const StructureTab& active = *tabs[activeTabIdx];
        std::string title = "AtomForge";
        if (active.state.hasContent())
            title = std::string(active.title) + (active.state.isModified() ? " \xe2\x80\xa2" : "") + " \xe2\x80\x94 AtomForge";
        if (title != shownTitle) {
            glfwSetWindowTitle(window, title.c_str());
            shownTitle = title;
        }
    }

    // Save in the close prompt: close once saved; give up when the save
    // dialog closes without saving.
    for (auto& tab : tabs) {
        if (tab->closeAfterSaveFrames <= 0) continue;
        if (!tab->state.isModified()) { tab->pendingClose = true; tab->closeAfterSaveFrames = 0; }
        else if (tab->state.fileBrowser.isAnyDialogOpen()) tab->closeAfterSaveFrames = 30;
        else --tab->closeAfterSaveFrames;
    }

    if (uiScript::take("app.quit") || tabs[activeTabIdx]->state.fileBrowser.consumeQuitRequest()) g_userCloseRequested = true;  // the UI tour's stand-in for the close button
    if (g_userCloseRequested) {
        g_userCloseRequested = false;
        bool modified = false;
        for (const auto& tab : tabs) modified |= tab->state.isModified();
        if (modified) quitPrompt = true;
        else glfwSetWindowShouldClose(window, GLFW_TRUE);
    }

    if (closePromptTab >= static_cast<int>(tabs.size())) closePromptTab = -1;
    const ImGuiWindowFlags flags = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings;
    // Each prompt opens once per request; closed any other way, it is cancelled.
    static bool closeShown = false, quitShown = false;
    if (closePromptTab >= 0 && !closeShown) { ImGui::OpenPopup("Unsaved Changes##close"); closeShown = true; }
    if (!ImGui::IsPopupOpen("Unsaved Changes##close") && closeShown) { closeShown = false; closePromptTab = -1; }
    if (ImGui::BeginPopupModal("Unsaved Changes##close", nullptr, flags)) {
        if (closePromptTab < 0) ImGui::CloseCurrentPopup();
        else {
            StructureTab& tab = *tabs[closePromptTab];
            ImGui::Text("Save the changes to %s before closing it?", tab.title);
            ImGui::TextDisabled("Changes that are not saved are lost.");
            ImGui::Spacing();
            switch (unsavedPromptButtons("Don't Save", "Save")) {
            case 1: tab.pendingClose = true; closePromptTab = -1; ImGui::CloseCurrentPopup(); break;
            case 2: closePromptTab = -1; ImGui::CloseCurrentPopup(); break;
            case 3:
                if (closePromptTab != activeTabIdx) pendingTabSwitch = closePromptTab;
                tab.state.fileBrowser.saveFile();  // Save As when it has no file yet
                tab.closeAfterSaveFrames = 30;
                closePromptTab = -1;
                ImGui::CloseCurrentPopup();
                break;
            default: break;
            }
        }
        ImGui::EndPopup();
    }

    if (quitPrompt && !quitShown) { ImGui::OpenPopup("Unsaved Changes##quit"); quitShown = true; }
    if (!ImGui::IsPopupOpen("Unsaved Changes##quit") && quitShown) { quitShown = false; quitPrompt = false; }
    if (ImGui::BeginPopupModal("Unsaved Changes##quit", nullptr, flags)) {
        std::string names;
        int count = 0;
        for (const auto& tab : tabs)
            if (tab->state.isModified()) { names += (count++ ? ", " : "") + std::string(tab->title); }
        ImGui::Text("%s unsaved changes: %s", count == 1 ? "One structure has" : (std::to_string(count) + " structures have").c_str(), names.c_str());
        ImGui::TextDisabled("Quit anyway and lose them, or cancel and save them first.");
        ImGui::Spacing();
        switch (unsavedPromptButtons("Quit Without Saving", nullptr)) {
        case 1: quitPrompt = false; glfwSetWindowShouldClose(window, GLFW_TRUE); ImGui::CloseCurrentPopup(); break;
        case 2: quitPrompt = false; ImGui::CloseCurrentPopup(); break;
        default: break;
        }
        ImGui::EndPopup();
    }
}

void setTabTitleFromPath(StructureTab& tab, const std::string& path)
{
    if (path.empty()) { std::snprintf(tab.title, sizeof(tab.title), "Untitled"); return; }
    tab.state.markSaved();
    const auto sep = path.find_last_of("/\\");
    const std::string fn = (sep == std::string::npos) ? path : path.substr(sep + 1);
    std::snprintf(tab.title, sizeof(tab.title), "%s", fn.c_str());
}

void initTabResources(StructureTab& tab,
                      const SphereMesh& sphere,
                      const LowPolyMesh& lowPolyMesh,
                      const BillboardMesh& billboardMesh,
                      const CylinderMesh& cylinder,
                      Renderer& renderer)
{
    tab.state.sceneBuffers.init(
        sphere.vbo,      sphere.ebo,      sphere.indexCount,
        lowPolyMesh.vbo, lowPolyMesh.ebo, lowPolyMesh.indexCount,
        billboardMesh.vbo, billboardMesh.ebo, billboardMesh.indexCount,
        cylinder.vbo,    cylinder.vertexCount);
    tab.state.fileBrowser.initNanoCrystalRenderResources(renderer);
    tab.state.fileBrowser.initCustomStructureRenderResources(renderer);
    tab.state.fileBrowser.initMergeStructuresRenderResources(renderer);
    tab.state.fileBrowser.initInterfaceBuilderRenderResources(renderer);
    tab.state.fileBrowser.initCSLGrainBoundaryRenderResources(renderer);
    tab.state.fileBrowser.initPolyCrystalRenderResources(renderer);
    tab.state.fileBrowser.initDislocationRenderResources(renderer);
#if ATOMFORGE_ENABLE_SFE_BUILDER
    tab.state.fileBrowser.initStackingFaultRenderResources(renderer);
#endif
#if ATOMFORGE_ENABLE_SSS_BUILDER
    tab.state.fileBrowser.initSubstitutionalSolidSolutionRenderResources(renderer);
#endif
    tab.state.fileBrowser.initCellSculptorRenderResources(renderer);
    tab.state.fileBrowser.initInterstitialAtomsRenderResources(renderer);
    updateBuffers(tab.state);
    tab.state.undoRedo.reset(captureSnapshot(tab.state));
}

void applyPendingDefaultViewReset(Camera& camera, EditorState& state, const FrameView& frame)
{
    if (!state.pendingDefaultViewReset)
        return;
    if (state.hasCloud())
    {
        // Large datasets fit their own box (there are no scene atoms).
        state.pendingCloudFit = true;
        state.pendingDefaultViewReset = false;
        return;
    }

    applyDefaultView(camera, state.sceneBuffers, frame.framebufferWidth, frame.framebufferHeight, true);
    state.pendingDefaultViewReset = false;
}

void mergeFileBrowserRequests(EditorState& state,
                             FrameActionRequests& requests,
                             ImageExportRequest& imageExportRequest,
                             bool& hasImageExportRequest)
{
    hasImageExportRequest = state.fileBrowser.consumeImageExportRequest(imageExportRequest);
    requests.requestUndo = requests.requestUndo || state.fileBrowser.consumeUndoRequest();
    requests.requestRedo = requests.requestRedo || state.fileBrowser.consumeRedoRequest();
    requests.requestStructureInfo = requests.requestStructureInfo || state.fileBrowser.consumeStructureInfoRequest();
    requests.requestMeasureDistance = requests.requestMeasureDistance || state.fileBrowser.consumeMeasureDistanceRequest();
    requests.requestMeasureAngle = requests.requestMeasureAngle || state.fileBrowser.consumeMeasureAngleRequest();
    requests.requestViewAxisX = requests.requestViewAxisX || state.fileBrowser.consumeViewAxisXRequest();
    requests.requestViewAxisY = requests.requestViewAxisY || state.fileBrowser.consumeViewAxisYRequest();
    requests.requestViewAxisZ = requests.requestViewAxisZ || state.fileBrowser.consumeViewAxisZRequest();
    requests.requestViewLatticeA = requests.requestViewLatticeA || state.fileBrowser.consumeViewLatticeARequest();
    requests.requestViewLatticeB = requests.requestViewLatticeB || state.fileBrowser.consumeViewLatticeBRequest();
    requests.requestViewLatticeC = requests.requestViewLatticeC || state.fileBrowser.consumeViewLatticeCRequest();
    requests.requestRotateCrystalX = requests.requestRotateCrystalX || state.fileBrowser.consumeRotateCrystalXRequest();
    requests.requestRotateCrystalY = requests.requestRotateCrystalY || state.fileBrowser.consumeRotateCrystalYRequest();
    requests.requestRotateCrystalZ = requests.requestRotateCrystalZ || state.fileBrowser.consumeRotateCrystalZRequest();
}

void handleUndoRedoRequest(EditorState& state, const FrameActionRequests& requests)
{
    if (const int step = state.fileBrowser.consumeHistoryJump(); step >= 0 && static_cast<std::size_t>(step) != state.undoRedo.index())
    {
        applySnapshot(state, state.undoRedo.jumpTo(static_cast<std::size_t>(step)));
        std::cout << "[Operation] History: step " << step << std::endl;
        return;
    }
    if (requests.requestUndo && state.undoRedo.canUndo())
    {
        applySnapshot(state, state.undoRedo.undo());
        std::cout << "[Operation] Undo" << std::endl;
        return;
    }

    if (requests.requestRedo && state.undoRedo.canRedo())
    {
        applySnapshot(state, state.undoRedo.redo());
        std::cout << "[Operation] Redo" << std::endl;
    }
}

void handleAxisViewRequest(Camera& camera,
                           const FrameActionRequests& requests,
                           const Structure& structure)
{
    auto setViewAlongVector = [&](int vectorIndex, const char* label) {
        if (!structure.hasUnitCell || structure.atoms.empty())
            return;

        const float vx = (float)structure.cellVectors[vectorIndex][0];
        const float vy = (float)structure.cellVectors[vectorIndex][1];
        const float vz = (float)structure.cellVectors[vectorIndex][2];
        const float len = std::sqrt(vx * vx + vy * vy + vz * vz);
        if (len < 1e-6f)
            return;

        const float nx = vx / len;
        const float ny = vy / len;
        const float nz = vz / len;

        float yClamped = ny;
        if (yClamped > 1.0f) yClamped = 1.0f;
        if (yClamped < -1.0f) yClamped = -1.0f;

        camera.pitch = std::asin(yClamped) * 57.2957795f;
        camera.yaw = std::atan2(nx, nz) * 57.2957795f;
        std::cout << "[Operation] View perpendicular to lattice vector " << label << std::endl;
    };

    if (requests.requestViewAxisX)
    {
        camera.yaw = 90.0f;
        camera.pitch = 0.0f;
            camera.roll = 0.0f;
        std::cout << "[Operation] View along X axis" << std::endl;
    }
    else if (requests.requestViewAxisY)
    {
        camera.yaw = 0.0f;
        camera.pitch = 90.0f;
            camera.roll = 0.0f;
        std::cout << "[Operation] View along Y axis" << std::endl;
    }
    else if (requests.requestViewAxisZ)
    {
        camera.yaw = 0.0f;
        camera.pitch = 0.0f;
        camera.roll = 0.0f;
        std::cout << "[Operation] View along Z axis" << std::endl;
    }
    else if (requests.requestViewLatticeA)
    {
        setViewAlongVector(0, "a");
    }
    else if (requests.requestViewLatticeB)
    {
        setViewAlongVector(1, "b");
    }
    else if (requests.requestViewLatticeC)
    {
        setViewAlongVector(2, "c");
    }
}

void mergeContextRequests(FrameActionRequests& requests, const AtomRequests& contextRequests)
{
    requests.doDeleteSelected = requests.doDeleteSelected || contextRequests.doDelete;
    requests.requestMeasureDistance = requests.requestMeasureDistance || contextRequests.measureDistance;
    requests.requestMeasureAngle = requests.requestMeasureAngle || contextRequests.measureAngle;
    requests.requestAtomInfo = requests.requestAtomInfo || contextRequests.atomInfo;
}

// The view's background and cell outline colours, and whether the background is light.
struct SceneColours
{
    glm::vec4 background;
    glm::vec3 cell;
    bool light;
};

SceneColours sceneColours(const EditorState& state)
{
    const atomforge::display::DisplaySettings& display = state.editMenuDialogs.display;
    const bool lightTheme = state.fileBrowser.isLightThemeEnabled();
    SceneColours colours;
    colours.background = display.customBackground ? glm::vec4(display.background, 1.0f)
                                                   : (lightTheme ? kLightBackground : kDarkBackground);
    colours.light = 0.2126f * colours.background.r + 0.7152f * colours.background.g + 0.0722f * colours.background.b > 0.5f;
    colours.cell = display.customCellColour ? display.cellColour : (colours.light ? glm::vec3(0.25f) : glm::vec3(0.85f));
    return colours;
}

void drawScene(Renderer& renderer,
               const FrameView& frame,
               const ShadowMap& shadow,
               const SceneBuffers& sceneBuffers,
               bool showBonds,
               bool showAtoms,
               bool showBoundingBox,
               bool showDislocationLines,
               const SceneColours& colours,
               CloudRenderer* cloud = nullptr,
               bool orthographic = false)
{
    // GPU-driven frustum cull: one compute dispatch per frame builds a compact list
    // of visible atom indices in sceneBuffers.visibleIndexSSBO and writes the
    // instanceCount into sceneBuffers.drawIndirectBuffer.
    // A second dispatch (cullAtomsForShadow) culls against the light frustum so
    // only shadow-casting atoms are submitted to the shadow pass.
    if (showAtoms && sceneBuffers.gpuCullingReady)
    {
        int indexCount = 0;
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:  indexCount = sceneBuffers.tabSphereIndexCount;    break;
            case RenderingMode::LowPolyInstancing:   indexCount = sceneBuffers.tabLowPolyIndexCount;   break;
            case RenderingMode::BillboardImposters:  indexCount = sceneBuffers.tabBillboardIndexCount; break;
        }
        renderer.cullAtoms(sceneBuffers, frame.projection, frame.view, indexCount);
        renderer.cullAtomsForShadow(sceneBuffers, frame.lightMVP, indexCount);
    }

    // Shadow passes (render into shadow FBO at shadow resolution).
    // When GPU culling is active, use indirect shadow draw (only light-frustum atoms).
    if (showAtoms)
    {
    if (sceneBuffers.gpuCullingReady)
    {
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:
                renderer.drawShadowPassIndirect(shadow, frame.lightMVP, sceneBuffers,
                    sceneBuffers.tabSphereVAO_indirect);
                break;
            case RenderingMode::LowPolyInstancing:
                renderer.drawShadowPassIndirect(shadow, frame.lightMVP, sceneBuffers,
                    sceneBuffers.tabLowPolyVAO_indirect);
                break;
            case RenderingMode::BillboardImposters:
                renderer.drawShadowPassBillboardIndirect(shadow, frame.lightMVP,
                    frame.cameraPosition, sceneBuffers,
                    sceneBuffers.tabBillboardVAO_indirect);
                break;
        }
    }
    else
    {
    switch (sceneBuffers.renderMode)
    {
        case RenderingMode::StandardInstancing:
            renderer.drawShadowPass(shadow,
                sceneBuffers.tabSphereVAO, sceneBuffers.tabSphereIndexCount,
                frame.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::LowPolyInstancing:
            renderer.drawShadowPassLowPoly(shadow,
                sceneBuffers.tabLowPolyVAO, sceneBuffers.tabLowPolyIndexCount,
                frame.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::BillboardImposters:
            renderer.drawShadowPassBillboard(shadow,
                sceneBuffers.tabBillboardVAO, sceneBuffers.tabBillboardIndexCount,
                frame.lightMVP, frame.view, sceneBuffers.atomCount);
            break;
    }
    }
    }

    if (showBonds && showAtoms)
        renderer.drawBondShadowPass(shadow,
            sceneBuffers.tabCylinderVAO, sceneBuffers.tabCylinderVertexCount,
            frame.lightMVP, sceneBuffers.bondCount);

    // Restore viewport to screen size after shadow passes
    glViewport(0, 0, frame.framebufferWidth, frame.framebufferHeight);
    const bool lightTheme = colours.light;
    const glm::vec4& bg = colours.background;
    glClearColor(bg.r, bg.g, bg.b, bg.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // A large dataset: the out-of-core renderer draws everything.
    if (cloud && cloud->active())
    {
        cloud->settings.showBox = showBoundingBox;
        cloud->draw(frame.projection, frame.view, frame.framebufferWidth, frame.framebufferHeight, orthographic, lightTheme);
        return;
    }

    if (showBonds && showAtoms)
    {
        renderer.drawBonds(
            frame.projection,
            frame.view,
            frame.lightPosition,
            frame.cameraPosition,
            sceneBuffers.tabCylinderVAO, sceneBuffers.tabCylinderVertexCount,
            sceneBuffers.bondCount);
    }

    // Depth prepass: write depth only (no color) for atoms so the color pass can use
    // GL_LEQUAL and skip all over-drawn fragments.  Billboard mode is excluded because
    // its FS computes per-fragment depth via gl_FragDepth, making a separate prepass
    // incorrect (the billboard quad depth ≠ the reconstructed sphere-surface depth).
    const bool doDepthPrepass = showAtoms &&
        sceneBuffers.renderMode != RenderingMode::BillboardImposters &&
        sceneBuffers.atomCount > 0;

    if (doDepthPrepass)
    {
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        if (sceneBuffers.gpuCullingReady)
        {
            GLuint prepassVAO = (sceneBuffers.renderMode == RenderingMode::StandardInstancing)
                                ? sceneBuffers.tabSphereVAO_indirect
                                : sceneBuffers.tabLowPolyVAO_indirect;
            renderer.drawDepthPrepassIndirect(frame.projection, frame.view,
                                              sceneBuffers, prepassVAO);
        }
        else
        {
            GLuint  prepassVAO    = (sceneBuffers.renderMode == RenderingMode::StandardInstancing)
                                    ? sceneBuffers.tabSphereVAO
                                    : sceneBuffers.tabLowPolyVAO;
            int     prepassCount  = (sceneBuffers.renderMode == RenderingMode::StandardInstancing)
                                    ? sceneBuffers.tabSphereIndexCount
                                    : sceneBuffers.tabLowPolyIndexCount;
            renderer.drawDepthPrepass(frame.projection, frame.view,
                                      prepassVAO, prepassCount, sceneBuffers.atomCount);
        }
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glDepthFunc(GL_LEQUAL);
    }

    // Draw atoms based on rendering mode.
    // When GPU culling is available, use the indirect-draw path: the compute has
    // already written only visible atom indices into visibleIndexSSBO so the vertex
    // shader is invoked only for in-frustum atoms.  Otherwise fall back to the
    // legacy vertex-shader-culling path (all instances submitted, clipped on GPU).
    if (showAtoms)
    {
    if (sceneBuffers.gpuCullingReady)
    {
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:
                renderer.drawAtomsIndirect(
                    frame.projection, frame.view, frame.lightMVP,
                    frame.lightPosition, frame.cameraPosition,
                    shadow, sceneBuffers);
                break;
            case RenderingMode::LowPolyInstancing:
                renderer.drawAtomsLowPolyIndirect(
                    frame.projection, frame.view, frame.lightMVP,
                    frame.lightPosition, frame.cameraPosition,
                    shadow, sceneBuffers);
                break;
            case RenderingMode::BillboardImposters:
                renderer.drawAtomsBillboardIndirect(
                    frame.projection, frame.view, frame.lightMVP,
                    frame.lightPosition, frame.cameraPosition,
                    shadow, sceneBuffers);
                break;
        }
    }
    else
    {
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:
                renderer.drawAtoms(
                    frame.projection,
                    frame.view,
                    frame.lightMVP,
                    frame.lightPosition,
                    frame.cameraPosition,
                    shadow,
                    sceneBuffers.tabSphereVAO, sceneBuffers.tabSphereIndexCount,
                    sceneBuffers.atomCount);
                break;
            case RenderingMode::LowPolyInstancing:
                renderer.drawAtomsLowPoly(
                    frame.projection,
                    frame.view,
                    frame.lightMVP,
                    frame.lightPosition,
                    frame.cameraPosition,
                    shadow,
                    sceneBuffers.tabLowPolyVAO, sceneBuffers.tabLowPolyIndexCount,
                    sceneBuffers.atomCount);
                break;
            case RenderingMode::BillboardImposters:
                renderer.drawAtomsBillboard(
                    frame.projection,
                    frame.view,
                    frame.lightMVP,
                    frame.lightPosition,
                    frame.cameraPosition,
                    shadow,
                    sceneBuffers.tabBillboardVAO, sceneBuffers.tabBillboardIndexCount,
                    sceneBuffers.atomCount);
                break;
        }
    }
    }

    if (doDepthPrepass)
        glDepthFunc(GL_LESS);

    if (showBoundingBox)
    {
    renderer.drawBoxLines(
        frame.projection,
        frame.view,
        sceneBuffers.lineVAO,
        sceneBuffers.boxLines.size(),
        colours.cell);
    }

    if (showDislocationLines
        && sceneBuffers.dislocationLineVAO != 0
        && sceneBuffers.dislocationLoopPointCount >= 2)
    {
        renderer.drawLineLoop(
            frame.projection,
            frame.view,
            sceneBuffers.dislocationLineVAO,
            sceneBuffers.dislocationLoopPointCount,
            glm::vec3(1.0f, 0.25f, 0.1f));  // red-orange
    }
}

#ifdef ATOMFORGE_VULKAN
// The scene in Vulkan: the shadow map (outside the window's pass), then the
// window's pass with the scene and the selection halos; the interface is
// added to the same pass before it ends.
void drawSceneVulkan(const Renderer& renderer,
                     const FrameView& frame,
                     const SceneBuffers& sceneBuffers,
                     bool showBonds,
                     bool showAtoms,
                     bool showBoundingBox,
                     bool showDislocationLines,
                     const SceneColours& colours,
                     const std::vector<glm::vec3>& haloPositions,
                     const std::vector<float>& haloRadii,
                     CloudRenderer* cloud,
                     bool orthographic)
{
    vkr::Context& ctx = vkr::context();
    vkr::SceneRenderer& scene = vkr::sceneRenderer();
    VkCommandBuffer cmd = ctx.commandBuffer();
    vkr::SceneView view;
    view.projection = frame.projection;
    view.view = frame.view;
    view.lightMVP = frame.lightMVP;
    view.lightPosition = frame.lightPosition;
    view.cameraPosition = frame.cameraPosition;
    vkr::SceneLighting lighting;
    lighting.ambient = renderer.lightAmbient;
    lighting.saturation = renderer.lightSaturation;
    lighting.contrast = renderer.lightContrast;
    lighting.shadowStrength = renderer.lightShadowStrength;
    lighting.specularIntensity = renderer.materialSpecularIntensity;
    lighting.shininessScale = renderer.materialShininessScale;
    lighting.shininessFloor = renderer.materialShininessFloor;
    vkr::SceneFlags flags;
    flags.atoms = showAtoms;
    flags.bonds = showBonds;
    flags.box = showBoundingBox;
    flags.dislocation = showDislocationLines;
    const bool lightTheme = colours.light;
    flags.lightTheme = lightTheme;
    const glm::vec4& bg = colours.background;
    const float clear[4] = {bg.r, bg.g, bg.b, bg.a};
    if (cloud && cloud->active())
    {
        // A large dataset: the out-of-core renderer draws everything.
        cloud->settings.showBox = showBoundingBox;
        cloud->prepareVulkan(cmd, frame.projection, frame.view, frame.framebufferWidth, frame.framebufferHeight, orthographic);
        ctx.beginMainPass(clear);
        cloud->drawVulkan(cmd, static_cast<int>(ctx.extent().width), static_cast<int>(ctx.extent().height), lightTheme);
        return;
    }
    const uint32_t offset = scene.prepare(view, lighting);
    const glm::vec3 boxColor = colours.cell;
    if (gpuDrivenScene())
    {
        scene.beginGpuScene(cmd, sceneBuffers, view, offset, flags, ctx.extent());
        ctx.beginMainPass(clear, scene.occlusionCulling());
        scene.drawGpuFirst(cmd, sceneBuffers, offset, flags, ctx.extent());
        ctx.pauseMainPass();
        scene.cullGpuSecond(cmd, sceneBuffers, flags, ctx.extent());
        ctx.resumeMainPass();
        scene.drawGpuSecond(cmd, sceneBuffers, offset, flags, ctx.extent(), boxColor);
    }
    else
    {
        scene.renderShadowMap(cmd, sceneBuffers, offset, flags);
        ctx.beginMainPass(clear);
        scene.drawScene(cmd, sceneBuffers, offset, flags, ctx.extent(), boxColor);
    }
    scene.drawHalos(cmd, offset, haloPositions, haloRadii, ctx.extent());
}
#endif

void handleImageExportIfRequested(bool hasImageExportRequest,
                                  const ImageExportRequest& imageExportRequest,
                                  const FrameView& frame,
                                  EditorState& state,
                                  Renderer& renderer,
                                  ShadowMap& shadow)
{
    if (!hasImageExportRequest)
        return;

    ImageExportView exportView;
    const int scale = std::max(1, imageExportRequest.resolutionScale);
    exportView.width = frame.framebufferWidth * scale;
    exportView.height = frame.framebufferHeight * scale;
    exportView.projection = frame.fullProjection;
    exportView.view = frame.view;
    exportView.lightMVP = frame.lightMVP;
    exportView.lightPosition = frame.lightPosition;
    exportView.cameraPosition = frame.cameraPosition;
    const SceneColours colours = sceneColours(state);
    exportView.boxColor = colours.cell;

    std::string exportError;
    const bool exportOk = exportStructureImage(
        imageExportRequest,
        exportView,
        colours.background,
        state.fileBrowser.isShowBondsEnabled(),
        state.fileBrowser.isShowAtomsEnabled(),
        state.fileBrowser.isShowBoundingBoxEnabled(),
        state.sceneBuffers,
        renderer,
        shadow,
        state.structure,
        exportError);

    if (exportOk)
    {
        state.fileBrowser.showNotification(std::string("Image exported: ") + imageExportRequest.outputPath);
        std::cout << "[Operation] Image exported: " << imageExportRequest.outputPath << std::endl;
        return;
    }

    const std::string message = exportError.empty()
        ? "Failed to export image."
        : ("Failed to export image: " + exportError);
    state.fileBrowser.showLoadError(message);
    std::cout << "[Operation] Image export failed: " << imageExportRequest.outputPath
              << " (" << (exportError.empty() ? "Unknown error" : exportError) << ")" << std::endl;
}

// Trajectory GIF export shows one frame per main-loop pass; capture it here
// with the current camera after the scene buffers have been updated.
void captureTrajectoryGifFrame(const FrameView& frame,
                               EditorState& state,
                               Renderer& renderer,
                               ShadowMap& shadow)
{
    if (!state.fileBrowser.trajectoryGifCapturePending())
        return;

    ImageExportView exportView;
    exportView.width = frame.framebufferWidth;
    exportView.height = frame.framebufferHeight;
    exportView.projection = frame.fullProjection;
    exportView.view = frame.view;
    exportView.lightMVP = frame.lightMVP;
    exportView.lightPosition = frame.lightPosition;
    exportView.cameraPosition = frame.cameraPosition;
    const SceneColours colours = sceneColours(state);
    exportView.boxColor = colours.cell;

    std::vector<unsigned char> pixels;
    std::string error;
    if (renderSceneToRgba(exportView,
                          colours.background,
                          state.fileBrowser.isShowBondsEnabled(),
                          state.fileBrowser.isShowAtomsEnabled(),
                          state.fileBrowser.isShowBoundingBoxEnabled(),
                          state.sceneBuffers,
                          renderer,
                          shadow,
                          pixels,
                          error))
        state.fileBrowser.addTrajectoryGifFrame(pixels, exportView.width, exportView.height);
    else
        state.fileBrowser.addTrajectoryGifFrame({}, 0, 0);  // stops the export
}

void loadStartupStructureIfRequested(StructureTab& tab, const std::string& startupStructurePath)
{
    EditorState& state = tab.state;
    if (startupStructurePath.empty())
    {
        state.fileBrowser.initFromPath("");
        return;
    }

    if (opensAsCloud(std::filesystem::u8path(startupStructurePath)))
    {
        state.fileBrowser.initFromPath("");
        state.pendingExternalLoadPaths.push_back(startupStructurePath);
        return;
    }

    Structure loadedStructure;
    std::string loadError;
    if (!loadStructureFromPath(startupStructurePath, loadedStructure, loadError))
    {
        state.fileBrowser.initFromPath(startupStructurePath);
        state.fileBrowser.showLoadError(loadError);
        std::cout << "[Operation] Startup load failed: " << startupStructurePath
                  << " (" << loadError << ")" << std::endl;
        return;
    }

    state.structure = std::move(loadedStructure);
    state.fileBrowser.initFromPath(startupStructurePath);
    state.fileBrowser.applyElementColorOverrides(state.structure);
    state.fileBrowser.showLoadInfo(loadedMessage(startupStructurePath, state.structure));
    updateBuffers(state);
    state.pendingDefaultViewReset = true;
    setTabTitleFromPath(tab, startupStructurePath);

    std::cout << "[Operation] Startup-loaded structure: " << startupStructurePath
              << " (atoms=" << state.structure.atoms.size() << ")" << std::endl;
}
} // namespace

int runAtomsEditor(const std::vector<std::string>& startupPaths)
{
    SplashScreen* splash = createSplashScreen();
    updateSplashScreen(splash, 0.15f, "Creating window");

    GLFWwindow* window = createMainWindow();
    graphics::setMainWindow(window);
    if (!window)
    {
        destroySplashScreen(splash);
        glfwTerminate();
        return -1;
    }

    Camera camera;
    configureCameraCallbacks(window, camera);
    glfwSetWindowCloseCallback(window, [](GLFWwindow* w) {
        g_userCloseRequested = true;
        glfwSetWindowShouldClose(w, GLFW_FALSE);
    });
    updateSplashScreen(splash, 0.28f, "Initializing UI");

    initImGui(window);
    // The user's theme, interface size and saved display settings; the review
    // and capture modes keep the defaults so their screenshots do not depend on them.
    atomforge::display::DisplaySettings savedDisplay;
    {
        bool harness = false;
        for (const char* mode : {"ATOMFORGE_UI_TOUR", "ATOMFORGE_TUTORIAL_CAPTURE", "ATOMFORGE_STEP_DIALOG_CHECK", "ATOMFORGE_CLOUD_BENCHMARK", "ATOMFORGE_SCENE_BENCHMARK"})
            if (const char* value = std::getenv(mode); value && *value) harness = true;
        if (!harness) {
            loadAppearancePreferences();
            savedDisplay = atomforge::display::loadDefaults();
        }
        // Modes taking screenshots keep every frame (the cloud benchmark only
        // with ATOMFORGE_CLOUD_SHOTS, as copying them would skew its timings).
        bool screenshots = false;
        for (const char* mode : {"ATOMFORGE_UI_TOUR", "ATOMFORGE_TUTORIAL_CAPTURE", "ATOMFORGE_STEP_DIALOG_CHECK", "ATOMFORGE_CLOUD_SHOTS", "ATOMFORGE_SCENE_SHOTS"})
            if (const char* value = std::getenv(mode); value && *value) screenshots = true;
        if (screenshots) graphics::enableFrameCapture();
    }

    SphereMesh sphere(40, 40);
    LowPolyMesh lowPolyMesh;
    BillboardMesh billboardMesh;
    CylinderMesh cylinder(32);
    updateSplashScreen(splash, 0.42f, "Preparing scene");

    Renderer renderer;
#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        if (!vkr::sceneRenderer().init(sphere, lowPolyMesh, billboardMesh, cylinder))
            std::cerr << "[Vulkan] The scene pipelines could not be created." << std::endl;
    }
    else
#endif
        renderer.init();
    updateSplashScreen(splash, 0.52f, startupPaths.empty() ? "Preparing workspace" : "Loading structure");

    // Multi-tab state — unique_ptr because EditorState holds non-movable atomics/threads
    std::vector<std::unique_ptr<StructureTab>> tabs;
    tabs.push_back(std::make_unique<StructureTab>());
    tabs[0]->state.editMenuDialogs.display = savedDisplay;
    int activeTabIdx     = 0;
    int pendingTabSwitch = -1;
    int closePromptTab = -1;   // tab whose close waits for the unsaved-changes prompt
    bool quitPrompt = false;   // the window close waits for it

    // Pending file drops (shared; routed to active tab each frame)
    std::vector<std::string> appPendingDrops;
    installDropFileCallback(window, appPendingDrops);

    // Initialise the first tab and load the first startup path (if any)
    initTabResources(*tabs[0], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
    loadStartupStructureIfRequested(*tabs[0], startupPaths.empty() ? std::string() : startupPaths[0]);
    restoreCameraFromTab(camera, *tabs[0]);
    if (!graphics::fallbackReason().empty())
        tabs[0]->state.fileBrowser.showLoadInfo("Vulkan could not start (" + graphics::fallbackReason() + "); drawing with OpenGL.");

    // Load any additional startup paths each into their own tab
    for (int si = 1; si < (int)startupPaths.size(); ++si)
    {
        tabs.push_back(std::make_unique<StructureTab>());
        int newIdx = (int)tabs.size() - 1;
        initTabResources(*tabs[newIdx], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
        loadStartupStructureIfRequested(*tabs[newIdx], startupPaths[si]);
    }

    ShadowMap shadow{};
    if (graphics::usingOpenGL())
        shadow = createShadowMap(2048, 2048);
    updateSplashScreen(splash, 0.92f, "Finalizing");

    showMainWindow(window);
    updateSplashScreen(splash, 1.0f, "Ready");
    destroySplashScreen(splash);

    // Propagate display settings (lighting, material, theme) from source tab to a new tab.
    auto copyDisplaySettings = [](const EditorState& src, EditorState& dst)
    {
        dst.editMenuDialogs.display = src.editMenuDialogs.display;
    };

    atomforge::BackgroundTask<bool> autosaveTask;
    double lastAutosave=glfwGetTime();
    std::string autosavePath;
    try {
        std::filesystem::path root;
#ifdef _WIN32
        if (const char* folder=std::getenv("LOCALAPPDATA")) root=std::filesystem::u8path(folder);
#else
        if (const char* folder=std::getenv("XDG_STATE_HOME")) root=std::filesystem::u8path(folder);
        else if (const char* userDirectory=std::getenv("HOME")) root=std::filesystem::u8path(userDirectory)/".local"/"state";
#endif
        if (!root.empty()) {
            root/="AtomForge"; std::filesystem::create_directories(root);
            autosavePath=(root/"recovery.afproject").u8string();
        }
    } catch (const std::exception& error) { std::cerr << "Autosave unavailable: " << error.what() << '\n'; }

    CloudConversion cloudConversion;
    std::unique_ptr<CloudBenchmark> cloudBenchmark = CloudBenchmark::fromEnvironment();
    std::unique_ptr<SceneBenchmark> sceneBenchmark = SceneBenchmark::fromEnvironment();
    double lastFrameTime = glfwGetTime();
    float frameMilliseconds = 0.0f;

    while (!glfwWindowShouldClose(window))
    {
        {
            const double now = glfwGetTime();
            frameMilliseconds = static_cast<float>(1000.0 * (now - lastFrameTime));
            lastFrameTime = now;
        }
        // --- Active-tab alias ---
        EditorState& state = tabs[activeTabIdx]->state;

        camera.allowPan    = !state.fileBrowser.isBoxSelectModeEnabled()
                           && !state.fileBrowser.isLassoSelectModeEnabled();
        camera.allowOrbit  = !state.grabState.active;

        // Idle: wait for input instead of redrawing continuously (at most
        // 0.1 s, so background results still appear). Full speed while the
        // user interacts (and a second after, for camera inertia), while
        // something animates, and always in the review and capture modes.
        {
            static double lastActivity = 0.0;
            static const bool harness = [] {
                for (const char* mode : {"ATOMFORGE_UI_TOUR", "ATOMFORGE_TUTORIAL_CAPTURE", "ATOMFORGE_STEP_DIALOG_CHECK", "ATOMFORGE_CLOUD_BENCHMARK", "ATOMFORGE_SCENE_BENCHMARK"})
                    if (const char* value = std::getenv(mode); value && *value) return true;
                return false;
            }();
            const ImGuiIO& io = ImGui::GetIO();
            bool keyDown = false;
            for (int key = ImGuiKey_NamedKey_BEGIN; key < ImGuiKey_NamedKey_END && !keyDown; ++key)
                keyDown = ImGui::IsKeyDown(static_cast<ImGuiKey>(key));
            const bool animating = FileBrowser::notificationsShowing() || state.fileBrowser.trajectoryPlaying() ||
                                   (state.hasCloud() && !state.cloud->stats().complete) || polyhedralOverlayFilling();
            const double now = glfwGetTime();
            if (io.MouseDelta.x != 0.0f || io.MouseDelta.y != 0.0f || io.MouseWheel != 0.0f || io.MouseWheelH != 0.0f ||
                ImGui::IsAnyMouseDown() || keyDown || io.WantTextInput || animating)
                lastActivity = now;
            if (harness || now - lastActivity < 1.0) glfwPollEvents();
            else
            {
                glfwWaitEventsTimeout(0.1);
                g_idleSeconds += glfwGetTime() - now;
            }
        }

        // Route app-level drops to the active tab
        for (auto& f : appPendingDrops)
            state.pendingDroppedFiles.push_back(std::move(f));
        appPendingDrops.clear();
        processDroppedFiles(state);

        // Helper: load a structure file into a target tab (creates new tab when
        // current tab already has atoms). Returns true on success.
        // Opens an atom cloud in the active tab when it is empty, else in a new tab.
        auto openCloudIntoTab = [&](const std::string& path) -> bool
        {
            auto cloudRenderer = std::make_unique<CloudRenderer>();
            try { cloudRenderer->open(std::filesystem::u8path(path)); }
            catch (const std::exception& error)
            {
                state.fileBrowser.showLoadError(error.what());
                std::cout << "[Operation] Cloud open failed: " << path << " (" << error.what() << ")" << std::endl;
                return false;
            }
            int targetIdx = activeTabIdx;
            if (tabs[activeTabIdx]->state.hasContent())
            {
                saveCameraToTab(camera, *tabs[activeTabIdx]);
                tabs.push_back(std::make_unique<StructureTab>());
                targetIdx = (int)tabs.size() - 1;
                initTabResources(*tabs[targetIdx], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
                copyDisplaySettings(tabs[activeTabIdx]->state, tabs[targetIdx]->state);
                activeTabIdx     = targetIdx;
                pendingTabSwitch = targetIdx;
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
            }
            EditorState& t = tabs[targetIdx]->state;
            t.structure = Structure();
            updateBuffers(t);
            t.cloud = std::move(cloudRenderer);
            t.sceneBuffers.orbitCenter = t.cloud->center();
            t.pendingDefaultViewReset = false;
            t.pendingCloudFit = true;
            t.fileBrowser.showLoadInfo("Large dataset: " + std::to_string(t.cloud->stats().totalAtoms) + " atoms (view only)");
            setTabTitleFromPath(*tabs[targetIdx], path);
            std::cout << "[Operation] Opened atom cloud: " << path << " (atoms=" << t.cloud->stats().totalAtoms << ")" << std::endl;
            return true;
        };

        auto loadPathIntoTab = [&](const std::string& path) -> bool
        {
            if (path.empty()) return false;
            if (opensAsCloud(std::filesystem::u8path(path)))
            {
                if (atomforge::cloud::isCloudFile(std::filesystem::u8path(path))) return openCloudIntoTab(path);
                // Converted in the background once, then opened.
                cloudConversion.start(std::filesystem::u8path(path));
                return true;
            }

            Structure newStructure;
            std::string loadError;
            if (!loadStructureFromPath(path, newStructure, loadError))
            {
                state.fileBrowser.showLoadError(loadError);
                std::cout << "[Operation] Load failed: " << path
                          << " (" << loadError << ")" << std::endl;
                return false;
            }

            // Use current tab only when it is empty; otherwise open a new tab.
            int targetIdx = activeTabIdx;
            if (tabs[activeTabIdx]->state.hasContent())
            {
                saveCameraToTab(camera, *tabs[activeTabIdx]);
                tabs.push_back(std::make_unique<StructureTab>());
                targetIdx = (int)tabs.size() - 1;
                initTabResources(*tabs[targetIdx], sphere, lowPolyMesh,
                                 billboardMesh, cylinder, renderer);
                copyDisplaySettings(tabs[activeTabIdx]->state, tabs[targetIdx]->state);
                activeTabIdx     = targetIdx;
                pendingTabSwitch = targetIdx;
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
            }

            EditorState& t = tabs[targetIdx]->state;
            t.structure = std::move(newStructure);
            t.fileBrowser.initFromPath(path);
            t.fileBrowser.applyElementColorOverrides(t.structure);
            t.fileBrowser.showLoadInfo(loadedMessage(path, t.structure));
            updateBuffers(t);
            t.pendingDefaultViewReset = true;
            setTabTitleFromPath(*tabs[targetIdx], path);
            std::cout << "[Operation] Loaded structure: " << path
                      << " (atoms=" << t.structure.atoms.size() << ")" << std::endl;
            return true;
        };

        // Handle file drop load requests (one or more files dropped at once)
        if (!state.pendingExternalLoadPaths.empty())
        {
            std::vector<std::string> dropPaths = std::move(state.pendingExternalLoadPaths);
            state.pendingExternalLoadPaths.clear();
            for (const auto& dropPath : dropPaths)
                loadPathIntoTab(dropPath);
        }

        // Drain builder-created structures into new tabs
        if (!state.pendingNewTabStructures.empty())
        {
            std::vector<Structure> newStructures = std::move(state.pendingNewTabStructures);
            state.pendingNewTabStructures.clear();
            for (auto& newStruct : newStructures)
            {
                // Reuse the current tab if it is empty; otherwise open a new one.
                // This prevents a stale empty tab from remaining after a build.
                int targetIdx = activeTabIdx;
                if (tabs[activeTabIdx]->state.hasContent())
                {
                    saveCameraToTab(camera, *tabs[activeTabIdx]);
                    tabs.push_back(std::make_unique<StructureTab>());
                    targetIdx = (int)tabs.size() - 1;
                    initTabResources(*tabs[targetIdx], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
                    copyDisplaySettings(tabs[activeTabIdx]->state, tabs[targetIdx]->state);
                    activeTabIdx     = targetIdx;
                    pendingTabSwitch = targetIdx;
                    restoreCameraFromTab(camera, *tabs[activeTabIdx]);
                }
                EditorState& t = tabs[targetIdx]->state;
                t.structure = std::move(newStruct);
                updateBuffers(t);
                t.pendingDefaultViewReset = true;
                std::snprintf(tabs[targetIdx]->title, sizeof(tabs[targetIdx]->title), "Untitled");
            }
        }

        // Update tab title if any initFromPath was called this frame
        {
            const std::string newPath = state.fileBrowser.consumeLastLoadedPath();
            if (!newPath.empty())
                setTabTitleFromPath(*tabs[activeTabIdx], newPath);
        }

        FrameView frame;
        if (!updateViewport(window, frame))
            continue;

        applyPendingDefaultViewReset(camera, state, frame);
        camera.applyScrollVelocity(glfwGetTime());
        if (state.hasCloud())
        {
            if (state.pendingCloudFit)
            {
                fitCameraToCloud(camera, *state.cloud, frame.framebufferWidth, frame.framebufferHeight);
                state.pendingCloudFit = false;
            }
            if (cloudBenchmark)
                cloudBenchmark->beforeFrame(camera, *state.cloud, frame.framebufferWidth, frame.framebufferHeight);
        }
        if (sceneBenchmark && !state.hasCloud())
            sceneBenchmark->beforeFrame(camera, state.sceneBuffers);
        buildFrameView(camera, state.sceneBuffers, state.fileBrowser.isOrthographicViewEnabled(), frame,
                       state.editMenuDialogs.display);
        if (state.hasCloud())
            adjustCloudProjection(frame, camera, *state.cloud, state.fileBrowser.isOrthographicViewEnabled());
        // The view fills what docked windows leave free (as of the last frame).
        for (int k = 0; k < 4; ++k) frame.area[k] = state.fileBrowser.viewArea()[static_cast<std::size_t>(k)];
        fitProjectionToArea(frame);

        if (state.grabState.active)
        {
            camera.pendingClick      = false;
            camera.pendingRightClick = false;
        }
        else
        {
            handlePendingAtomPick(camera, state, frame.cameraPosition,
                                  frame.windowWidth, frame.windowHeight,
                                  frame.projection, frame.view);
            handleRightClick(camera, state);
        }

#ifdef ATOMFORGE_VULKAN
        // The frame's command buffer is open from here: previews in the
        // interface record their passes into it.
        if (graphics::usingVulkan())
        {
            vkr::sceneRenderer().newFrame();
            if (!vkr::context().beginFrame())
                continue;
        }
#endif
        updateImGuiScale(window);
        updateSystemTheme();
        newImGuiRendererFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        drawAtomHoverTooltip(state, frame.cameraPosition, frame.windowWidth, frame.windowHeight, frame.projection, frame.view);

        FrameActionRequests requests = beginFrameActionRequests(state);
        applyKeyboardShortcuts(state, requests);

        handleGrabMode(state, camera, frame.projection, frame.view,
                       frame.windowWidth, frame.windowHeight);

        // Menu operations see the view's selection as atom indices.
        {
            std::vector<int> atoms;
            for (int instance : state.selectedInstanceIndices)
                if (instance >= 0 && instance < static_cast<int>(state.sceneBuffers.atomIndices.size()))
                    atoms.push_back(state.sceneBuffers.atomIndices[static_cast<std::size_t>(instance)]);
            std::sort(atoms.begin(), atoms.end());
            atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
            state.fileBrowser.setViewSelection(std::move(atoms));
        }
        state.fileBrowser.setLargeDatasetShown(state.hasCloud());
        state.fileBrowser.setUndoDescriptions(state.undoRedo.undoDescription(), state.undoRedo.redoDescription());
        {
            std::vector<std::string> history;
            history.reserve(state.undoRedo.size());
            for (std::size_t i = 0; i < state.undoRedo.size(); ++i) history.push_back(state.undoRedo.description(i));
            state.fileBrowser.setHistory(std::move(history), static_cast<int>(state.undoRedo.index()));
        }
        state.fileBrowser.draw(
            state.structure,
            state.editMenuDialogs,
            [&](Structure& structure) {
                const bool previous=state.suppressHistoryCommit;
                state.suppressHistoryCommit=previous || state.fileBrowser.trajectoryPlaying() || state.fileBrowser.pipelineUpdating();
                updateBuffers(state, structure);
                state.suppressHistoryCommit=previous;
            },
            [&](Structure s) { state.pendingNewTabStructures.push_back(std::move(s)); },
            state.undoRedo.canUndo(),
            state.undoRedo.canRedo());

        // A selection made by a menu operation selects every image of those atoms.
        if (std::vector<int> atoms; state.fileBrowser.consumeSelectionRequest(atoms)) {
            std::vector<char> wanted(state.structure.atoms.size(), 0);
            for (int a : atoms) if (a >= 0 && static_cast<std::size_t>(a) < wanted.size()) wanted[static_cast<std::size_t>(a)] = 1;
            state.selectedInstanceIndices.clear();
            for (std::size_t instance = 0; instance < state.sceneBuffers.atomIndices.size(); ++instance) {
                const int atom = state.sceneBuffers.atomIndices[instance];
                if (atom >= 0 && static_cast<std::size_t>(atom) < wanted.size() && wanted[static_cast<std::size_t>(atom)])
                    state.selectedInstanceIndices.push_back(static_cast<int>(instance));
            }
        }
        const auto captureProject=[&]() {
            saveCameraToTab(camera,*tabs[activeTabIdx]);
            std::vector<atomforge::Workspace> saved;
            for (const auto& tab:tabs) {
                if (tab->state.fileBrowser.workspaceBusy()) throw std::runtime_error("Wait for or cancel calculations before saving the project");
                auto item=tab->state.fileBrowser.workspace(tab->state.structure);
                const auto& c=tab->cameraState;
                item.camera={c.yaw,c.pitch,c.roll,c.distance,c.panOffset.x,c.panOffset.y,c.panOffset.z};
                item.title=tab->title; saved.push_back(std::move(item));
            }
            return saved;
        };
        if (autosaveTask.poll() && !autosaveTask.error().empty())
            state.fileBrowser.showLoadError("Autosave failed: "+autosaveTask.error());
        const auto projectRequest=state.fileBrowser.drawProjectPicker();
        try {
            if (projectRequest.action==2) {
                atomforge::saveWorkspace(captureProject(),projectRequest.path);
                for (auto& tab:tabs) tab->state.markSaved();
                state.fileBrowser.showNotification("Project saved: "+projectRequest.path);
            } else if (projectRequest.action==1 || projectRequest.action==3) {
                const auto path=projectRequest.action==3 ? autosavePath : projectRequest.path;
                auto saved=atomforge::loadWorkspace(path);
                // Build all tabs before publishing any of them. Existing work remains open.
                std::vector<std::unique_ptr<StructureTab>> restored;
                for (const auto& item:saved) {
                    auto tab=std::make_unique<StructureTab>();
                    initTabResources(*tab,sphere,lowPolyMesh,billboardMesh,cylinder,renderer);
                    tab->state.fileBrowser.restoreWorkspace(item,tab->state.structure);
                    tab->state.editMenuDialogs.display=tabs[activeTabIdx]->state.editMenuDialogs.display;
                    if (item.camera.size()==7) {
                        auto& c=tab->cameraState;
                        c.yaw=item.camera[0]; c.pitch=item.camera[1]; c.roll=item.camera[2];
                        c.distance=std::max(.01, item.camera[3]);
                        c.panOffset={item.camera[4],item.camera[5],item.camera[6]};
                    }
                    std::snprintf(tab->title,sizeof(tab->title),"%s",item.title.empty() ? "Restored project" : item.title.c_str());
                    updateBuffers(tab->state); tab->state.pendingDefaultViewReset=false;
                    tab->state.undoRedo.reset(captureSnapshot(tab->state));
                    tab->state.markSaved();
                    restored.push_back(std::move(tab));
                }
                saveCameraToTab(camera,*tabs[activeTabIdx]);
                activeTabIdx=static_cast<int>(tabs.size()); pendingTabSwitch=activeTabIdx;
                for (auto& tab:restored) tabs.push_back(std::move(tab));
                restoreCameraFromTab(camera,*tabs[activeTabIdx]);
            }
            if (state.fileBrowser.autosaveEnabled && !autosavePath.empty() && !autosaveTask.running() &&
                glfwGetTime()-lastAutosave>=60) {
                lastAutosave=glfwGetTime();
                bool anyContent=false, busy=false;
                for (const auto& tab:tabs) {
                    busy |= tab->state.fileBrowser.workspaceBusy();
                    anyContent |= !tab->state.structure.atoms.empty();
                }
                if (!busy) {
                    auto saved=captureProject();
                    for (const auto& item:saved) anyContent |= !item.volume.fields.empty();
                    if (anyContent) autosaveTask.start([saved=std::move(saved),autosavePath] {
                        atomforge::saveWorkspace(saved,autosavePath); return true;
                    });
                }
            }
        } catch (const std::exception& error) { state.fileBrowser.showLoadError(error.what()); }

        // A finished large-data conversion opens its cloud.
        {
            std::filesystem::path converted;
            std::string conversionError;
            if (cloudConversion.finished(converted, conversionError))
            {
                if (!converted.empty()) openCloudIntoTab(converted.u8string());
                else if (conversionError != "Cancelled") state.fileBrowser.showLoadError("Cannot convert the large dataset: " + conversionError);
            }
        }
        drawConversionProgress(cloudConversion);
        if (state.hasCloud())
            drawCloudPanel(*state.cloud, frameMilliseconds);

        // Handle File > Open requests: load into current tab (if empty) or new tab
        {
            const std::string openPath = state.fileBrowser.consumePendingOpenPath();
            if (!openPath.empty())
                loadPathIntoTab(openPath);
        }

        // --- Tab bar window (positioned below toolbar) ---
        // Only show the tab bar when at least one structure is loaded, or multiple tabs are open.
        {
            const bool anyTabHasStructure = [&]() {
                for (const auto& t : tabs)
                    if (t->state.hasContent()) return true;
                return false;
            }();
            const bool showTabBar = anyTabHasStructure || (int)tabs.size() > 1;

            float tabBarY = ImGui::GetFrameHeight();
            if (ImGuiWindow* tbWin = ImGui::FindWindowByName("##ViewToolbar"))
                tabBarY = tbWin->Pos.y + tbWin->Size.y;

            const ImGuiIO& io = ImGui::GetIO();
            if (showTabBar)
            {
            ImGui::SetNextWindowPos(ImVec2(0.0f, tabBarY));
            ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, 0.0f));
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(4.0f, 3.0f));
            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing,   ImVec2(4.0f, 2.0f));
            constexpr ImGuiWindowFlags kTabBarFlags =
                ImGuiWindowFlags_NoTitleBar    | ImGuiWindowFlags_NoMove    |
                ImGuiWindowFlags_NoResize      | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoScrollbar |
                ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoDocking;

            if (ImGui::Begin("##StructureTabBar", nullptr, kTabBarFlags))
            {
                constexpr ImGuiTabBarFlags kTBFlags =
                    ImGuiTabBarFlags_Reorderable          |
                    ImGuiTabBarFlags_AutoSelectNewTabs    |
                    ImGuiTabBarFlags_FittingPolicyResizeDown;

                if (ImGui::BeginTabBar("##StructureTabs", kTBFlags))
                {
                    for (int i = 0; i < (int)tabs.size(); ++i)
                    {
                        bool tabOpen = true;
                        ImGuiTabItemFlags itemFlags = 0;
                        if (tabs[i]->state.isModified())
                            itemFlags |= ImGuiTabItemFlags_UnsavedDocument;
                        if (pendingTabSwitch == i)
                        {
                            itemFlags |= ImGuiTabItemFlags_SetSelected;
                            pendingTabSwitch = -1;
                        }

                        // Append ##<index> so ImGui IDs are unique even when
                        // multiple tabs share the same display title ("Untitled").
                        char tabLabel[272];
                        std::snprintf(tabLabel, sizeof(tabLabel), "%s##tab%d",
                                      tabs[i]->title, i);

                        if (ImGui::BeginTabItem(tabLabel, &tabOpen, itemFlags))
                        {
                            if (activeTabIdx != i)
                            {
                                saveCameraToTab(camera, *tabs[activeTabIdx]);
                                activeTabIdx = i;
                                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
                            }
                            ImGui::EndTabItem();
                        }

                        if (!tabOpen)
                        {
                            if (tabs[i]->state.isModified()) closePromptTab = i;
                            else tabs[i]->pendingClose = true;
                        }
                    }

                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
            ImGui::PopStyleVar(2);
            } // if (showTabBar)
        }

        // Unsaved changes: the window title, the close-tab and quit prompts.
        drawUnsavedChangesPrompts(window, tabs, activeTabIdx, pendingTabSwitch, closePromptTab, quitPrompt);

        // Remove closed tabs
        {
            const StructureTab* previousTab = tabs[activeTabIdx].get();
            for (int i = 0; i < (int)tabs.size(); )
            {
                if (tabs[i]->pendingClose)
                {
                    tabs[i]->state.sceneBuffers.destroy();
                    tabs.erase(tabs.begin() + i);
                    if (i < activeTabIdx) --activeTabIdx;
                }
                else { ++i; }
            }
            if (tabs.empty())
            {
                // Last tab closed: keep one empty tab instead of no tabs
                tabs.push_back(std::make_unique<StructureTab>());
                tabs[0]->state.editMenuDialogs.display = savedDisplay;
                initTabResources(*tabs[0], sphere, lowPolyMesh,
                                 billboardMesh, cylinder, renderer);
                restoreCameraFromTab(camera, *tabs[0]);
            }
            activeTabIdx = std::max(0, std::min(activeTabIdx, (int)tabs.size() - 1));
            if (tabs[activeTabIdx].get() != previousTab)
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
        }

        // Re-alias after possible tab vector changes
        EditorState& activeState = tabs[activeTabIdx]->state;

        ImageExportRequest imageExportRequest;
        bool hasImageExportRequest = false;
        mergeFileBrowserRequests(activeState, requests, imageExportRequest, hasImageExportRequest);

        // Close-structure: close tab if multiple, else clear structure
        // (unsaved changes are confirmed first).
        bool closeStructure = activeState.fileBrowser.consumeCloseStructureRequest();
        if (closeStructure && activeState.isModified())
        {
            closePromptTab = activeTabIdx;
            closeStructure = false;
        }
        if (closeStructure)
        {
            if ((int)tabs.size() > 1)
            {
                activeState.sceneBuffers.destroy();
                tabs.erase(tabs.begin() + activeTabIdx);
                activeTabIdx = std::max(0, activeTabIdx - 1);
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
                // activeState refers to the destroyed tab; rebind next frame.
                ImGui::EndFrame();
#ifdef ATOMFORGE_VULKAN
                if (graphics::usingVulkan())
                {
                    // The frame was begun: present it with the background only.
                    const glm::vec4 bg = sceneColours(tabs[activeTabIdx]->state).background;
                    const float clear[4] = {bg.r, bg.g, bg.b, bg.a};
                    vkr::context().beginMainPass(clear);
                    vkr::context().endMainPass();
                    vkr::context().endFrame();
                }
#endif
                continue;
            }
            else
            {
                clearSelection(activeState);
                activeState.structure = Structure();
                activeState.fileBrowser.clearLatticePlanes();
                activeState.fileBrowser.clearMillerDirections();
                updateBuffers(activeState);
                activeState.pendingDefaultViewReset = true;
                std::snprintf(tabs[activeTabIdx]->title, sizeof(tabs[activeTabIdx]->title), "Untitled");
                std::cout << "[Operation] Structure unloaded" << std::endl;
            }
        }
        if (activeState.fileBrowser.consumeResetDefaultViewRequest())
            activeState.pendingDefaultViewReset = true;

        handleUndoRedoRequest(activeState, requests);
        handleAxisViewRequest(camera, requests, activeState.structure);

        if (requests.requestRotateCrystalX || requests.requestRotateCrystalY || requests.requestRotateCrystalZ)
        {
            const double angleDeg = (double)activeState.fileBrowser.getRotateCrystalAngle();
            if (requests.requestRotateCrystalX)      rotateCrystalAroundAxis(camera, 0, angleDeg);
            else if (requests.requestRotateCrystalY) rotateCrystalAroundAxis(camera, 1, angleDeg);
            else                                     rotateCrystalAroundAxis(camera, 2, angleDeg);
        }


        ImDrawList* drawList = ImGui::GetBackgroundDrawList();

        handleBoxSelection(activeState, frame.windowWidth, frame.windowHeight,
                           frame.projection, frame.view, drawList);
        handleLassoSelection(activeState, frame.windowWidth, frame.windowHeight,
                             frame.projection, frame.view, drawList);

        if (uiScript::take("atom.context"))
            activeState.contextMenu.openAt(ImGui::GetIO().DisplaySize.x * 0.55f, ImGui::GetIO().DisplaySize.y * 0.62f);
        if (uiScript::take("atom.substitute"))
            activeState.contextMenu.openSubstitute();
        AtomRequests contextRequests;
        activeState.contextMenu.draw(
            activeState.structure, activeState.sceneBuffers,
            activeState.editMenuDialogs.elementColors,
            activeState.selectedInstanceIndices, contextRequests,
            [&](Structure& structure) { updateBuffers(activeState, structure); });
        mergeContextRequests(requests, contextRequests);

        processMeasurementRequests(
            activeState.measurementState,
            requests.requestMeasureDistance, requests.requestMeasureAngle,
            requests.requestAtomInfo,
            activeState.selectedInstanceIndices, activeState.sceneBuffers,
            activeState.structure, activeState.editMenuDialogs.elementRadii);
        {
            // Results refer to the structure they were measured on: an edit of
            // its content (not a display change) closes them.
            auto& measurement = activeState.measurementState;
            if (requests.requestMeasureDistance || requests.requestMeasureAngle)
                measurement.revision = activeState.undoRedo.currentId();
            else if (measurement.revision != activeState.undoRedo.currentId()) {
                measurement.distanceCardOpen = measurement.angleCardOpen = false;
                measurement.clearVisuals();
            }
        }
        drawMeasurementPopups(activeState.measurementState);
        drawStructureInfoDialog(activeState.structureInfoDialog,
                                requests.requestStructureInfo, activeState.structure);

        if (requests.doDeleteSelected)
            deleteSelectedAtoms(activeState);

        if (requests.requestInsertMidpoint && activeState.selectedInstanceIndices.size() >= 2)
            activeState.contextMenu.openInsertMidpoint();

        refreshSelectionHighlights(activeState);

        drawSelectionOverlay(activeState, drawList, frame.projection, frame.view,
                             frame.windowWidth, frame.windowHeight);
        drawMeasurementOverlays(activeState.measurementState, drawList,
                                frame.projection, frame.view,
                                frame.framebufferWidth, frame.framebufferHeight,
                                activeState.sceneBuffers);
        drawLatticePlanesOverlay(drawList, frame.projection, frame.view,
                                 frame.framebufferWidth, frame.framebufferHeight,
                                 activeState.structure,
                                 activeState.fileBrowser.getLatticePlanes(),
                                 activeState.fileBrowser.isShowLatticePlanesEnabled());
        drawMillerDirectionsOverlay(drawList, frame.projection, frame.view,
                                    frame.framebufferWidth, frame.framebufferHeight,
                                    activeState.structure,
                                    activeState.fileBrowser.getMillerDirections(),
                                    activeState.fileBrowser.isShowMillerDirectionsEnabled());

        if (activeState.fileBrowser.isShowVoronoiEnabled())
        {
            if (activeState.voronoiDirty)
            {
                activeState.voronoiDiagram = computeVoronoi(activeState.structure);
                activeState.voronoiDirty   = false;
            }
            drawVoronoiOverlay(drawList, frame.projection, frame.view,
                               frame.framebufferWidth, frame.framebufferHeight,
                               activeState.voronoiDiagram,
                               activeState.selectedInstanceIndices, true);
        }

        drawPolyhedralOverlay(drawList, frame.projection, frame.view,
                              frame.framebufferWidth, frame.framebufferHeight,
                              activeState.structure,
                              activeState.selectedInstanceIndices,
                              activeState.sceneBuffers.atomIndices,
                              activeState.fileBrowser.getPolyhedralOverlaySettings(),
                              activeState.editMenuDialogs.elementColors,
                              activeState.fileBrowser.isShowPolyhedralViewerEnabled());

        drawSurfaceMeshOverlay(drawList, frame.projection, frame.view, frame.framebufferWidth, frame.framebufferHeight);
        drawTrajectoryLinesOverlay(drawList, frame.projection, frame.view, frame.framebufferWidth, frame.framebufferHeight);

        activeState.fileBrowser.drawInterstitialVoidOverlay(drawList,
                                                            frame.projection,
                                                            frame.view,
                                                            frame.framebufferWidth,
                                                            frame.framebufferHeight,
                                                            activeState.structure);

        if (activeState.fileBrowser.isShowElementEnabled())
            drawElementLabelsOverlay(drawList, frame.projection, frame.view,
                                     frame.framebufferWidth, frame.framebufferHeight,
                                     activeState.sceneBuffers, activeState.structure,
                                     static_cast<AtomLabelContent>(activeState.fileBrowser.getAtomLabelContent()));

        if (!activeState.structure.atoms.empty())
            drawOrientationAxesOverlay(drawList, frame.view,
                                       frame.framebufferWidth, frame.framebufferHeight,
                                       activeState.structure);

        if (activeState.fileBrowser.isShowElementLegendEnabled() && !activeState.hasCloud() &&
            activeState.fileBrowser.getAtomColorMode() == AtomColorMode::ElementType)
        {
            // Below the toolbar and the tab bar.
            float top = ImGui::GetFrameHeight();
            for (const char* name : {"##ViewToolbar", "##StructureTabBar"})
                if (ImGuiWindow* bar = ImGui::FindWindowByName(name); bar && bar->WasActive)
                    top = std::max(top, bar->Pos.y + bar->Size.y);
            const ElementLegendClick clicked = drawElementLegend(top, activeState.structure, activeState.undoRedo.currentId());
            if (!clicked.symbol.empty())
            {
                // Select that element's atoms (with Ctrl, in addition to the selection).
                std::vector<char> chosen(activeState.structure.atoms.size(), 0);
                if (clicked.additive)
                    for (const int instance : activeState.selectedInstanceIndices)
                        if (instance >= 0 && static_cast<std::size_t>(instance) < activeState.sceneBuffers.atomIndices.size()) {
                            const int atom = activeState.sceneBuffers.atomIndices[static_cast<std::size_t>(instance)];
                            if (atom >= 0 && static_cast<std::size_t>(atom) < chosen.size()) chosen[static_cast<std::size_t>(atom)] = 1;
                        }
                std::vector<int> atoms;
                for (std::size_t i = 0; i < activeState.structure.atoms.size(); ++i)
                    if (chosen[i] || activeState.structure.atoms[i].symbol == clicked.symbol) atoms.push_back(static_cast<int>(i));
                activeState.fileBrowser.requestSelection(std::move(atoms));
            }
        }

        if (activeState.fileBrowser.getAtomColorMode() == AtomColorMode::CrystalOrientation)
            drawIPFTriangleLegend(drawList, frame.framebufferWidth, frame.framebufferHeight);
        else if (activeState.fileBrowser.getAtomColorMode() == AtomColorMode::AtomProperty &&
                 activeState.structure.atomProperty.size() == activeState.structure.atoms.size())
        {
            const auto& display = activeState.fileBrowser.getPropertyDisplay();
            const auto range = display.autoRange ? atomforge::science::finiteRange(activeState.structure.atomProperty) : display.range;
            drawAtomPropertyLegend(drawList, frame.framebufferWidth, frame.framebufferHeight,
                                   activeState.structure.atomPropertyName, range.low, range.high,
                                   activeState.structure.atomProperty);
        }

        drawGrabOverlay(activeState, drawList, frame.projection, frame.view,
                        frame.windowWidth, frame.windowHeight);

        ImGui::Render();

        // Apply per-tab lighting/material parameters to the renderer before drawing.
        {
            const atomforge::display::DisplaySettings& display = activeState.editMenuDialogs.display;
            renderer.lightAmbient              = display.ambient;
            renderer.lightSaturation           = display.saturation;
            renderer.lightContrast             = display.contrast;
            renderer.lightShadowStrength       = display.shadowStrength;
            renderer.materialSpecularIntensity = display.specularIntensity;
            renderer.materialShininessScale    = display.shininessScale;
            renderer.materialShininessFloor    = display.shininessFloor;
        }

        // Selection halos: drawn after atoms so nearer atoms hide them.
        std::vector<glm::vec3> selPos;
        std::vector<float>     selRad;
        if (!activeState.selectedInstanceIndices.empty() &&
            !activeState.sceneBuffers.cpuCachesDisabled)
        {
            selPos.reserve(activeState.selectedInstanceIndices.size());
            selRad.reserve(activeState.selectedInstanceIndices.size());
            for (int idx : activeState.selectedInstanceIndices)
            {
                if (idx >= 0 && idx < (int)activeState.sceneBuffers.atomPositions.size())
                {
                    selPos.push_back(activeState.sceneBuffers.atomPositions[idx]);
                    selRad.push_back(idx < (int)activeState.sceneBuffers.atomRadii.size()
                                         ? activeState.sceneBuffers.atomRadii[idx]
                                         : 1.0f);
                }
            }
        }

        // --- 3D scene rendering ---
        if (gpuDrivenScene())
            activeState.sceneBuffers.renderMode = RenderingMode::StandardInstancing;
        else
            tabs[activeTabIdx]->adaptiveRender.update(activeState.sceneBuffers);
#ifdef ATOMFORGE_VULKAN
        if (graphics::usingVulkan())
        {
            // Off-screen work (image export, GIF frames) is recorded before the window's pass.
            setImageExportCloud(activeState.cloud.get());
            handleImageExportIfRequested(hasImageExportRequest, imageExportRequest,
                                         frame, activeState, renderer, shadow);
            captureTrajectoryGifFrame(frame, activeState, renderer, shadow);
            drawSceneVulkan(renderer, frame, activeState.sceneBuffers,
                            activeState.fileBrowser.isShowBondsEnabled(),
                            activeState.fileBrowser.isShowAtomsEnabled(),
                            activeState.fileBrowser.isShowBoundingBoxEnabled(),
                            activeState.fileBrowser.isShowDislocationLinesEnabled(),
                            sceneColours(activeState),
                            selPos, selRad, activeState.cloud.get(),
                            activeState.fileBrowser.isOrthographicViewEnabled());
            saveCameraToTab(camera, *tabs[activeTabIdx]);
            renderImGuiDrawData();
            vkr::context().endMainPass();
            vkr::context().endFrame();
            if (cloudBenchmark && activeState.hasCloud())
            {
                cloudBenchmark->afterFrame(*activeState.cloud, frame.framebufferWidth, frame.framebufferHeight);
                if (cloudBenchmark->done()) glfwSetWindowShouldClose(window, true);
            }
            if (sceneBenchmark && !activeState.hasCloud())
            {
                sceneBenchmark->afterFrame(activeState.sceneBuffers);
                if (sceneBenchmark->done()) glfwSetWindowShouldClose(window, true);
            }
            continue;
        }
#endif
        drawScene(renderer, frame, shadow,
                  activeState.sceneBuffers,
                  activeState.fileBrowser.isShowBondsEnabled(),
                  activeState.fileBrowser.isShowAtomsEnabled(),
                  activeState.fileBrowser.isShowBoundingBoxEnabled(),
                  activeState.fileBrowser.isShowDislocationLinesEnabled(),
                  sceneColours(activeState),
                  activeState.cloud.get(),
                  activeState.fileBrowser.isOrthographicViewEnabled());
        if (!selPos.empty())
            renderer.drawSelectionHalos(frame.projection, frame.view, selPos, selRad);

        setImageExportCloud(activeState.cloud.get());
        handleImageExportIfRequested(hasImageExportRequest, imageExportRequest,
                                     frame, activeState, renderer, shadow);
        captureTrajectoryGifFrame(frame, activeState, renderer, shadow);

        // Save active tab camera back
        saveCameraToTab(camera, *tabs[activeTabIdx]);

        renderImGuiDrawData();
        if (cloudBenchmark && activeState.hasCloud())
        {
            cloudBenchmark->afterFrame(*activeState.cloud, frame.framebufferWidth, frame.framebufferHeight);
            if (cloudBenchmark->done()) glfwSetWindowShouldClose(window, true);
        }
        glfwSwapBuffers(window);
        if (sceneBenchmark && !activeState.hasCloud())
        {
            sceneBenchmark->afterFrame(activeState.sceneBuffers);
            if (sceneBenchmark->done()) glfwSetWindowShouldClose(window, true);
        }
    }

#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        // Everything holding Vulkan resources goes before the device.
        vkDeviceWaitIdle(vkr::context().device());
        tabs.clear();
        shutdownImGui();
        vkr::sceneRenderer().shutdown();
        vkr::context().shutdown();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 0;
    }
#endif
    shutdownImGui();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

