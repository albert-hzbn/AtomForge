#include "app/EditorApplication.h"

#include "app/FileDropHandler.h"
#include "app/EditorOps.h"
#include "app/SceneView.h"
#include "app/EditorState.h"
#include "app/ImageExport.h"
#include "app/InteractionHandlers.h"
#include "app/StructureFileService.h"
#include "app/WindowSetup.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#include "Camera.h"
#include "CylinderMesh.h"
#include "Renderer.h"
#include "ShadowMap.h"
#include "SphereMesh.h"
#include "LowPolyMesh.h"
#include "BillboardMesh.h"
#include "ui/ImGuiSetup.h"
#include "ui/LatticePlaneOverlay.h"
#include "ui/VoronoiOverlay.h"

#include <glm/glm.hpp>

#include <cmath>
#include <iostream>
#include <string>

namespace
{
const glm::vec4 kDarkBackground(0.09f, 0.11f, 0.15f, 1.0f);
const glm::vec4 kLightBackground(0.92f, 0.93f, 0.94f, 1.0f);

// Dynamically selects rendering mode based on measured frame performance.
// Starts from an initial estimate (set by SceneBuffers::upload based on atom
// count) and adjusts up or down if the GPU can handle a better mode or is
// struggling with the current one.
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

        double dtMs = (now - lastTime) * 1000.0;
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

void applyPendingDefaultViewReset(Camera& camera, EditorState& state, const FrameView& frame)
{
    if (!state.pendingDefaultViewReset)
        return;

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

void handleStructureResetRequests(EditorState& state)
{
    if (state.fileBrowser.consumeCloseStructureRequest())
    {
        clearSelection(state);
        state.structure = Structure();
        state.fileBrowser.clearLatticePlanes();
        updateBuffers(state);
        state.pendingDefaultViewReset = true;

        std::cout << "[Operation] Structure unloaded" << std::endl;
    }

    if (state.fileBrowser.consumeResetDefaultViewRequest())
        state.pendingDefaultViewReset = true;
}

void handleUndoRedoRequest(EditorState& state, const FrameActionRequests& requests)
{
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

void drawScene(Renderer& renderer,
               const FrameView& frame,
               const ShadowMap& shadow,
               const SphereMesh& sphere,
               const LowPolyMesh& lowPolyMesh,
               const BillboardMesh& billboardMesh,
               const CylinderMesh& cylinder,
               const SceneBuffers& sceneBuffers,
               bool showBonds,
               bool lightTheme)
{
    // Shadow passes first (render into shadow FBO at shadow resolution)
    switch (sceneBuffers.renderMode)
    {
        case RenderingMode::StandardInstancing:
            renderer.drawShadowPass(shadow, sphere, frame.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::LowPolyInstancing:
            renderer.drawShadowPassLowPoly(shadow, lowPolyMesh, frame.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::BillboardImposters:
            renderer.drawShadowPassBillboard(shadow, billboardMesh, frame.lightMVP, frame.view, sceneBuffers.atomCount);
            break;
    }

    if (showBonds)
        renderer.drawBondShadowPass(shadow, cylinder, frame.lightMVP, sceneBuffers.bondCount);

    // Restore viewport to screen size after shadow passes
    glViewport(0, 0, frame.framebufferWidth, frame.framebufferHeight);
    const glm::vec4& bg = lightTheme ? kLightBackground : kDarkBackground;
    glClearColor(bg.r, bg.g, bg.b, bg.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (showBonds)
    {
        renderer.drawBonds(
            frame.projection,
            frame.view,
            frame.lightPosition,
            frame.cameraPosition,
            cylinder,
            sceneBuffers.bondCount);
    }

    // Draw atoms based on rendering mode
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
                sphere,
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
                lowPolyMesh,
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
                billboardMesh,
                sceneBuffers.atomCount);
            break;
    }

    renderer.drawBoxLines(
        frame.projection,
        frame.view,
        sceneBuffers.lineVAO,
        sceneBuffers.boxLines.size(),
        lightTheme ? glm::vec3(0.25f) : glm::vec3(0.85f));
}

void handleImageExportIfRequested(bool hasImageExportRequest,
                                  const ImageExportRequest& imageExportRequest,
                                  const FrameView& frame,
                                  EditorState& state,
                                  Renderer& renderer,
                                  ShadowMap& shadow,
                                  SphereMesh& sphere,
                                  LowPolyMesh& lowPolyMesh,
                                  BillboardMesh& billboardMesh,
                                  CylinderMesh& cylinder)
{
    if (!hasImageExportRequest)
        return;

    ImageExportView exportView;
    const int scale = std::max(1, imageExportRequest.resolutionScale);
    exportView.width = frame.framebufferWidth * scale;
    exportView.height = frame.framebufferHeight * scale;
    exportView.projection = frame.projection;
    exportView.view = frame.view;
    exportView.lightMVP = frame.lightMVP;
    exportView.lightPosition = frame.lightPosition;
    exportView.cameraPosition = frame.cameraPosition;

    std::string exportError;
    const bool exportOk = exportStructureImage(
        imageExportRequest,
        exportView,
        state.fileBrowser.isLightThemeEnabled() ? kLightBackground : kDarkBackground,
        state.fileBrowser.isShowBondsEnabled(),
        state.sceneBuffers,
        renderer,
        shadow,
        sphere,
        lowPolyMesh,
        billboardMesh,
        cylinder,
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

void loadStartupStructureIfRequested(EditorState& state, const std::string& startupStructurePath)
{
    if (startupStructurePath.empty())
    {
        state.fileBrowser.initFromPath("");
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
    state.fileBrowser.showLoadInfo(std::string("Structure loaded. ") + state.structure.ipfLoadStatus);
    state.pendingDefaultViewReset = true;

    std::cout << "[Operation] Startup-loaded structure: " << startupStructurePath
              << " (atoms=" << state.structure.atoms.size() << ")" << std::endl;
}
} // namespace

int runAtomsEditor(const std::string& startupStructurePath)
{
    GLFWwindow* window = createMainWindow();
    if (!window)
        return -1;

    Camera camera;
    configureCameraCallbacks(window, camera);

    initImGui(window);

    SphereMesh sphere(40, 40);
    LowPolyMesh lowPolyMesh;      // 12-facet icosahedron for 100k-10M atoms
    BillboardMesh billboardMesh;  // Quad billboard for 10M+ atoms
    CylinderMesh cylinder(32);

    EditorState state;
    loadStartupStructureIfRequested(state, startupStructurePath);
    installDropFileCallback(window, state);

    state.sceneBuffers.init(sphere.vao, lowPolyMesh.vao, billboardMesh.vao, cylinder.vao);

    Renderer renderer;
    renderer.init();
    state.fileBrowser.initNanoCrystalRenderResources(renderer);
    state.fileBrowser.initCustomStructureRenderResources(renderer);
    state.fileBrowser.initMergeStructuresRenderResources(renderer);
    state.fileBrowser.initInterfaceBuilderRenderResources(renderer);
    state.fileBrowser.initCSLGrainBoundaryRenderResources(renderer);
    state.fileBrowser.initPolyCrystalRenderResources(renderer);
    state.fileBrowser.initSubstitutionalSolidSolutionRenderResources(renderer);

    ShadowMap shadow = createShadowMap(1024, 1024);

    updateBuffers(state);
    state.undoRedo.reset(captureSnapshot(state));

    AdaptiveRenderState adaptiveRender;

    while (!glfwWindowShouldClose(window))
    {
        camera.allowPan = !state.fileBrowser.isBoxSelectModeEnabled();
        camera.allowOrbit = !state.grabState.active;

        glfwPollEvents();
        processDroppedFiles(state);

        FrameView frame;
        if (!updateViewport(window, frame))
            continue;

        applyPendingDefaultViewReset(camera, state, frame);

        buildFrameView(camera,
                   state.sceneBuffers,
                   state.fileBrowser.isOrthographicViewEnabled(),
                   frame);

        // Suppress camera orbit and atom picking during grab mode
        if (state.grabState.active)
        {
            camera.pendingClick = false;
            camera.pendingRightClick = false;
        }
        else
        {
            handlePendingAtomPick(
                camera,
                state,
                frame.cameraPosition,
                frame.windowWidth,
                frame.windowHeight,
                frame.projection,
                frame.view);
            handleRightClick(camera, state);
        }

        updateImGuiScale(window);
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        FrameActionRequests requests = beginFrameActionRequests(state);
        applyKeyboardShortcuts(state, requests);

        handleGrabMode(
            state,
            camera,
            frame.projection,
            frame.view,
            frame.windowWidth,
            frame.windowHeight);

        state.fileBrowser.draw(
            state.structure,
            state.editMenuDialogs,
            [&](Structure& structure) { updateBuffers(state, structure); },
            state.undoRedo.canUndo(),
            state.undoRedo.canRedo());

        ImageExportRequest imageExportRequest;
        bool hasImageExportRequest = false;
        mergeFileBrowserRequests(state, requests, imageExportRequest, hasImageExportRequest);
        handleStructureResetRequests(state);
        handleUndoRedoRequest(state, requests);
        handleAxisViewRequest(camera, requests, state.structure);

        if (requests.requestRotateCrystalX || requests.requestRotateCrystalY || requests.requestRotateCrystalZ)
        {
            const double angleDeg = (double)state.fileBrowser.getRotateCrystalAngle();
            if (requests.requestRotateCrystalX)
                rotateCrystalAroundAxis(camera, 0, angleDeg);
            else if (requests.requestRotateCrystalY)
                rotateCrystalAroundAxis(camera, 1, angleDeg);
            else
                rotateCrystalAroundAxis(camera, 2, angleDeg);
        }

        // Keep scene overlays visible over the viewport but underneath UI popups/dialogs.
        ImDrawList* drawList = ImGui::GetBackgroundDrawList();
        handleBoxSelection(
            state,
            frame.windowWidth,
            frame.windowHeight,
            frame.projection,
            frame.view,
            drawList);

        AtomRequests contextRequests;
        state.contextMenu.draw(
            state.structure,
            state.sceneBuffers,
            state.editMenuDialogs.elementColors,
            state.selectedInstanceIndices,
            contextRequests,
            [&](Structure& structure) { updateBuffers(state, structure); });
        mergeContextRequests(requests, contextRequests);

        processMeasurementRequests(
            state.measurementState,
            requests.requestMeasureDistance,
            requests.requestMeasureAngle,
            requests.requestAtomInfo,
            state.selectedInstanceIndices,
            state.sceneBuffers,
            state.structure,
            state.editMenuDialogs.elementRadii);
        drawMeasurementPopups(state.measurementState);
        drawStructureInfoDialog(state.structureInfoDialog, requests.requestStructureInfo, state.structure);

        if (requests.doDeleteSelected)
            deleteSelectedAtoms(state);

        refreshSelectionHighlights(state);

        drawMeasurementOverlays(
            state.measurementState,
            drawList,
            frame.projection,
            frame.view,
            frame.framebufferWidth,
            frame.framebufferHeight,
            state.sceneBuffers);

        drawLatticePlanesOverlay(
            drawList,
            frame.projection,
            frame.view,
            frame.framebufferWidth,
            frame.framebufferHeight,
            state.structure,
            state.fileBrowser.getLatticePlanes(),
            state.fileBrowser.isShowLatticePlanesEnabled());

        if (state.fileBrowser.isShowVoronoiEnabled())
        {
            if (state.voronoiDirty)
            {
                state.voronoiDiagram = computeVoronoi(state.structure);
                state.voronoiDirty = false;
            }
            drawVoronoiOverlay(
                drawList,
                frame.projection,
                frame.view,
                frame.framebufferWidth,
                frame.framebufferHeight,
                state.voronoiDiagram,
                state.selectedInstanceIndices,
                true);
        }

        if (state.fileBrowser.isShowElementEnabled())
        {
            drawElementLabelsOverlay(
                drawList,
                frame.projection,
                frame.view,
                frame.framebufferWidth,
                frame.framebufferHeight,
                state.sceneBuffers,
                state.structure);
        }

        drawOrientationAxesOverlay(
            drawList,
            frame.view,
            frame.framebufferWidth,
            frame.framebufferHeight);

        // Draw IPF triangle legend when Crystal Orientation coloring is active
        if (state.fileBrowser.getAtomColorMode() == AtomColorMode::CrystalOrientation)
        {
            drawIPFTriangleLegend(
                drawList,
                frame.framebufferWidth,
                frame.framebufferHeight);
        }

        // Draw grab mode overlay with real-time atom coordinates
        drawGrabOverlay(
            state,
            drawList,
            frame.projection,
            frame.view,
            frame.windowWidth,
            frame.windowHeight);

        ImGui::Render();

        // Adaptive rendering: measure frame performance and adjust mode.
        adaptiveRender.update(state.sceneBuffers);

        drawScene(
            renderer,
            frame,
            shadow,
            sphere,
            lowPolyMesh,
            billboardMesh,
            cylinder,
            state.sceneBuffers,
            state.fileBrowser.isShowBondsEnabled(),
            state.fileBrowser.isLightThemeEnabled());

        handleImageExportIfRequested(
            hasImageExportRequest,
            imageExportRequest,
            frame,
            state,
            renderer,
            shadow,
            sphere,
            lowPolyMesh,
            billboardMesh,
            cylinder);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    shutdownImGui();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}