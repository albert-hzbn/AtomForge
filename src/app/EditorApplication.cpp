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
#include "imgui_internal.h"

#include "Camera.h"
#include "CylinderMesh.h"
#include "Renderer.h"
#include "ShadowMap.h"
#include "SphereMesh.h"
#include "LowPolyMesh.h"
#include "BillboardMesh.h"
#include "ui/ImGuiSetup.h"
#include "ui/LatticePlaneOverlay.h"
#include "ui/MillerDirectionOverlay.h"
#include "ui/PolyhedralOverlay.h"
#include "ui/VoronoiOverlay.h"

#include <glm/glm.hpp>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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
};

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

void setTabTitleFromPath(StructureTab& tab, const std::string& path)
{
    if (path.empty()) { std::snprintf(tab.title, sizeof(tab.title), "Untitled"); return; }
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
    tab.state.fileBrowser.initStackingFaultRenderResources(renderer);
    tab.state.fileBrowser.initSubstitutionalSolidSolutionRenderResources(renderer);
    updateBuffers(tab.state);
    tab.state.undoRedo.reset(captureSnapshot(tab.state));
}

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
               const SceneBuffers& sceneBuffers,
               bool showBonds,
               bool showAtoms,
               bool showBoundingBox,
               bool lightTheme)
{
    // Shadow passes first (render into shadow FBO at shadow resolution)
    if (showAtoms)
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

    if (showBonds && showAtoms)
        renderer.drawBondShadowPass(shadow,
            sceneBuffers.tabCylinderVAO, sceneBuffers.tabCylinderVertexCount,
            frame.lightMVP, sceneBuffers.bondCount);

    // Restore viewport to screen size after shadow passes
    glViewport(0, 0, frame.framebufferWidth, frame.framebufferHeight);
    const glm::vec4& bg = lightTheme ? kLightBackground : kDarkBackground;
    glClearColor(bg.r, bg.g, bg.b, bg.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

    // Draw atoms based on rendering mode
    if (showAtoms)
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

    if (showBoundingBox)
    {
    renderer.drawBoxLines(
        frame.projection,
        frame.view,
        sceneBuffers.lineVAO,
        sceneBuffers.boxLines.size(),
        lightTheme ? glm::vec3(0.25f) : glm::vec3(0.85f));
    }
}

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
        state.fileBrowser.isShowAtomsEnabled(),
        state.fileBrowser.isShowBoundingBoxEnabled(),
        state.sceneBuffers,
        renderer,
        shadow,
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

void loadStartupStructureIfRequested(StructureTab& tab, const std::string& startupStructurePath)
{
    EditorState& state = tab.state;
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
    if (!window)
    {
        destroySplashScreen(splash);
        glfwTerminate();
        return -1;
    }

    Camera camera;
    configureCameraCallbacks(window, camera);
    updateSplashScreen(splash, 0.28f, "Initializing UI");

    initImGui(window);

    SphereMesh sphere(40, 40);
    LowPolyMesh lowPolyMesh;
    BillboardMesh billboardMesh;
    CylinderMesh cylinder(32);
    updateSplashScreen(splash, 0.42f, "Preparing scene");

    Renderer renderer;
    renderer.init();
    updateSplashScreen(splash, 0.52f, startupPaths.empty() ? "Preparing workspace" : "Loading structure");

    // Multi-tab state — unique_ptr because EditorState holds non-movable atomics/threads
    std::vector<std::unique_ptr<StructureTab>> tabs;
    tabs.push_back(std::make_unique<StructureTab>());
    int activeTabIdx     = 0;
    int pendingTabSwitch = -1;

    // Pending file drops (shared; routed to active tab each frame)
    std::vector<std::string> appPendingDrops;
    installDropFileCallback(window, appPendingDrops);

    // Initialise the first tab and load the first startup path (if any)
    initTabResources(*tabs[0], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
    loadStartupStructureIfRequested(*tabs[0], startupPaths.empty() ? std::string() : startupPaths[0]);
    restoreCameraFromTab(camera, *tabs[0]);

    // Load any additional startup paths each into their own tab
    for (int si = 1; si < (int)startupPaths.size(); ++si)
    {
        tabs.push_back(std::make_unique<StructureTab>());
        int newIdx = (int)tabs.size() - 1;
        initTabResources(*tabs[newIdx], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
        loadStartupStructureIfRequested(*tabs[newIdx], startupPaths[si]);
    }

    ShadowMap shadow = createShadowMap(1024, 1024);
    updateSplashScreen(splash, 0.92f, "Finalizing");

    showMainWindow(window);
    updateSplashScreen(splash, 1.0f, "Ready");
    destroySplashScreen(splash);

    while (!glfwWindowShouldClose(window))
    {
        // --- Active-tab alias ---
        EditorState& state = tabs[activeTabIdx]->state;

        camera.allowPan    = !state.fileBrowser.isBoxSelectModeEnabled()
                           && !state.fileBrowser.isLassoSelectModeEnabled();
        camera.allowOrbit  = !state.grabState.active;

        glfwPollEvents();

        // Route app-level drops to the active tab
        for (auto& f : appPendingDrops)
            state.pendingDroppedFiles.push_back(std::move(f));
        appPendingDrops.clear();
        processDroppedFiles(state);

        // Helper: load a structure file into a target tab (creates new tab when
        // current tab already has atoms). Returns true on success.
        auto loadPathIntoTab = [&](const std::string& path) -> bool
        {
            if (path.empty()) return false;

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
            if (!tabs[activeTabIdx]->state.structure.atoms.empty())
            {
                saveCameraToTab(camera, *tabs[activeTabIdx]);
                tabs.push_back(std::make_unique<StructureTab>());
                targetIdx = (int)tabs.size() - 1;
                initTabResources(*tabs[targetIdx], sphere, lowPolyMesh,
                                 billboardMesh, cylinder, renderer);
                activeTabIdx     = targetIdx;
                pendingTabSwitch = targetIdx;
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
            }

            EditorState& t = tabs[targetIdx]->state;
            t.structure = std::move(newStructure);
            t.fileBrowser.initFromPath(path);
            t.fileBrowser.applyElementColorOverrides(t.structure);
            t.fileBrowser.showLoadInfo(std::string("Structure loaded. ") + t.structure.ipfLoadStatus);
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
                saveCameraToTab(camera, *tabs[activeTabIdx]);
                tabs.push_back(std::make_unique<StructureTab>());
                int newIdx = (int)tabs.size() - 1;
                initTabResources(*tabs[newIdx], sphere, lowPolyMesh, billboardMesh, cylinder, renderer);
                EditorState& t = tabs[newIdx]->state;
                t.structure = std::move(newStruct);
                updateBuffers(t);
                t.pendingDefaultViewReset = true;
                std::snprintf(tabs[newIdx]->title, sizeof(tabs[newIdx]->title), "Untitled");
                activeTabIdx = newIdx;
                pendingTabSwitch = newIdx;
                restoreCameraFromTab(camera, *tabs[newIdx]);
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
        buildFrameView(camera, state.sceneBuffers, state.fileBrowser.isOrthographicViewEnabled(), frame);

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

        updateImGuiScale(window);
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        FrameActionRequests requests = beginFrameActionRequests(state);
        applyKeyboardShortcuts(state, requests);

        handleGrabMode(state, camera, frame.projection, frame.view,
                       frame.windowWidth, frame.windowHeight);

        state.fileBrowser.draw(
            state.structure,
            state.editMenuDialogs,
            [&](Structure& structure) { updateBuffers(state, structure); },
            [&](Structure s) { state.pendingNewTabStructures.push_back(std::move(s)); },
            state.undoRedo.canUndo(),
            state.undoRedo.canRedo());

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
                    if (!t->state.structure.atoms.empty()) return true;
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
                ImGuiWindowFlags_NoScrollWithMouse;

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
                        if (pendingTabSwitch == i)
                        {
                            itemFlags |= ImGuiTabItemFlags_SetSelected;
                            pendingTabSwitch = -1;
                        }

                        if (ImGui::BeginTabItem(tabs[i]->title, &tabOpen, itemFlags))
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
                            tabs[i]->pendingClose = true;
                    }

                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
            ImGui::PopStyleVar(2);
            } // if (showTabBar)
        }

        // Remove closed tabs
        {
            const int prevActive = activeTabIdx;
            int removed = 0;
            for (int i = 0; i < (int)tabs.size(); )
            {
                if (tabs[i]->pendingClose)
                {
                    tabs[i]->state.sceneBuffers.destroy();
                    tabs.erase(tabs.begin() + i);
                    if (i < prevActive) ++removed;
                }
                else { ++i; }
            }
            if (tabs.empty())
            {
                // Last tab closed: keep one empty tab instead of no tabs
                tabs.push_back(std::make_unique<StructureTab>());
                initTabResources(*tabs[0], sphere, lowPolyMesh,
                                 billboardMesh, cylinder, renderer);
                restoreCameraFromTab(camera, *tabs[0]);
            }
            activeTabIdx = std::max(0, std::min(prevActive - removed, (int)tabs.size() - 1));
        }

        // Re-alias after possible tab vector changes
        EditorState& activeState = tabs[activeTabIdx]->state;

        ImageExportRequest imageExportRequest;
        bool hasImageExportRequest = false;
        mergeFileBrowserRequests(activeState, requests, imageExportRequest, hasImageExportRequest);

        // Close-structure: close tab if multiple, else clear structure
        if (activeState.fileBrowser.consumeCloseStructureRequest())
        {
            if ((int)tabs.size() > 1)
            {
                tabs.erase(tabs.begin() + activeTabIdx);
                activeTabIdx = std::max(0, activeTabIdx - 1);
                restoreCameraFromTab(camera, *tabs[activeTabIdx]);
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
        drawMeasurementPopups(activeState.measurementState);
        drawStructureInfoDialog(activeState.structureInfoDialog,
                                requests.requestStructureInfo, activeState.structure);

        if (requests.doDeleteSelected)
            deleteSelectedAtoms(activeState);

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

        if (activeState.fileBrowser.isShowElementEnabled())
            drawElementLabelsOverlay(drawList, frame.projection, frame.view,
                                     frame.framebufferWidth, frame.framebufferHeight,
                                     activeState.sceneBuffers, activeState.structure);

        drawOrientationAxesOverlay(drawList, frame.view,
                                   frame.framebufferWidth, frame.framebufferHeight);

        if (activeState.fileBrowser.getAtomColorMode() == AtomColorMode::CrystalOrientation)
            drawIPFTriangleLegend(drawList, frame.framebufferWidth, frame.framebufferHeight);

        drawGrabOverlay(activeState, drawList, frame.projection, frame.view,
                        frame.windowWidth, frame.windowHeight);

        ImGui::Render();

        // --- 3D scene rendering ---
        tabs[activeTabIdx]->adaptiveRender.update(activeState.sceneBuffers);
        drawScene(renderer, frame, shadow,
                  activeState.sceneBuffers,
                  activeState.fileBrowser.isShowBondsEnabled(),
                  activeState.fileBrowser.isShowAtomsEnabled(),
                  activeState.fileBrowser.isShowBoundingBoxEnabled(),
                  activeState.fileBrowser.isLightThemeEnabled());

        handleImageExportIfRequested(hasImageExportRequest, imageExportRequest,
                                     frame, activeState, renderer, shadow);

        // Save active tab camera back
        saveCameraToTab(camera, *tabs[activeTabIdx]);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    shutdownImGui();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

