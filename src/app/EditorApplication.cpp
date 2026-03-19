#include "app/EditorApplication.h"

#include "app/EditorOps.h"
#include "app/EditorState.h"
#include "app/InteractionHandlers.h"
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
#include "ui/ImGuiSetup.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cmath>
#include <string>

namespace {

float estimateSceneRadius(const SceneBuffers& sceneBuffers)
{
    float maxRadius = 0.0f;

    if (!sceneBuffers.atomPositions.empty())
    {
        for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
        {
            const float radius = (i < sceneBuffers.atomRadii.size()) ? sceneBuffers.atomRadii[i] : 0.0f;
            const float distance = glm::length(sceneBuffers.atomPositions[i] - sceneBuffers.orbitCenter) + radius;
            maxRadius = std::max(maxRadius, distance);
        }
    }
    else if (!sceneBuffers.boxLines.empty())
    {
        for (const glm::vec3& point : sceneBuffers.boxLines)
        {
            const float distance = glm::length(point - sceneBuffers.orbitCenter);
            maxRadius = std::max(maxRadius, distance);
        }
    }

    return std::max(maxRadius, 1.0f);
}

struct FrameView
{
    int framebufferWidth = 0;
    int framebufferHeight = 0;
    int windowWidth = 0;
    int windowHeight = 0;

    glm::mat4 projection = glm::mat4(1.0f);
    glm::mat4 view = glm::mat4(1.0f);
    glm::mat4 lightMVP = glm::mat4(1.0f);
    glm::vec3 lightPosition = glm::vec3(0.0f);
    glm::vec3 cameraPosition = glm::vec3(0.0f);
};

bool updateViewport(GLFWwindow* window, FrameView& frame)
{
    glfwGetFramebufferSize(window, &frame.framebufferWidth, &frame.framebufferHeight);
    if (frame.framebufferWidth == 0 || frame.framebufferHeight == 0)
    {
        glfwSwapBuffers(window);
        return false;
    }

    glfwGetWindowSize(window, &frame.windowWidth, &frame.windowHeight);
    if (frame.windowWidth == 0)
        frame.windowWidth = frame.framebufferWidth;
    if (frame.windowHeight == 0)
        frame.windowHeight = frame.framebufferHeight;

    return true;
}

void buildFrameView(Camera& camera,
                    const SceneBuffers& sceneBuffers,
                    bool useOrthographicView,
                    FrameView& frame)
{
    const float aspect = (float)frame.framebufferWidth / (float)frame.framebufferHeight;
    const float verticalFov = glm::radians(45.0f);
    const float sceneRadius = estimateSceneRadius(sceneBuffers);
    const float depthPadding = std::max(10.0f, sceneRadius * 0.25f);

    if (useOrthographicView)
    {
        const float halfHeight = std::max(0.1f, camera.distance * std::tan(verticalFov * 0.5f));
        const float halfWidth = halfHeight * aspect;
        const float depthRange = std::max(1000.0f, camera.distance + sceneRadius + depthPadding);
        frame.projection = glm::ortho(
            -halfWidth,
            halfWidth,
            -halfHeight,
            halfHeight,
            -depthRange,
            depthRange);
    }
    else
    {
        const float nearestSurface = camera.distance - sceneRadius;
        const float nearClip = (nearestSurface > 0.0f)
            ? std::max(0.01f, nearestSurface * 0.25f)
            : 0.01f;
        const float farClip = std::max(nearClip + 100.0f, camera.distance + sceneRadius + depthPadding);

        frame.projection = glm::perspective(
            verticalFov,
            aspect,
            nearClip,
            farClip);
    }

    float yaw = glm::radians(camera.yaw);
    float pitch = glm::radians(camera.pitch);

    glm::vec3 cameraOffset(
        camera.distance * std::cos(pitch) * std::sin(yaw),
        camera.distance * std::sin(pitch),
        camera.distance * std::cos(pitch) * std::cos(yaw));

    frame.cameraPosition = sceneBuffers.orbitCenter + cameraOffset;
    frame.view = glm::lookAt(frame.cameraPosition, sceneBuffers.orbitCenter, glm::vec3(0, 1, 0));

    glm::mat4 lightProjection = glm::perspective(glm::radians(45.0f), 1.0f, 0.1f, 1000.0f);
    frame.lightPosition = sceneBuffers.orbitCenter + glm::vec3(40.0f, 40.0f, 40.0f);
    frame.lightMVP = lightProjection * glm::lookAt(
        frame.lightPosition,
        sceneBuffers.orbitCenter,
        glm::vec3(0, 1, 0));
}

void dropFileCallback(GLFWwindow* window, int count, const char** paths)
{
    if (count <= 0 || paths == nullptr)
        return;

    EditorState* state = static_cast<EditorState*>(glfwGetWindowUserPointer(window));
    if (!state)
        return;

    for (int i = 0; i < count; ++i)
    {
        if (paths[i] == nullptr || paths[i][0] == '\0')
            continue;
        state->pendingDroppedFiles.push_back(paths[i]);
    }
}

void processDroppedFiles(EditorState& state)
{
    if (state.pendingDroppedFiles.empty())
        return;

    const std::string droppedFile = state.pendingDroppedFiles.back();
    state.pendingDroppedFiles.clear();

    Structure loadedStructure;
    std::string loadError;
    if (!loadStructureFromFile(droppedFile, loadedStructure, loadError))
    {
        state.fileBrowser.showLoadError(loadError);
        return;
    }

    state.structure = std::move(loadedStructure);
    state.fileBrowser.initFromPath(droppedFile);
    state.fileBrowser.applyElementColorOverrides(state.structure);
    updateBuffers(state);
    state.pendingDefaultViewReset = true;
}

} // namespace

int runAtomsEditor()
{
    GLFWwindow* window = createMainWindow();
    if (!window)
        return -1;

    Camera camera;
    configureCameraCallbacks(window, camera);

    initImGui(window);

    SphereMesh sphere(40, 40);
    CylinderMesh cylinder(32);

    std::string filename;
    EditorState state;
    state.structure = loadStructure(filename);
    state.fileBrowser.initFromPath(filename);

    glfwSetWindowUserPointer(window, &state);
    glfwSetDropCallback(window, dropFileCallback);

    state.sceneBuffers.init(sphere.vao, cylinder.vao);

    Renderer renderer;
    renderer.init();

    ShadowMap shadow = createShadowMap(1024, 1024);

    updateBuffers(state);
    state.undoRedo.reset(captureSnapshot(state));

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        processDroppedFiles(state);

        FrameView frame;
        if (!updateViewport(window, frame))
            continue;

        if (state.pendingDefaultViewReset)
        {
            applyDefaultView(camera, state.sceneBuffers, frame.framebufferWidth, frame.framebufferHeight, true);
            state.pendingDefaultViewReset = false;
        }

        buildFrameView(camera,
                   state.sceneBuffers,
                   state.fileBrowser.isOrthographicViewEnabled(),
                   frame);

        handlePendingAtomPick(
            camera,
            state,
            frame.cameraPosition,
            frame.windowWidth,
            frame.windowHeight,
            frame.projection,
            frame.view);
        handleRightClick(camera, state);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        FrameActionRequests requests = beginFrameActionRequests(state);
        applyKeyboardShortcuts(state, requests);

        state.fileBrowser.draw(
            state.structure,
            state.editMenuDialogs,
            [&](Structure& structure) { updateBuffers(state, structure); },
            state.undoRedo.canUndo(),
            state.undoRedo.canRedo());

        requests.requestUndo = requests.requestUndo || state.fileBrowser.consumeUndoRequest();
        requests.requestRedo = requests.requestRedo || state.fileBrowser.consumeRedoRequest();
        requests.requestStructureInfo = requests.requestStructureInfo || state.fileBrowser.consumeStructureInfoRequest();

        if (state.fileBrowser.consumeResetDefaultViewRequest())
            state.pendingDefaultViewReset = true;

        if (requests.requestUndo && state.undoRedo.canUndo())
            applySnapshot(state, state.undoRedo.undo());
        else if (requests.requestRedo && state.undoRedo.canRedo())
            applySnapshot(state, state.undoRedo.redo());

        ImDrawList* drawList = ImGui::GetForegroundDrawList();
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
        requests.doDeleteSelected = requests.doDeleteSelected || contextRequests.doDelete;
        requests.requestMeasureDistance = requests.requestMeasureDistance || contextRequests.measureDistance;
        requests.requestMeasureAngle = requests.requestMeasureAngle || contextRequests.measureAngle;
        requests.requestAtomInfo = requests.requestAtomInfo || contextRequests.atomInfo;

        processMeasurementRequests(
            state.measurementState,
            requests.requestMeasureDistance,
            requests.requestMeasureAngle,
            requests.requestAtomInfo,
            state.selectedInstanceIndices,
            state.sceneBuffers,
            state.structure);
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

        ImGui::Render();

        glViewport(0, 0, frame.framebufferWidth, frame.framebufferHeight);
        glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (state.fileBrowser.isShowBondsEnabled())
        {
            renderer.drawBonds(
                frame.projection,
                frame.view,
                cylinder,
                state.sceneBuffers.bondCount);
        }

        renderer.drawAtoms(
            frame.projection,
            frame.view,
            frame.lightMVP,
            frame.lightPosition,
            frame.cameraPosition,
            shadow,
            sphere,
            state.sceneBuffers.atomCount);

        renderer.drawBoxLines(
            frame.projection,
            frame.view,
            state.sceneBuffers.lineVAO,
            state.sceneBuffers.boxLines.size());

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    shutdownImGui();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}