#include "app/EditorApplication.h"

#include "app/FileDropHandler.h"
#include "app/EditorOps.h"
#include "app/SceneView.h"
#include "app/EditorState.h"
#include "app/ImageExport.h"
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

#include <iostream>
#include <string>

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
    installDropFileCallback(window, state);

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

        ImageExportRequest imageExportRequest;
        const bool hasImageExportRequest = state.fileBrowser.consumeImageExportRequest(imageExportRequest);

        requests.requestUndo = requests.requestUndo || state.fileBrowser.consumeUndoRequest();
        requests.requestRedo = requests.requestRedo || state.fileBrowser.consumeRedoRequest();
        requests.requestStructureInfo = requests.requestStructureInfo || state.fileBrowser.consumeStructureInfoRequest();

        if (state.fileBrowser.consumeCloseStructureRequest())
        {
            clearSelection(state);
            state.structure = Structure();
            updateBuffers(state);
            state.pendingDefaultViewReset = true;

            std::cout << "[Operation] Structure unloaded" << std::endl;
        }

        if (state.fileBrowser.consumeResetDefaultViewRequest())
            state.pendingDefaultViewReset = true;

        if (requests.requestUndo && state.undoRedo.canUndo())
        {
            applySnapshot(state, state.undoRedo.undo());
            std::cout << "[Operation] Undo" << std::endl;
        }
        else if (requests.requestRedo && state.undoRedo.canRedo())
        {
            applySnapshot(state, state.undoRedo.redo());
            std::cout << "[Operation] Redo" << std::endl;
        }

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

        const glm::vec4 sceneBackgroundColor(0.09f, 0.11f, 0.15f, 1.0f);
        glViewport(0, 0, frame.framebufferWidth, frame.framebufferHeight);
        glClearColor(sceneBackgroundColor.r, sceneBackgroundColor.g, sceneBackgroundColor.b, sceneBackgroundColor.a);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (state.fileBrowser.isShowBondsEnabled())
        {
            renderer.drawBonds(
                frame.projection,
                frame.view,
                frame.lightPosition,
                frame.cameraPosition,
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

        if (hasImageExportRequest)
        {
            ImageExportView exportView;
            exportView.width = frame.framebufferWidth;
            exportView.height = frame.framebufferHeight;
            exportView.projection = frame.projection;
            exportView.view = frame.view;
            exportView.lightMVP = frame.lightMVP;
            exportView.lightPosition = frame.lightPosition;
            exportView.cameraPosition = frame.cameraPosition;

            std::string exportError;
            const bool exportOk = exportStructureImage(
                imageExportRequest,
                exportView,
                sceneBackgroundColor,
                state.fileBrowser.isShowBondsEnabled(),
                state.sceneBuffers,
                renderer,
                shadow,
                sphere,
                cylinder,
                exportError);

            if (exportOk)
            {
                std::cout << "[Operation] Image exported: " << imageExportRequest.outputPath << std::endl;
            }
            else
            {
                const std::string message = exportError.empty()
                    ? "Failed to export image."
                    : ("Failed to export image: " + exportError);
                state.fileBrowser.showLoadError(message);
                std::cout << "[Operation] Image export failed: " << imageExportRequest.outputPath
                          << " (" << (exportError.empty() ? "Unknown error" : exportError) << ")" << std::endl;
            }
        }

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    shutdownImGui();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}