#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#include "Camera.h"
#include "CylinderMesh.h"
#include "SphereMesh.h"
#include "ShadowMap.h"
#include "StructureLoader.h"
#include "StructureInstanceBuilder.h"
#include "SceneBuffers.h"
#include "Renderer.h"
#include "graphics/Picking.h"
#include "ui/FileBrowser.h"
#include "ui/ImGuiSetup.h"
#include "ui/EditMenuDialogs.h"
#include "ui/AtomContextMenu.h"
#include "ui/MeasurementOverlay.h"
#include "ui/StructureInfoDialog.h"
#include "UndoRedo.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cmath>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
int main()
{
    // ----------------------------------------------------------------
    // Window
    // ----------------------------------------------------------------

    if (!glfwInit())
        return -1;

    GLFWwindow* window =
        glfwCreateWindow(1280, 800, "Atoms Editor", nullptr, nullptr);

    if (!window)
        return -1;

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
        return -1;

    glEnable(GL_DEPTH_TEST);

    // ----------------------------------------------------------------
    // Camera
    // ----------------------------------------------------------------

    Camera camera;
    Camera::instance = &camera;

    glfwSetMouseButtonCallback(window, Camera::mouseButton);
    glfwSetCursorPosCallback(window,   Camera::cursor);
    glfwSetScrollCallback(window,      Camera::scroll);

    // ----------------------------------------------------------------
    // ImGui
    // ----------------------------------------------------------------

    initImGui(window);

    // ----------------------------------------------------------------
    // Geometry, structure, UI
    // ----------------------------------------------------------------

    SphereMesh sphere(40, 40);
    CylinderMesh cylinder(32);

    std::string filename = "";
    Structure structure = loadStructure(filename);

    FileBrowser fileBrowser;
    fileBrowser.initFromPath(filename);

    // ----------------------------------------------------------------
    // GPU resources
    // ----------------------------------------------------------------

    SceneBuffers sceneBuffers;
    sceneBuffers.init(sphere.vao, cylinder.vao);

    Renderer renderer;
    renderer.init();

    ShadowMap shadow = createShadowMap(1024, 1024);

    // ----------------------------------------------------------------
    // Picking / selection state
    // ----------------------------------------------------------------

    std::vector<int> selectedInstanceIndices;
    EditMenuDialogs  editMenuDialogs;
    AtomContextMenu  contextMenu;
    MeasurementOverlayState measurementState;
    StructureInfoDialogState structureInfoDialog;
    UndoRedoManager undoRedo;
    bool suppressHistoryCommit = false;

    // ----------------------------------------------------------------
    // Buffer update helper
    // ----------------------------------------------------------------

    auto captureSnapshot = [&]() -> EditorSnapshot {
        EditorSnapshot snapshot;
        snapshot.structure = structure;
        snapshot.elementRadii = editMenuDialogs.elementRadii;
        snapshot.elementColors = editMenuDialogs.elementColors;
        snapshot.elementShininess = editMenuDialogs.elementShininess;
        return snapshot;
    };

    auto updateBuffers = [&](Structure& s) {
        // Treat Transform Structure as an explicit supercell build step.
        // This prevents subsequent atom edits (insert/substitute/delete)
        // from being echoed to symmetry-equivalent render replicas.
        if (fileBrowser.isTransformMatrixEnabled() && s.hasUnitCell)
        {
            s = buildSupercell(s, fileBrowser.getTransformMatrix());
            fileBrowser.clearTransformMatrix();
        }

        for (auto& atom : s.atoms)
        {
            int z = atom.atomicNumber;
            if (z >= 0 && z < (int)editMenuDialogs.elementColors.size())
            {
                atom.r = editMenuDialogs.elementColors[z].r;
                atom.g = editMenuDialogs.elementColors[z].g;
                atom.b = editMenuDialogs.elementColors[z].b;
            }
        }

        StructureInstanceData data = buildStructureInstanceData(
            s,
            fileBrowser.isTransformMatrixEnabled(),
            fileBrowser.getTransformMatrix(),
            editMenuDialogs.elementRadii,
            editMenuDialogs.elementShininess);

        sceneBuffers.upload(data);
        selectedInstanceIndices.clear();

        if (!suppressHistoryCommit)
            undoRedo.commit(captureSnapshot());
    };

    updateBuffers(structure);
    undoRedo.reset(captureSnapshot());

    bool pendingDefaultViewReset = true;

    auto applyDefaultView = [&](int viewportW, int viewportH, bool fitToStructure) {
        constexpr float kIsoYawDeg = 45.0f;
        constexpr float kIsoPitchDeg = 35.2643897f;

        camera.yaw = kIsoYawDeg;
        camera.pitch = kIsoPitchDeg;

        if (!fitToStructure || sceneBuffers.atomPositions.empty())
        {
            camera.distance = 10.0f;
            return;
        }

        float maxRadius = 0.0f;
        for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
        {
            float radius = (i < sceneBuffers.atomRadii.size()) ? sceneBuffers.atomRadii[i] : 0.0f;
            float dist = glm::length(sceneBuffers.atomPositions[i] - sceneBuffers.orbitCenter) + radius;
            maxRadius = std::max(maxRadius, dist);
        }

        maxRadius = std::max(maxRadius, 1.0f);

        float aspect = (viewportH > 0) ? (float)viewportW / (float)viewportH : 1.0f;
        float vFov = glm::radians(45.0f);
        float hFov = 2.0f * std::atan(std::tan(vFov * 0.5f) * aspect);
        float halfFov = 0.5f * std::min(vFov, hFov);
        halfFov = std::max(halfFov, glm::radians(10.0f));

        float framedDistance = (maxRadius / std::sin(halfFov)) * 1.15f;
        camera.distance = std::max(2.0f, std::min(500.0f, framedDistance));
    };

    // ----------------------------------------------------------------
    // Render loop
    // ----------------------------------------------------------------

    // Clears the atom selection and resets all measurement overlays.
    auto clearSelection = [&]() {
        for (int idx : selectedInstanceIndices)
            sceneBuffers.restoreAtomColor(idx);
        selectedInstanceIndices.clear();
        measurementState.clearVisuals();
    };

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        if (w == 0 || h == 0) { glfwSwapBuffers(window); continue; }

        if (pendingDefaultViewReset)
        {
            applyDefaultView(w, h, true);
            pendingDefaultViewReset = false;
        }

        // ------------------------------------------------------------
        // Matrices  (computed early so picking can use them)
        // ------------------------------------------------------------

        glm::mat4 projection =
            glm::perspective(glm::radians(45.0f), (float)w / h, 0.1f, 1000.0f);

        float yaw   = glm::radians(camera.yaw);
        float pitch = glm::radians(camera.pitch);

        glm::vec3 camOffset(
            camera.distance * std::cos(pitch) * std::sin(yaw),
            camera.distance * std::sin(pitch),
            camera.distance * std::cos(pitch) * std::cos(yaw));

        glm::vec3 camPos = sceneBuffers.orbitCenter + camOffset;

        glm::mat4 view =
            glm::lookAt(camPos, sceneBuffers.orbitCenter, glm::vec3(0, 1, 0));

        // use window size (logical pixels) to match GLFW cursor coordinates
        int winW, winH;
        glfwGetWindowSize(window, &winW, &winH);
        if (winW == 0) winW = w;
        if (winH == 0) winH = h;

        glm::mat4 lightProj =
            glm::perspective(glm::radians(45.0f), 1.0f, 0.1f, 1000.0f);

        glm::vec3 lightPos = sceneBuffers.orbitCenter + glm::vec3(40.0f, 40.0f, 40.0f);
        glm::mat4 lightMVP =
            lightProj * glm::lookAt(lightPos, sceneBuffers.orbitCenter, glm::vec3(0, 1, 0));

        // ------------------------------------------------------------
        // Atom picking  (left click without drag)
        // ------------------------------------------------------------

        if (camera.pendingClick)
        {
            camera.pendingClick = false;
            glm::vec3 ray = pickRayDir(camera.clickX, camera.clickY,
                                       winW, winH, projection, view);
            int newIdx = pickAtom(camPos, ray,
                                  sceneBuffers.atomPositions,
                                  sceneBuffers.atomRadii,
                                  1.0f);
            
            if (newIdx >= 0)
            {
                // Check if Ctrl is held for multi-select
                bool ctrlHeld = ImGui::GetIO().KeyCtrl;
                
                if (ctrlHeld)
                {
                    // Multi-select: toggle this atom in selection
                    auto it = std::find(selectedInstanceIndices.begin(), 
                                       selectedInstanceIndices.end(), newIdx);
                    if (it != selectedInstanceIndices.end())
                    {
                        // Already selected, deselect it
                        sceneBuffers.restoreAtomColor(newIdx);
                        selectedInstanceIndices.erase(it);
                    }
                    else
                    {
                        // Not selected, add it
                        selectedInstanceIndices.push_back(newIdx);
                    }
                }
                else
                {
                    // Regular click: replace selection
                    // Restore colors of previously selected atoms
                    clearSelection();
                    selectedInstanceIndices.push_back(newIdx);
                }
            }
            else
            {
                // Clicked on empty space: clear selection
                clearSelection();
            }
        }

        // ------------------------------------------------------------
        // Right-click  → open context menu if an atom is selected
        // ------------------------------------------------------------

        if (camera.pendingRightClick)
        {
            camera.pendingRightClick = false;
            if (!selectedInstanceIndices.empty())
                contextMenu.open();
        }

        // ------------------------------------------------------------
        // ImGui frame
        // ------------------------------------------------------------

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ------------------------------------------------------------
        // Keyboard shortcuts
        // ------------------------------------------------------------

        bool doDeleteSelected = false;
        bool requestMeasureDistance = fileBrowser.consumeMeasureDistanceRequest();
        bool requestMeasureAngle    = fileBrowser.consumeMeasureAngleRequest();
        bool requestAtomInfo        = fileBrowser.consumeAtomInfoRequest();
        bool requestStructureInfo   = fileBrowser.consumeStructureInfoRequest();

        if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !selectedInstanceIndices.empty())
            doDeleteSelected = true;

        if (ImGui::IsKeyPressed(ImGuiKey_D) && ImGui::GetIO().KeyCtrl &&
            !selectedInstanceIndices.empty())
        {
            clearSelection();
        }

        if (ImGui::IsKeyPressed(ImGuiKey_Escape) && !selectedInstanceIndices.empty())
        {
            clearSelection();
        }

        if (ImGui::IsKeyPressed(ImGuiKey_A) && ImGui::GetIO().KeyCtrl &&
            !structure.atoms.empty())
        {
            selectedInstanceIndices.clear();
            for (int i = 0; i < (int)sceneBuffers.atomIndices.size(); ++i)
                selectedInstanceIndices.push_back(i);
        }

        if (ImGui::IsKeyPressed(ImGuiKey_O) && ImGui::GetIO().KeyCtrl)
            fileBrowser.openFileDialog();

        if (ImGui::IsKeyPressed(ImGuiKey_S) && ImGui::GetIO().KeyCtrl)
            fileBrowser.saveFileDialog();

        // ------------------------------------------------------------
        // UI modules
        // ------------------------------------------------------------

        bool requestUndo = false;
        bool requestRedo = false;

        if (ImGui::IsKeyPressed(ImGuiKey_Z) && ImGui::GetIO().KeyCtrl && !ImGui::GetIO().KeyShift)
            requestUndo = true;
        if ((ImGui::IsKeyPressed(ImGuiKey_Y) && ImGui::GetIO().KeyCtrl) ||
            (ImGui::IsKeyPressed(ImGuiKey_Z) && ImGui::GetIO().KeyCtrl && ImGui::GetIO().KeyShift))
            requestRedo = true;

        fileBrowser.draw(structure, editMenuDialogs, updateBuffers,
                         undoRedo.canUndo(), undoRedo.canRedo());

        requestUndo = requestUndo || fileBrowser.consumeUndoRequest();
        requestRedo = requestRedo || fileBrowser.consumeRedoRequest();
        requestStructureInfo = requestStructureInfo || fileBrowser.consumeStructureInfoRequest();

        if (fileBrowser.consumeResetDefaultViewRequest())
            pendingDefaultViewReset = true;

        auto applySnapshot = [&](const EditorSnapshot& snapshot) {
            structure = snapshot.structure;
            editMenuDialogs.elementRadii = snapshot.elementRadii;
            editMenuDialogs.elementColors = snapshot.elementColors;
            editMenuDialogs.elementShininess = snapshot.elementShininess;
            suppressHistoryCommit = true;
            updateBuffers(structure);
            suppressHistoryCommit = false;
            measurementState.clearVisuals();
        };

        if (requestUndo && undoRedo.canUndo())
        {
            applySnapshot(undoRedo.undo());
        }
        else if (requestRedo && undoRedo.canRedo())
        {
            applySnapshot(undoRedo.redo());
        }

        AtomRequests ctxReq;
        contextMenu.draw(structure, sceneBuffers,
                 editMenuDialogs.elementColors,
                 selectedInstanceIndices,
                 ctxReq,
                 updateBuffers);
        doDeleteSelected       = doDeleteSelected       || ctxReq.doDelete;
        requestMeasureDistance = requestMeasureDistance || ctxReq.measureDistance;
        requestMeasureAngle    = requestMeasureAngle    || ctxReq.measureAngle;
        requestAtomInfo        = requestAtomInfo        || ctxReq.atomInfo;

        processMeasurementRequests(measurementState,
                                   requestMeasureDistance,
                                   requestMeasureAngle,
                                   requestAtomInfo,
                                   selectedInstanceIndices,
                                   sceneBuffers,
                                   structure);
        drawMeasurementPopups(measurementState);
        drawStructureInfoDialog(structureInfoDialog, requestStructureInfo, structure);

        // Handle deletion of selected atoms
        if (doDeleteSelected && !selectedInstanceIndices.empty())
        {
            // Collect base indices and sort in reverse order for safe deletion
            std::vector<int> baseIndicesToDelete;
            for (int selectedIdx : selectedInstanceIndices)
            {
                if (selectedIdx >= 0 && selectedIdx < (int)sceneBuffers.atomIndices.size())
                {
                    int baseIdx = sceneBuffers.atomIndices[selectedIdx];
                    if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                    {
                        baseIndicesToDelete.push_back(baseIdx);
                    }
                }
            }
            
            // Sort in reverse order to delete from end to start (avoid index shifting)
            std::sort(baseIndicesToDelete.begin(), baseIndicesToDelete.end(), std::greater<int>());
            
            // Remove duplicates
            baseIndicesToDelete.erase(
                std::unique(baseIndicesToDelete.begin(), baseIndicesToDelete.end()),
                baseIndicesToDelete.end()
            );

            // Delete atoms in reverse order
            for (int baseIdx : baseIndicesToDelete)
            {
                if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                {
                    structure.atoms.erase(structure.atoms.begin() + baseIdx);
                }
            }

            // Rebuild the scene with updated atom list
            updateBuffers(structure);
        }

        // Apply / maintain highlight for all selected atoms
        for (int& idx : selectedInstanceIndices)
        {
            if (idx >= (int)sceneBuffers.atomCount)
            {
                sceneBuffers.restoreAtomColor(idx);
                idx = -1;
            }
            else if (idx >= 0)
            {
                sceneBuffers.highlightAtom(idx, glm::vec3(1.0f, 1.0f, 0.0f));
            }
        }
        // Remove -1 entries (failed atoms)
        selectedInstanceIndices.erase(
            std::remove(selectedInstanceIndices.begin(), selectedInstanceIndices.end(), -1),
            selectedInstanceIndices.end()
        );

        ImDrawList* drawList = ImGui::GetForegroundDrawList();

        drawMeasurementOverlays(measurementState,
                                drawList,
                                projection,
                                view,
                                w,
                                h,
                                sceneBuffers);

        if (fileBrowser.isShowElementEnabled())
        {
            drawElementLabelsOverlay(drawList,
                                     projection,
                                     view,
                                     w,
                                     h,
                                     sceneBuffers,
                                     structure);
        }

        ImGui::Render();

        // ------------------------------------------------------------
        // Draw
        // ------------------------------------------------------------
        glViewport(0, 0, w, h);
        glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (fileBrowser.isShowBondsEnabled())
        {
            renderer.drawBonds(projection, view, lightMVP,
                               shadow, cylinder, sceneBuffers.bondCount);
        }

        renderer.drawAtoms(projection, view, lightMVP,
                           lightPos, camPos,
                           shadow, sphere, sceneBuffers.atomCount);

        renderer.drawBoxLines(projection, view,
                              sceneBuffers.lineVAO,
                              sceneBuffers.boxLines.size());

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // ----------------------------------------------------------------
    // Cleanup
    // ----------------------------------------------------------------

    shutdownImGui();
    glfwTerminate();
}

