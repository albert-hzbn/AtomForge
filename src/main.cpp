#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#include "Camera.h"
#include "SphereMesh.h"
#include "ShadowMap.h"
#include "StructureLoader.h"
#include "StructureInstanceBuilder.h"
#include "SceneBuffers.h"
#include "Renderer.h"
#include "ElementData.h"
#include "graphics/Picking.h"
#include "ui/FileBrowser.h"
#include "ui/ImGuiSetup.h"
#include "ui/EditMenuDialogs.h"
#include "ui/AtomContextMenu.h"

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

    std::string filename = "";
    Structure structure = loadStructure(filename);

    FileBrowser fileBrowser;
    fileBrowser.initFromPath(filename);

    // ----------------------------------------------------------------
    // GPU resources
    // ----------------------------------------------------------------

    SceneBuffers sceneBuffers;
    sceneBuffers.init(sphere.vao);

    Renderer renderer;
    renderer.init();

    ShadowMap shadow = createShadowMap(1024, 1024);

    // ----------------------------------------------------------------
    // Picking / selection state
    // ----------------------------------------------------------------

    std::vector<int> selectedInstanceIndices;
    EditMenuDialogs  editMenuDialogs;
    AtomContextMenu  contextMenu;

    // ----------------------------------------------------------------
    // Buffer update helper
    // ----------------------------------------------------------------

    auto updateBuffers = [&](Structure& s) {
        // Treat Transform Atoms as an explicit supercell build step.
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
            editMenuDialogs.elementRadii);

        sceneBuffers.upload(data);
        selectedInstanceIndices.clear();
    };

    updateBuffers(structure);

    // ----------------------------------------------------------------
    // Render loop
    // ----------------------------------------------------------------

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        if (w == 0 || h == 0) { glfwSwapBuffers(window); continue; }

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
                    for (int idx : selectedInstanceIndices)
                        sceneBuffers.restoreAtomColor(idx);
                    
                    selectedInstanceIndices.clear();
                    selectedInstanceIndices.push_back(newIdx);
                }
            }
            else
            {
                // Clicked on empty space: clear selection
                for (int idx : selectedInstanceIndices)
                    sceneBuffers.restoreAtomColor(idx);
                selectedInstanceIndices.clear();
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

        if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !selectedInstanceIndices.empty())
            doDeleteSelected = true;

        if (ImGui::IsKeyPressed(ImGuiKey_D) && ImGui::GetIO().KeyCtrl &&
            !selectedInstanceIndices.empty())
        {
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
        }

        if (ImGui::IsKeyPressed(ImGuiKey_Escape) && !selectedInstanceIndices.empty())
        {
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
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

        fileBrowser.draw(structure, editMenuDialogs, updateBuffers);

        contextMenu.draw(structure, sceneBuffers,
                         editMenuDialogs.elementColors,
                         selectedInstanceIndices,
                         doDeleteSelected,
                         updateBuffers);

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

        if (fileBrowser.isShowElementEnabled())
        {
            ImDrawList* drawList = ImGui::GetForegroundDrawList();
            for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
            {
                int baseIdx = (i < sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[i] : -1;
                if (baseIdx < 0 || baseIdx >= (int)structure.atoms.size())
                    continue;

                const glm::vec3& p = sceneBuffers.atomPositions[i];
                glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
                if (clip.w <= 0.0f)
                    continue;

                float invW = 1.0f / clip.w;
                float ndcX = clip.x * invW;
                float ndcY = clip.y * invW;
                float ndcZ = clip.z * invW;
                if (ndcX < -1.0f || ndcX > 1.0f || ndcY < -1.0f || ndcY > 1.0f || ndcZ < -1.0f || ndcZ > 1.0f)
                    continue;

                float sx = (ndcX * 0.5f + 0.5f) * (float)w;
                float sy = (1.0f - (ndcY * 0.5f + 0.5f)) * (float)h;

                const std::string& label = structure.atoms[baseIdx].symbol;
                ImVec2 textSize = ImGui::CalcTextSize(label.c_str());
                drawList->AddText(ImVec2(sx - textSize.x * 0.5f, sy - textSize.y * 0.5f),
                                  IM_COL32(255, 255, 255, 255),
                                  label.c_str());
            }
        }

        ImGui::Render();

        // ------------------------------------------------------------
        // Draw
        // ------------------------------------------------------------
        renderer.drawShadowPass(shadow, sphere, lightMVP, sceneBuffers.atomCount);

        glViewport(0, 0, w, h);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.drawAtoms(projection, view, lightMVP,
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

