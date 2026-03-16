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
#include "ui/FileBrowser.h"
#include "ui/ImGuiSetup.h"
#include "ui/PeriodicTableDialog.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>

#include <cmath>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Ray-picking helpers
// ---------------------------------------------------------------------------

static glm::vec3 pickRayDir(double mx, double my, int w, int h,
                             const glm::mat4& projection,
                             const glm::mat4& view)
{
    float ndcX =  2.0f * (float)mx / (float)w - 1.0f;
    float ndcY = -2.0f * (float)my / (float)h + 1.0f;

    glm::mat4 invVP = glm::inverse(projection * view);

    glm::vec4 nearP = invVP * glm::vec4(ndcX, ndcY, -1.0f, 1.0f);
    glm::vec4 farP  = invVP * glm::vec4(ndcX, ndcY,  1.0f, 1.0f);

    nearP /= nearP.w;
    farP  /= farP.w;

    return glm::normalize(glm::vec3(farP - nearP));
}

static int pickAtom(const glm::vec3& origin, const glm::vec3& dir,
                    const std::vector<glm::vec3>& positions, float radius)
{
    int   best  = -1;
    float bestT = 1e30f;

    for (int i = 0; i < (int)positions.size(); ++i)
    {
        glm::vec3 oc = positions[i] - origin;
        float t  = glm::dot(oc, dir);
        if (t < 0.0f) continue;
        float d2 = glm::dot(oc, oc) - t * t;
        if (d2 < radius * radius && t < bestT)
        {
            best  = i;
            bestT = t;
        }
    }
    return best;
}

enum class PeriodicAction
{
    None,
    Substitute,
    InsertMidpoint,
};

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
    // Picking / selection state (declared before lambda so it is captured)
    // ----------------------------------------------------------------

    std::vector<int> selectedInstanceIndices;  // Multiple selected atoms
    bool openContextMenu     = false;
    PeriodicAction pendingPeriodicAction = PeriodicAction::None;

    // ----------------------------------------------------------------
    // Buffer update helper
    // ----------------------------------------------------------------

    auto updateBuffers = [&](const Structure& s) {
        StructureInstanceData data = buildStructureInstanceData(
            s,
            fileBrowser.isTransformMatrixEnabled(),
            fileBrowser.getTransformMatrix());

        sceneBuffers.upload(data);
        selectedInstanceIndices.clear();  // clear selection whenever structure rebuilds
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
                                  sceneBuffers.atomPositions, 1.0f);
            
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
                openContextMenu = true;
        }

        // ------------------------------------------------------------
        // ImGui frame
        // ------------------------------------------------------------

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ------------------------------------------------------------
        bool doDeleteSelected = false;
        bool doOpenPeriodicTable = false;
        // Keyboard shortcuts
        // ------------------------------------------------------------

        if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !selectedInstanceIndices.empty())
        {
            doDeleteSelected = true;
        }

        if (ImGui::IsKeyPressed(ImGuiKey_D) && ImGui::GetIO().KeyCtrl && !selectedInstanceIndices.empty())
        {
            // Ctrl+D: Deselect all
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
        }

        if (ImGui::IsKeyPressed(ImGuiKey_Escape) && !selectedInstanceIndices.empty())
        {
            // Escape: Deselect all
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
        }

        // Ctrl+A: Select all atoms
        if (ImGui::IsKeyPressed(ImGuiKey_A) && ImGui::GetIO().KeyCtrl && structure.atoms.size() > 0)
        {
            selectedInstanceIndices.clear();
            for (int i = 0; i < (int)sceneBuffers.atomIndices.size(); ++i)
            {
                selectedInstanceIndices.push_back(i);
            }
        }

        // Ctrl+O: Open file dialog
        if (ImGui::IsKeyPressed(ImGuiKey_O) && ImGui::GetIO().KeyCtrl)
        {
            fileBrowser.openFileDialog();
        }

        fileBrowser.draw(structure, updateBuffers);

        // Context menu (must be outside any Begin/End window block)
        if (openContextMenu)
        {
            ImGui::OpenPopup("##atomCtx");
            openContextMenu = false;
        }
        if (ImGui::BeginPopup("##atomCtx"))
        {
            if (ImGui::MenuItem("Substitute Atom..."))
            {
                pendingPeriodicAction = PeriodicAction::Substitute;
                doOpenPeriodicTable = true;
            }

            bool canInsertAtMidpoint = selectedInstanceIndices.size() >= 2;
            if (!canInsertAtMidpoint)
                ImGui::BeginDisabled();
            if (ImGui::MenuItem("Insert Atom at Midpoint..."))
            {
                pendingPeriodicAction = PeriodicAction::InsertMidpoint;
                doOpenPeriodicTable = true;
            }
            if (!canInsertAtMidpoint)
                ImGui::EndDisabled();

            if (ImGui::MenuItem("Delete"))
                doDeleteSelected = true;
            if (ImGui::MenuItem("Deselect"))
            {
                for (int idx : selectedInstanceIndices)
                    sceneBuffers.restoreAtomColor(idx);
                selectedInstanceIndices.clear();
            }
            ImGui::EndPopup();
        }
        // Open periodic table at top-level (outside any popup scope)
        if (doOpenPeriodicTable)
            openPeriodicTable();

        // Periodic table picker + apply selected action
        {
            std::vector<ElementSelection> selections;
            if (drawPeriodicTable(selections))
            {
                if (!selections.empty())
                {
                    const auto& sel = selections[0];

                    if (pendingPeriodicAction == PeriodicAction::Substitute && !selectedInstanceIndices.empty())
                    {
                        // Apply selected element to all selected atoms
                        for (int selectedIdx : selectedInstanceIndices)
                        {
                            if (selectedIdx >= 0 && selectedIdx < (int)sceneBuffers.atomIndices.size())
                            {
                                int baseIdx = sceneBuffers.atomIndices[selectedIdx];
                                if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                                {
                                    structure.atoms[baseIdx].symbol      = sel.symbol;
                                    structure.atoms[baseIdx].atomicNumber = sel.atomicNumber;
                                    float r, g, b;
                                    getDefaultElementColor(sel.atomicNumber, r, g, b);
                                    structure.atoms[baseIdx].r = r;
                                    structure.atoms[baseIdx].g = g;
                                    structure.atoms[baseIdx].b = b;
                                }
                            }
                        }
                        updateBuffers(structure);
                    }

                    if (pendingPeriodicAction == PeriodicAction::InsertMidpoint && selectedInstanceIndices.size() >= 2)
                    {
                        // Insert one new atom at the centroid of all selected atoms.
                        double sumX = 0.0;
                        double sumY = 0.0;
                        double sumZ = 0.0;
                        int validCount = 0;

                        for (int selectedIdx : selectedInstanceIndices)
                        {
                            if (selectedIdx >= 0 && selectedIdx < (int)sceneBuffers.atomIndices.size())
                            {
                                int baseIdx = sceneBuffers.atomIndices[selectedIdx];
                                if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                                {
                                    sumX += structure.atoms[baseIdx].x;
                                    sumY += structure.atoms[baseIdx].y;
                                    sumZ += structure.atoms[baseIdx].z;
                                    ++validCount;
                                }
                            }
                        }

                        if (validCount >= 2)
                        {
                            AtomSite newAtom;
                            newAtom.symbol = sel.symbol;
                            newAtom.atomicNumber = sel.atomicNumber;
                            newAtom.x = sumX / (double)validCount;
                            newAtom.y = sumY / (double)validCount;
                            newAtom.z = sumZ / (double)validCount;
                            getDefaultElementColor(sel.atomicNumber, newAtom.r, newAtom.g, newAtom.b);
                            structure.atoms.push_back(newAtom);
                            updateBuffers(structure);
                        }
                    }
                }

                pendingPeriodicAction = PeriodicAction::None;
            }
        }

        if (pendingPeriodicAction != PeriodicAction::None &&
            !ImGui::IsPopupOpen("Periodic Table##picker"))
        {
            pendingPeriodicAction = PeriodicAction::None;
        }

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

        ImGui::Render();

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

