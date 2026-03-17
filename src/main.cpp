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
#include <cstdio>
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
    bool             showDistancePopup = false;
    char             distanceMessage[256] = {0};
    bool             showDistanceLine  = false;
    int              distanceLineIdx0  = -1;
    int              distanceLineIdx1  = -1;
    bool             showAnglePopup    = false;
    char             angleMessage[256] = {0};
    bool             showAngleLines    = false;
    int              angleLineIdx0     = -1;  // first atom
    int              angleLineIdx1     = -1;  // vertex atom
    int              angleLineIdx2     = -1;  // third atom
    bool             showAtomInfoPopup = false;
    char             atomInfoMessage[512] = {0};

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
            showDistanceLine = false;
            showAngleLines = false;
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
        bool requestMeasureDistance = false;
        bool requestMeasureAngle = false;
        bool requestAtomInfo = false;

        if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !selectedInstanceIndices.empty())
            doDeleteSelected = true;

        if (ImGui::IsKeyPressed(ImGuiKey_D) && ImGui::GetIO().KeyCtrl &&
            !selectedInstanceIndices.empty())
        {
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
            showDistanceLine = false;
            showAngleLines = false;
        }

        if (ImGui::IsKeyPressed(ImGuiKey_Escape) && !selectedInstanceIndices.empty())
        {
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
            showDistanceLine = false;
            showAngleLines = false;
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
        requestMeasureDistance = requestMeasureDistance || fileBrowser.consumeMeasureDistanceRequest();
        requestMeasureAngle    = requestMeasureAngle    || fileBrowser.consumeMeasureAngleRequest();
        requestAtomInfo        = requestAtomInfo        || fileBrowser.consumeAtomInfoRequest();

        contextMenu.draw(structure, sceneBuffers,
                         editMenuDialogs.elementColors,
                         selectedInstanceIndices,
                         doDeleteSelected,
                         requestMeasureDistance,
                         requestMeasureAngle,
                         requestAtomInfo,
                         updateBuffers);

        if (requestMeasureDistance)
        {
            if (selectedInstanceIndices.size() != 2)
            {
                std::snprintf(distanceMessage, sizeof(distanceMessage),
                              "Select exactly 2 atoms to measure distance.");
            }
            else
            {
                int idxA = selectedInstanceIndices[0];
                int idxB = selectedInstanceIndices[1];

                bool validA = idxA >= 0 && idxA < (int)sceneBuffers.atomPositions.size();
                bool validB = idxB >= 0 && idxB < (int)sceneBuffers.atomPositions.size();

                if (!validA || !validB)
                {
                    std::snprintf(distanceMessage, sizeof(distanceMessage),
                                  "Unable to measure distance for current selection.");
                    showDistanceLine = false;
                }
                else
                {
                    glm::vec3 pA = sceneBuffers.atomPositions[idxA];
                    glm::vec3 pB = sceneBuffers.atomPositions[idxB];
                    float distance = glm::length(pA - pB);

                    int baseA = (idxA < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idxA] : -1;
                    int baseB = (idxB < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idxB] : -1;
                    const char* symA = (baseA >= 0 && baseA < (int)structure.atoms.size()) ? structure.atoms[baseA].symbol.c_str() : "A";
                    const char* symB = (baseB >= 0 && baseB < (int)structure.atoms.size()) ? structure.atoms[baseB].symbol.c_str() : "B";

                    std::snprintf(distanceMessage, sizeof(distanceMessage),
                                  "Distance %s-%s: %.4f", symA, symB, distance);
                    showDistanceLine = true;
                    distanceLineIdx0 = idxA;
                    distanceLineIdx1 = idxB;
                }
            }
            showDistancePopup = true;
        }

        // ---- Measure Angle ----
        if (requestMeasureAngle)
        {
            if (selectedInstanceIndices.size() != 3)
            {
                std::snprintf(angleMessage, sizeof(angleMessage),
                              "Select exactly 3 atoms to measure angle.");
                showAngleLines = false;
            }
            else
            {
                int idx0 = selectedInstanceIndices[0];
                int idx1 = selectedInstanceIndices[1]; // vertex
                int idx2 = selectedInstanceIndices[2];

                bool ok = idx0 >= 0 && idx0 < (int)sceneBuffers.atomPositions.size() &&
                          idx1 >= 0 && idx1 < (int)sceneBuffers.atomPositions.size() &&
                          idx2 >= 0 && idx2 < (int)sceneBuffers.atomPositions.size();
                if (!ok)
                {
                    std::snprintf(angleMessage, sizeof(angleMessage),
                                  "Unable to measure angle for current selection.");
                    showAngleLines = false;
                }
                else
                {
                    glm::vec3 p0 = sceneBuffers.atomPositions[idx0];
                    glm::vec3 p1 = sceneBuffers.atomPositions[idx1]; // vertex
                    glm::vec3 p2 = sceneBuffers.atomPositions[idx2];

                    glm::vec3 v0 = glm::normalize(p0 - p1);
                    glm::vec3 v2 = glm::normalize(p2 - p1);
                    float cosA = glm::clamp(glm::dot(v0, v2), -1.0f, 1.0f);
                    float angleDeg = glm::degrees(std::acos(cosA));

                    int b0 = (idx0 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx0] : -1;
                    int b1 = (idx1 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx1] : -1;
                    int b2 = (idx2 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx2] : -1;
                    const char* s0 = (b0 >= 0 && b0 < (int)structure.atoms.size()) ? structure.atoms[b0].symbol.c_str() : "A";
                    const char* s1 = (b1 >= 0 && b1 < (int)structure.atoms.size()) ? structure.atoms[b1].symbol.c_str() : "B";
                    const char* s2 = (b2 >= 0 && b2 < (int)structure.atoms.size()) ? structure.atoms[b2].symbol.c_str() : "C";

                    std::snprintf(angleMessage, sizeof(angleMessage),
                                  "Angle %s-%s-%s: %.2f deg", s0, s1, s2, angleDeg);
                    showAngleLines = true;
                    angleLineIdx0  = idx0;
                    angleLineIdx1  = idx1;
                    angleLineIdx2  = idx2;
                }
            }
            showAnglePopup = true;
        }

        if (showAnglePopup)
        {
            ImGui::OpenPopup("Measure Angle");
            showAnglePopup = false;
        }

        ImGui::SetNextWindowSize(ImVec2(540.0f, 0.0f), ImGuiCond_Appearing);
        if (ImGui::BeginPopupModal("Measure Angle", nullptr))
        {
            ImGui::TextWrapped("%s", angleMessage);
            if (ImGui::Button("OK"))
            {
                ImGui::CloseCurrentPopup();
                showAngleLines = false;
            }
            ImGui::EndPopup();
        }

        if (showDistancePopup)
        {
            ImGui::OpenPopup("Measure Distance");
            showDistancePopup = false;
        }

        ImGui::SetNextWindowSize(ImVec2(540.0f, 0.0f), ImGuiCond_Appearing);
        if (ImGui::BeginPopupModal("Measure Distance", nullptr))
        {
            ImGui::TextWrapped("%s", distanceMessage);
            if (ImGui::Button("OK"))
            {
                ImGui::CloseCurrentPopup();
                showDistanceLine = false;
            }
            ImGui::EndPopup();
        }

        // ---- Atom Info ----
        if (requestAtomInfo)
        {
            if (selectedInstanceIndices.size() != 1)
            {
                std::snprintf(atomInfoMessage, sizeof(atomInfoMessage),
                              "Select exactly 1 atom to view info.");
            }
            else
            {
                int idx = selectedInstanceIndices[0];
                bool valid = idx >= 0 && idx < (int)sceneBuffers.atomPositions.size()
                          && idx < (int)sceneBuffers.atomIndices.size();
                if (!valid)
                {
                    std::snprintf(atomInfoMessage, sizeof(atomInfoMessage),
                                  "Unable to retrieve atom info.");
                }
                else
                {
                    int baseIdx = sceneBuffers.atomIndices[idx];
                    glm::vec3 pos = sceneBuffers.atomPositions[idx];
                    int len = 0;

                    if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                    {
                        const AtomSite& atom = structure.atoms[baseIdx];
                        len += std::snprintf(atomInfoMessage + len, sizeof(atomInfoMessage) - len,
                                            "Element:  %s (%s)\n",
                                            elementName(atom.atomicNumber), atom.symbol.c_str());
                        len += std::snprintf(atomInfoMessage + len, sizeof(atomInfoMessage) - len,
                                            "Atomic number:  %d\n", atom.atomicNumber);
                    }

                    len += std::snprintf(atomInfoMessage + len, sizeof(atomInfoMessage) - len,
                                        "Cartesian:  (%.6f, %.6f, %.6f) \xc3\x85\n",
                                        pos.x, pos.y, pos.z);

                    if (structure.hasUnitCell)
                    {
                        glm::mat3 cellMat(
                            glm::vec3((float)structure.cellVectors[0][0],
                                      (float)structure.cellVectors[0][1],
                                      (float)structure.cellVectors[0][2]),
                            glm::vec3((float)structure.cellVectors[1][0],
                                      (float)structure.cellVectors[1][1],
                                      (float)structure.cellVectors[1][2]),
                            glm::vec3((float)structure.cellVectors[2][0],
                                      (float)structure.cellVectors[2][1],
                                      (float)structure.cellVectors[2][2]));
                        glm::vec3 origin((float)structure.cellOffset[0],
                                         (float)structure.cellOffset[1],
                                         (float)structure.cellOffset[2]);
                        glm::vec3 frac = glm::inverse(cellMat) * (pos - origin);
                        std::snprintf(atomInfoMessage + len, sizeof(atomInfoMessage) - len,
                                      "Direct:  (%.6f, %.6f, %.6f)",
                                      frac.x, frac.y, frac.z);
                    }
                }
            }
            showAtomInfoPopup = true;
        }

        if (showAtomInfoPopup)
        {
            ImGui::OpenPopup("Atom Info");
            showAtomInfoPopup = false;
        }

        ImGui::SetNextWindowSize(ImVec2(480.0f, 0.0f), ImGuiCond_Appearing);
        if (ImGui::BeginPopupModal("Atom Info", nullptr))
        {
            ImGui::TextWrapped("%s", atomInfoMessage);
            if (ImGui::Button("OK"))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
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

        // Draw dotted line between the two atoms measured last (only when active).
        if (showDistanceLine)
        {
            int idxA = distanceLineIdx0;
            int idxB = distanceLineIdx1;
            bool validA = idxA >= 0 && idxA < (int)sceneBuffers.atomPositions.size();
            bool validB = idxB >= 0 && idxB < (int)sceneBuffers.atomPositions.size();

            if (validA && validB)
            {
                auto projectToScreen = [&](const glm::vec3& p, float& sx, float& sy) -> bool {
                    glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
                    if (clip.w <= 0.0f)
                        return false;

                    float invW = 1.0f / clip.w;
                    float ndcX = clip.x * invW;
                    float ndcY = clip.y * invW;
                    float ndcZ = clip.z * invW;
                    if (ndcX < -1.0f || ndcX > 1.0f || ndcY < -1.0f || ndcY > 1.0f || ndcZ < -1.0f || ndcZ > 1.0f)
                        return false;

                    sx = (ndcX * 0.5f + 0.5f) * (float)w;
                    sy = (1.0f - (ndcY * 0.5f + 0.5f)) * (float)h;
                    return true;
                };

                float ax = 0.0f, ay = 0.0f, bx = 0.0f, by = 0.0f;
                if (projectToScreen(sceneBuffers.atomPositions[idxA], ax, ay) &&
                    projectToScreen(sceneBuffers.atomPositions[idxB], bx, by))
                {
                    float dx = bx - ax;
                    float dy = by - ay;
                    float len = std::sqrt(dx * dx + dy * dy);
                    if (len > 1.0f)
                    {
                        const float dashLen = 7.0f;
                        const float gapLen = 5.0f;
                        float ux = dx / len;
                        float uy = dy / len;

                        for (float t = 0.0f; t < len; t += (dashLen + gapLen))
                        {
                            float t2 = std::min(t + dashLen, len);
                            ImVec2 p0(ax + ux * t,  ay + uy * t);
                            ImVec2 p1(ax + ux * t2, ay + uy * t2);
                            drawList->AddLine(p0, p1, IM_COL32(255, 255, 80, 230), 2.0f);
                        }
                    }
                }
            }
        }

        // Draw dotted angle legs (atom0-vertex and atom2-vertex) when active.
        if (showAngleLines)
        {
            bool ok = angleLineIdx0 >= 0 && angleLineIdx0 < (int)sceneBuffers.atomPositions.size() &&
                      angleLineIdx1 >= 0 && angleLineIdx1 < (int)sceneBuffers.atomPositions.size() &&
                      angleLineIdx2 >= 0 && angleLineIdx2 < (int)sceneBuffers.atomPositions.size();
            if (ok)
            {
                auto projectPt = [&](const glm::vec3& p, float& sx, float& sy) -> bool {
                    glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
                    if (clip.w <= 0.0f) return false;
                    float invW = 1.0f / clip.w;
                    float nx = clip.x * invW, ny = clip.y * invW, nz = clip.z * invW;
                    if (nx < -1.0f || nx > 1.0f || ny < -1.0f || ny > 1.0f || nz < -1.0f || nz > 1.0f) return false;
                    sx = (nx * 0.5f + 0.5f) * (float)w;
                    sy = (1.0f - (ny * 0.5f + 0.5f)) * (float)h;
                    return true;
                };

                auto drawDashedLine = [&](ImVec2 a, ImVec2 b, ImU32 col) {
                    float dx = b.x - a.x, dy = b.y - a.y;
                    float len = std::sqrt(dx*dx + dy*dy);
                    if (len < 1.0f) return;
                    float ux = dx/len, uy = dy/len;
                    for (float t = 0.0f; t < len; t += 12.0f)
                    {
                        float t2 = std::min(t + 7.0f, len);
                        drawList->AddLine(ImVec2(a.x+ux*t, a.y+uy*t),
                                          ImVec2(a.x+ux*t2, a.y+uy*t2), col, 2.0f);
                    }
                };

                float x0,y0, x1,y1, x2,y2;
                bool v0ok = projectPt(sceneBuffers.atomPositions[angleLineIdx0], x0, y0);
                bool v1ok = projectPt(sceneBuffers.atomPositions[angleLineIdx1], x1, y1);
                bool v2ok = projectPt(sceneBuffers.atomPositions[angleLineIdx2], x2, y2);

                ImU32 col = IM_COL32(80, 220, 255, 230);
                if (v0ok && v1ok) drawDashedLine({x0,y0}, {x1,y1}, col);
                if (v2ok && v1ok) drawDashedLine({x2,y2}, {x1,y1}, col);

                // Draw short arc at the vertex in screen space.
                if (v0ok && v1ok && v2ok)
                {
                    float d0x = x0-x1, d0y = y0-y1;
                    float d2x = x2-x1, d2y = y2-y1;
                    float len0 = std::sqrt(d0x*d0x + d0y*d0y);
                    float len2 = std::sqrt(d2x*d2x + d2y*d2y);
                    float r = std::min({len0, len2, 30.0f}) * 0.35f;
                    if (r > 4.0f)
                    {
                        float a0 = std::atan2(d0y, d0x);
                        float a2 = std::atan2(d2y, d2x);
                        // Normalise so we sweep the short way.
                        float da = a2 - a0;
                        while (da >  (float)M_PI) da -= 2.0f*(float)M_PI;
                        while (da < -(float)M_PI) da += 2.0f*(float)M_PI;
                        const int segs = 20;
                        for (int s = 0; s < segs; ++s)
                        {
                            float t0 = a0 + da * ((float)s       / segs);
                            float t1a = a0 + da * ((float)(s + 1) / segs);
                            drawList->AddLine(
                                ImVec2(x1 + r*std::cos(t0),  y1 + r*std::sin(t0)),
                                ImVec2(x1 + r*std::cos(t1a), y1 + r*std::sin(t1a)),
                                col, 1.5f);
                        }
                    }
                }
            }
        }

        if (fileBrowser.isShowElementEnabled())
        {
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

