#include "app/InteractionHandlers.h"

#include "Camera.h"
#include "app/EditorOps.h"
#include "graphics/Picking.h"
#include "math/StructureMath.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_set>

namespace
{
bool isInsideSelectionRect(const ImVec2& p, const ImVec2& a, const ImVec2& b)
{
    const float minX = std::min(a.x, b.x);
    const float maxX = std::max(a.x, b.x);
    const float minY = std::min(a.y, b.y);
    const float maxY = std::max(a.y, b.y);
    return p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY;
}

bool isCtrlHeld()
{
    return ImGui::GetIO().KeyCtrl;
}

void selectAllInstances(EditorState& state)
{
    state.selectedInstanceIndices.clear();
    for (int i = 0; i < (int)state.sceneBuffers.atomIndices.size(); ++i)
        state.selectedInstanceIndices.push_back(i);
}

void toggleSelectedInstance(EditorState& state, int pickedIndex)
{
    auto it = std::find(
        state.selectedInstanceIndices.begin(),
        state.selectedInstanceIndices.end(),
        pickedIndex);
    if (it != state.selectedInstanceIndices.end())
    {
        state.sceneBuffers.restoreAtomColor(pickedIndex);
        state.selectedInstanceIndices.erase(it);
        return;
    }

    state.selectedInstanceIndices.push_back(pickedIndex);
}

void selectSingleInstance(EditorState& state, int pickedIndex)
{
    clearSelection(state);
    state.selectedInstanceIndices.push_back(pickedIndex);
}

bool isDragLargeEnough(const ImVec2& start, const ImVec2& end)
{
    const float width = std::abs(end.x - start.x);
    const float height = std::abs(end.y - start.y);
    return width >= 4.0f && height >= 4.0f;
}

void drawSelectionRect(ImDrawList* drawList, const ImVec2& start, const ImVec2& end)
{
    const ImVec2 minCorner(std::min(start.x, end.x), std::min(start.y, end.y));
    const ImVec2 maxCorner(std::max(start.x, end.x), std::max(start.y, end.y));
    drawList->AddRectFilled(minCorner, maxCorner, IM_COL32(80, 160, 255, 45));
    drawList->AddRect(minCorner, maxCorner, IM_COL32(80, 160, 255, 220), 0.0f, 0, 1.5f);
}
}

FrameActionRequests beginFrameActionRequests(EditorState& state)
{
    FrameActionRequests requests;
    requests.requestMeasureDistance = state.fileBrowser.consumeMeasureDistanceRequest();
    requests.requestMeasureAngle = state.fileBrowser.consumeMeasureAngleRequest();
    requests.requestAtomInfo = state.fileBrowser.consumeAtomInfoRequest();
    requests.requestStructureInfo = state.fileBrowser.consumeStructureInfoRequest();
    return requests;
}

void applyKeyboardShortcuts(EditorState& state, FrameActionRequests& requests)
{
    const bool ctrlHeld = isCtrlHeld();

    // When grab mode is active, only grab-related keys should be processed.
    if (state.grabState.active)
        return;

    if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !state.selectedInstanceIndices.empty())
        requests.doDeleteSelected = true;

    if (ImGui::IsKeyPressed(ImGuiKey_D) && ctrlHeld && !state.selectedInstanceIndices.empty())
        clearSelection(state);

    if (ImGui::IsKeyPressed(ImGuiKey_Escape) && !state.selectedInstanceIndices.empty())
        clearSelection(state);

    if (ImGui::IsKeyPressed(ImGuiKey_A) && ctrlHeld && !state.structure.atoms.empty())
        selectAllInstances(state);

    if (ImGui::IsKeyPressed(ImGuiKey_O) && ctrlHeld)
        state.fileBrowser.openFileDialog();

    if (ImGui::IsKeyPressed(ImGuiKey_S) && ctrlHeld && ImGui::GetIO().KeyShift)
        state.fileBrowser.exportImageDialog();

    if (ImGui::IsKeyPressed(ImGuiKey_S) && ctrlHeld && !ImGui::GetIO().KeyShift)
        state.fileBrowser.saveFileDialog();

    if (ImGui::IsKeyPressed(ImGuiKey_S) && !ctrlHeld && !state.selectedInstanceIndices.empty()
        && !state.sceneBuffers.cpuCachesDisabled)
        state.contextMenu.openSubstitute();

    if (ImGui::IsKeyPressed(ImGuiKey_W) && ctrlHeld)
        state.fileBrowser.closeStructure();

    if (ImGui::IsKeyPressed(ImGuiKey_R)
        && !ctrlHeld
        && !ImGui::GetIO().WantTextInput
        && !ImGui::GetIO().KeyAlt
        && !ImGui::GetIO().KeySuper)
    {
        state.pendingDefaultViewReset = true;
    }

    if (ImGui::IsKeyPressed(ImGuiKey_Z) && ctrlHeld && !ImGui::GetIO().KeyShift)
        requests.requestUndo = true;

    if ((ImGui::IsKeyPressed(ImGuiKey_Y) && ctrlHeld) ||
        (ImGui::IsKeyPressed(ImGuiKey_Z) && ctrlHeld && ImGui::GetIO().KeyShift))
    {
        requests.requestRedo = true;
    }
}

void handlePendingAtomPick(
    Camera& camera,
    EditorState& state,
    const glm::vec3& cameraPosition,
    int windowWidth,
    int windowHeight,
    const glm::mat4& projection,
    const glm::mat4& view)
{
    if (!camera.pendingClick)
        return;

    // Picking disabled for large structures (CPU caches not available)
    if (state.sceneBuffers.cpuCachesDisabled)
    {
        camera.pendingClick = false;
        return;
    }

    camera.pendingClick = false;
    const bool isOrthographicProjection = std::abs(projection[3][3] - 1.0f) <= 1e-4f;
    const glm::vec3 rayOrigin = isOrthographicProjection
        ? pickRayOrigin(camera.clickX, camera.clickY, windowWidth, windowHeight, projection, view)
        : cameraPosition;
    glm::vec3 ray = pickRayDir(camera.clickX, camera.clickY, windowWidth, windowHeight, projection, view);
    int pickedIndex = pickAtom(
        rayOrigin,
        ray,
        state.sceneBuffers.atomPositions,
        state.sceneBuffers.atomRadii,
        1.0f);

    if (pickedIndex >= 0)
    {
        if (isCtrlHeld())
            toggleSelectedInstance(state, pickedIndex);
        else
            selectSingleInstance(state, pickedIndex);
    }
    else
    {
        clearSelection(state);
    }
}

void handleRightClick(Camera& camera, EditorState& state)
{
    if (!camera.pendingRightClick)
        return;

    if (state.fileBrowser.isBoxSelectModeEnabled())
    {
        camera.pendingRightClick = false;
        return;
    }

    camera.pendingRightClick = false;
    if (!state.selectedInstanceIndices.empty())
        state.contextMenu.open();
}

void handleBoxSelection(
    EditorState& state,
    int windowWidth,
    int windowHeight,
    const glm::mat4& projection,
    const glm::mat4& view,
    ImDrawList* drawList)
{
    if (!state.fileBrowser.isBoxSelectModeEnabled())
        return;

    if (state.sceneBuffers.cpuCachesDisabled)
        return;

    static bool dragging = false;
    static ImVec2 dragStart(0.0f, 0.0f);
    static ImVec2 dragEnd(0.0f, 0.0f);

    ImGuiIO& io = ImGui::GetIO();
    const ImVec2 mousePos = io.MousePos;

    if (ImGui::IsMouseClicked(ImGuiMouseButton_Right) && !io.WantCaptureMouse)
    {
        dragging = true;
        dragStart = mousePos;
        dragEnd = mousePos;
    }

    if (dragging && ImGui::IsMouseDown(ImGuiMouseButton_Right))
        dragEnd = mousePos;

    if (dragging)
        drawSelectionRect(drawList, dragStart, dragEnd);

    if (dragging && ImGui::IsMouseReleased(ImGuiMouseButton_Right))
    {
        dragging = false;

        if (!isDragLargeEnough(dragStart, dragEnd))
            return;

        const bool additive = io.KeyCtrl;
        if (!additive)
            clearSelection(state);

        for (int i = 0; i < (int)state.sceneBuffers.atomPositions.size(); ++i)
        {
            const glm::vec3& worldPos = state.sceneBuffers.atomPositions[i];
            const glm::vec4 clip = projection * view * glm::vec4(worldPos, 1.0f);
            if (clip.w <= 1e-6f)
                continue;

            const glm::vec3 ndc = glm::vec3(clip) / clip.w;
            if (ndc.z < -1.0f || ndc.z > 1.0f)
                continue;

            const ImVec2 screenPos(
                (ndc.x * 0.5f + 0.5f) * (float)windowWidth,
                (1.0f - (ndc.y * 0.5f + 0.5f)) * (float)windowHeight);

            if (!isInsideSelectionRect(screenPos, dragStart, dragEnd))
                continue;

            auto it = std::find(state.selectedInstanceIndices.begin(), state.selectedInstanceIndices.end(), i);
            if (it == state.selectedInstanceIndices.end())
                state.selectedInstanceIndices.push_back(i);
        }
    }
}

// --- Blender-style Grab Mode ---

namespace
{
constexpr float kGrabPbcTol = 1e-4f;

float grabWrapFrac(float value)
{
    value -= std::floor(value);
    if (std::abs(value) <= kGrabPbcTol)
        return 0.0f;
    if (std::abs(1.0f - value) <= kGrabPbcTol)
        return 0.0f;
    return value;
}

void startGrab(EditorState& state)
{
    GrabState& grab = state.grabState;
    grab.active = true;
    grab.axisConstraint = GrabAxisConstraint::None;

    const ImVec2 mousePos = ImGui::GetIO().MousePos;
    grab.startMousePos = glm::vec2(mousePos.x, mousePos.y);

    grab.originalPositions.clear();
    for (int instanceIdx : state.selectedInstanceIndices)
    {
        if (instanceIdx >= 0 && instanceIdx < (int)state.sceneBuffers.atomPositions.size())
            grab.originalPositions.push_back(state.sceneBuffers.atomPositions[instanceIdx]);
        else
            grab.originalPositions.push_back(glm::vec3(0.0f));
    }

    // Compute cell info for periodic handling
    grab.hasCellInfo = false;
    if (state.structure.hasUnitCell)
    {
        grab.cellA = glm::vec3((float)state.structure.cellVectors[0][0],
                               (float)state.structure.cellVectors[0][1],
                               (float)state.structure.cellVectors[0][2]);
        grab.cellB = glm::vec3((float)state.structure.cellVectors[1][0],
                               (float)state.structure.cellVectors[1][1],
                               (float)state.structure.cellVectors[1][2]);
        grab.cellC = glm::vec3((float)state.structure.cellVectors[2][0],
                               (float)state.structure.cellVectors[2][1],
                               (float)state.structure.cellVectors[2][2]);
        grab.cellMatrix = glm::mat3(grab.cellA, grab.cellB, grab.cellC);
        const float det = glm::determinant(grab.cellMatrix);
        if (std::abs(det) > 1e-8f)
        {
            grab.invCellMatrix = glm::inverse(grab.cellMatrix);
            grab.cellOrigin = glm::vec3((float)state.structure.cellOffset[0],
                                        (float)state.structure.cellOffset[1],
                                        (float)state.structure.cellOffset[2]);
            grab.pbcTolerance = state.structure.pbcBoundaryTol > 0.0f
                ? state.structure.pbcBoundaryTol : kGrabPbcTol;
            grab.hasCellInfo = true;
        }
    }

    // Build periodic siblings map with cell shifts
    grab.periodicSiblings.clear();
    grab.originalBasePositions.clear();

    // Collect the set of base atom indices involved in the grab
    std::unordered_set<int> grabbedBaseIndices;
    for (int instanceIdx : state.selectedInstanceIndices)
    {
        if (instanceIdx >= 0 && instanceIdx < (int)state.sceneBuffers.atomIndices.size())
            grabbedBaseIndices.insert(state.sceneBuffers.atomIndices[instanceIdx]);
    }

    // Store original base atom (canonical) positions from structure
    for (int baseIdx : grabbedBaseIndices)
    {
        if (baseIdx >= 0 && baseIdx < (int)state.structure.atoms.size())
        {
            const auto& atom = state.structure.atoms[baseIdx];
            grab.originalBasePositions[baseIdx] = glm::vec3((float)atom.x, (float)atom.y, (float)atom.z);
        }
    }

    // Scan all instances and group by base atom index, computing cell shifts
    for (int i = 0; i < (int)state.sceneBuffers.atomIndices.size(); ++i)
    {
        const int baseIdx = state.sceneBuffers.atomIndices[i];
        if (!grabbedBaseIndices.count(baseIdx))
            continue;

        PeriodicSibling sib;
        sib.instanceIdx = i;
        sib.originalPos = state.sceneBuffers.atomPositions[i];
        sib.cellShift = glm::ivec3(0);

        // Compute cell shift: how many cell vectors offset is this instance from canonical?
        if (grab.hasCellInfo)
        {
            const glm::vec3 frac = grab.invCellMatrix * (sib.originalPos - grab.cellOrigin);
            // Cell shift is the integer part: images at frac >= 1 in some dimension
            for (int d = 0; d < 3; ++d)
                sib.cellShift[d] = (frac[d] > 0.5f) ? 1 : 0;
        }

        grab.periodicSiblings[baseIdx].push_back(sib);
    }
}

void cancelGrab(EditorState& state)
{
    GrabState& grab = state.grabState;

    // Restore all instances (including periodic images) to original positions
    for (auto& pair : grab.periodicSiblings)
    {
        for (const PeriodicSibling& sib : pair.second)
            state.sceneBuffers.updateAtomPosition(sib.instanceIdx, sib.originalPos);
    }

    // Restore base atom positions from stored originals (not from selected instance
    // positions, which could be periodic images with cell-vector offsets).
    for (auto& pair : grab.originalBasePositions)
    {
        const int baseIdx = pair.first;
        const glm::vec3& origPos = pair.second;
        if (baseIdx >= 0 && baseIdx < (int)state.structure.atoms.size())
        {
            state.structure.atoms[baseIdx].x = (double)origPos.x;
            state.structure.atoms[baseIdx].y = (double)origPos.y;
            state.structure.atoms[baseIdx].z = (double)origPos.z;
        }
    }

    grab.active = false;
    grab.originalPositions.clear();
    grab.periodicSiblings.clear();
    grab.originalBasePositions.clear();
}

void confirmGrab(EditorState& state)
{
    GrabState& grab = state.grabState;

    // Wrap moved atoms into the unit cell before rebuilding
    if (grab.hasCellInfo)
    {
        for (auto& pair : grab.originalBasePositions)
        {
            const int baseIdx = pair.first;
            if (baseIdx < 0 || baseIdx >= (int)state.structure.atoms.size())
                continue;

            AtomSite& atom = state.structure.atoms[baseIdx];
            const glm::vec3 pos((float)atom.x, (float)atom.y, (float)atom.z);
            glm::vec3 frac = grab.invCellMatrix * (pos - grab.cellOrigin);
            frac.x = grabWrapFrac(frac.x);
            frac.y = grabWrapFrac(frac.y);
            frac.z = grabWrapFrac(frac.z);
            const glm::vec3 wrapped = grab.cellOrigin + grab.cellMatrix * frac;
            atom.x = (double)wrapped.x;
            atom.y = (double)wrapped.y;
            atom.z = (double)wrapped.z;
        }
    }

    grab.active = false;
    grab.originalPositions.clear();
    grab.periodicSiblings.clear();
    grab.originalBasePositions.clear();

    // Full rebuild to get correct periodic images at the new position
    updateBuffers(state);
    state.voronoiDirty = true;
}

glm::vec3 screenToWorldDelta(
    const glm::vec2& screenDelta,
    const glm::vec3& refWorldPos,
    const glm::mat4& projection,
    const glm::mat4& view,
    int windowWidth,
    int windowHeight)
{
    // Project reference point to screen
    const glm::vec4 clip = projection * view * glm::vec4(refWorldPos, 1.0f);
    if (clip.w <= 1e-6f)
        return glm::vec3(0.0f);

    const glm::vec3 ndc = glm::vec3(clip) / clip.w;

    // Convert screen delta to NDC delta
    const float ndcDx = (screenDelta.x / (float)windowWidth) * 2.0f;
    const float ndcDy = -(screenDelta.y / (float)windowHeight) * 2.0f;

    // Unproject the original and offset NDC positions at the same depth
    const glm::mat4 invVP = glm::inverse(projection * view);

    const glm::vec4 origNDC(ndc.x, ndc.y, ndc.z, 1.0f);
    const glm::vec4 offsetNDC(ndc.x + ndcDx, ndc.y + ndcDy, ndc.z, 1.0f);

    glm::vec4 origWorld = invVP * origNDC;
    glm::vec4 offsetWorld = invVP * offsetNDC;

    if (std::abs(origWorld.w) < 1e-8f || std::abs(offsetWorld.w) < 1e-8f)
        return glm::vec3(0.0f);

    origWorld /= origWorld.w;
    offsetWorld /= offsetWorld.w;

    return glm::vec3(offsetWorld) - glm::vec3(origWorld);
}
} // namespace

void handleGrabMode(
    EditorState& state,
    Camera& camera,
    const glm::mat4& projection,
    const glm::mat4& view,
    int windowWidth,
    int windowHeight)
{
    GrabState& grab = state.grabState;

    // Start grab with G key when atoms are selected and grab is not active
    if (!grab.active)
    {
        if (ImGui::IsKeyPressed(ImGuiKey_G) && !state.selectedInstanceIndices.empty()
            && !state.sceneBuffers.cpuCachesDisabled && !isCtrlHeld())
        {
            startGrab(state);
        }
        return;
    }

    // --- Grab mode is active ---

    // Cancel with Escape or right-click
    if (ImGui::IsKeyPressed(ImGuiKey_Escape) || ImGui::IsMouseClicked(ImGuiMouseButton_Right))
    {
        cancelGrab(state);
        return;
    }

    // Axis constraint keys
    if (ImGui::IsKeyPressed(ImGuiKey_X))
        grab.axisConstraint = (grab.axisConstraint == GrabAxisConstraint::X)
            ? GrabAxisConstraint::None : GrabAxisConstraint::X;
    if (ImGui::IsKeyPressed(ImGuiKey_Y))
        grab.axisConstraint = (grab.axisConstraint == GrabAxisConstraint::Y)
            ? GrabAxisConstraint::None : GrabAxisConstraint::Y;
    if (ImGui::IsKeyPressed(ImGuiKey_Z))
        grab.axisConstraint = (grab.axisConstraint == GrabAxisConstraint::Z)
            ? GrabAxisConstraint::None : GrabAxisConstraint::Z;

    // Confirm with left-click or Enter
    if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsKeyPressed(ImGuiKey_Enter))
    {
        confirmGrab(state);
        // Consume the pending click so it doesn't trigger atom picking
        camera.pendingClick = false;
        return;
    }

    // Compute mouse delta from start
    const ImVec2 mousePos = ImGui::GetIO().MousePos;
    const glm::vec2 currentMouse(mousePos.x, mousePos.y);
    const glm::vec2 screenDelta = currentMouse - grab.startMousePos;

    // Use the centroid of selected atoms as the reference point for depth
    glm::vec3 centroid(0.0f);
    for (const glm::vec3& p : grab.originalPositions)
        centroid += p;
    if (!grab.originalPositions.empty())
        centroid /= (float)grab.originalPositions.size();

    // Convert screen delta to world-space delta
    glm::vec3 worldDelta = screenToWorldDelta(
        screenDelta, centroid, projection, view, windowWidth, windowHeight);

    // Apply axis constraint
    switch (grab.axisConstraint)
    {
        case GrabAxisConstraint::X: worldDelta = glm::vec3(worldDelta.x, 0.0f, 0.0f); break;
        case GrabAxisConstraint::Y: worldDelta = glm::vec3(0.0f, worldDelta.y, 0.0f); break;
        case GrabAxisConstraint::Z: worldDelta = glm::vec3(0.0f, 0.0f, worldDelta.z); break;
        default: break;
    }

    // Apply delta to all selected atoms and their periodic images
    std::unordered_set<int> updatedBaseIndices;
    for (size_t i = 0; i < state.selectedInstanceIndices.size(); ++i)
    {
        const int instanceIdx = state.selectedInstanceIndices[i];

        if (instanceIdx < 0 || instanceIdx >= (int)state.sceneBuffers.atomIndices.size())
            continue;
        const int baseIdx = state.sceneBuffers.atomIndices[instanceIdx];
        if (!updatedBaseIndices.insert(baseIdx).second)
            continue; // already processed this base atom

        auto sibIt = grab.periodicSiblings.find(baseIdx);
        if (sibIt == grab.periodicSiblings.end())
            continue;
        auto baseIt = grab.originalBasePositions.find(baseIdx);
        if (baseIt == grab.originalBasePositions.end())
            continue;

        // Compute new canonical position from original base atom + world delta.
        // This avoids depending on which sibling instance was selected.
        const glm::vec3 canonNew = baseIt->second + worldDelta;

        // For non-periodic structures, just move all instances by the same delta
        if (!grab.hasCellInfo)
        {
            for (const PeriodicSibling& sib : sibIt->second)
                state.sceneBuffers.updateAtomPosition(sib.instanceIdx, sib.originalPos + worldDelta);
            if (baseIdx >= 0 && baseIdx < (int)state.structure.atoms.size())
            {
                state.structure.atoms[baseIdx].x = (double)canonNew.x;
                state.structure.atoms[baseIdx].y = (double)canonNew.y;
                state.structure.atoms[baseIdx].z = (double)canonNew.z;
            }
            continue;
        }

        // Wrap canonical position into the unit cell [0,1) in fractional space
        glm::vec3 frac = grab.invCellMatrix * (canonNew - grab.cellOrigin);
        frac.x = grabWrapFrac(frac.x);
        frac.y = grabWrapFrac(frac.y);
        frac.z = grabWrapFrac(frac.z);
        const glm::vec3 wrappedCanon = grab.cellOrigin + grab.cellMatrix * frac;

        // Update structure data with wrapped canonical position
        if (baseIdx >= 0 && baseIdx < (int)state.structure.atoms.size())
        {
            state.structure.atoms[baseIdx].x = (double)wrappedCanon.x;
            state.structure.atoms[baseIdx].y = (double)wrappedCanon.y;
            state.structure.atoms[baseIdx].z = (double)wrappedCanon.z;
        }

        // Determine which boundary images are needed from wrapped fractional coords
        const float tol = grab.pbcTolerance;
        std::vector<int> shiftsX = {0}, shiftsY = {0}, shiftsZ = {0};
        if (std::abs(frac.x) <= tol) shiftsX.push_back(1);
        if (std::abs(frac.y) <= tol) shiftsY.push_back(1);
        if (std::abs(frac.z) <= tol) shiftsZ.push_back(1);

        // Build list of needed image positions (primary + boundary images)
        std::vector<glm::vec3> neededPositions;
        for (int sx : shiftsX)
            for (int sy : shiftsY)
                for (int sz : shiftsZ)
                    neededPositions.push_back(
                        wrappedCanon
                        + (float)sx * grab.cellA
                        + (float)sy * grab.cellB
                        + (float)sz * grab.cellC);

        // Assign positions to ALL sibling instances (including the selected one).
        // First N get unique visible positions; the rest are hidden at wrappedCanon.
        auto& siblings = sibIt->second;
        for (size_t j = 0; j < siblings.size(); ++j)
        {
            if (j < neededPositions.size())
                state.sceneBuffers.updateAtomPosition(siblings[j].instanceIdx, neededPositions[j]);
            else
                state.sceneBuffers.updateAtomPosition(siblings[j].instanceIdx, wrappedCanon);
        }
    }
}

namespace
{
ImVec2 projectToScreen(const glm::vec3& worldPos,
                       const glm::mat4& projection,
                       const glm::mat4& view,
                       int windowWidth,
                       int windowHeight,
                       bool& visible)
{
    const glm::vec4 clip = projection * view * glm::vec4(worldPos, 1.0f);
    if (clip.w <= 1e-6f)
    {
        visible = false;
        return ImVec2(0.0f, 0.0f);
    }
    const glm::vec3 ndc = glm::vec3(clip) / clip.w;
    if (ndc.z < -1.0f || ndc.z > 1.0f)
    {
        visible = false;
        return ImVec2(0.0f, 0.0f);
    }
    visible = true;
    return ImVec2(
        (ndc.x * 0.5f + 0.5f) * (float)windowWidth,
        (1.0f - (ndc.y * 0.5f + 0.5f)) * (float)windowHeight);
}

void drawArrowhead(ImDrawList* drawList, const ImVec2& tip, const ImVec2& tail, ImU32 color, float size)
{
    const float dx = tip.x - tail.x;
    const float dy = tip.y - tail.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1.0f)
        return;

    const float nx = dx / len;
    const float ny = dy / len;

    // Two points forming the arrowhead wings
    const float px = -ny;
    const float py = nx;

    const ImVec2 wing1(tip.x - nx * size + px * size * 0.4f,
                       tip.y - ny * size + py * size * 0.4f);
    const ImVec2 wing2(tip.x - nx * size - px * size * 0.4f,
                       tip.y - ny * size - py * size * 0.4f);

    drawList->AddTriangleFilled(tip, wing1, wing2, color);
}
} // namespace

void drawGrabOverlay(
    const EditorState& state,
    ImDrawList* drawList,
    const glm::mat4& projection,
    const glm::mat4& view,
    int windowWidth,
    int windowHeight)
{
    const GrabState& grab = state.grabState;
    if (!grab.active || state.selectedInstanceIndices.empty())
        return;

    // Header overlay showing grab mode status
    const char* axisLabel = "";
    ImU32 axisColor = isLightTheme() ? IM_COL32(30, 30, 30, 220) : IM_COL32(255, 255, 255, 220);
    switch (grab.axisConstraint)
    {
        case GrabAxisConstraint::X: axisLabel = " [X axis]"; axisColor = IM_COL32(255, 80, 80, 220); break;
        case GrabAxisConstraint::Y: axisLabel = " [Y axis]"; axisColor = IM_COL32(80, 255, 80, 220); break;
        case GrabAxisConstraint::Z: axisLabel = " [Z axis]"; axisColor = IM_COL32(80, 130, 255, 220); break;
        default: break;
    }

    char headerBuf[128];
    std::snprintf(headerBuf, sizeof(headerBuf), "Grab Mode%s  |  LMB/Enter: confirm  |  Esc/RMB: cancel  |  X/Y/Z: constrain", axisLabel);

    const ImVec2 headerSize = ImGui::CalcTextSize(headerBuf);
    const float headerX = ((float)windowWidth - headerSize.x) * 0.5f;
    const float headerY = 8.0f;

    // Background rectangle for header
    drawList->AddRectFilled(
        ImVec2(headerX - 8.0f, headerY - 4.0f),
        ImVec2(headerX + headerSize.x + 8.0f, headerY + headerSize.y + 4.0f),
        isLightTheme() ? IM_COL32(240, 240, 240, 210) : IM_COL32(20, 20, 20, 200), 4.0f);
    drawList->AddText(ImVec2(headerX, headerY), axisColor, headerBuf);

    // Draw 3D arrows and coordinate info for each selected atom (up to 8 to avoid clutter)
    const int maxLabels = std::min((int)state.selectedInstanceIndices.size(), 8);
    for (int i = 0; i < maxLabels; ++i)
    {
        const int instanceIdx = state.selectedInstanceIndices[i];
        if (instanceIdx < 0 || instanceIdx >= (int)state.sceneBuffers.atomPositions.size())
            continue;
        if (i >= (int)grab.originalPositions.size())
            continue;

        const glm::vec3& worldPos = state.sceneBuffers.atomPositions[instanceIdx];
        const glm::vec3& origPos = grab.originalPositions[i];

        // --- Draw 3D arrow from original position to current position ---
        bool origVisible = false;
        bool curVisible = false;
        const ImVec2 origScreen = projectToScreen(origPos, projection, view, windowWidth, windowHeight, origVisible);
        const ImVec2 curScreen = projectToScreen(worldPos, projection, view, windowWidth, windowHeight, curVisible);

        if (origVisible && curVisible)
        {
            const float arrowDx = curScreen.x - origScreen.x;
            const float arrowDy = curScreen.y - origScreen.y;
            const float arrowLen = std::sqrt(arrowDx * arrowDx + arrowDy * arrowDy);

            if (arrowLen > 3.0f)
            {
                // Pick color based on axis constraint
                ImU32 arrowColor;
                switch (grab.axisConstraint)
                {
                    case GrabAxisConstraint::X: arrowColor = IM_COL32(255, 80, 80, 200); break;
                    case GrabAxisConstraint::Y: arrowColor = IM_COL32(80, 255, 80, 200); break;
                    case GrabAxisConstraint::Z: arrowColor = IM_COL32(80, 130, 255, 200); break;
                    default:                    arrowColor = IM_COL32(255, 200, 50, 200); break;
                }

                // Draw dashed origin marker (small circle at original position)
                drawList->AddCircle(origScreen, 5.0f,
                                    isLightTheme() ? IM_COL32(100, 100, 100, 150) : IM_COL32(180, 180, 180, 150), 12, 1.5f);

                // Draw arrow shaft
                drawList->AddLine(origScreen, curScreen, arrowColor, 2.5f);

                // Draw arrowhead at current position
                drawArrowhead(drawList, curScreen, origScreen, arrowColor, 12.0f);
            }
        }

        // --- Draw coordinate label next to atom ---
        if (!curVisible)
            continue;

        // Build coordinate text
        char cartBuf[128];
        std::snprintf(cartBuf, sizeof(cartBuf), "Cart: (%.4f, %.4f, %.4f)",
                      worldPos.x, worldPos.y, worldPos.z);

        char fracBuf[128];
        fracBuf[0] = '\0';
        glm::vec3 frac;
        if (tryCartesianToFractional(state.structure, worldPos, frac))
        {
            std::snprintf(fracBuf, sizeof(fracBuf), "Frac: (%.4f, %.4f, %.4f)",
                          frac.x, frac.y, frac.z);
        }

        const float labelX = curScreen.x + 12.0f;
        float labelY = curScreen.y - 10.0f;

        const ImVec2 cartSize = ImGui::CalcTextSize(cartBuf);
        const ImVec2 fracSize = fracBuf[0] ? ImGui::CalcTextSize(fracBuf) : ImVec2(0, 0);
        const float bgWidth = std::max(cartSize.x, fracSize.x) + 12.0f;
        const float bgHeight = cartSize.y + (fracBuf[0] ? fracSize.y + 2.0f : 0.0f) + 8.0f;

        drawList->AddRectFilled(
            ImVec2(labelX - 4.0f, labelY - 4.0f),
            ImVec2(labelX + bgWidth, labelY + bgHeight),
            isLightTheme() ? IM_COL32(240, 240, 240, 210) : IM_COL32(10, 10, 10, 200), 3.0f);

        drawList->AddText(ImVec2(labelX, labelY),
                          isLightTheme() ? IM_COL32(120, 100, 20, 255) : IM_COL32(255, 255, 180, 255), cartBuf);
        labelY += cartSize.y + 2.0f;

        if (fracBuf[0])
            drawList->AddText(ImVec2(labelX, labelY),
                              isLightTheme() ? IM_COL32(30, 80, 150, 255) : IM_COL32(180, 220, 255, 255), fracBuf);
    }

    // Draw arrows for remaining atoms beyond maxLabels (no coordinate labels)
    for (int i = maxLabels; i < (int)state.selectedInstanceIndices.size(); ++i)
    {
        const int instanceIdx = state.selectedInstanceIndices[i];
        if (instanceIdx < 0 || instanceIdx >= (int)state.sceneBuffers.atomPositions.size())
            continue;
        if (i >= (int)grab.originalPositions.size())
            continue;

        const glm::vec3& worldPos = state.sceneBuffers.atomPositions[instanceIdx];
        const glm::vec3& origPos = grab.originalPositions[i];

        bool origVisible = false;
        bool curVisible = false;
        const ImVec2 origScreen = projectToScreen(origPos, projection, view, windowWidth, windowHeight, origVisible);
        const ImVec2 curScreen = projectToScreen(worldPos, projection, view, windowWidth, windowHeight, curVisible);

        if (origVisible && curVisible)
        {
            const float arrowDx = curScreen.x - origScreen.x;
            const float arrowDy = curScreen.y - origScreen.y;
            const float arrowLen = std::sqrt(arrowDx * arrowDx + arrowDy * arrowDy);

            if (arrowLen > 3.0f)
            {
                ImU32 arrowColor;
                switch (grab.axisConstraint)
                {
                    case GrabAxisConstraint::X: arrowColor = IM_COL32(255, 80, 80, 140); break;
                    case GrabAxisConstraint::Y: arrowColor = IM_COL32(80, 255, 80, 140); break;
                    case GrabAxisConstraint::Z: arrowColor = IM_COL32(80, 130, 255, 140); break;
                    default:                    arrowColor = IM_COL32(255, 200, 50, 140); break;
                }
                drawList->AddCircle(origScreen, 4.0f,
                                    isLightTheme() ? IM_COL32(100, 100, 100, 100) : IM_COL32(150, 150, 150, 100), 12, 1.0f);
                drawList->AddLine(origScreen, curScreen, arrowColor, 1.5f);
                drawArrowhead(drawList, curScreen, origScreen, arrowColor, 9.0f);
            }
        }
    }

    if ((int)state.selectedInstanceIndices.size() > maxLabels)
    {
        char moreBuf[64];
        std::snprintf(moreBuf, sizeof(moreBuf), "... and %d more atoms",
                      (int)state.selectedInstanceIndices.size() - maxLabels);
        const ImVec2 moreSize = ImGui::CalcTextSize(moreBuf);
        const float moreX = ((float)windowWidth - moreSize.x) * 0.5f;
        const float moreY = (float)windowHeight - 30.0f;
        drawList->AddRectFilled(
            ImVec2(moreX - 4.0f, moreY - 2.0f),
            ImVec2(moreX + moreSize.x + 4.0f, moreY + moreSize.y + 2.0f),
            isLightTheme() ? IM_COL32(240, 240, 240, 200) : IM_COL32(10, 10, 10, 180), 3.0f);
        drawList->AddText(ImVec2(moreX, moreY),
                          isLightTheme() ? IM_COL32(60, 60, 60, 200) : IM_COL32(200, 200, 200, 200), moreBuf);
    }
}