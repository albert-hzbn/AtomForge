#include "app/EditorOps.h"

#include "Camera.h"
#include "graphics/StructureInstanceBuilder.h"

#include <algorithm>
#include <cmath>

namespace
{
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
}

EditorSnapshot captureSnapshot(const EditorState& state)
{
    EditorSnapshot snapshot;
    snapshot.structure = state.structure;
    snapshot.elementRadii = state.editMenuDialogs.elementRadii;
    snapshot.elementColors = state.editMenuDialogs.elementColors;
    snapshot.elementShininess = state.editMenuDialogs.elementShininess;
    return snapshot;
}

void updateBuffers(EditorState& state)
{
    if (state.fileBrowser.isTransformMatrixEnabled() && state.structure.hasUnitCell)
    {
        state.structure = buildSupercell(state.structure, state.fileBrowser.getTransformMatrix());
        state.fileBrowser.clearTransformMatrix();
    }

    for (auto& atom : state.structure.atoms)
    {
        int atomicNumber = atom.atomicNumber;
        if (atomicNumber >= 0 && atomicNumber < (int)state.editMenuDialogs.elementColors.size())
        {
            atom.r = state.editMenuDialogs.elementColors[atomicNumber].r;
            atom.g = state.editMenuDialogs.elementColors[atomicNumber].g;
            atom.b = state.editMenuDialogs.elementColors[atomicNumber].b;
        }
    }

    StructureInstanceData data = buildStructureInstanceData(
        state.structure,
        state.fileBrowser.isTransformMatrixEnabled(),
        state.fileBrowser.getTransformMatrix(),
        state.editMenuDialogs.elementRadii,
        state.editMenuDialogs.elementShininess);

    state.sceneBuffers.upload(
        data,
        state.fileBrowser.isBondElementFilterEnabled(),
        state.fileBrowser.getBondElementFilterMask());
    state.selectedInstanceIndices.clear();

    if (!state.suppressHistoryCommit)
        state.undoRedo.commit(captureSnapshot(state));
}

void updateBuffers(EditorState& state, Structure& structure)
{
    state.structure = structure;
    updateBuffers(state);
    structure = state.structure;
}

void applySnapshot(EditorState& state, const EditorSnapshot& snapshot)
{
    state.structure = snapshot.structure;
    state.editMenuDialogs.elementRadii = snapshot.elementRadii;
    state.editMenuDialogs.elementColors = snapshot.elementColors;
    state.editMenuDialogs.elementShininess = snapshot.elementShininess;
    state.suppressHistoryCommit = true;
    updateBuffers(state);
    state.suppressHistoryCommit = false;
    state.measurementState.clearVisuals();
}

void applyDefaultView(
    Camera& camera,
    const SceneBuffers& sceneBuffers,
    int viewportWidth,
    int viewportHeight,
    bool fitToStructure)
{
    constexpr float kIsoYawDeg = 45.0f;
    constexpr float kIsoPitchDeg = 35.2643897f;

    camera.yaw = kIsoYawDeg;
    camera.pitch = kIsoPitchDeg;

    if (!fitToStructure || sceneBuffers.atomCount == 0)
    {
        camera.distance = 10.0f;
        return;
    }

    float maxRadius = estimateSceneRadius(sceneBuffers);

    float aspect = (viewportHeight > 0) ? (float)viewportWidth / (float)viewportHeight : 1.0f;
    float verticalFov = glm::radians(45.0f);
    float horizontalFov = 2.0f * std::atan(std::tan(verticalFov * 0.5f) * aspect);
    float halfFov = 0.5f * std::min(verticalFov, horizontalFov);
    halfFov = std::max(halfFov, glm::radians(10.0f));

    float framedDistance = (maxRadius / std::sin(halfFov)) * 1.15f;
    camera.distance = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, framedDistance));
}

void clearSelection(EditorState& state)
{
    for (int index : state.selectedInstanceIndices)
        state.sceneBuffers.restoreAtomColor(index);

    state.selectedInstanceIndices.clear();
    state.measurementState.clearVisuals();
}

void deleteSelectedAtoms(EditorState& state)
{
    if (state.selectedInstanceIndices.empty())
        return;

    std::vector<int> baseIndicesToDelete;
    for (int selectedIndex : state.selectedInstanceIndices)
    {
        if (selectedIndex >= 0 && selectedIndex < (int)state.sceneBuffers.atomIndices.size())
        {
            int baseIndex = state.sceneBuffers.atomIndices[selectedIndex];
            if (baseIndex >= 0 && baseIndex < (int)state.structure.atoms.size())
                baseIndicesToDelete.push_back(baseIndex);
        }
    }

    std::sort(baseIndicesToDelete.begin(), baseIndicesToDelete.end(), std::greater<int>());
    baseIndicesToDelete.erase(
        std::unique(baseIndicesToDelete.begin(), baseIndicesToDelete.end()),
        baseIndicesToDelete.end());

    for (int baseIndex : baseIndicesToDelete)
    {
        if (baseIndex >= 0 && baseIndex < (int)state.structure.atoms.size())
            state.structure.atoms.erase(state.structure.atoms.begin() + baseIndex);
    }

    updateBuffers(state);
}

void refreshSelectionHighlights(EditorState& state)
{
    for (int& index : state.selectedInstanceIndices)
    {
        if (index >= (int)state.sceneBuffers.atomCount)
        {
            state.sceneBuffers.restoreAtomColor(index);
            index = -1;
        }
        else if (index >= 0)
        {
            state.sceneBuffers.highlightAtom(index, glm::vec3(1.0f, 1.0f, 0.0f));
        }
    }

    state.selectedInstanceIndices.erase(
        std::remove(state.selectedInstanceIndices.begin(), state.selectedInstanceIndices.end(), -1),
        state.selectedInstanceIndices.end());
}