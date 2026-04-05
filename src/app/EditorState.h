#pragma once

#include "UndoRedo.h"
#include "algorithms/VoronoiComputation.h"
#include "graphics/SceneBuffers.h"
#include "io/StructureLoader.h"
#include "ui/AtomContextMenu.h"
#include "ui/EditMenuDialogs.h"
#include "ui/FileBrowser.h"
#include "ui/MeasurementOverlay.h"
#include "ui/StructureInfoDialog.h"

#include <glm/glm.hpp>

#include <string>
#include <unordered_map>
#include <vector>

enum class GrabAxisConstraint
{
    None,
    X,
    Y,
    Z
};

struct PeriodicSibling
{
    int instanceIdx;
    glm::vec3 originalPos;
    glm::ivec3 cellShift;  // integer cell vector shift (e.g., (1,0,0) = image at +a)
};

struct GrabState
{
    bool active = false;
    GrabAxisConstraint axisConstraint = GrabAxisConstraint::None;
    glm::vec2 startMousePos = glm::vec2(0.0f);
    // Original positions of grabbed atoms (indexed same as selectedInstanceIndices)
    std::vector<glm::vec3> originalPositions;
    // Maps base atom index -> list of ALL instances (including periodic images) with their original positions.
    // Built at grab start so periodic siblings move together.
    std::unordered_map<int, std::vector<PeriodicSibling>> periodicSiblings;
    // Original base atom positions (canonical cartesian) for correct cancel.
    std::unordered_map<int, glm::vec3> originalBasePositions;

    // Unit cell info for periodic image handling during grab
    bool hasCellInfo = false;
    glm::mat3 cellMatrix = glm::mat3(1.0f);
    glm::mat3 invCellMatrix = glm::mat3(1.0f);
    glm::vec3 cellOrigin = glm::vec3(0.0f);
    glm::vec3 cellA = glm::vec3(0.0f);
    glm::vec3 cellB = glm::vec3(0.0f);
    glm::vec3 cellC = glm::vec3(0.0f);
    float pbcTolerance = 1e-4f;
};

struct EditorState
{
    Structure structure;
    FileBrowser fileBrowser;
    EditMenuDialogs editMenuDialogs;
    SceneBuffers sceneBuffers;
    std::vector<int> selectedInstanceIndices;
    AtomContextMenu contextMenu;
    MeasurementOverlayState measurementState;
    StructureInfoDialogState structureInfoDialog;
    UndoRedoManager undoRedo;
    std::vector<std::string> pendingDroppedFiles;
    bool suppressHistoryCommit = false;
    bool pendingDefaultViewReset = true;
    VoronoiDiagram voronoiDiagram;
    bool voronoiDirty = true;  // recompute when structure changes
    GrabState grabState;
};