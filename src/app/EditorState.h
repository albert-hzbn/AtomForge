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
#include <vector>

enum class GrabAxisConstraint
{
    None,
    X,
    Y,
    Z
};

struct GrabState
{
    bool active = false;
    GrabAxisConstraint axisConstraint = GrabAxisConstraint::None;
    glm::vec2 startMousePos = glm::vec2(0.0f);
    // Original positions of grabbed atoms (indexed same as selectedInstanceIndices)
    std::vector<glm::vec3> originalPositions;
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