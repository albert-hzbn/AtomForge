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

#include <string>
#include <vector>

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
};