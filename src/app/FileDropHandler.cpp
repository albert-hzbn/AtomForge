#include "app/FileDropHandler.h"

#include "app/EditorOps.h"
#include "app/EditorState.h"

#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <utility>

namespace
{
void dropFileCallback(GLFWwindow* window, int count, const char** paths)
{
    if (count <= 0 || paths == nullptr)
        return;

    EditorState* state = static_cast<EditorState*>(glfwGetWindowUserPointer(window));
    if (!state)
        return;

    for (int i = 0; i < count; ++i)
    {
        if (paths[i] == nullptr || paths[i][0] == '\0')
            continue;
        state->pendingDroppedFiles.push_back(paths[i]);
    }
}
}

void installDropFileCallback(GLFWwindow* window, EditorState& state)
{
    glfwSetWindowUserPointer(window, &state);
    glfwSetDropCallback(window, dropFileCallback);
}

void processDroppedFiles(EditorState& state)
{
    if (state.pendingDroppedFiles.empty())
        return;

    // Custom Structure dialog accepts multiple files (structure + model).
    if (state.fileBrowser.isCustomStructureDialogOpen())
    {
        for (const auto& f : state.pendingDroppedFiles)
            state.fileBrowser.feedDropToCustomStructureDialog(f);
        state.pendingDroppedFiles.clear();
        return;
    }

    if (state.fileBrowser.isMergeStructuresDialogOpen())
    {
        for (const auto& f : state.pendingDroppedFiles)
            state.fileBrowser.feedDropToMergeStructuresDialog(f);
        state.pendingDroppedFiles.clear();
        return;
    }

    const std::string droppedFile = state.pendingDroppedFiles.back();
    state.pendingDroppedFiles.clear();

    if (state.fileBrowser.isNanoCrystalDialogOpen())
    {
        // Nano dialog accepts both structure references and OBJ/STL model volumes.
        state.fileBrowser.feedDropToNanoCrystalDialog(droppedFile);
        return;
    }

    if (state.fileBrowser.isCSLGrainBoundaryDialogOpen())
    {
        state.fileBrowser.feedDropToCSLGrainBoundaryDialog(droppedFile);
        return;
    }

    if (state.fileBrowser.isInterfaceBuilderDialogOpen())
    {
        state.fileBrowser.feedDropToInterfaceBuilderDialog(droppedFile);
        return;
    }

    if (state.fileBrowser.isPolyCrystalDialogOpen())
    {
        state.fileBrowser.feedDropToPolyCrystalDialog(droppedFile);
        return;
    }

    if (state.fileBrowser.isSubstitutionalSolidSolutionDialogOpen())
    {
        state.fileBrowser.feedDropToSubstitutionalSolidSolutionDialog(droppedFile);
        return;
    }

    Structure loadedStructure;
    std::string loadError;
    if (!loadStructureFromFile(droppedFile, loadedStructure, loadError))
    {
        std::cout << "[Operation] Drop-load failed: " << droppedFile
                  << " (" << loadError << ")" << std::endl;
        state.fileBrowser.showLoadError(loadError);
        return;
    }

    state.structure = std::move(loadedStructure);
    state.fileBrowser.initFromPath(droppedFile);
    state.fileBrowser.applyElementColorOverrides(state.structure);
    state.fileBrowser.showLoadInfo(std::string("Structure loaded. ") + state.structure.ipfLoadStatus);
    updateBuffers(state);
    state.pendingDefaultViewReset = true;

    std::cout << "[Operation] Drop-loaded structure: " << droppedFile
              << " (atoms=" << state.structure.atoms.size() << ")" << std::endl;
}
