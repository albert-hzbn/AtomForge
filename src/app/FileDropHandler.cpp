#include "app/FileDropHandler.h"

#include "app/EditorOps.h"
#include "app/EditorState.h"

#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <utility>

namespace
{
static std::vector<std::string>* s_pendingDrops = nullptr;

void dropFileCallback(GLFWwindow*, int count, const char** paths)
{
    if (count <= 0 || paths == nullptr || !s_pendingDrops)
        return;

    for (int i = 0; i < count; ++i)
    {
        if (paths[i] != nullptr && paths[i][0] != '\0')
            s_pendingDrops->push_back(paths[i]);
    }
}
}

void installDropFileCallback(GLFWwindow* window, std::vector<std::string>& pendingDrops)
{
    s_pendingDrops = &pendingDrops;
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

    // No special dialog is open — process every dropped file.
    std::vector<std::string> pendingFiles = std::move(state.pendingDroppedFiles);
    state.pendingDroppedFiles.clear();

    for (const auto& droppedFile : pendingFiles)
    {
        if (state.fileBrowser.isNanoCrystalDialogOpen())
        {
            // Nano dialog accepts both structure references and OBJ/STL model volumes.
            state.fileBrowser.feedDropToNanoCrystalDialog(droppedFile);
            continue;
        }

        if (state.fileBrowser.isCSLGrainBoundaryDialogOpen())
        {
            state.fileBrowser.feedDropToCSLGrainBoundaryDialog(droppedFile);
            continue;
        }

        if (state.fileBrowser.isInterfaceBuilderDialogOpen())
        {
            state.fileBrowser.feedDropToInterfaceBuilderDialog(droppedFile);
            continue;
        }

        if (state.fileBrowser.isPolyCrystalDialogOpen())
        {
            state.fileBrowser.feedDropToPolyCrystalDialog(droppedFile);
            continue;
        }

        if (state.fileBrowser.isStackingFaultDialogOpen())
        {
            state.fileBrowser.feedDropToStackingFaultDialog(droppedFile);
            continue;
        }

        if (state.fileBrowser.isSubstitutionalSolidSolutionDialogOpen())
        {
            state.fileBrowser.feedDropToSubstitutionalSolidSolutionDialog(droppedFile);
            continue;
        }

        Structure loadedStructure;
        std::string loadError;
        if (!loadStructureFromFile(droppedFile, loadedStructure, loadError))
        {
            std::cout << "[Operation] Drop-load failed: " << droppedFile
                      << " (" << loadError << ")" << std::endl;
            state.fileBrowser.showLoadError(loadError);
            continue;
        }

        // Signal the main loop to load this into a (possibly new) tab.
        state.pendingExternalLoadPaths.push_back(droppedFile);
        std::cout << "[Operation] Drop-queued for load: " << droppedFile << std::endl;
    }
}
