#pragma once

#include <vector>
#include <functional>
#include <glm/glm.hpp>

#include "ElementData.h"
#include "io/StructureLoader.h"
#include "graphics/SceneBuffers.h"

// Manages the right-click context menu and the periodic-table picker that
// drives Substitute and Insert-at-Midpoint actions.

// Requests produced by the context menu (and forwarded from the View menu).
// All fields start false; main.cpp ORs them with the per-frame booleans.
struct AtomRequests {
    bool doDelete        = false;
    bool measureDistance = false;
    bool measureAngle    = false;
    bool atomInfo        = false;
};

struct AtomContextMenu
{
    // Signal that the menu should open on the next draw() call.
    void open();

    // Call once per ImGui frame.
    // doDelete is set to true when the user picks "Delete" from the menu.
    // selectedInstanceIndices may be cleared by the "Deselect" action.
    void draw(Structure& structure,
              SceneBuffers& sceneBuffers,
              const std::vector<glm::vec3>& elementColors,
              std::vector<int>& selectedInstanceIndices,
              AtomRequests& requests,
              const std::function<void(Structure&)>& updateBuffers);

private:
    bool           m_openRequested = false;
    bool           m_ownsPeriodicPopup = false;
    PeriodicAction m_pendingAction = PeriodicAction::None;
};
