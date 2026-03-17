#pragma once

#include <vector>
#include <functional>
#include <glm/glm.hpp>

#include "io/StructureLoader.h"

// Owns the per-element radius and colour tables, and draws the Edit menu
// items + the "Atomic Sizes" and "Element Colors" modal dialogs.
struct EditMenuDialogs
{
    // Public state so the updateBuffers lambda in main can read these.
    std::vector<float>     elementRadii;
    std::vector<glm::vec3> elementColors;
    std::vector<float>     elementShininess;

    // Initialise radii and colours to literature / CPK defaults.
    EditMenuDialogs();

    // Call inside an already-open "Edit" ImGui menu to add the Edit items.
    void drawMenuItems();

    // Call once per ImGui frame (outside any menu) to service the modal popups.
    // Calls updateBuffers(structure) when the user changes a value.
    void drawPopups(Structure& structure,
                    const std::function<void(Structure&)>& updateBuffers);

private:
    bool m_openAtomicSize   = false;
    bool m_openElementColor = false;
    bool m_openEditStructure = false;
    bool m_showEditStructureElementPicker = false;
    bool m_restoreEditStructureAfterElementPicker = false;
    bool m_useDirectCoords = false;
    int  m_selectedRadiusElement = 6;
    int  m_selectedColorElement  = 6;
    int  m_selectedEditAtom = 0;
    int  m_selectedEditElement = 6;
    int  m_editStructureElementTargetAtom = -1;
    bool m_scrollEditRowsToBottom = false;
};
