#pragma once

#include <vector>
#include <functional>
#include <string>
#include <glm/glm.hpp>

#include "model/DisplaySettings.h"
#include "model/Structure.h"
#include "ui/StepEditing.h"

// Owns the per-element radius and colour tables, and draws the Edit menu
// items + the "Atomic Sizes" and "Display Settings" modal dialogs.
struct EditMenuDialogs
{
    // Public state so the updateBuffers lambda in main can read these.
    std::vector<float>     elementRadii;
    std::vector<glm::vec3> elementColors;
    std::vector<float>     elementShininess;

    // Lighting, surface, sizes and colours of the view (Display Settings).
    atomforge::display::DisplaySettings display;

    // Initialise radii and colours to literature / CPK defaults.
    EditMenuDialogs();

    // Call inside an already-open "Edit" ImGui menu to add the Edit items.
    void drawMenuItems();

    // Call inside an already-open "Settings" ImGui menu to add Atomic Sizes
    // and Display Settings items (without Edit Structure).
    void drawSettingsMenuItems();
    void openAtomicSizes() { m_openAtomicSize = true; }
    void openDisplaySettings() { m_openElementColor = true; }
    void openEditStructure() { m_openEditStructure = true; }

    // Call once per ImGui frame (outside any menu) to service the modal popups.
    // Calls updateBuffers(structure) when the user changes a value.
    void drawPopups(Structure& structure,
                    const std::function<void(Structure&)>& updateBuffers);

    // Opens Edit Structure on a pipeline step (set-cell or add-atom): filled
    // from the step's settings, with "Update step" writing them back instead of
    // editing the structure.
    bool editStep(StepEdit edit);
    // The dialog's current settings as the step's parameters.
    atomforge::pipeline::Json stepParameters() const;

private:
    bool m_openAtomicSize   = false;
    bool m_openElementColor = false;
    bool m_openEditStructure = false;
    bool m_showEditStructureElementPicker = false;
    bool m_restoreEditStructureAfterElementPicker = false;
    bool m_useDirectCoords = false;
    bool m_displayRebuildPending = false;
    std::string m_displayStatus;
    int  m_selectedRadiusElement = 6;
    int  m_selectedColorElement  = 6;
    int  m_selectedEditAtom = 0;
    int  m_selectedEditElement = 6;
    int  m_editStructureElementTargetAtom = -1;
    bool m_scrollEditRowsToBottom = false;

    // Pipeline step being edited (set-cell: lattice vectors; add-atom: one atom).
    StepEdit m_step;
    double m_stepCell[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    bool   m_stepScaleAtoms = true;
    int    m_stepElement = 29;
    double m_stepPosition[3] = {0, 0, 0};
    bool   m_stepFractional = false;
    bool   m_stepSelect = true;
    void drawStepEditor();
    void drawDisplaySettings(Structure& structure, const std::function<void(Structure&)>& updateBuffers);
    void drawRadiusEditor(Structure& structure, const std::function<void(Structure&)>& updateBuffers);
};
