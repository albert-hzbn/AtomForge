#pragma once

#include "pipeline/PipelineEditor.h"
#include "ui/PathPicker.h"
#include "ui/StepEditing.h"

#include <array>
#include <functional>
#include <string>
#include <utility>
#include <vector>

// The desktop "Structure pipeline" panel: a movable window listing the input
// structure and the modifiers it flows through. Modifiers are added from a
// menu, edited in place, enabled or disabled, reordered by dragging (or with
// the arrows) and deleted; the viewport shows the output, the input is kept.
class PipelineDialog
{
public:
    void drawMenuItem(bool hasStructure);
    // update shows a structure in the viewport; showAtomProperty switches the
    // view to atom-property colouring (used by "Colour by selection").
    void draw(Structure& structure, const std::function<void(Structure&)>& update,
              const std::function<void()>& showAtomProperty = {});
    // True while the pipeline writes its output into the structure, so that the
    // change stays out of undo history.
    bool isUpdating() const { return m_updating; }
    bool isActive() const { return m_editor.active(); }
    // Appends steps (from a menu operation), starting the pipeline on `structure`
    // when none is active, and opens the panel.
    void addSteps(const std::vector<atomforge::pipeline::Modifier>& steps, const Structure& structure);

    // "Edit in dialog" on a step (or a double click on it) asks for the Build or
    // Edit dialog of that step; the application opens it with this request,
    // whose commit writes the dialog's settings back into the step.
    bool consumeDialogRequest(StepEdit& request);
    void requestDialog(std::size_t index);
    // Writes parameters into step `index` if it is still a `step` step.
    void updateStep(std::size_t index, const std::string& step, const atomforge::pipeline::Json& parameters);

    atomforge::pipeline::PipelineEditor& editor() { return m_editor; }
    const atomforge::pipeline::PipelineEditor& editor() const { return m_editor; }
    // Project state: the pipeline (JSON text) and options; the input structure
    // is stored by the project separately.
    std::string snapshot() const;
    void restore(const std::string& state, const Structure* input);
    // Screen rectangles (x0, y0, x1, y1) of the step controls drawn in the last
    // frame, by name ("step 2", "up 2", "down 2", "remove 2", "move up", ...),
    // so interface tests can click and drag them.
    const std::vector<std::pair<std::string, std::array<float, 4>>>& controls() const { return m_controls; }

private:
    void show(Structure& structure, const std::function<void(Structure&)>& update, const std::function<void()>& showAtomProperty);
    void drawList();
    void drawEditor();
    void drawText();
    void drawAddMenu(const Structure& structure);
    void record(const std::string& name);
    void moveStep(std::size_t from, std::size_t to);
    void removeStep(std::size_t index);
    void duplicateStep(std::size_t index);
    bool m_open = false;
    bool m_updating = false;
    bool m_autoUpdate = true;
    bool m_colourSelection = false;
    bool m_showText = false;
    int m_selected = -1;
    int m_pickerAction = 0;  // 1 load, 2 save
    std::array<char, 4096> m_textBuffer{};
    std::string m_message;
    std::vector<std::pair<std::string, std::array<float, 4>>> m_controls;
    StepEdit m_dialogRequest;
    bool m_dialogRequested = false;
    atomforge::pipeline::PipelineEditor m_editor;
    PathPicker m_picker;
};
