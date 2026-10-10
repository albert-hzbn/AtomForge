#pragma once

#include "pipeline/Pipeline.h"
#include "ui/StepEditing.h"

#include <array>
#include <functional>
#include <string>
#include <vector>

// One-step operation dialog: a pipeline step applied directly to the active
// structure from the Build or Edit menu (see MenuParity). Selection steps set
// the selection in the view; other steps change the structure as one undoable
// edit. "Add to pipeline" hands the same step to the Structure pipeline panel.
class OperationDialog
{
public:
    // Menu items of `menu` ("Build") or of an Edit submenu ("Structure Operations", "Select Atoms").
    void drawMenuItems(const char* menu, const char* submenu, bool hasStructure);
    bool open(const std::string& step);
    // Opens the dialog on a pipeline step: filled from its parameters, with
    // "Update step" handing them back instead of changing the structure.
    bool editStep(StepEdit edit);
    bool editingStep() const { return m_open && m_edit.active(); }
    void close() { m_open = false; m_edit.finish(); }
    bool isOpen() const { return m_open; }
    // Screen rectangle (x0, y0, x1, y1) of "Update step" in the last frame, for interface tests.
    const std::array<float, 4>& updateButton() const { return m_updateButton; }
    const std::string& step() const { return m_step; }
    const atomforge::pipeline::Json& parameters() const { return m_parameters; }

    struct Callbacks
    {
        std::function<void(Structure&)> update;                                   // structure edited (undoable)
        std::function<void(const std::vector<int>&)> select;                      // atoms to select in the view
        std::function<void(const std::vector<atomforge::pipeline::Modifier>&)> addToPipeline;
        std::function<void()> showAtomProperty;
    };
    // viewSelection: atom indices selected in the view.
    void draw(Structure& structure, const std::vector<int>& viewSelection, const Callbacks& callbacks);

    // The step's work, without drawing: applies the dialog's step to `structure`
    // with the chosen atoms; returns the new selection for selection steps.
    std::vector<int> apply(Structure& structure, const std::vector<int>& viewSelection);
    // Which atoms an operation acts on.
    enum Target { AllAtoms = 0, ViewSelection = 1, Condition = 2 };
    void setTarget(Target target, const std::string& condition = "");
    void setParameter(const std::string& name, const atomforge::pipeline::Json& value) { m_parameters[name] = value; }
    // Steps that act on a choice of atoms (delete, assign element, displace, ...).
    static bool actsOnAtoms(const std::string& step);
    static bool isSelectionStep(const std::string& step);

private:
    bool m_open = false;
    std::string m_step;
    atomforge::pipeline::Json m_parameters = atomforge::pipeline::Json::object();
    int m_target = AllAtoms;
    std::array<char, 512> m_condition{};
    std::string m_message;
    bool m_error = false;
    StepEdit m_edit;
    std::array<float, 4> m_updateButton{};
};
