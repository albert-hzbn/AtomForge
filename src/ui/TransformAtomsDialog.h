#pragma once

#include "ui/StepEditing.h"

#include <functional>

struct TransformAtomsDialog
{
    bool isEnabled() const { return useTransformMatrix; }
    const int (&getMatrix() const)[3][3] { return transformMatrix; }
    void clearTransform();

    void drawMenuItem(bool hasUnitCell);
    void drawDialog(const std::function<void()>& onApply);

    // Opens the dialog on a pipeline step (supercell): filled from the step's
    // settings, with "Update step" writing them back instead of applying.
    bool editStep(StepEdit edit);
    // The dialog's matrix as the supercell step's parameters ("matrix": nine
    // integers row by row).
    atomforge::pipeline::Json stepParameters() const;

private:
    bool showDialog = false;
    StepEdit m_step;
    bool useTransformMatrix = false;

    int transformMatrix[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    int pendingMatrix[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
};
