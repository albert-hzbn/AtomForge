#pragma once

#include "algorithms/CommonNeighbourAnalysis.h"
#include "util/BackgroundTask.h"

struct CommonNeighbourAnalysisDialog
{
    void drawMenuItem(bool enabled);
    void drawDialog(const Structure& structure);

private:
    bool m_openRequested = false;
    bool m_computeCompleted = false;
    atomforge::analysis::CnaParams m_params;
    atomforge::BackgroundTask<atomforge::analysis::CnaResult> m_task;
};
