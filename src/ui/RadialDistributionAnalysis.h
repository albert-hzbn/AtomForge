#pragma once

#include "algorithms/RadialDistributionAnalysis.h"
#include "util/BackgroundTask.h"

struct RadialDistributionAnalysisDialog
{
    void drawMenuItem(bool enabled);
    void drawDialog(const Structure& structure);

private:
    bool m_openRequested = false;
    bool m_computeCompleted = false;
    atomforge::analysis::RdfParams m_params;
    atomforge::BackgroundTask<atomforge::analysis::RdfResult> m_task;
    bool m_showRawCounts = false;
    bool m_showCumulative = false;
    int m_refSpeciesIndex = 0;
    int m_targetSpeciesIndex = 0;
};
