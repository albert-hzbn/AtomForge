#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Time-averaged structure of a trajectory (positions unwrapped across periodic
// boundaries before averaging) with each atom's RMS fluctuation.
ToolOutput trajectoryAverage(const Parameters& p);

// Unwrapped (continuous) trajectory frames, optionally smoothed by a centred
// moving average over a window of frames.
ToolOutput unwrapTrajectory(const Parameters& p);

// A per-atom analysis (structure types, centrosymmetry, entropy, ...) run on
// every frame: type counts or the mean property versus time, plus the cell
// volume and number density.
Json timeSeries(const Parameters& p);
}
