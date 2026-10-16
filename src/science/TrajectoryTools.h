#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Time-averaged structure of a trajectory (positions unwrapped across periodic
// boundaries before averaging) with each atom's RMS fluctuation.
ToolOutput trajectoryAverage(const Parameters& p);

// Unwrapped paths of the given atoms through every `stride`-th frame, each
// sampled frame continuing the previous one by the minimum-image step (the
// sampling must resolve the motion); `frames` receives the frame indices.
std::vector<std::vector<Vec3>> unwrappedPaths(const TrajectoryStream& stream, const std::vector<std::size_t>& atoms,
                                              std::size_t stride, std::vector<std::size_t>& frames);

// Unwrapped (continuous) trajectory frames, optionally smoothed by a centred
// moving average over a window of frames.
ToolOutput unwrapTrajectory(const Parameters& p);

// A per-atom analysis (structure types, centrosymmetry, entropy, ...) run on
// every frame: type counts or the mean property versus time, plus the cell
// volume and number density.
Json timeSeries(const Parameters& p);

// The change of a per-atom quantity between two frames: the displacement
// (minimum image) or any per-atom analysis (after minus reference).
Json frameDifference(const Parameters& p);
}
