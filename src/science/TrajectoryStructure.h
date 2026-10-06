#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Time-averaged structure of a trajectory file (streamed frame by frame):
// total and partial radial distribution functions g_ab(r) (b around a),
// coordination numbers within `cutoff_A` (per element and per element pair, with
// their distribution) and the bond-angle distribution of neighbour pairs.
Json trajectoryStructure(const Parameters& p);
}
