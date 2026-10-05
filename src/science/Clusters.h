#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Connected clusters of selected atoms (bond length below a cutoff, periodic
// images included), with sizes, unwrapped centroids and radii of gyration.
Json clusterAnalysis(const Parameters& p);

// Empty space accessible to a probe sphere: grid points farther than
// atomic radius + probe radius from every atom, grouped into connected voids.
// With open boundaries, regions touching the bounding box are exterior.
Json voidAnalysis(const Parameters& p);
}
