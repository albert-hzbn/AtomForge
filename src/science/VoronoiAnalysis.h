#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Voronoi tessellation of each atom's environment: atomic volume, coordination
// (faces above the area threshold) and Voronoi index (faces by edge count,
// short edges below the edge threshold ignored). Cells are cut from all
// neighbours (periodic images included) until no further atom can cut them.
Json voronoiAnalysis(const Parameters& p);

// Adaptive common-neighbour analysis (Stukowski, 2012): each atom's cutoff
// follows its own 12 (fcc, hcp, icosahedral) or 14 (bcc) nearest neighbours;
// fixed-cutoff CNA when cutoff_A is given.
Json adaptiveCna(const Parameters& p);
}
