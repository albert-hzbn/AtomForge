#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Voronoi tessellation of each atom's environment: atomic volume, coordination
// (faces above the area threshold) and Voronoi index (faces by edge count,
// short edges below the edge threshold ignored). Cells are cut from all
// neighbours (periodic images included) until no further atom can cut them.
Json voronoiAnalysis(const Parameters& p);

// VoroTop (Lazar, Han and Srolovitz, 2015): each atom's structure from the
// topology of its Voronoi cell (canonical Weinberg code) looked up in a
// filter: a VoroTop filter file, or one sampled from thermally perturbed fcc,
// bcc and hcp crystals.
Json voroTop(const Parameters& p);

// Adaptive common-neighbour analysis (Stukowski, 2012): each atom's cutoff
// follows its own 12 (fcc, hcp, icosahedral) or 14 (bcc) nearest neighbours;
// fixed-cutoff CNA when cutoff_A is given.
Json adaptiveCna(const Parameters& p);
}
