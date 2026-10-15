#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// A per-atom value aggregated on a 1D, 2D or 3D grid of bins along the cell
// vectors (mean, sum, density = sum per bin volume, min, max or atom count).
Json spatialBinning(const Parameters& p);

// Histogram and reductions (count, sum, mean, standard deviation, min,
// median, max) of a per-atom value, optionally of selected atoms only.
Json propertyHistogram(const Parameters& p);

// Spatial correlation <a_i b_j> of two per-atom values as a function of the
// distance between atoms i and j, computed on an FFT grid in a periodic cell.
Json spatialCorrelation(const Parameters& p);
}
