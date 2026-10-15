#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Surface of a solid by the probe-sphere (alpha-shape-like) construction:
// space a probe sphere of the given radius can reach without enclosing an
// atom centre is empty, the rest is solid. Computed on a grid with exact
// distance transforms; the surface is the isosurface between solid and
// empty. Reports surface area, solid volume and fraction, the separate solid
// regions and empty pores, surface atoms, and the mesh as an OBJ file.
Json surfaceMesh(const Parameters& p);
}
