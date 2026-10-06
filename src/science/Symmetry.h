#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Space group, Wyckoff positions and standardized cells of a periodic structure
// (spglib). Frames: the conventional standard cell, the primitive standard cell
// and the input cell with symmetrized positions; the one named by `output_cell`
// comes last, so "Open final structure" loads it.
ToolOutput symmetryAnalysis(const StructureInput& structure, const Parameters& p);
}
