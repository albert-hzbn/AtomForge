#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Space group, Wyckoff positions and standardized cells of a periodic structure
// (spglib). Frames: the conventional standard cell, the primitive standard cell
// and the input cell with symmetrized positions; the one named by `output_cell`
// comes last, so "Open final structure" loads it.
ToolOutput symmetryAnalysis(const StructureInput& structure, const Parameters& p);

// The matrix taking Miller indices of a cell to those of its conventional
// standard cell (spglib), h_conventional = h * transform with h a row: a
// supercell's (2 1 1) is the conventional (1 1 1). False without spglib or
// when the symmetry cannot be found.
bool conventionalIndexTransform(const Mat3& cell, const std::vector<Vec3>& fractions, const std::vector<int>& numbers,
                                double symprec, Mat3& transform);
}
