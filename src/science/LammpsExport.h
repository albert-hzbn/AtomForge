#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// LAMMPS restricted-triclinic box for a cell: rows (lx,0,0), (xy,ly,0),
// (xz,yz,lz), lattice-reduced so |xy| <= lx/2, |xz| <= lx/2, |yz| <= ly/2.
// Fractional coordinates are preserved; returns the new cell rows.
Mat3 lammpsCell(const Mat3& cell);

// The "lammps-export" tool: data file and input script for minimise, NVT, NPT
// or NEB with the selected potential, in result["files"].
ToolOutput lammpsExport(const Parameters& p);
}
