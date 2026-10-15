#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Bond orders from bond lengths (Pauling, J. Am. Chem. Soc. 69, 542 (1947)):
// atoms closer than tolerance x the sum of their covalent radii are bonded,
// with order n = 10^((r1 - r) / 0.6 A), r1 the single-bond length (sum of the
// covalent radii); n is also rounded to 1, 1.5, 2 or 3. Reports every bond and
// each atom's bond-valence sum.
Json bondOrders(const Parameters& p);

// LAMMPS regions (block, sphere, cylinder, cone, plane, prism; side in/out;
// units box or lattice) given as their region commands: which atoms each
// region holds, and the regions as a surface mesh (OBJ) within the cell.
Json lammpsRegions(const Parameters& p);
}
