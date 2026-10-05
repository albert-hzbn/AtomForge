#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// HPKOT high-symmetry path (SeeK-path convention) in the standardized primitive
// cell, from spglib symmetry. Requires a build with spglib.
ToolOutput reciprocalPath(const StructureInput& structure, double spacingInvA, double symprecA, bool timeReversal);

// Band-structure inputs on the HPKOT path: VASP KPOINTS (line mode) and POSCAR
// of the standardized primitive cell, and Quantum ESPRESSO K_POINTS crystal_b,
// CELL_PARAMETERS, ATOMIC_SPECIES and ATOMIC_POSITIONS blocks. File texts are
// returned in result["files"].
ToolOutput dftInputs(const StructureInput& structure, int pointsPerSegment, double symprecA, bool timeReversal);
}
