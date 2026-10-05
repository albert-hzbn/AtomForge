#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Dislocation lines and Burgers vectors of a deformed structure relative to
// its perfect reference (same atom order), in the spirit of the dislocation
// extraction algorithm (Stukowski and Albe 2010). Core atoms are those whose
// Ackland-Jones structure type changed (slip by a lattice vector leaves atoms
// crystalline; free surfaces are excluded because they are already
// non-crystalline in the reference). Core atoms are clustered into lines; the
// direction is the lattice period of a cluster wrapping the cell, else its
// principal axis. The Burgers vector is the closure failure of a right-handed
// Burgers circuit about the line, with each displacement step reduced modulo
// the reference lattice's nearest-neighbour vectors.
ToolOutput dislocationLines(const StructureInput& reference, const StructureInput& deformed, const Parameters& p);
}
