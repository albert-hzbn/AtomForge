#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// HPKOT high-symmetry path (SeeK-path convention) in the standardized primitive
// cell, from spglib symmetry. Requires a build with spglib.
ToolOutput reciprocalPath(const StructureInput& structure, double spacingInvA, double symprecA, bool timeReversal);
}
