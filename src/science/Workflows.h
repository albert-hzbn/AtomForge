#pragma once

#include "science/Potentials.h"
#include "science/ScienceData.h"

namespace atomforge::science
{
// Elastic constants from a structure and a potential: every Voigt strain
// component at +/- each amplitude, ions relaxed at each fixed strained cell, and
// the stress-strain fit of the elastic-tensor tool (stiffness, moduli, stability).
ToolOutput elasticConstants(const StructureInput& structure, const Potential& potential, const Parameters& p);

// Energy-volume scan: the cell is scaled isotropically over +/- volume_range,
// ions are relaxed at each volume and a Birch-Murnaghan equation of state is fitted.
ToolOutput equationOfStateScan(const StructureInput& structure, const Potential& potential, const Parameters& p);

// Point-defect formation energy E_f = E_defect - E_bulk - sum_i dn_i mu_i, with both
// cells relaxed (ions; optionally the defect cell's shape and volume too).
ToolOutput formationEnergy(const StructureInput& defect, const StructureInput& bulk, const Potential& potential, const Parameters& p);

// Energy per area of planar defects (free surfaces, grain boundaries, stacking
// faults): gamma = (E_cell - (N / N_bulk) E_bulk) / (interfaces * A), A spanned by
// the two cell vectors other than normal_axis.
ToolOutput planarDefectEnergy(const StructureInput& cell, const StructureInput& bulk, const Potential& potential, const Parameters& p);
}
