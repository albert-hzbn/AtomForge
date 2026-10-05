#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Trajectory statistics. Positions: Angstrom; velocities: Angstrom/fs; time: fs.
Json meanSquareDisplacement(const Parameters& p);
Json diffusionCoefficient(const Parameters& p);
Json velocityAutocorrelation(const Parameters& p);
Json vibrationalSpectrum(const Parameters& p);

// Local environment and defect analyses with explicit periodic geometry.
Json localStrain(const Parameters& p);
Json centrosymmetry(const Parameters& p);
// Ackland-Jones bond-angle structure type (other, fcc, hcp, bcc, icosahedral).
Json structureType(const Parameters& p);
// Classification of one environment from its neighbour vectors (any order).
int acklandJonesType(std::vector<Vec3> bonds);
Json bondOrder(const Parameters& p);
Json wignerSeitz(const Parameters& p);
Json staticStructureFactor(const Parameters& p);

// Band edges, local effective masses and vacuum-referenced work functions.
Json bandGap(const Parameters& p);
Json effectiveMass(const Parameters& p);
Json workFunction(const Parameters& p);

// Equation of state, linear elastic response and harmonic thermal properties.
Json equationOfState(const Parameters& p);
Json elasticTensor(const Parameters& p);
Json phononDos(const Parameters& p);
Json harmonicThermodynamics(const Parameters& p);
}
