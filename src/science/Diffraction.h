#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Kinematic powder X-ray diffraction: structure factors from tabulated
// atomic scattering factors, optional isotropic Debye-Waller factor,
// Lorentz-polarization correction and Gaussian-broadened profile.
ToolOutput powderXrd(const StructureInput& structure, const Parameters& p);

// Kinematic zone-axis electron diffraction spot pattern (zero-order Laue
// zone) with Mott-Bethe electron scattering factors and a relativistic
// electron wavelength.
ToolOutput electronDiffraction(const StructureInput& structure, const Parameters& p);

// X-ray atomic scattering factor f(s) and electron factor f_e(s) (Angstrom)
// for s = sin(theta)/lambda in 1/Angstrom.
double xrayScatteringFactor(const std::string& symbol, double s);
double electronScatteringFactor(const std::string& symbol, double s);
}
