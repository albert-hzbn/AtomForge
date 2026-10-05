#pragma once

#include "science/Potentials.h"
#include "science/ScienceData.h"

#include <array>
#include <vector>

namespace atomforge::science
{
// Real-space harmonic force constants from central finite displacements.
struct ForceConstants
{
    Configuration unit;            // unit cell whose phonons are computed
    std::array<int, 3> supercell{};
    Configuration super;           // supercell, unit-cell images first
    // phi[a][b'] = d2E / du_a du_b' (eV/A^2) for unit atom a and supercell atom b'.
    std::vector<std::vector<Mat3>> phi;
};

ForceConstants forceConstants(const Configuration& unit, const Potential& potential,
                              const std::array<int, 3>& supercell, double displacement);

// Supercell repetitions so each periodic image direction is at least minimumLength long.
std::array<int, 3> supercellFor(const Configuration& unit, double minimumLength);

// Phonon frequencies (THz, negative for imaginary modes) at Cartesian q (rad/A),
// using minimum-image weighting of the supercell force constants.
std::vector<double> phononFrequencies(const ForceConstants& constants, const Vec3& q);

constexpr double kThzToEv = 4.135667696e-3;  // h * 1 THz in eV

// The "phonons" tool: band structure along a k-path, mesh DOS and harmonic thermodynamics.
ToolOutput phononCalculation(const StructureInput& structure, const Potential& potential, const Parameters& p);
}
