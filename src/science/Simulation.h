#pragma once

#include "science/Potentials.h"
#include "science/ScienceData.h"

namespace atomforge::science
{
struct NebOptions
{
    long long images = 7;      // including both endpoints
    double fmax = 0.03;        // eV/Angstrom
    long long steps = 300;
    double spring = 0.1;       // eV/Angstrom^2
    bool climb = true;
    bool mic = false;
};

// Fixed-cell climbing-image NEB with the improved tangent (Henkelman and
// Jonsson 2000) and FIRE optimisation, with ASE's default FIRE parameters.
ToolOutput migrationPath(const Configuration& initial, const Configuration& final,
                         const PotentialFactory& factory, const NebOptions& options);

struct DynamicsOptions
{
    long long steps = 1000;
    double timestepFs = 1.0;
    double temperatureK = 300.0;
    double thermostatFs = 100.0;
    double pressureGPa = 0.0;
    double barostatFs = 1000.0;
    unsigned long long seed = 0;
    long long sampleInterval = 10;
};

// Langevin NVT (BAOAB splitting) with Maxwell-Boltzmann initial velocities.
ToolOutput nvtDynamics(const Configuration& start, const Potential& potential, const DynamicsOptions& options);
// Isotropic Martyna-Tobias-Klein NPT with Nose-Hoover chains (Tuckerman et al. 2006 splitting).
ToolOutput nptDynamics(const Configuration& start, const Potential& potential, const DynamicsOptions& options);
}
