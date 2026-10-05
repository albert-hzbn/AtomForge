#pragma once

#include "model/Structure.h"
#include "science/ScienceCore.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace atomforge::science
{
// Atoms for force evaluation. Units: Angstrom, amu.
struct Configuration
{
    std::vector<std::string> symbols;
    std::vector<int> numbers;
    std::vector<Vec3> positions;
    std::vector<double> masses;
    Mat3 cell{};
    Pbc pbc{};

    std::size_t size() const { return positions.size(); }
    bool fullyPeriodic() const { return pbc[0] && pbc[1] && pbc[2]; }
};

Configuration configurationFrom(const Structure& structure, const Pbc& pbc);
Structure structureFrom(const Configuration& configuration);

struct PotentialResult
{
    double energy = 0.0;       // eV
    std::vector<Vec3> forces;  // eV/Angstrom
    Mat3 stress{};             // eV/Angstrom^3, tension positive (ASE sign convention)
    bool hasStress = false;
};

// Native interatomic potential; implementations must be stateless between calls.
class Potential
{
public:
    virtual ~Potential() = default;
    virtual PotentialResult compute(const Configuration& configuration, bool stress) const = 0;
    virtual std::string description() const = 0;
};

using PotentialFactory = std::function<std::unique_ptr<Potential>()>;

// Effective-medium theory with ASE's parameters (Jacobsen, Stoltze, Norskov 1996)
// for Al, Cu, Ag, Au, Ni, Pd, Pt and the illustrative H, C, N, O set.
std::unique_ptr<Potential> makeEmt();
// 12-6 Lennard-Jones, energy shifted to zero at the cutoff like ASE's LennardJones.
std::unique_ptr<Potential> makeLennardJones(double epsilon, double sigma, double cutoff);

// Parses {"potential": "EMT"} or {"potential": "LennardJones", "epsilon", "sigma", "cutoff"};
// ASE-style {"module", "attribute", "kwargs"} descriptions of EMT/LennardJones are accepted.
PotentialFactory potentialFactory(const Json& specification);
}
