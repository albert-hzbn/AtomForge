#pragma once

#include "model/Structure.h"
#include "science/ScienceCore.h"

#include <filesystem>
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

// Embedded-atom method from LAMMPS-format tables: "setfl" (eam/alloy), "fs"
// (eam/fs, Finnis-Sinclair) or "funcfl" (single-element eam). Tables are
// interpolated with LAMMPS's cubic scheme; pair terms are stored as r*phi.
std::unique_ptr<Potential> makeEam(const std::filesystem::path& file, const std::string& format = "auto");

// Tersoff bond-order potential (LAMMPS pair_style tersoff file format); an empty
// path gives the built-in Si parameters (Tersoff 1988).
std::unique_ptr<Potential> makeTersoff(const std::filesystem::path& file);
// Stillinger-Weber (LAMMPS pair_style sw file format); empty path: built-in Si (1985).
std::unique_ptr<Potential> makeStillingerWeber(const std::filesystem::path& file);
// Buckingham A exp(-r/rho) - C/r^6 pairs plus Coulomb interactions of fixed
// charges, Ewald-summed for periodic cells (LAMMPS buck/coul/long):
// {"pairs": [{"elements": [a, b], "A", "rho", "C"}], "charges": {element: q},
//  "cutoff": 10, "ewald_accuracy": 1e-6}.
std::unique_ptr<Potential> makeBuckingham(const Json& options);

// Parses {"potential": "EMT"}, {"potential": "LennardJones", "epsilon", "sigma", "cutoff"} or
// {"potential": "EAM", "file": PATH, "format": "auto|setfl|fs|funcfl"},
// {"potential": "Tersoff" | "StillingerWeber", "file": PATH (optional)} or
// {"potential": "Buckingham", ...}; relative files are
// resolved against base. ASE-style {"module", "attribute", "kwargs"} descriptions of
// EMT/LennardJones are accepted.
PotentialFactory potentialFactory(const Json& specification, const std::filesystem::path& base = ".");
}
