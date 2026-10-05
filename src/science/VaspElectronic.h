#pragma once

#include "science/ScienceData.h"

#include <filesystem>
#include <string>
#include <vector>

namespace atomforge::science
{
struct Eigenvalues
{
    std::vector<Vec3> kpoints;  // fractional
    std::vector<double> weights;
    // energies[spin][k][band] (eV)
    std::vector<std::vector<std::vector<double>>> energies;
};
Eigenvalues readEigenval(const std::filesystem::path& path);

struct DensityOfStates
{
    double fermi = 0;
    std::vector<double> energy;               // eV
    std::vector<std::vector<double>> total;   // [spin][point]
    // Site-projected DOS summed over orbitals: [atom][spin][point]; empty when absent.
    std::vector<std::vector<std::vector<double>>> projected;
};
DensityOfStates readDoscar(const std::filesystem::path& path);

// Labels of a line-mode KPOINTS file: one entry per segment end point,
// in order (start, end, start, end, ...).
std::vector<std::string> readLineModeLabels(const std::filesystem::path& path, int& pointsPerSegment);

struct Projections
{
    // weight[spin][k][band][ion][orbital]; orbital names in `orbitals`.
    std::vector<std::string> orbitals;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> weight;
};
Projections readProcar(const std::filesystem::path& path);

// The "vasp-electronic" tool: band structure, DOS and fat-band weights.
ToolOutput vaspElectronic(const Parameters& p);
}
