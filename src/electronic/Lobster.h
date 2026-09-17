#pragma once

#include <string>
#include <vector>

// Readers for chemical-bonding analysis output from LOBSTER (Local Orbital
// Basis Suite Towards Electronic-Structure Reconstruction; Nelson et al.,
// J. Comput. Chem. 2020; Dronskowski and Blochl, J. Phys. Chem. 1993 for the
// original COHP method). COHPCAR/COOPCAR/COBICAR give pCOHP/pCOOP/pCOBI
// curves and their energy integral (IpCOHP etc.) per bond; ICOHPLIST/
// ICOOPLIST/ICOBILIST give the same integral evaluated at the Fermi level,
// one row per bond, without the full curve. Column layouts follow LOBSTER's
// documented conventions exactly (https://www.cohp.de). When orbital-resolved
// sub-entries are present (LOBSTER's default output style, interleaving each
// bond's total with per-orbital-pair contributions), only the per-bond
// totals are kept; orbital-resolved data itself is not exposed. Multi-center
// COBI files are not supported and raise a clear error rather than being
// silently mis-parsed.
namespace atomforge::electronic
{
struct CohpBond { int atom1 = -1; int atom2 = -1; double length = 0.0; };  // 0-based site indices

struct CohpData
{
    bool spinPolarized = false;
    double fermiEnergy = 0.0;
    std::vector<double> energies;                        // size numEnergies, Fermi level at 0 eV
    std::vector<CohpBond> bonds;                          // size numBonds (excludes the "average" entry)
    std::vector<std::vector<double>> averageCohp;         // [spin][energy]; size 1 (nonpolarized) or 2 (up, down)
    std::vector<std::vector<double>> averageIcohp;        // [spin][energy]
    std::vector<std::vector<std::vector<double>>> cohp;   // [spin][bond][energy]
    std::vector<std::vector<std::vector<double>>> icohp;  // [spin][bond][energy]
};

struct IcohpEntry { int atom1 = -1; int atom2 = -1; double length = 0.0; int numBonds = 1; double icohp[2]{0.0, 0.0}; };

struct IcohpList { bool spinPolarized = false; std::vector<IcohpEntry> entries; };

CohpData readCohpcar(const std::string& path);
IcohpList readIcohplist(const std::string& path);
}
