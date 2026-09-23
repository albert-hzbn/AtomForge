// Validates detectDislocationLattice/detectStackingFaultFamily's common-neighbor
// classification against ideal FCC, BCC and HCP lattices built directly from their
// textbook fractional coordinates (not loaded from any file), so a wrong answer here
// points at the classification algorithm itself rather than a structure-builder bug.
//
// Reference common-neighbor-analysis (Honeycutt-Andersen) pair signatures used by
// the implementation were hand- and numerically-verified against the ideal lattices'
// nearest-neighbor geometry (see the fix that accompanies this test): FCC is pure
// 4-2-1; HCP is an even 4-2-1 / 4-2-2 split; BCC (14-neighbor coordination under the
// same covalent-radius cutoff) is 6-6-3 (first shell) plus 4-4-2 (second shell).
#include "algorithms/DislocationBuilder.h"
#include "model/Structure.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
Structure buildFcc(double a, int n)
{
    const double frac[4][3] = {{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
    Structure s;
    s.hasUnitCell = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            s.cellVectors[i][j] = (i == j) ? a * n : 0.0;
    for (int ix = 0; ix < n; ++ix)
    for (int iy = 0; iy < n; ++iy)
    for (int iz = 0; iz < n; ++iz)
    for (int k = 0; k < 4; ++k)
    {
        AtomSite atom;
        atom.symbol = "Cu"; atom.atomicNumber = 29;
        atom.x = (ix + frac[k][0]) * a;
        atom.y = (iy + frac[k][1]) * a;
        atom.z = (iz + frac[k][2]) * a;
        s.atoms.push_back(atom);
    }
    return s;
}

Structure buildBcc(double a, int n)
{
    const double frac[2][3] = {{0,0,0},{0.5,0.5,0.5}};
    Structure s;
    s.hasUnitCell = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            s.cellVectors[i][j] = (i == j) ? a * n : 0.0;
    for (int ix = 0; ix < n; ++ix)
    for (int iy = 0; iy < n; ++iy)
    for (int iz = 0; iz < n; ++iz)
    for (int k = 0; k < 2; ++k)
    {
        AtomSite atom;
        atom.symbol = "Fe"; atom.atomicNumber = 26;
        atom.x = (ix + frac[k][0]) * a;
        atom.y = (iy + frac[k][1]) * a;
        atom.z = (iz + frac[k][2]) * a;
        s.atoms.push_back(atom);
    }
    return s;
}

// Ideal c/a = sqrt(8/3), where all 12 nearest neighbors are exactly equidistant.
Structure buildHcp(double a, int n, int nz)
{
    const double c = a * std::sqrt(8.0 / 3.0);
    const double a1[3] = {a, 0, 0};
    const double a2[3] = {-0.5 * a, 0.5 * a * std::sqrt(3.0), 0};
    const double a3[3] = {0, 0, c};
    const double basisFrac[2][3] = {{0,0,0}, {2.0/3.0, 1.0/3.0, 0.5}};

    Structure s;
    s.hasUnitCell = true;
    for (int j = 0; j < 3; ++j) s.cellVectors[0][j] = a1[j] * n;
    for (int j = 0; j < 3; ++j) s.cellVectors[1][j] = a2[j] * n;
    for (int j = 0; j < 3; ++j) s.cellVectors[2][j] = a3[j] * nz;

    for (int ix = 0; ix < n; ++ix)
    for (int iy = 0; iy < n; ++iy)
    for (int iz = 0; iz < nz; ++iz)
    for (int k = 0; k < 2; ++k)
    {
        const double fx = ix + basisFrac[k][0];
        const double fy = iy + basisFrac[k][1];
        const double fz = iz + basisFrac[k][2];
        AtomSite atom;
        atom.symbol = "Mg"; atom.atomicNumber = 12;
        atom.x = fx * a1[0] + fy * a2[0] + fz * a3[0];
        atom.y = fx * a1[1] + fy * a2[1] + fz * a3[1];
        atom.z = fx * a1[2] + fy * a2[2] + fz * a3[2];
        s.atoms.push_back(atom);
    }
    return s;
}

void expectFamily(const char* label, const Structure& s, DislocationLatticeFamily expected,
                  double minFractionPbc, double minFractionNoPbc)
{
    for (bool pbc : {true, false})
    {
        // The non-PBC path has no periodic images, so atoms near any of the
        // finite slab's surfaces genuinely lack a full coordination shell and
        // are correctly left unrecognized -- a real, expected effect of
        // boundary truncation, not a classification bug. Only the PBC path
        // (full bulk periodicity, no surfaces) should recognize essentially
        // every atom.
        const double minFraction = pbc ? minFractionPbc : minFractionNoPbc;
        const DislocationDetectionResult r = detectDislocationLattice(s, pbc);
        if (r.family != expected)
        {
            throw std::runtime_error(std::string(label) + " (pbc=" + (pbc ? "true" : "false")
                + "): expected " + dislocationLatticeFamilyName(expected)
                + ", got " + dislocationLatticeFamilyName(r.family));
        }
        const double fraction = (double)r.recognizedCount / (double)s.atoms.size();
        if (fraction < minFraction)
        {
            throw std::runtime_error(std::string(label) + " (pbc=" + (pbc ? "true" : "false")
                + "): only " + std::to_string(fraction) + " of atoms recognized, expected >= "
                + std::to_string(minFraction));
        }
        // A correctly classified ideal lattice should have essentially none of its
        // recognized atoms attributed to the other two families.
        int other = r.fccCount + r.hcpCount + r.bccCount
                  - (expected == DislocationLatticeFamily::Fcc ? r.fccCount
                     : expected == DislocationLatticeFamily::Hcp ? r.hcpCount : r.bccCount);
        if (other != 0)
        {
            throw std::runtime_error(std::string(label) + " (pbc=" + (pbc ? "true" : "false")
                + "): " + std::to_string(other) + " atoms misclassified into another family");
        }
    }
}
}

int main()
{
    try
    {
        // With full periodic images (box >> interaction cutoff in every direction),
        // every single atom in an ideal lattice should be classified correctly.
        // Without PBC, a lower bound just confirms the majority of bulk atoms
        // (away from the finite slab's surfaces) are still recognized correctly.
        expectFamily("FCC", buildFcc(3.61, 10), DislocationLatticeFamily::Fcc, 0.99, 0.90);
        expectFamily("BCC", buildBcc(2.87, 12), DislocationLatticeFamily::Bcc, 0.99, 0.75);
        expectFamily("HCP", buildHcp(3.21, 14, 6), DislocationLatticeFamily::Hcp, 0.99, 0.65);

        std::cout << "Lattice family regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
