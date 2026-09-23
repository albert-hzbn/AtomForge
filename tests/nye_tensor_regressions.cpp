// Validates computeNyeTensor (src/algorithms/NyeTensor.h), a native
// reimplementation of BABEL's documented Hartley & Mishin (2005) Nye-tensor
// algorithm (see NyeTensor.h for the exact source references). BABEL's own
// worked example for this feature (Examples/BccScrewDipole/4_NyeMap) starts
// from an MD-relaxed structure this environment has no way to reproduce, so
// an exact numeric match to that example isn't attempted here (unlike
// babel_reference_regressions.cpp, which does reproduce BABEL numerically
// for the pure-elasticity dislocation displacement field). Instead this
// checks two properties that must hold for ANY correct implementation of
// the algorithm, independent of BABEL's specific example data:
//   1. A uniform (spatially constant) displacement gradient has zero curl
//      by construction, so the Nye tensor must be exactly zero everywhere
//      -- this is a strict analytic identity, not an approximation.
//   2. Around a real edge dislocation (built with AtomForge's own,
//      already-validated isotropic solver), the Nye tensor must be
//      concentrated near the dislocation line and small far from it.
#include "algorithms/DislocationBuilder.h"
#include "algorithms/NyeTensor.h"
#include "model/Structure.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge;

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
}

int main()
{
    try
    {
        const double a = 3.61;
        Structure reference = buildFcc(a, 13);

        // --- 1. Uniform strain: Nye tensor must be exactly (near machine
        // precision) zero everywhere, since the displacement gradient is
        // spatially constant.
        {
            Structure strained = reference;
            const double eps[3][3] = {{0.01, 0.002, 0.0}, {0.0, -0.005, 0.001}, {0.0, 0.0, 0.003}};
            auto applyStrain = [&](double x, double y, double z, double& outX, double& outY, double& outZ) {
                outX = x + eps[0][0]*x + eps[0][1]*y + eps[0][2]*z;
                outY = y + eps[1][0]*x + eps[1][1]*y + eps[1][2]*z;
                outZ = z + eps[2][0]*x + eps[2][1]*y + eps[2][2]*z;
            };
            for (auto& atom : strained.atoms)
                applyStrain(atom.x, atom.y, atom.z, atom.x, atom.y, atom.z);
            // A uniform strain applied to a PERIODIC crystal must also strain
            // the periodicity vectors themselves, or atoms near the periodic
            // boundary see a spurious discontinuity under minimum-image
            // wrapping against the unstrained cell.
            for (int i = 0; i < 3; ++i)
                applyStrain(reference.cellVectors[i][0], reference.cellVectors[i][1], reference.cellVectors[i][2],
                           strained.cellVectors[i][0], strained.cellVectors[i][1], strained.cellVectors[i][2]);

            const NyeTensorResult result = computeNyeTensor(reference, strained, 3.0f, true);
            if (!result.success)
                throw std::runtime_error("Uniform-strain Nye tensor computation failed: " + result.message);

            double maxNorm = 0.0;
            for (double norm : result.alphaNorm)
                maxNorm = std::max(maxNorm, norm);
            std::cout << "Uniform strain: max |alpha| = " << maxNorm << " (expect ~0)\n";
            if (maxNorm > 1e-8)
                throw std::runtime_error("Uniform strain produced a non-zero Nye tensor (max |alpha| = "
                    + std::to_string(maxNorm) + "); the curl computation has a bug");
        }

        // --- 2. Real edge dislocation, checked two ways:
        //   (a) atoms entirely outside the cylinder the displacement field
        //       was masked to must have EXACTLY zero Nye tensor -- they and
        //       every one of their neighbors received zero displacement, so
        //       the local displacement gradient (and hence its curl) is
        //       identically zero there. This is as exact as check 1 above.
        //   (b) atoms actually inside the displaced region, close enough to
        //       the line to be affected by real dislocation-core strain,
        //       must show a clearly non-zero Nye tensor (confirms the
        //       algorithm detects real signal, not just correctly reporting
        //       zero everywhere).
        {
            DislocationParams params;
            params.character = DislocationCharacter::Edge;
            params.shape = DislocationShape::Cylinder;
            params.cylinderRadius = 15.0f;
            params.coreRadius = 1.2f;
            params.cutoffRadius = 0.0f; // no attenuation, matches a "textbook" field out to the cylinder edge

            const DislocationResult dislo = buildDislocation(reference, params);
            if (!dislo.success)
                throw std::runtime_error("Failed to build reference edge dislocation: " + dislo.message);

            // A single dislocation has a non-zero net Burgers vector, so its
            // displacement field is not periodicity-compatible (same reason
            // DislocationBuilder.cpp's own post-build validation avoids PBC
            // neighbor search for this case) -- using PBC minimum-image
            // neighbors here would manufacture spurious wraparound
            // "neighbors" near every box face, not just near the line.
            const NyeTensorResult result = computeNyeTensor(reference, dislo.output, 3.0f, false);
            if (!result.success)
                throw std::runtime_error("Dislocation Nye tensor computation failed: " + result.message);

            const glm::dvec3 linePoint(dislo.linePoint);
            const glm::dvec3 lineDir(dislo.lineDirection);

            double outsideMax = 0.0;  // beyond the cylinder, plus a margin for the neighbor cutoff
            double insideMax = 0.0;   // inside the displaced cylinder
            for (size_t i = 0; i < reference.atoms.size(); ++i)
            {
                const glm::dvec3 pos(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);
                const glm::dvec3 rel = pos - linePoint;
                const double alongLine = glm::dot(rel, lineDir);
                const double radial = glm::length(rel - alongLine * lineDir);
                // A neighbor exactly at the cylinder edge is still displaced,
                // so "outside" needs a margin comfortably past the neighbor
                // cutoff (3 A) to guarantee every neighbor of a sampled atom
                // is also strictly outside the displaced region.
                if (radial > 15.0f + 6.0)
                    outsideMax = std::max(outsideMax, result.alphaNorm[i]);
                else if (radial < 15.0f - 3.0)  // inside cylinder, away from its edge
                    insideMax = std::max(insideMax, result.alphaNorm[i]);
            }

            std::cout << "Edge dislocation: outside-cylinder max |alpha| = " << outsideMax
                      << ", inside max |alpha| = " << insideMax << "\n";
            if (outsideMax > 1e-8)
                throw std::runtime_error("Atoms entirely outside the displaced region have a non-zero "
                    "Nye tensor (max = " + std::to_string(outsideMax) + "); expected exactly zero");
            if (!(insideMax > 0.01))
                throw std::runtime_error("No clearly non-zero Nye tensor was found inside the displaced "
                    "region: max = " + std::to_string(insideMax));
        }

        std::cout << "Nye tensor regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
