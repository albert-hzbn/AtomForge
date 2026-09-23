// Validates computeVitekMap (src/algorithms/VitekMap.h), a native
// reimplementation of BABEL's documented differential-displacement
// algorithm (Vitek, Phil. Mag. 18 (1968) 773; see BABEL's vitek.f90).
#include "algorithms/DislocationBuilder.h"
#include "algorithms/VitekMap.h"
#include "model/Structure.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge;

namespace
{
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
}

int main()
{
    try
    {
        const double a = 2.87;
        Structure reference = buildBcc(a, 12);

        // --- 1. Identical structures: every pair's differential
        // displacement must be exactly zero.
        {
            const VitekResult result = computeVitekMap(reference, reference, glm::dvec3(0, 0, 1), 0.0, 0.0f, false);
            if (!result.success)
                throw std::runtime_error("Zero-displacement Vitek map failed: " + result.message);
            if (result.pairs.empty())
                throw std::runtime_error("No neighbor pairs found for the auto-detected cutoff");
            double maxAbs = 0.0;
            for (const auto& p : result.pairs)
                maxAbs = std::max({maxAbs, std::abs(p.screwComponent), std::abs(p.edgeComponent.x), std::abs(p.edgeComponent.y)});
            std::cout << "Identical structures: " << result.pairs.size() << " pairs, max |dU| = " << maxAbs << "\n";
            if (maxAbs > 1e-9)
                throw std::runtime_error("Identical structures produced non-zero differential displacement: " + std::to_string(maxAbs));
        }

        // --- 2. Real 1/2<111> BCC screw dislocation: near-core pairs must
        // show a clearly larger differential displacement than far-field
        // pairs (the classic screw-core signature).
        {
            DislocationParams params;
            params.character = DislocationCharacter::Screw;
            params.shape = DislocationShape::Cylinder;
            params.cylinderRadius = 20.0f;
            params.coreRadius = 1.0f;
            params.cutoffRadius = 0.0f;
            params.autoDirections = false;
            params.lineUvw = glm::ivec3(1, 1, 1);
            params.burgersUvw = glm::ivec3(1, 1, 1);

            const DislocationResult dislo = buildDislocation(reference, params);
            if (!dislo.success)
                throw std::runtime_error("Failed to build reference screw dislocation: " + dislo.message);

            const glm::dvec3 lineDir(dislo.lineDirection);
            const double bz = dislo.burgersMagnitude;

            const VitekResult result = computeVitekMap(reference, dislo.output, lineDir, bz, 0.0f, false);
            if (!result.success)
                throw std::runtime_error("Dislocation Vitek map failed: " + result.message);

            // The branch cut of a Volterra dislocation's displacement field
            // extends as a ray from the core out to infinity (the same cut
            // this session's earlier isotropic-formula work already
            // established), so a handful of neighbor pairs crossing that
            // ray show a large differential displacement no matter how far
            // from the core they are -- a real, expected feature of the
            // construction, not noise. The screw field's gradient magnitude
            // scales as 1/r away from the cut, so the MEAN over a full
            // angular ring (dominated by the many near-zero, non-cut pairs,
            // unlike the median which can land exactly on a similarly-sized
            // cluster of near-threshold values in both rings) is what
            // actually captures the 1/r falloff here.
            const glm::dvec3 linePoint(dislo.linePoint);
            double nearCoreSum = 0.0, farFieldSum = 0.0;
            int nearCoreCount = 0, farFieldCount = 0;
            for (const auto& p : result.pairs)
            {
                const glm::dvec3 mid = 0.5 * (p.positionI + p.positionJ);
                const glm::dvec3 rel = mid - linePoint;
                const double radial = glm::length(rel - lineDir * glm::dot(rel, lineDir));
                const double screwFrac = std::abs(p.screwComponent) / bz;
                if (radial < 3.0)
                    { nearCoreSum += screwFrac; ++nearCoreCount; }
                else if (radial > 10.0 && radial < 15.0)
                    { farFieldSum += screwFrac; ++farFieldCount; }
            }
            const double nearCoreMean = nearCoreCount ? nearCoreSum / nearCoreCount : 0.0;
            const double farFieldMean = farFieldCount ? farFieldSum / farFieldCount : 0.0;

            std::cout << "BCC screw dislocation: near-core mean |dU_screw/b| = " << nearCoreMean
                      << " (n=" << nearCoreCount << "), far-field mean = " << farFieldMean
                      << " (n=" << farFieldCount << ", bz=" << bz << ")\n";
            if (!(nearCoreMean > 0.02))
                throw std::runtime_error("No significant differential displacement found near the core: " + std::to_string(nearCoreMean));
            if (!(farFieldMean < 0.5 * nearCoreMean))
                throw std::runtime_error("Differential displacement is not localized near the core (far field mean "
                    + std::to_string(farFieldMean) + " vs near-core mean " + std::to_string(nearCoreMean) + ")");
        }

        std::cout << "Vitek map regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
