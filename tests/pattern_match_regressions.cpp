// Validates PatternMatch (src/algorithms/PatternMatch.h), a native
// reimplementation of BABEL's documented pattern-based defect detection
// (see patternInit.f90/patternDetect.f90), generalized to auto-extract the
// reference pattern from any perfect structure instead of BABEL's hardcoded
// bcc/fcc/hcp direction tables.
#include "algorithms/PatternMatch.h"
#include "model/Structure.h"

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
        Structure reference = buildFcc(a, 8);
        const float cutoff = 3.2f; // between 1st (a/sqrt2=2.55) and 2nd (a=3.61) shells

        const CrystalPattern pattern = buildPattern(reference, cutoff, true);
        if (!pattern.success)
            throw std::runtime_error("buildPattern failed: " + pattern.message);
        if (pattern.directions.size() != 12)
            throw std::runtime_error("Expected 12 first-neighbor directions for FCC, got "
                + std::to_string(pattern.directions.size()));

        // --- 1. A perfect (undisturbed) crystal must match everywhere,
        // with zero angular deviation (it IS the pattern, just translated).
        {
            const PatternDetectResult result = detectPattern(reference, pattern, cutoff, 5.0, true);
            if (!result.success)
                throw std::runtime_error("detectPattern on perfect crystal failed: " + result.message);
            for (const auto& row : result.rows)
            {
                if (!row.matched)
                    throw std::runtime_error("Perfect crystal atom " + std::to_string(row.index) + " did not match the pattern");
                if (row.maxAngleDeviationDeg > 1e-3)
                    throw std::runtime_error("Perfect crystal atom " + std::to_string(row.index)
                        + " matched with non-zero deviation: " + std::to_string(row.maxAngleDeviationDeg));
            }
            std::cout << "Perfect crystal: all " << result.rows.size() << " atoms matched exactly\n";
        }

        // --- 2. A single displaced atom must break the pattern match for
        // itself AND for its former neighbors (their neighbor set changed),
        // while atoms far from the disturbance still match.
        {
            Structure disturbed = reference;
            const glm::dvec3 disturbedAtomPos(disturbed.atoms[0].x, disturbed.atoms[0].y, disturbed.atoms[0].z);
            disturbed.atoms[0].x += 5.0;

            const PatternDetectResult result = detectPattern(disturbed, pattern, cutoff, 5.0, true);
            if (!result.success)
                throw std::runtime_error("detectPattern on disturbed crystal failed: " + result.message);

            if (result.rows[0].matched)
                throw std::runtime_error("Displaced atom itself should not match the pattern");

            int unmatchedNearby = 0, matchedFar = 0;
            for (size_t i = 1; i < disturbed.atoms.size(); ++i)
            {
                const glm::dvec3 pos(disturbed.atoms[i].x, disturbed.atoms[i].y, disturbed.atoms[i].z);
                const double distToOrigin = glm::length(pos - disturbedAtomPos);
                if (distToOrigin < cutoff && !result.rows[i].matched)
                    ++unmatchedNearby;
                if (distToOrigin > 3.0 * cutoff && result.rows[i].matched)
                    ++matchedFar;
            }
            std::cout << "Disturbed crystal: " << unmatchedNearby << " former neighbors unmatched, "
                      << matchedFar << " far atoms still matched\n";
            if (unmatchedNearby == 0)
                throw std::runtime_error("No former neighbors of the displaced atom were flagged as unmatched");
            if (matchedFar == 0)
                throw std::runtime_error("No far-away atoms matched; detection is not localized");
        }

        std::cout << "Pattern match regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
