#include "algorithms/AmorphousBuilder.h"

#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <sstream>
#include <unordered_map>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace
{

// Encode three grid-cell indices into a single 64-bit hash key.
uint64_t encodeCell(int ix, int iy, int iz)
{
    // Large primes to mix the three signed integer coordinates.
    return  (uint64_t)(unsigned int)ix * 100003ULL
          ^ (uint64_t)(unsigned int)iy * 200003ULL
          ^ (uint64_t)(unsigned int)iz * 300007ULL;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// buildAmorphousStructure
// ---------------------------------------------------------------------------

AmorphousResult buildAmorphousStructure(
    const AmorphousParams&        params,
    const std::vector<float>&     covalentRadii,
    const std::vector<glm::vec3>& elementColors)
{
    AmorphousResult result;

    // ------------------------------------------------------------------
    // 1. Validate composition
    // ------------------------------------------------------------------
    if (params.elements.empty())
    {
        result.message = "No element species specified.";
        return result;
    }

    int totalAtoms = 0;
    for (const auto& e : params.elements)
    {
        if (e.count < 0)
        {
            result.message = "Negative atom count for element Z=" + std::to_string(e.atomicNumber);
            return result;
        }
        totalAtoms += e.count;
    }
    if (totalAtoms == 0)
    {
        result.message = "Total atom count is zero.";
        return result;
    }
    result.requestedAtoms = totalAtoms;

    // ------------------------------------------------------------------
    // 2. Compute box dimensions
    // ------------------------------------------------------------------
    float boxA, boxB, boxC;

    if (params.boxMode == AmorphousBoxMode::Manual)
    {
        boxA = params.boxA;
        boxB = params.boxB;
        boxC = params.boxC;
    }
    else // AutoFromDensity
    {
        // Total mass in atomic mass units (g/mol = u).
        double totalMass = 0.0;
        for (const auto& e : params.elements)
        {
            if (e.atomicNumber >= 1 && e.atomicNumber <= 118)
                totalMass += (double)e.count * elementAtomicMass(e.atomicNumber);
        }

        // density [g/cm³] = mass[u] * 1.66054e-24 [g/u] / vol[Å³] * 1e24 [Å³/cm³]
        // => volume [Å³]  = mass[u] * 1.66054 / density [g/cm³]
        const double targetDensity = std::max(1e-3, (double)params.targetDensity);
        const double volume        = totalMass * 1.66054 / targetDensity;
        const float  side          = (float)std::cbrt(volume);
        boxA = side;
        boxB = side;
        boxC = side;
    }

    // Apply user-specified scale factor.
    boxA *= params.cellScaleFactor;
    boxB *= params.cellScaleFactor;
    boxC *= params.cellScaleFactor;

    if (boxA < 0.1f || boxB < 0.1f || boxC < 0.1f)
    {
        result.message = "Computed box dimensions are too small (< 0.1 Å).";
        return result;
    }

    result.boxA = boxA;
    result.boxB = boxB;
    result.boxC = boxC;

    // ------------------------------------------------------------------
    // 3. Build element-pair minimum-distance lookup
    // ------------------------------------------------------------------
    std::unordered_map<int, float> pairDistMap;
    pairDistMap.reserve(params.pairDistances.size() * 2);
    for (const auto& pd : params.pairDistances)
    {
        const int za = std::min(pd.z1, pd.z2);
        const int zb = std::max(pd.z1, pd.z2);
        pairDistMap[za * 200 + zb] = pd.minDist;
    }

    auto getMinDist = [&](int z1, int z2) -> float
    {
        const int za = std::min(z1, z2);
        const int zb = std::max(z1, z2);
        const auto it = pairDistMap.find(za * 200 + zb);
        if (it != pairDistMap.end())
            return it->second;
        // Default: (r_cov1 + r_cov2) * tolerance
        const float r1 = (z1 >= 1 && z1 < (int)covalentRadii.size()) ? covalentRadii[z1] : 1.5f;
        const float r2 = (z2 >= 1 && z2 < (int)covalentRadii.size()) ? covalentRadii[z2] : 1.5f;
        return (r1 + r2) * params.covalentTolerance;
    };

    // Global maximum min-distance across all requested element pairs –
    // used to size the spatial hash grid cells.
    float maxMinDist = 0.01f;
    for (const auto& e1 : params.elements)
        for (const auto& e2 : params.elements)
            maxMinDist = std::max(maxMinDist, getMinDist(e1.atomicNumber, e2.atomicNumber));

    // ------------------------------------------------------------------
    // 4. Build and shuffle the flat atom list
    // ------------------------------------------------------------------
    std::vector<int> atomList; // atomic numbers, one entry per atom to place
    atomList.reserve((size_t)totalAtoms);
    for (const auto& e : params.elements)
        for (int k = 0; k < e.count; ++k)
            atomList.push_back(e.atomicNumber);

    unsigned int seed = params.seed;
    if (seed == 0)
        seed = (unsigned int)std::time(nullptr);
    std::mt19937 rng(seed);
    std::shuffle(atomList.begin(), atomList.end(), rng);

    // ------------------------------------------------------------------
    // 5. Initialise spatial grid
    // ------------------------------------------------------------------
    // Grid cell size = maxMinDist ensures that when checking the 3×3×3
    // neighbourhood of a candidate cell we cannot miss any conflict.
    const float cellSize = maxMinDist;
    const int numA = std::max(1, (int)std::floor(boxA / cellSize));
    const int numB = std::max(1, (int)std::floor(boxB / cellSize));
    const int numC = std::max(1, (int)std::floor(boxC / cellSize));

    // grid[encodeCell(ix,iy,iz)] = list of placed-atom indices into `positions`.
    std::unordered_map<uint64_t, std::vector<int>> grid;
    grid.reserve((size_t)totalAtoms);

    std::vector<glm::vec3> positions;
    std::vector<int>       placedZ;
    positions.reserve((size_t)totalAtoms);
    placedZ.reserve((size_t)totalAtoms);

    std::uniform_real_distribution<float> randA(0.0f, boxA);
    std::uniform_real_distribution<float> randB(0.0f, boxB);
    std::uniform_real_distribution<float> randC(0.0f, boxC);

    // ------------------------------------------------------------------
    // 6. Random Sequential Addition loop
    // ------------------------------------------------------------------
    int placed  = 0;
    int skipped = 0;

    for (int ai = 0; ai < totalAtoms; ++ai)
    {
        const int z   = atomList[(size_t)ai];
        bool success  = false;

        for (int attempt = 0; attempt < params.maxAttempts; ++attempt)
        {
            const float cx = randA(rng);
            const float cy = randB(rng);
            const float cz = randC(rng);

            // Grid cell of candidate position.
            const int gx = ((int)std::floor(cx / cellSize) % numA + numA) % numA;
            const int gy = ((int)std::floor(cy / cellSize) % numB + numB) % numB;
            const int gz = ((int)std::floor(cz / cellSize) % numC + numC) % numC;

            bool conflict = false;

            // Check 3×3×3 neighbourhood with PBC wrapping.
            for (int dx = -1; dx <= 1 && !conflict; ++dx)
            for (int dy = -1; dy <= 1 && !conflict; ++dy)
            for (int dz = -1; dz <= 1 && !conflict; ++dz)
            {
                const int nx = ((gx + dx) % numA + numA) % numA;
                const int ny = ((gy + dy) % numB + numB) % numB;
                const int nz = ((gz + dz) % numC + numC) % numC;

                const auto it = grid.find(encodeCell(nx, ny, nz));
                if (it == grid.end()) continue;

                for (const int idx : it->second)
                {
                    // Minimum-image-convention distance.
                    float ddx = cx - positions[(size_t)idx].x;
                    float ddy = cy - positions[(size_t)idx].y;
                    float ddz = cz - positions[(size_t)idx].z;

                    ddx -= std::round(ddx / boxA) * boxA;
                    ddy -= std::round(ddy / boxB) * boxB;
                    ddz -= std::round(ddz / boxC) * boxC;

                    const float dist2 = ddx*ddx + ddy*ddy + ddz*ddz;
                    const float minD  = getMinDist(z, placedZ[(size_t)idx]);
                    if (dist2 < minD * minD)
                    {
                        conflict = true;
                        break;
                    }
                }
            }

            if (!conflict)
            {
                const int atomIdx = (int)positions.size();
                positions.push_back({ cx, cy, cz });
                placedZ.push_back(z);

                const int igx = ((int)std::floor(cx / cellSize) % numA + numA) % numA;
                const int igy = ((int)std::floor(cy / cellSize) % numB + numB) % numB;
                const int igz = ((int)std::floor(cz / cellSize) % numC + numC) % numC;
                grid[encodeCell(igx, igy, igz)].push_back(atomIdx);

                placed++;
                success = true;
                break;
            }
        }

        if (!success)
            skipped++;
    }

    // ------------------------------------------------------------------
    // 7. Build output Structure
    // ------------------------------------------------------------------
    Structure output;
    for (int ai = 0; ai < (int)positions.size(); ++ai)
    {
        AtomSite site;
        site.atomicNumber = placedZ[(size_t)ai];
        site.symbol       = elementSymbol(placedZ[(size_t)ai]);
        site.x            = (double)positions[(size_t)ai].x;
        site.y            = (double)positions[(size_t)ai].y;
        site.z            = (double)positions[(size_t)ai].z;

        const int z = placedZ[(size_t)ai];
        if (z >= 0 && z < (int)elementColors.size())
        {
            site.r = elementColors[(size_t)z].r;
            site.g = elementColors[(size_t)z].g;
            site.b = elementColors[(size_t)z].b;
        }

        output.atoms.push_back(std::move(site));
    }

    if (params.periodic)
    {
        output.hasUnitCell      = true;
        output.cellVectors[0][0] = (double)boxA; output.cellVectors[0][1] = 0.0; output.cellVectors[0][2] = 0.0;
        output.cellVectors[1][0] = 0.0; output.cellVectors[1][1] = (double)boxB; output.cellVectors[1][2] = 0.0;
        output.cellVectors[2][0] = 0.0; output.cellVectors[2][1] = 0.0; output.cellVectors[2][2] = (double)boxC;
        output.cellOffset[0]     = 0.0;
        output.cellOffset[1]     = 0.0;
        output.cellOffset[2]     = 0.0;
    }

    // ------------------------------------------------------------------
    // 8. Compute actual density and assemble result
    // ------------------------------------------------------------------
    double totalMassPlaced = 0.0;
    for (const int z : placedZ)
        if (z >= 1 && z <= 118)
            totalMassPlaced += elementAtomicMass(z);

    const double volumeAng3  = (double)boxA * (double)boxB * (double)boxC;
    result.actualDensity     = (float)(totalMassPlaced * 1.66054 / volumeAng3);
    result.output            = std::move(output);
    result.placedAtoms       = placed;

    if (skipped == 0)
    {
        result.success = true;
        result.message = "Successfully placed " + std::to_string(placed) + " atoms.";
    }
    else
    {
        result.success = placed > 0;
        std::ostringstream msg;
        msg << "Placed " << placed << " / " << totalAtoms << " atoms. "
            << skipped << " atom(s) could not be placed — try increasing "
            << "box size, scale factor, or max attempts.";
        result.message = msg.str();
    }

    return result;
}
