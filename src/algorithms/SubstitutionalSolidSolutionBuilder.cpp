#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <numeric>
#include <random>
#include <sstream>

// ---------------------------------------------------------------------------
// buildSubstitutionalSolidSolution
// ---------------------------------------------------------------------------
SSSResult buildSubstitutionalSolidSolution(const Structure& base,
                                            const SSSParams& params)
{
    SSSResult result;

    if (base.atoms.empty())
    {
        result.message = "Base structure has no atoms.";
        return result;
    }

    if (params.composition.empty())
    {
        result.message = "No composition entries specified.";
        return result;
    }

    // Validate atomic numbers.
    for (const SSSElementFraction& ef : params.composition)
    {
        if (ef.atomicNumber < 1 || ef.atomicNumber > 118)
        {
            result.message = "Invalid atomic number in composition.";
            return result;
        }
        if (ef.fraction < 0.0f)
        {
            result.message = "Negative fraction in composition.";
            return result;
        }
    }

    // Normalise fractions so they sum to exactly 1.
    float totalFrac = 0.0f;
    for (const SSSElementFraction& ef : params.composition)
        totalFrac += ef.fraction;

    if (totalFrac < 1e-6f)
    {
        result.message = "All composition fractions are zero.";
        return result;
    }

    const int N = static_cast<int>(base.atoms.size());

    // Compute per-element target counts using the "largest remainder" method
    // so that the counts always sum exactly to N.
    const int numEl = static_cast<int>(params.composition.size());
    std::vector<int> counts(numEl, 0);
    std::vector<float> remainders(numEl, 0.0f);

    int assigned = 0;
    for (int i = 0; i < numEl; ++i)
    {
        const float exact = (params.composition[i].fraction / totalFrac) * (float)N;
        counts[i]     = static_cast<int>(std::floor(exact));
        remainders[i] = exact - (float)counts[i];
        assigned     += counts[i];
    }

    // Distribute leftover seats by largest remainder.
    int leftover = N - assigned;
    if (leftover > 0)
    {
        std::vector<int> order(numEl);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b)
        {
            return remainders[a] > remainders[b];
        });
        for (int k = 0; k < leftover; ++k)
            counts[order[k % numEl]]++;
    }

    // Build a shuffled index list so substitutions are random.
    const unsigned int rngSeed = (params.seed != 0)
                                ? params.seed
                                : static_cast<unsigned int>(std::time(nullptr));
    std::mt19937 rng(rngSeed);

    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    // Assign element types.
    std::vector<int> atomicNumbers(N, 0);
    int pos = 0;
    for (int i = 0; i < numEl; ++i)
    {
        for (int c = 0; c < counts[i]; ++c)
            atomicNumbers[indices[pos++]] = params.composition[i].atomicNumber;
    }

    // Build output structure – copy geometry from base, overwrite element info.
    result.output = base;
    result.output.grainColors.clear();
    result.output.grainRegionIds.clear();

    // Pre-compute element colours from ElementData defaults.
    const auto defaultColors = makeDefaultElementColors();

    for (int i = 0; i < N; ++i)
    {
        AtomSite& atom = result.output.atoms[i];
        const int z    = atomicNumbers[i];

        atom.atomicNumber = z;
        atom.symbol       = elementSymbol(z);

        if (z >= 1 && z < static_cast<int>(defaultColors.size()))
        {
            atom.r = defaultColors[z].r;
            atom.g = defaultColors[z].g;
            atom.b = defaultColors[z].b;
        }
    }

    // Fill result metadata.
    for (int i = 0; i < numEl; ++i)
        result.elementCounts.emplace_back(params.composition[i].atomicNumber, counts[i]);

    std::ostringstream msg;
    msg << "Built substitutional solid solution: " << N << " atoms";
    for (int i = 0; i < numEl; ++i)
    {
        msg << ",  " << elementSymbol(params.composition[i].atomicNumber)
            << " " << counts[i];
    }
    result.message = msg.str();
    result.success = true;
    return result;
}
