#include "algorithms/SQSBuilder.h"
#include "algorithms/ShortRangeOrderAnalysis.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace atomforge
{

namespace
{
int resolveAtomicNumber(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    return 0;
}

void setSpecies(AtomSite& atom, int z)
{
    atom.atomicNumber = z;
    atom.symbol = elementSymbol(z);
    getDefaultElementColor(z, atom.r, atom.g, atom.b);
}

double objective(const Structure& structure, int shells, double tolerance)
{
    const SroReport report = analyzeWarrenCowleySRO(structure, shells, tolerance);
    double sum = 0.0;
    for (const auto& entry : report.entries) sum += entry.alpha * entry.alpha;
    return sum;
}
}

SQSResult buildSQS(const Structure& source, const SQSParams& params)
{
    SQSResult result;
    if (source.atoms.empty())
    {
        result.message = "Structure has no atoms.";
        return result;
    }
    if (params.composition.size() < 2)
    {
        result.message = "composition needs at least two elements.";
        return result;
    }

    double total = 0.0;
    std::vector<std::pair<int, double>> fractions;
    for (const auto& [symbol, fraction] : params.composition)
    {
        const int z = resolveAtomicNumber(symbol);
        if (z == 0)
        {
            result.message = "Unknown element symbol '" + symbol + "'.";
            return result;
        }
        if (fraction < 0.0)
        {
            result.message = "Composition fractions must be non-negative.";
            return result;
        }
        fractions.push_back({z, fraction});
        total += fraction;
    }
    if (std::abs(total - 1.0) > 1e-2)
    {
        result.message = "Composition fractions must sum to 1.0 (got " + std::to_string(total) + ").";
        return result;
    }
    if (params.steps < 0)
    {
        result.message = "steps must be non-negative.";
        return result;
    }

    const int n = (int)source.atoms.size();
    std::mt19937 rng(params.seed);

    // Assign target atom counts per species (largest-remainder rounding so
    // the counts sum exactly to n), then shuffle species labels onto the
    // fixed lattice positions -- geometry never changes after this point,
    // only which species occupies which site.
    std::vector<int> counts(fractions.size());
    int assigned = 0;
    std::vector<double> remainders(fractions.size());
    for (size_t i = 0; i < fractions.size(); ++i)
    {
        const double exact = fractions[i].second * n;
        counts[i] = (int)std::floor(exact);
        remainders[i] = exact - counts[i];
        assigned += counts[i];
    }
    std::vector<size_t> order(fractions.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return remainders[a] > remainders[b]; });
    for (size_t k = 0; k < order.size() && assigned < n; ++k, ++assigned)
        counts[order[k]] += 1;

    std::vector<int> species;
    species.reserve((size_t)n);
    for (size_t i = 0; i < fractions.size(); ++i)
        for (int c = 0; c < counts[i]; ++c) species.push_back(fractions[i].first);
    std::shuffle(species.begin(), species.end(), rng);

    Structure structure = source;
    for (int i = 0; i < n; ++i) setSpecies(structure.atoms[(size_t)i], species[(size_t)i]);

    result.initialObjective = objective(structure, params.shells, params.shellTolerance);
    double current = result.initialObjective;

    std::uniform_int_distribution<int> pick(0, n - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int step = 0; step < params.steps; ++step)
    {
        const int i = pick(rng);
        int j = pick(rng);
        int guard = 0;
        while ((j == i || structure.atoms[(size_t)j].atomicNumber == structure.atoms[(size_t)i].atomicNumber)
               && guard++ < 64)
            j = pick(rng);
        if (structure.atoms[(size_t)j].atomicNumber == structure.atoms[(size_t)i].atomicNumber) continue;

        const int zi = structure.atoms[(size_t)i].atomicNumber;
        const int zj = structure.atoms[(size_t)j].atomicNumber;
        setSpecies(structure.atoms[(size_t)i], zj);
        setSpecies(structure.atoms[(size_t)j], zi);

        const double candidate = objective(structure, params.shells, params.shellTolerance);
        const double delta = candidate - current;
        const double t = params.startTemperature
            + (params.endTemperature - params.startTemperature) * (params.steps > 1 ? (double)step / (double)(params.steps - 1) : 1.0);
        const bool accept = delta <= 0.0 || (t > 1e-9 && std::exp(-delta / t) > unit(rng));
        if (accept)
        {
            current = candidate;
            ++result.acceptedSwaps;
        }
        else
        {
            setSpecies(structure.atoms[(size_t)i], zi);
            setSpecies(structure.atoms[(size_t)j], zj);
        }
    }

    result.structure = structure;
    result.finalObjective = current;
    result.success = true;
    result.message = "SQS optimization: sum(alpha^2) over " + std::to_string(params.shells) + " shells went from "
        + std::to_string(result.initialObjective) + " to " + std::to_string(result.finalObjective) + " over "
        + std::to_string(params.steps) + " steps (" + std::to_string(result.acceptedSwaps) + " accepted swaps).";
    return result;
}

}
