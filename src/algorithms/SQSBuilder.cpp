#include "algorithms/SQSBuilder.h"
#include "algorithms/ShortRangeOrderAnalysis.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>
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

// Sum of squared Warren-Cowley parameters over the neighbour pairs, as
// analyzeWarrenCowleySRO reports them; `kind` is each site's species index.
double objective(const std::vector<SroPair>& pairs, int shells, const std::vector<int>& kind, const std::vector<double>& concentration)
{
    const std::size_t m = concentration.size();
    std::vector<double> n(static_cast<std::size_t>(shells) * m * m, 0.0);
    for (const SroPair& p : pairs)
        n[(static_cast<std::size_t>(p.shell) * m + static_cast<std::size_t>(kind[static_cast<std::size_t>(p.i)])) * m
          + static_cast<std::size_t>(kind[static_cast<std::size_t>(p.j)])] += 1.0;
    double sum = 0.0;
    for (std::size_t s = 0; s < static_cast<std::size_t>(shells); ++s)
        for (std::size_t a = 0; a < m; ++a) {
            double z = 0.0;
            for (std::size_t b = 0; b < m; ++b) z += n[(s * m + a) * m + b];
            if (z == 0.0) continue;
            for (std::size_t b = 0; b < m; ++b) {
                const double alpha = concentration[b] > 1e-12 ? 1.0 - n[(s * m + a) * m + b] / z / concentration[b] : 0.0;
                sum += alpha * alpha;
            }
        }
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

    // The geometry never changes: find the neighbour pairs once, then each
    // swap only recounts species over them.
    const std::vector<double> shellDistances = buildSroShellDistances(structure, params.shells, params.shellTolerance);
    const std::vector<SroPair> pairs = buildSroPairs(structure, shellDistances, params.shellTolerance);
    const int shells = static_cast<int>(shellDistances.size());
    std::vector<double> concentration(fractions.size());
    for (size_t k = 0; k < fractions.size(); ++k) concentration[k] = static_cast<double>(counts[k]) / n;
    std::vector<int> kind((size_t)n);
    for (int i = 0; i < n; ++i)
        for (size_t k = 0; k < fractions.size(); ++k)
            if (fractions[k].first == species[(size_t)i]) kind[(size_t)i] = static_cast<int>(k);

    result.initialObjective = objective(pairs, shells, kind, concentration);
    double current = result.initialObjective;
    // Uphill moves at a finite temperature can end the walk above its best
    // point: keep and return the best arrangement seen.
    Structure best = structure;
    double bestObjective = current;

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
        std::swap(kind[(size_t)i], kind[(size_t)j]);

        const double candidate = objective(pairs, shells, kind, concentration);
        const double delta = candidate - current;
        const double t = params.startTemperature
            + (params.endTemperature - params.startTemperature) * (params.steps > 1 ? (double)step / (double)(params.steps - 1) : 1.0);
        const bool accept = delta <= 0.0 || (t > 1e-9 && std::exp(-delta / t) > unit(rng));
        if (accept)
        {
            current = candidate;
            ++result.acceptedSwaps;
            if (current < bestObjective)
            {
                bestObjective = current;
                best = structure;
            }
        }
        else
        {
            setSpecies(structure.atoms[(size_t)i], zi);
            setSpecies(structure.atoms[(size_t)j], zj);
            std::swap(kind[(size_t)i], kind[(size_t)j]);
        }
    }

    result.structure = std::move(best);
    result.finalObjective = bestObjective;
    result.success = true;
    result.message = "SQS optimization: sum(alpha^2) over " + std::to_string(params.shells) + " shells went from "
        + std::to_string(result.initialObjective) + " to " + std::to_string(result.finalObjective) + " over "
        + std::to_string(params.steps) + " steps (" + std::to_string(result.acceptedSwaps) + " accepted swaps).";
    return result;
}

}
