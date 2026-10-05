#pragma once

#include "model/Structure.h"

#include <map>
#include <string>

// Alloy tool: builds a special-quasirandom-structure-style atom arrangement
// on a fixed lattice by simulated annealing over atom-species swaps,
// minimizing Warren-Cowley short-range-order parameters (see
// ShortRangeOrderAnalysis.h, Cowley, Phys. Rev. 138 (1965) A1384) toward
// zero across the first `shells` neighbor shells -- the same target most
// practical SQS generation (e.g. ATAT's mcsqs, icet) uses when no
// DFT-derived cluster-expansion target is supplied. This is a lighter-weight
// local optimizer, not a full cluster-expansion correlation-function match.
namespace atomforge
{

struct SQSParams
{
    // Target composition as element symbol -> fraction (must sum to ~1).
    std::map<std::string, double> composition;

    int shells = 2;
    double shellTolerance = 0.2;

    int steps = 3000;
    double startTemperature = 1.0;
    double endTemperature = 0.02;
    unsigned seed = 1;
};

struct SQSResult
{
    bool success = false;
    std::string message;
    Structure structure;

    double initialObjective = 0.0;
    double finalObjective = 0.0;
    int acceptedSwaps = 0;
};

SQSResult buildSQS(const Structure& source, const SQSParams& params);

}
