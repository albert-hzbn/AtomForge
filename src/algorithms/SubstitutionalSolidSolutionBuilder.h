#pragma once

#include "io/StructureLoader.h"

#include <string>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// SubstitutionalSolidSolutionBuilder
//
// Randomly replaces atom species in a host structure to match a user-supplied
// composition.  The total number of atoms is preserved; only element types are
// changed.
// ---------------------------------------------------------------------------

struct SSSElementFraction
{
    int   atomicNumber = 0;
    float fraction     = 0.0f; // target fraction in [0, 1]; fractions must sum to 1
};

struct SSSParams
{
    // Desired composition.  The fractions should sum to 1 (clamped internally).
    std::vector<SSSElementFraction> composition;

    // RNG seed – set to 0 to use a time-based seed.
    unsigned int seed = 12345;
};

struct SSSResult
{
    bool        success      = false;
    std::string message;
    Structure   output;

    // Per-element atom counts in the output structure, in the same order as
    // SSSParams::composition.
    std::vector<std::pair<int, int>> elementCounts; // (atomicNumber, count)
};

// Build a substitutional solid solution from `base`.
// The function randomises which lattice sites carry which element species
// according to the target fractions supplied in `params`.
SSSResult buildSubstitutionalSolidSolution(const Structure& base,
                                            const SSSParams& params);
