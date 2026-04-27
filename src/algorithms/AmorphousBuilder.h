#pragma once

#include "io/StructureLoader.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

// How the simulation box dimensions are determined.
enum class AmorphousBoxMode
{
    Manual = 0,       // User specifies a, b, c directly in Angstroms
    AutoFromDensity,  // Box volume computed from total mass and target density
};

// One element species in the amorphous mixture.
struct AmorphousElementSpec
{
    int atomicNumber = 14;  // Si by default
    int count        = 100;
};

// Explicit minimum centre-to-centre separation for one element pair (z1 <= z2).
// When provided, overrides the covalent-radii default for that pair.
struct AmorphousPairDist
{
    int   z1      = 0;
    int   z2      = 0;
    float minDist = 0.0f;  // Angstroms
};

struct AmorphousParams
{
    // -----------------------------------------------------------------------
    // Composition
    // -----------------------------------------------------------------------
    std::vector<AmorphousElementSpec> elements;

    // -----------------------------------------------------------------------
    // Box geometry
    // -----------------------------------------------------------------------
    AmorphousBoxMode boxMode = AmorphousBoxMode::AutoFromDensity;

    // Manual box (used when boxMode == Manual).
    float boxA = 20.0f;
    float boxB = 20.0f;
    float boxC = 20.0f;

    // Target density in g/cm³ (used when boxMode == AutoFromDensity).
    // A cubic box is derived whose volume gives exactly this density.
    float targetDensity = 2.0f;

    // Multiplicative scale applied to every box vector before packing begins.
    // Values > 1 reduce effective packing density and make placement easier
    // (analogous to the 1.1× cell-expansion trick used in melt-quench MD).
    float cellScaleFactor = 1.0f;

    // -----------------------------------------------------------------------
    // Minimum-distance constraints
    // -----------------------------------------------------------------------
    // Per-pair overrides; any (z1, z2) pair not listed uses the covalent default.
    std::vector<AmorphousPairDist> pairDistances;

    // Fraction of (r_cov1 + r_cov2) used as the default minimum separation
    // when no explicit override is supplied.
    float covalentTolerance = 0.75f;

    // -----------------------------------------------------------------------
    // Placement parameters
    // -----------------------------------------------------------------------
    // RNG seed (0 = time-based).
    unsigned int seed = 42;

    // Maximum random placement attempts per atom before skipping it.
    int maxAttempts = 1000;

    // -----------------------------------------------------------------------
    // Output options
    // -----------------------------------------------------------------------
    // Attach orthogonal PBC unit cell vectors to the output structure.
    bool periodic = true;
};

struct AmorphousResult
{
    bool        success       = false;
    std::string message;
    Structure   output;

    int   requestedAtoms = 0;
    int   placedAtoms    = 0;
    float actualDensity  = 0.0f;  // g/cm³ of the placed atoms in the box
    float boxA           = 0.0f;
    float boxB           = 0.0f;
    float boxC           = 0.0f;
};

// Build an amorphous structure using Random Sequential Addition (RSA).
//
// Atoms are placed one by one at uniformly-random positions inside an
// orthogonal periodic box.  A candidate position is accepted if and only if
// every already-placed atom satisfies the element-pair minimum-distance
// constraint.  A spatial hash grid makes each acceptance check O(1).
//
// covalentRadii  – size-119 vector indexed by atomic number (from makeLiteratureCovalentRadii).
// elementColors  – size-119 vector indexed by atomic number (from makeDefaultElementColors).
[[nodiscard]] AmorphousResult buildAmorphousStructure(
    const AmorphousParams&        params,
    const std::vector<float>&     covalentRadii,
    const std::vector<glm::vec3>& elementColors);
