#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <string>
#include <vector>

// Prepares input for a constrained-minimization ("drag") calculation of a
// migration barrier (e.g. the Peierls barrier of a dislocation moving
// between two neighboring equilibrium positions): linearly interpolate
// between an initial and a final atomic configuration at a reaction
// coordinate zeta in [0,1], and report the reaction-coordinate direction
// per atom (the constraint an external relaxation/MD code would fix the
// projection onto, while letting every perpendicular degree of freedom
// relax). AtomForge does not itself perform the constrained minimization --
// like BABEL's prepareDrag program (mainPrepareDrag.f90), this only
// prepares the interpolated structure and constraint direction for an
// external code to consume.
namespace atomforge
{

struct DragPrepResult
{
    bool success = false;
    std::string message;

    Structure interpolated;                        // initial + zeta * (final - initial)
    std::vector<glm::dvec3> constraintDirection;    // final[i] - initial[i], one per atom
};

// initial and final must have the same atom count and ordering. If
// clipDisplacement, each atom's final-minus-initial vector is wrapped to
// its minimum-image representation under `initial`'s periodicity before
// interpolating (needed when an atom's migration path was stored as having
// crossed a periodic boundary the "long way round").
DragPrepResult prepareDrag(const Structure& initial,
                          const Structure& final_,
                          double zeta,
                          bool clipDisplacement,
                          bool usePbc);

}
