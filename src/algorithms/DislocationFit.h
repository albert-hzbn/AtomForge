#pragma once

#include "algorithms/NyeTensor.h"

#include <glm/glm.hpp>

#include <string>
#include <vector>

// Determines a dislocation's position and Burgers vector from a measured
// Nye tensor field, the same purpose BABEL's (self-described as
// "not documented yet") displacementFit/vitekFit programs serve --
// extracting dislocation properties by fitting, rather than the forward
// direction (DislocationBuilder.h) of inserting a known dislocation. Unlike
// BABEL's nonlinear least-squares fit against a fixed elastic solution,
// this uses the Nye tensor's own defining property directly: it is a
// density whose weighted centroid IS the best-fit dislocation line
// position, and whose integral over the cross-section IS the Burgers
// vector (Nye, Acta Metall. 1 (1953) 153) -- a simpler, non-iterative
// estimator built on the same physical quantity BABEL's tools consume.
// Position and Burgers-direction recovery are validated precisely
// (tests/dislocation_fit_regressions.cpp) for an edge dislocation, where
// Burgers is perpendicular to the line; a pure screw dislocation (Burgers
// parallel to the line) has not been separately validated here and may be
// less reliable -- flagged rather than silently assumed to generalize.
namespace atomforge
{

struct DislocationFitResult
{
    bool success = false;
    std::string message;

    // Point on the fitted dislocation line, in the plane perpendicular to
    // lineDirection through the |alpha|-weighted centroid of the input
    // positions.
    glm::dvec3 linePosition{0.0};

    // Burgers vector estimate: the Nye tensor alpha integrated over the
    // supplied atoms (each weighted by an equal representative area/atom),
    // contracted with lineDirection -- Nye's theorem, sum_i alpha_i . dA_i
    // = b for a small loop threaded once by the dislocation line.
    glm::dvec3 burgersVector{0.0};

    double totalWeight = 0.0; // sum of |alpha| used in the fit, for diagnostics
};

// positions and the Nye tensor result must correspond atom-for-atom (i.e.
// pass the same `reference` structure and its NyeTensorResult here).
// lineDirection is the assumed dislocation line direction (fit only
// resolves position in the plane perpendicular to it, and the Burgers
// component along it); areaPerAtom is the cross-sectional area (Angstrom^2)
// represented by each atom, used to convert the per-atom Nye tensor into an
// integrated Burgers vector (e.g. the reference lattice's atomic volume
// divided by the neighbor-shell thickness along lineDirection -- pass 1.0
// to get a per-atom-count-normalized rather than physically-integrated
// Burgers estimate if that spacing isn't known).
DislocationFitResult fitDislocationFromNye(const std::vector<glm::dvec3>& positions,
                                          const NyeTensorResult& nye,
                                          const glm::dvec3& lineDirection,
                                          double areaPerAtom);

}
