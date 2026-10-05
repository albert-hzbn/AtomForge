#pragma once

#include "model/Structure.h"

#include <string>

// Surface tool: cleaves a single bulk crystal along a Miller plane (hkl),
// producing a 2D-periodic slab of the requested layer count with vacuum
// padding along the surface normal -- the counterpart to the existing
// two-grain tools (CSL Grain Boundary, Interface Builder), which both
// already build vacuum-padded, plane-oriented cells but only for a pair of
// grains. This finds the smallest integer lattice basis whose first two
// vectors lie in the (hkl) plane family via a bounded search (same bounded
// "nmax" convention InterfaceBuilder.h's generateUniqueSupercells uses),
// then reuses the CSL module's supercell/orthogonalization primitives
// (CSLComputation.h) to build and vacuum-pad the resulting cell.
namespace atomforge
{

struct SurfaceParams
{
    int h = 1, k = 1, l = 1;
    int layers = 4;
    double vacuum = 15.0;

    // Bound on the integer search for the in-plane and out-of-plane basis
    // vectors; raise it for high-index or long-period planes that need
    // larger coefficients.
    int nmax = 8;

    // Reduce the input to its primitive cell before searching (fewer atoms
    // per layer, matching the CSL Grain Boundary dialog's own default).
    bool primitiveInput = true;
    double primitiveSymprec = 1e-3;
};

struct SurfaceResult
{
    bool success = false;
    std::string message;
    Structure structure;
    int atomsPerLayer = 0;
};

SurfaceResult buildSurface(const Structure& source, const SurfaceParams& params);

}
