#pragma once

#include "model/Structure.h"

#include <array>
#include <string>
#include <vector>

// Nye (dislocation density) tensor: measures the local lattice curvature
// needed to describe a crystal containing dislocations without reference to
// individual defect lines (Nye, Acta Metall. 1 (1953) 153). Follows the
// same finite-difference-on-neighbors algorithm and Hartley & Mishin
// (Acta Mater. 53 (2005) 1313) lattice-correspondence convention BABEL's
// own nyeTensor.f90/gradDisplacement.f90 implement (including the documented
// "David Olmsted" index-order correction to the original H&M paper), so
// this is a native reimplementation of that specific published method
// rather than a new algorithm.
namespace atomforge
{

struct NyeTensorResult
{
    bool success = false;
    std::string message;

    // alpha[n] is the 3x3 Nye tensor at atom n (reference-structure indexing);
    // alpha[n][i][j] has units of 1/length (dislocation density x Burgers
    // vector component), Angstrom^-1 when input coordinates are Angstrom.
    std::vector<std::array<std::array<double, 3>, 3>> alpha;

    // Frobenius norm of alpha[n], a convenient scalar field for
    // visualization/thresholding ("how much dislocation content is here").
    std::vector<double> alphaNorm;
};

// reference and deformed must have the same atom count and ordering (atom n
// in `deformed` is the same physical atom as atom n in `reference`, e.g. a
// perfect crystal and the same crystal after buildDislocation()). Only
// atoms with at least 4 non-coplanar neighbors within cutoffRadius (in the
// reference structure) get a defined Nye tensor; others are left at zero
// with a note in `message`. usePbc requires `reference.hasUnitCell`.
NyeTensorResult computeNyeTensor(const Structure& reference,
                                 const Structure& deformed,
                                 float cutoffRadius,
                                 bool usePbc);

}
