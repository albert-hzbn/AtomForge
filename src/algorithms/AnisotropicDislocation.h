#pragma once

#include "algorithms/ElasticConstants.h"

#include <glm/glm.hpp>

#include <complex>
#include <string>

// Straight-dislocation displacement field for a fully anisotropic elastic
// medium, via the Stroh sextic formalism (Stroh, Phil. Mag. 3 (1958) 625;
// J. Math. Phys. 41 (1962) 77; Barnett & Lothe, Physica Norvegica 7 (1973)
// 13; see also Hirth & Lothe, "Theory of Dislocations", 2nd ed., ch. 13).
// This is the same physical model and reference literature used by the
// BABEL anisotropic-dislocation code (Douin, Clouet et al.), reimplemented
// natively here rather than ported, and validated against AtomForge's own
// isotropic closed-form solution in the near-isotropic limit (see
// tests/anisotropic_dislocation_regressions.cpp) since no BABEL binary was
// available in this environment for direct output comparison.
namespace atomforge::dislocation
{

struct StrohSolution
{
    bool valid = false;
    std::string error;

    // Three roots (Im > 0 by convention) of the sextic characteristic
    // equation, and their associated displacement (A) and traction (L)
    // eigenvectors, plus the Burgers-vector-fixing coefficients D.
    std::complex<double> p[3];
    std::complex<double> A[3][3]; // A[alpha][component]
    std::complex<double> L[3][3]; // L[alpha][component]
    std::complex<double> D[3];
};

// Solve the sextic eigenproblem for elastic constants already expressed in
// the dislocation's local frame (x = e1, y = e2, z = line direction), then
// fix the free coefficients D from the Burgers vector (also expressed in
// that same local frame) via the standard two conditions: the displacement
// jump on circling the line equals b, and the "stress function" jump is
// zero (single-valued stresses away from the cut).
//
// If the input is elastically isotropic (or too close to it), the sextic
// roots become degenerate and the solve is refused (matching BABEL's own
// documented limitation); pass a small relative perturbation via
// `noiseSeed`/`noiseAmplitude` to break the degeneracy the same way BABEL's
// CVoigt_noise option does, or use the isotropic closed form directly for
// exactly isotropic material.
StrohSolution solveStroh(const ElasticTensor& elasticLocal,
                         const glm::dvec3& burgersLocal,
                         unsigned noiseSeed = 0,
                         double noiseAmplitude = 0.0);

// Displacement u(x,y) in the same local frame the solution was built in
// (x = e1, y = e2); the z (line-direction) component is always zero for a
// pure Volterra straight dislocation.
glm::dvec3 anisotropicDisplacement(const StrohSolution& solution, double x, double y);

}
