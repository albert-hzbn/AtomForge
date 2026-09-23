#pragma once

#include <glm/glm.hpp>

// Full anisotropic elastic stiffness tensor C_ijkl (GPa), 3x3x3x3, always
// fully symmetric (major + minor symmetries already applied by the
// constructors below). Values are expected in GPa; the dislocation solver
// only ever uses ratios of these constants so the unit is not otherwise
// load-bearing, but GPa matches the convention used by reference codes
// (e.g. BABEL, Clouet et al.) and by the cubic/hexagonal literature values
// a user is likely to type in directly.
struct ElasticTensor
{
    double c[3][3][3][3] = {};
};

// Cubic symmetry: 3 independent constants (C11, C12, C44), crystal axes
// assumed aligned with x,y,z.
ElasticTensor makeCubicElasticTensor(double c11, double c12, double c44);

// Hexagonal symmetry: 5 independent constants; the 6-fold axis (c-axis)
// must be along z (C66 = (C11-C12)/2 is derived automatically), matching
// the convention BABEL documents for hexagonal_elasticity.
ElasticTensor makeHexagonalElasticTensor(double c11, double c12, double c13, double c33, double c44);

// General (triclinic-or-lower) symmetry from the 6x6 Voigt stiffness matrix
// (GPa), crystal axes assumed aligned with x,y,z. voigt[a][b] for
// a,b in 0..5 using the standard Voigt index order xx,yy,zz,yz,xz,xy. Only
// the upper triangle needs to be filled in; symmetry is enforced.
ElasticTensor makeVoigtElasticTensor(const double voigt[6][6]);

// Rotate a tensor from the crystal frame into a new orthonormal frame whose
// x/y/z axes are given (as unit vectors, expressed in the crystal frame) by
// the columns of `axes` (axes[0]=new x, axes[1]=new y, axes[2]=new z).
ElasticTensor rotateElasticTensor(const ElasticTensor& c, const glm::dmat3& axes);

// True if the tensor is (numerically) elastically isotropic, i.e.
// C11-C12 == 2*C44 in every equivalent cubic-like orientation. The Stroh
// sextic formalism is singular at exact isotropy (repeated roots), so
// callers should perturb the input slightly before solving when this is
// true -- the same workaround BABEL's own documentation recommends
// (CVoigt_noise) for elastically degenerate/high-symmetry orientations.
bool isNearlyIsotropic(const ElasticTensor& c, double relativeTolerance = 1e-6);
