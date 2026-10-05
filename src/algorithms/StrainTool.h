#pragma once

#include "model/Structure.h"

#include <string>

// Lattice tool: applies a homogeneous deformation to a periodic structure's
// cell, holding fractional atomic coordinates fixed (the standard way to
// impose a macroscopic elastic/engineering strain on a crystal -- the same
// convention VASP/LAMMPS deformation-gradient workflows use). Distinct from
// the integer-matrix supercell/basis change already available via Edit >
// Transform: this applies a general, not-necessarily-integer deformation.
namespace atomforge
{

struct StrainParams
{
    // Deformation gradient F applied to the cell: new lattice vector i is
    // F * (old lattice vector i). Row-major; identity means no change.
    double f[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
};

struct StrainResult
{
    bool success = false;
    std::string message;
    Structure structure;
};

// Builds F = I + strain from engineering strain components (exx, eyy, ezz on
// the diagonal; exy, exz, eyz as the symmetric off-diagonal shear terms).
StrainParams engineeringStrain(double exx, double eyy, double ezz,
                               double exy, double exz, double eyz);

StrainResult applyStrain(const Structure& source, const StrainParams& params);

}
