#include "algorithms/StrainTool.h"

#include <cmath>

namespace atomforge
{

StrainParams engineeringStrain(double exx, double eyy, double ezz,
                               double exy, double exz, double eyz)
{
    StrainParams p;
    p.f[0][0] = 1.0 + exx; p.f[0][1] = exy;       p.f[0][2] = exz;
    p.f[1][0] = exy;       p.f[1][1] = 1.0 + eyy; p.f[1][2] = eyz;
    p.f[2][0] = exz;       p.f[2][1] = eyz;       p.f[2][2] = 1.0 + ezz;
    return p;
}

StrainResult applyStrain(const Structure& source, const StrainParams& params)
{
    StrainResult result;
    if (!source.hasUnitCell)
    {
        result.message = "Structure has no unit cell to strain.";
        return result;
    }

    const auto& oldCell = source.cellVectors;
    double det = oldCell[0][0] * (oldCell[1][1] * oldCell[2][2] - oldCell[1][2] * oldCell[2][1])
               - oldCell[0][1] * (oldCell[1][0] * oldCell[2][2] - oldCell[1][2] * oldCell[2][0])
               + oldCell[0][2] * (oldCell[1][0] * oldCell[2][1] - oldCell[1][1] * oldCell[2][0]);
    if (std::abs(det) < 1e-12)
    {
        result.message = "Structure's cell is degenerate (zero volume).";
        return result;
    }

    // Inverse of the (row-vector) cell, so frac = cart * cellInv.
    double inv[3][3];
    const double invDet = 1.0 / det;
    inv[0][0] = (oldCell[1][1]*oldCell[2][2]-oldCell[1][2]*oldCell[2][1]) * invDet;
    inv[0][1] = (oldCell[0][2]*oldCell[2][1]-oldCell[0][1]*oldCell[2][2]) * invDet;
    inv[0][2] = (oldCell[0][1]*oldCell[1][2]-oldCell[0][2]*oldCell[1][1]) * invDet;
    inv[1][0] = (oldCell[1][2]*oldCell[2][0]-oldCell[1][0]*oldCell[2][2]) * invDet;
    inv[1][1] = (oldCell[0][0]*oldCell[2][2]-oldCell[0][2]*oldCell[2][0]) * invDet;
    inv[1][2] = (oldCell[0][2]*oldCell[1][0]-oldCell[0][0]*oldCell[1][2]) * invDet;
    inv[2][0] = (oldCell[1][0]*oldCell[2][1]-oldCell[1][1]*oldCell[2][0]) * invDet;
    inv[2][1] = (oldCell[0][1]*oldCell[2][0]-oldCell[0][0]*oldCell[2][1]) * invDet;
    inv[2][2] = (oldCell[0][0]*oldCell[1][1]-oldCell[0][1]*oldCell[1][0]) * invDet;

    Structure out = source;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double v = 0.0;
            for (int k = 0; k < 3; ++k) v += params.f[i][k] * oldCell[k][j];
            out.cellVectors[i][j] = v;
        }

    for (auto& atom : out.atoms)
    {
        const double cart[3] = {atom.x, atom.y, atom.z};
        double frac[3] = {0.0, 0.0, 0.0};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                frac[j] += cart[i] * inv[i][j];
        double newCart[3] = {0.0, 0.0, 0.0};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                newCart[j] += frac[i] * out.cellVectors[i][j];
        atom.x = newCart[0];
        atom.y = newCart[1];
        atom.z = newCart[2];
    }

    result.structure = out;
    result.success = true;
    result.message = "Applied deformation gradient to cell and " + std::to_string(out.atoms.size()) + " atoms.";
    return result;
}

}
