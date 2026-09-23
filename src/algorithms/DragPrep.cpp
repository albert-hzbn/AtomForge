#include "algorithms/DragPrep.h"

#include <Eigen/Dense>

#include <cmath>

namespace atomforge
{

namespace
{
using Mat3 = Eigen::Matrix3d;
using Vec3 = Eigen::Vector3d;

Vec3 toEigen(const glm::dvec3& v) { return Vec3(v.x, v.y, v.z); }
glm::dvec3 toGlm(const Vec3& v) { return glm::dvec3(v.x(), v.y(), v.z()); }

bool buildCellMatrices(const Structure& s, Mat3& cell, Mat3& invCell)
{
    if (!s.hasUnitCell)
        return false;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            cell(j, i) = s.cellVectors[i][j];
    if (std::abs(cell.determinant()) < 1e-12)
        return false;
    invCell = cell.inverse();
    return true;
}
}

DragPrepResult prepareDrag(const Structure& initial,
                          const Structure& final_,
                          double zeta,
                          bool clipDisplacement,
                          bool usePbc)
{
    DragPrepResult result;
    const size_t n = initial.atoms.size();
    if (n != final_.atoms.size())
    {
        result.message = "Initial and final structures must have the same atom count and ordering.";
        return result;
    }
    if (n == 0)
    {
        result.message = "Structures have no atoms.";
        return result;
    }
    if (zeta < 0.0 || zeta > 1.0)
    {
        result.message = "zeta must be in [0, 1].";
        return result;
    }

    Mat3 cell = Mat3::Identity(), invCell = Mat3::Identity();
    const bool havePbc = clipDisplacement && usePbc && buildCellMatrices(initial, cell, invCell);

    result.interpolated = initial;
    result.constraintDirection.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        glm::dvec3 delta(final_.atoms[i].x - initial.atoms[i].x,
                         final_.atoms[i].y - initial.atoms[i].y,
                         final_.atoms[i].z - initial.atoms[i].z);
        if (havePbc)
        {
            Vec3 frac = invCell * toEigen(delta);
            frac -= frac.array().round().matrix();
            delta = toGlm(cell * frac);
        }
        result.constraintDirection[i] = delta;
        result.interpolated.atoms[i].x = initial.atoms[i].x + zeta * delta.x;
        result.interpolated.atoms[i].y = initial.atoms[i].y + zeta * delta.y;
        result.interpolated.atoms[i].z = initial.atoms[i].z + zeta * delta.z;
    }

    result.success = true;
    result.message = "Interpolated " + std::to_string(n) + " atoms at zeta=" + std::to_string(zeta) + ".";
    return result;
}

}
