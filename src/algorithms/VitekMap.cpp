#include "algorithms/VitekMap.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>

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

// Nearest lattice-vector wrap in the cell's first two lattice-vector
// directions only (matching BABEL's dsu(1)/dsu(2) check, which never
// adjusts the along-line component here).
glm::dvec3 wrapInPlane(const glm::dvec3& v, bool havePbc, const Mat3& cell, const Mat3& invCell)
{
    if (!havePbc)
        return v;
    Vec3 frac = invCell * toEigen(v);
    const double n0 = std::round(frac(0));
    const double n1 = std::round(frac(1));
    if (n0 == 0.0 && n1 == 0.0)
        return v;
    Vec3 shift = n0 * cell.col(0) + n1 * cell.col(1);
    return v - toGlm(shift);
}
}

VitekResult computeVitekMap(const Structure& reference,
                            const Structure& deformed,
                            const glm::dvec3& lineDirectionIn,
                            double burgersMagnitude,
                            float neighborCutoff,
                            bool usePbcXY)
{
    VitekResult result;
    const size_t n = reference.atoms.size();
    if (n != deformed.atoms.size())
    {
        result.message = "Reference and deformed structures must have the same atom count and ordering.";
        return result;
    }
    if (n < 2)
    {
        result.message = "Need at least two atoms.";
        return result;
    }

    const glm::dvec3 lineDir = glm::normalize(lineDirectionIn);
    glm::dvec3 e1 = glm::abs(glm::dot(lineDir, glm::dvec3(1, 0, 0))) < 0.9
        ? glm::normalize(glm::cross(lineDir, glm::dvec3(1, 0, 0)))
        : glm::normalize(glm::cross(lineDir, glm::dvec3(0, 1, 0)));
    const glm::dvec3 e2 = glm::normalize(glm::cross(lineDir, e1));

    std::vector<glm::dvec3> pos0(n), u(n);
    for (size_t i = 0; i < n; ++i)
    {
        pos0[i] = glm::dvec3(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);
        u[i] = glm::dvec3(deformed.atoms[i].x, deformed.atoms[i].y, deformed.atoms[i].z) - pos0[i];
    }

    Mat3 refCell = Mat3::Identity(), refInv = Mat3::Identity();
    Mat3 defCell = Mat3::Identity(), defInv = Mat3::Identity();
    const bool haveRefPbc = usePbcXY && buildCellMatrices(reference, refCell, refInv);
    const bool haveDefPbc = buildCellMatrices(deformed, defCell, defInv);

    // Auto-detect a first/second-neighbor-midpoint cutoff (in the plane
    // perpendicular to the line) from a handful of sample atoms, matching
    // the intent of BABEL's per-crystal-type Rcut table without needing to
    // hardcode bcc/fcc/hcp cases.
    if (neighborCutoff <= 0.0f)
    {
        std::vector<double> shellDistances;
        const size_t sampleCount = std::min<size_t>(n, 8);
        for (size_t s = 0; s < sampleCount; ++s)
        {
            std::vector<double> distances;
            for (size_t k = 0; k < n; ++k)
            {
                if (k == s) continue;
                glm::dvec3 delta = pos0[k] - pos0[s];
                delta -= lineDir * glm::dot(delta, lineDir);
                const double d = glm::length(delta);
                if (d > 1e-6) distances.push_back(d);
            }
            std::sort(distances.begin(), distances.end());
            for (double d : distances)
            {
                if (shellDistances.empty() || d > shellDistances.back() * 1.05)
                    shellDistances.push_back(d);
                if (shellDistances.size() >= 2) break;
            }
        }
        if (shellDistances.size() < 2)
        {
            result.message = "Could not auto-detect a neighbor cutoff; pass neighborCutoff explicitly.";
            return result;
        }
        std::sort(shellDistances.begin(), shellDistances.end());
        neighborCutoff = (float)(0.5 * (shellDistances.front() + shellDistances[1]));
    }
    const double cutoff2 = (double)neighborCutoff * (double)neighborCutoff;

    for (size_t i = 0; i < n; ++i)
    {
        glm::dvec3 ui = wrapInPlane(u[i], haveDefPbc, defCell, defInv);

        for (size_t j = i + 1; j < n; ++j)
        {
            glm::dvec3 posJ = pos0[j];
            glm::dvec3 dR0 = posJ - pos0[i];
            if (haveRefPbc)
            {
                Vec3 frac = refInv * toEigen(dR0);
                const double n0 = std::round(frac(0));
                const double n1 = std::round(frac(1));
                if (n0 != 0.0 || n1 != 0.0)
                {
                    const glm::dvec3 shift = -toGlm(n0 * refCell.col(0) + n1 * refCell.col(1));
                    dR0 += shift;
                    posJ += shift;
                }
            }

            glm::dvec3 perp = dR0 - lineDir * glm::dot(dR0, lineDir);
            const double r2 = glm::dot(perp, perp);
            if (r2 > cutoff2 || r2 < 1e-8)
                continue;

            glm::dvec3 uj = wrapInPlane(u[j], haveDefPbc, defCell, defInv);

            glm::dvec3 dU = uj - ui;
            double screw = glm::dot(dU, lineDir);
            if (burgersMagnitude > 0.0)
                screw -= burgersMagnitude * std::round(screw / burgersMagnitude);

            VitekPair pair;
            pair.indexI = (int)i;
            pair.indexJ = (int)j;
            pair.positionI = pos0[i];
            pair.positionJ = posJ;
            pair.screwComponent = screw;
            pair.edgeComponent = glm::dvec2(glm::dot(dU, e1), glm::dot(dU, e2));
            result.pairs.push_back(pair);
        }
    }

    result.success = true;
    result.message = "Found " + std::to_string(result.pairs.size()) + " neighbor pairs within cutoff "
        + std::to_string(neighborCutoff) + " A.";
    return result;
}

}
