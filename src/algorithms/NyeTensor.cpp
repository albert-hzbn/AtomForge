#include "algorithms/NyeTensor.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>

namespace atomforge
{

namespace
{
using Mat3 = Eigen::Matrix3d;
using Vec3 = Eigen::Vector3d;

Vec3 toEigen(const glm::dvec3& v) { return Vec3(v.x, v.y, v.z); }

// Minimum-image displacement under periodicity `cell`/`invCell` (columns of
// `cell` are the lattice vectors); identity when usePbc is false.
glm::dvec3 minimumImage(const glm::dvec3& delta, bool usePbc, const Eigen::Matrix3d& cell, const Eigen::Matrix3d& invCell)
{
    if (!usePbc)
        return delta;
    Vec3 frac = invCell * toEigen(delta);
    frac -= frac.array().round().matrix();
    const Vec3 wrapped = cell * frac;
    return glm::dvec3(wrapped.x(), wrapped.y(), wrapped.z());
}
}

NyeTensorResult computeNyeTensor(const Structure& reference,
                                 const Structure& deformed,
                                 float cutoffRadius,
                                 bool usePbc)
{
    NyeTensorResult result;
    const size_t n = reference.atoms.size();
    if (n != deformed.atoms.size())
    {
        result.message = "Reference and deformed structures must have the same atom count and ordering.";
        return result;
    }
    result.alpha.assign(n, {});
    result.alphaNorm.assign(n, 0.0);
    if (n == 0)
    {
        result.success = true;
        return result;
    }
    if (cutoffRadius <= 0.0f)
    {
        result.message = "cutoffRadius must be positive.";
        return result;
    }

    usePbc = usePbc && reference.hasUnitCell;
    Mat3 cell0 = Mat3::Identity(), invCell0 = Mat3::Identity();
    Mat3 cellDelta = Mat3::Zero(); // deformed cell - reference cell, for the BABEL-matching strain correction
    if (usePbc)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                cell0(j, i) = reference.cellVectors[i][j]; // column i = lattice vector i
        const double det = cell0.determinant();
        if (std::abs(det) < 1e-12)
        {
            result.message = "Reference unit cell is singular.";
            return result;
        }
        invCell0 = cell0.inverse();
        if (deformed.hasUnitCell)
        {
            Mat3 cellD = Mat3::Identity();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    cellD(j, i) = deformed.cellVectors[i][j];
            cellDelta = cellD - cell0;
        }
    }

    std::vector<glm::dvec3> pos0(n), u(n);
    for (size_t i = 0; i < n; ++i)
    {
        pos0[i] = glm::dvec3(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);
        u[i] = glm::dvec3(deformed.atoms[i].x, deformed.atoms[i].y, deformed.atoms[i].z) - pos0[i];
    }

    // O(n^2) neighbor search: simple and always correct; fine for the
    // few-thousand-atom scale this analysis (and BABEL's own reference
    // examples) is used at. A grid-accelerated search would only matter for
    // much larger systems than this feature is intended for.
    const double cutoff2 = (double)cutoffRadius * (double)cutoffRadius;
    std::vector<std::vector<int>> neighbors(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            if (i == j) continue;
            const glm::dvec3 delta = minimumImage(pos0[j] - pos0[i], usePbc, cell0, invCell0);
            if (glm::dot(delta, delta) <= cutoff2)
                neighbors[i].push_back((int)j);
        }
    }

    // Step 1: per-atom displacement gradient du[i][a][b] = d U_a / d X_b,
    // fit by least squares from neighbor displacement differences (BABEL's
    // gradDisplacement.f90 BuildGradDisplacement).
    std::vector<Mat3> du(n, Mat3::Zero());
    std::vector<bool> defined(n, false);
    int skipped = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const auto& nbrs = neighbors[i];
        if (nbrs.size() < 4) { ++skipped; continue; }

        Mat3 amat = Mat3::Zero(); // sum R0_n (outer) dU_n
        Mat3 bmat = Mat3::Zero(); // sum R0_n (outer) R0_n
        for (int j : nbrs)
        {
            const glm::dvec3 r0 = minimumImage(pos0[j] - pos0[i], usePbc, cell0, invCell0);
            glm::dvec3 du_n = u[j] - u[i];
            if (usePbc)
            {
                // Integer periodic-image count (NOT the wrapped remainder):
                // a neighbor found through the k-th periodic image needs its
                // displacement corrected by that many copies of how much the
                // cell itself changed, matching BABEL's `ds0 = aNInt(ds0)`.
                const Vec3 ds0 = (invCell0 * toEigen(pos0[j] - pos0[i])).array().round().matrix();
                const Vec3 correction = cellDelta * ds0;
                du_n -= glm::dvec3(correction.x(), correction.y(), correction.z());
            }
            const Vec3 r0e = toEigen(r0), due = toEigen(du_n);
            amat += r0e * due.transpose();
            bmat += r0e * r0e.transpose();
        }

        if (std::abs(bmat.determinant()) < 1e-12 * std::pow(bmat.norm(), 3) + 1e-300)
        {
            ++skipped; // neighbors are (near-)coplanar; underdetermined fit
            continue;
        }
        const Mat3 invBmat = bmat.inverse();
        const Mat3 ftotal = amat * invBmat; // ftotal(a,b) = d U_b / d X_a
        du[i] = ftotal.transpose();          // du(a,b) = d U_a / d X_b
        defined[i] = true;
    }

    // Step 2: lattice correspondence tensor G (Hartley & Mishin 2005),
    // G = -[(I+du)^-1 * du]^T.
    std::vector<Mat3> G(n, Mat3::Zero());
    for (size_t i = 0; i < n; ++i)
    {
        if (!defined[i]) continue;
        const Mat3 F = Mat3::Identity() + du[i];
        const Mat3 Finv = F.inverse();
        G[i] = -(Finv * du[i]).transpose();
    }

    // Step 3: fit the spatial gradient of G over the same neighbor
    // geometry, then take its curl to get the Nye tensor (nyeTensor.f90
    // BuildNyeTensor, including the documented "David Olmsted" index-order
    // correction to the original Hartley & Mishin formula).
    for (size_t i = 0; i < n; ++i)
    {
        const auto& nbrs = neighbors[i];
        if (!defined[i] || nbrs.size() < 4) continue;

        Mat3 bmat = Mat3::Zero();
        std::vector<Vec3> r0s;
        std::vector<Mat3> dGs;
        r0s.reserve(nbrs.size());
        dGs.reserve(nbrs.size());
        bool anyUndefinedNeighbor = false;
        for (int j : nbrs)
        {
            if (!defined[j]) { anyUndefinedNeighbor = true; continue; }
            const glm::dvec3 r0 = minimumImage(pos0[j] - pos0[i], usePbc, cell0, invCell0);
            const Vec3 r0e = toEigen(r0);
            r0s.push_back(r0e);
            dGs.push_back(G[j] - G[i]);
            bmat += r0e * r0e.transpose();
        }
        if (r0s.size() < 4 || std::abs(bmat.determinant()) < 1e-12 * std::pow(bmat.norm(), 3) + 1e-300)
        {
            if (anyUndefinedNeighbor || r0s.size() < 4) ++skipped;
            continue;
        }
        const Mat3 invBmat = bmat.inverse();

        Mat3 dGdx[3]; // dGdx[l](a,b) = d G(b,l) / d x(a)
        for (int l = 0; l < 3; ++l)
        {
            Mat3 vmat = Mat3::Zero();
            for (size_t k = 0; k < r0s.size(); ++k)
                vmat += r0s[k] * dGs[k].col(l).transpose();
            dGdx[l] = invBmat * vmat;
        }

        std::array<std::array<double, 3>, 3> alpha{};
        for (int l = 0; l < 3; ++l)
        {
            alpha[l][0] = -dGdx[l](1, 2) + dGdx[l](2, 1);
            alpha[l][1] = -dGdx[l](2, 0) + dGdx[l](0, 2);
            alpha[l][2] = -dGdx[l](0, 1) + dGdx[l](1, 0);
        }
        result.alpha[i] = alpha;

        double normSq = 0.0;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                normSq += alpha[a][b] * alpha[a][b];
        result.alphaNorm[i] = std::sqrt(normSq);
    }

    result.success = true;
    result.message = (skipped > 0)
        ? std::to_string(skipped) + " of " + std::to_string(n) + " atoms had too few/degenerate "
          "neighbors for a defined Nye tensor (left at zero); increase cutoffRadius if this is unexpected."
        : "Nye tensor computed for all atoms.";
    return result;
}

}
