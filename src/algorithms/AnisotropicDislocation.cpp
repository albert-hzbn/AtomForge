#include "algorithms/AnisotropicDislocation.h"

#include <Eigen/Dense>

#include <cmath>
#include <random>

namespace atomforge::dislocation
{

namespace
{
constexpr double kPi = 3.14159265358979323846;
using Mat3 = Eigen::Matrix3d;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Cplx = std::complex<double>;

ElasticTensor addVoigtNoise(const ElasticTensor& c, unsigned seed, double amplitude)
{
    if (amplitude <= 0.0)
        return c;

    double norm = 0.0;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
        norm += c.c[i][j][k][l] * c.c[i][j][k][l];
    norm = std::sqrt(norm);
    const double scale = amplitude * norm;

    // Perturb the 21 independent Voigt entries consistently (same random
    // value applied to every index permutation mapping to that entry), the
    // same workaround BABEL documents (CVoigt_noise) for breaking the
    // sextic formalism's degeneracy at elastically isotropic/high-symmetry
    // orientations.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    double delta[6][6] = {};
    for (int a = 0; a < 6; ++a)
        for (int b = a; b < 6; ++b)
            delta[a][b] = delta[b][a] = dist(rng) * scale;

    auto voigtIndex = [](int i, int j) {
        if (i == j) return i;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
        return 5;
    };

    ElasticTensor out = c;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
        out.c[i][j][k][l] += delta[voigtIndex(i, j)][voigtIndex(k, l)];
    return out;
}
}

StrohSolution solveStroh(const ElasticTensor& elasticLocalIn,
                         const glm::dvec3& burgersLocal,
                         unsigned noiseSeed,
                         double noiseAmplitude)
{
    StrohSolution sol;

    if (noiseAmplitude <= 0.0 && isNearlyIsotropic(elasticLocalIn))
    {
        sol.error = "Elastic constants are (nearly) isotropic; the Stroh sextic formalism is "
                    "singular there (repeated roots). Use the isotropic model instead, or pass "
                    "a small noiseAmplitude (e.g. 1e-4) to perturb the constants, matching "
                    "BABEL's own documented CVoigt_noise workaround for this case.";
        return sol;
    }

    const ElasticTensor elasticLocal = addVoigtNoise(elasticLocalIn, noiseSeed, noiseAmplitude);

    // Q_ik = C_i1k1, R_ik = C_i1k2, T_ik = C_i2k2 (0-indexed: 0=x=e1, 1=y=e2,
    // 2=z=line direction) -- Hirth & Lothe / Barnett-Lothe notation.
    Mat3 Q, R, T;
    for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 3; ++k)
    {
        Q(i, k) = elasticLocal.c[i][0][k][0];
        R(i, k) = elasticLocal.c[i][0][k][1];
        T(i, k) = elasticLocal.c[i][1][k][1];
    }

    if (std::abs(T.determinant()) < 1e-9 * std::pow(T.cwiseAbs().maxCoeff(), 3))
    {
        sol.error = "Elastic tensor is singular in the plane perpendicular to the line direction.";
        return sol;
    }
    const Mat3 Tinv = T.inverse();

    // Barnett-Lothe sextic (Stroh) matrix: eigenvalues are the roots p_alpha
    // of the characteristic equation, eigenvectors are [A; L] stacked.
    const Mat3 N1 = -Tinv * R.transpose();
    const Mat3 N2 = Tinv;
    const Mat3 N3 = R * Tinv * R.transpose() - Q;

    Mat6 N;
    N.block<3, 3>(0, 0) = N1;
    N.block<3, 3>(0, 3) = N2;
    N.block<3, 3>(3, 0) = N3;
    N.block<3, 3>(3, 3) = N1.transpose();

    const Eigen::EigenSolver<Mat6> solver(N);
    if (solver.info() != Eigen::Success)
    {
        sol.error = "Sextic eigenvalue solve did not converge.";
        return sol;
    }

    const auto& values = solver.eigenvalues();
    const auto& vectors = solver.eigenvectors();

    int found = 0;
    for (int i = 0; i < 6 && found < 3; ++i)
    {
        if (values(i).imag() > 1e-8)
        {
            sol.p[found] = values(i);
            for (int c = 0; c < 3; ++c)
            {
                sol.A[found][c] = vectors(c, i);
                sol.L[found][c] = vectors(c + 3, i);
            }
            ++found;
        }
    }

    if (found != 3)
    {
        sol.error = "Sextic roots are degenerate (elastically isotropic or a high-symmetry "
                    "orientation); retry with a small noiseAmplitude to break the degeneracy.";
        return sol;
    }

    // Fix D from: b = 2*Re[A_mat * D] (Burgers-vector circuit condition) and
    // 0 = 2*Re[L_mat * D] (single-valued stress / stress-function circuit
    // condition), stacked as one real 6x6 linear system for Re(D), Im(D).
    Eigen::Matrix3cd Amat, Lmat;
    for (int a = 0; a < 3; ++a)
        for (int c = 0; c < 3; ++c)
        {
            Amat(c, a) = sol.A[a][c];
            Lmat(c, a) = sol.L[a][c];
        }

    Mat6 system;
    system.block<3, 3>(0, 0) = Amat.real();
    system.block<3, 3>(0, 3) = -Amat.imag();
    system.block<3, 3>(3, 0) = Lmat.real();
    system.block<3, 3>(3, 3) = -Lmat.imag();

    Eigen::Matrix<double, 6, 1> rhs;
    rhs << 0.5 * burgersLocal.x, 0.5 * burgersLocal.y, 0.5 * burgersLocal.z, 0.0, 0.0, 0.0;

    const Eigen::Matrix<double, 6, 1> solution = system.fullPivLu().solve(rhs);
    const double residual = (system * solution - rhs).norm();
    if (!(residual < 1e-6 * std::max(1.0, rhs.norm())))
    {
        sol.error = "Burgers-vector linear system was singular or ill-conditioned.";
        return sol;
    }

    for (int a = 0; a < 3; ++a)
        sol.D[a] = Cplx(solution(a), solution(a + 3));

    sol.valid = true;
    return sol;
}

glm::dvec3 anisotropicDisplacement(const StrohSolution& solution, double x, double y)
{
    if (!solution.valid)
        return glm::dvec3(0.0);

    // Each term is singular only where x + p_alpha*y == 0, which (since
    // Im(p_alpha) > 0 strictly) is exactly the origin for a well-separated
    // root -- but a root whose imaginary part is only just above the
    // degeneracy threshold can make |x + p_alpha*y| collapse to near zero
    // along almost the whole line x + Re(p_alpha)*y = 0, far from the
    // origin, blowing log() up to -infinity there too. Floor |z| at a small
    // fraction of the characteristic in-plane scale (here, max(|x|,|y|,1))
    // so the field stays finite everywhere except a true near-core region,
    // matching the isotropic branch's own core-radius regularization.
    const double floor = 1e-6 * std::max({std::abs(x), std::abs(y), 1.0});

    glm::dvec3 u(0.0);
    for (int a = 0; a < 3; ++a)
    {
        Cplx z = Cplx(x, 0.0) + solution.p[a] * y;
        if (std::abs(z) < floor)
            z = z / std::max(std::abs(z), 1e-300) * floor;
        const Cplx term = solution.A[a][0] * solution.D[a] * std::log(z);
        const Cplx termY = solution.A[a][1] * solution.D[a] * std::log(z);
        const Cplx termZ = solution.A[a][2] * solution.D[a] * std::log(z);
        u.x += term.imag();
        u.y += termY.imag();
        u.z += termZ.imag();
    }
    return u / kPi;
}

}
