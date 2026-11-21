// Reference SelfAdjointEigenSolver<Matrix3d> decompositions for
// af_core::eigen::self_adjoint_eigen3 (tests/golden/self_adjoint_eigen.txt):
// each line is the matrix | eigenvalues, eigenvectors (row-major) as C hex
// floats. Random, near-degenerate, diagonal, rank-deficient and covariance
// matrices. Usage: golden_selfadjoint > rust/tests/golden/self_adjoint_eigen.txt
#include <Eigen/Dense>
#include <cstdio>
#include <random>

int main()
{
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int c = 0; c < 600; ++c) {
        Eigen::Matrix3d m;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) m(i, j) = u(rng);
        m = (m + m.transpose()).eval();
        switch (c % 6) {
        case 1: m = Eigen::Matrix3d::Identity() * 3.0 + m * 1e-9; break;
        case 2: m(2, 0) = m(0, 2) = 0; break;
        case 3: m = Eigen::Vector3d(u(rng), u(rng), u(rng)).asDiagonal(); break;
        case 4: {
            Eigen::Vector3d v(u(rng), u(rng), u(rng));
            m = v * v.transpose() * 50.0;
            break;
        }
        case 5: {
            Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
            const Eigen::Vector3d axis = Eigen::Vector3d(u(rng), u(rng), u(rng)).normalized();
            for (int k = 0; k < 12; ++k) {
                const Eigen::Vector3d d = axis * (k - 5.5) * 2.5 + Eigen::Vector3d(u(rng), u(rng), u(rng)) * 0.3;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) cov(i, j) += d[i] * d[j];
            }
            m = cov;
            break;
        }
        default: break;
        }
        if (c == 0) m.setZero();
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> e(m);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) std::printf("%a ", m(i, j));
        std::printf("|");
        for (int i = 0; i < 3; ++i) std::printf(" %a", e.eigenvalues()(i));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) std::printf(" %a", e.eigenvectors()(i, j));
        std::printf("\n");
    }
}
