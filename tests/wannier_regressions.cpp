#include "electronic/Wannier.h"
#include "electronic/Volume.h"  // pi

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace atomforge::electronic;
using Complex = std::complex<double>;

namespace
{
void close(double actual, double expected, double tolerance = 1e-8)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error("Expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}
template<class Work> void rejects(Work work)
{
    try { work(); } catch (const std::exception&) { return; }
    throw std::runtime_error("Invalid input was accepted");
}

// Writes the Qi-Wu-Zhang two-band Chern-insulator model,
// H(k) = sin(2*pi*kx)*sx + sin(2*pi*ky)*sy + (m+cos(2*pi*kx)+cos(2*pi*ky))*sz,
// as a nearest-neighbour Wannier90 seedname_hr.dat fixture.
std::filesystem::path writeQwzFixture(double m)
{
    const std::vector<std::pair<glm::ivec3, std::array<Complex, 4>>> blocks{
        {{0, 0, 0}, {Complex(m, 0), Complex(0, 0), Complex(0, 0), Complex(-m, 0)}},
        {{1, 0, 0}, {Complex(0.5, 0), Complex(0, -0.5), Complex(0, -0.5), Complex(-0.5, 0)}},
        {{-1, 0, 0}, {Complex(0.5, 0), Complex(0, 0.5), Complex(0, 0.5), Complex(-0.5, 0)}},
        {{0, 1, 0}, {Complex(0.5, 0), Complex(-0.5, 0), Complex(0.5, 0), Complex(-0.5, 0)}},
        {{0, -1, 0}, {Complex(0.5, 0), Complex(0.5, 0), Complex(-0.5, 0), Complex(-0.5, 0)}},
    };
    const auto path = std::filesystem::temp_directory_path() / "atomforge_wannier_qwz_fixture_hr.dat";
    std::ofstream out(path);
    out << "fixture\n2\n" << blocks.size() << "\n";
    for (std::size_t i = 0; i < blocks.size(); ++i) out << "1" << (i + 1 == blocks.size() ? "\n" : " ");
    for (const auto& [R, matrix] : blocks)
        for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col)
            {
                const auto& v = matrix[static_cast<std::size_t>(row) * 2 + col];
                out << R.x << " " << R.y << " " << R.z << " " << (row + 1) << " " << (col + 1)
                    << " " << v.real() << " " << v.imag() << "\n";
            }
    return path;
}

Eigen::Vector3d dVector(double m, double kx, double ky)
{
    return {std::sin(2 * pi * kx), std::sin(2 * pi * ky), m + std::cos(2 * pi * kx) + std::cos(2 * pi * ky)};
}
}

int main()
{
    try
    {
        const double m = -1.0;
        const auto path = writeQwzFixture(m);
        const auto hr = readWannierHr(path.string());
        if (hr.numWann != 2 || hr.vectors.size() != 5) throw std::runtime_error("Unexpected Wannier90 Hamiltonian dimensions");

        const double kx = 0.13, ky = -0.27;
        const glm::dvec3 kpoint(kx, ky, 0.0);
        const auto d = dVector(m, kx, ky);
        const auto h = blochHamiltonian(hr, kpoint);
        const Complex sx[4] = {0, 1, 1, 0}, sy[4] = {0, Complex(0, -1), Complex(0, 1), 0}, sz[4] = {1, 0, 0, -1};
        for (int i = 0; i < 4; ++i)
        {
            const Complex expected = d(0) * sx[i] + d(1) * sy[i] + d(2) * sz[i];
            close(h[static_cast<std::size_t>(i)].real(), expected.real(), 1e-10);
            close(h[static_cast<std::size_t>(i)].imag(), expected.imag(), 1e-10);
        }

        const auto bands = wannierBands(hr, {kpoint});
        const double norm = d.norm();
        close(std::min(bands[0], bands[1]), -norm, 1e-10);
        close(std::max(bands[0], bands[1]), norm, 1e-10);

        // Independent numerical cross-check: a small gauge-invariant plaquette's
        // discretized Berry phase / area, sharing no formula with berryCurvature().
        auto plaquetteCurvature = [&](int band, double x0, double y0, double dk = 1e-4)
        {
            const double xs[4] = {x0, x0 + dk, x0 + dk, x0}, ys[4] = {y0, y0, y0 + dk, y0 + dk};
            Eigen::VectorXcd states[4];
            for (int c = 0; c < 4; ++c)
            {
                const auto hk = blochHamiltonian(hr, {xs[c], ys[c], 0});
                Eigen::Matrix2cd matrix; matrix << hk[0], hk[1], hk[2], hk[3];
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd> solver(matrix);
                states[c] = solver.eigenvectors().col(band);
            }
            Complex product(1, 0);
            for (int c = 0; c < 4; ++c)
            {
                const Complex link = states[c].dot(states[(c + 1) % 4]);
                product *= link / std::abs(link);
            }
            return std::atan2(product.imag(), product.real()) / (dk * dk);
        };
        const auto curvature = berryCurvature(hr, kpoint, 0, 1);
        close(curvature[0], plaquetteCurvature(0, kx, ky), 5e-3);
        close(curvature[0] + curvature[1], 0, 1e-8);  // 2-band identity

        // Integrating the Kubo curvature over the fractional Brillouin zone
        // must reproduce the same integer as the independent, non-perturbative
        // Fukui-Hatsugai-Suzuki method (chernNumber()).
        const int gridN = 40;
        double total = 0.0;
        for (int i = 0; i < gridN; ++i)
            for (int j = 0; j < gridN; ++j)
                total += berryCurvature(hr, {double(i) / gridN, double(j) / gridN, 0}, 0, 1)[0];
        const double integratedChern = total / (2 * pi * gridN * gridN);
        const double fhsChern = chernNumber(hr, 0, 0, 1, gridN, 0.0);
        close(integratedChern, fhsChern, 1e-3);
        close(std::abs(fhsChern), 1, 1e-6);

        const auto flippedHr = readWannierHr(writeQwzFixture(-m).string());
        close(chernNumber(flippedHr, 0, 0, 1, 24, 0.0), -fhsChern, 1e-6);

        const auto trivialHr = readWannierHr(writeQwzFixture(3.0).string());
        close(chernNumber(trivialHr, 0, 0, 1, 24, 0.0), 0, 1e-6);

        // At m=-2 the gap closes exactly at k=0: Kubo curvature must diverge.
        const auto criticalHr = readWannierHr(writeQwzFixture(-2.0).string());
        rejects([&] { berryCurvature(criticalHr, {0, 0, 0}, 0, 1); });
        rejects([&] { berryCurvature(hr, {0, 0, 0}, 0, 0); });
        rejects([&] { chernNumber(hr, 5, 0, 1, 24, 0.0); });
        rejects([&] { readWannierHr((std::filesystem::temp_directory_path() / "does_not_exist_hr.dat").string()); });

        std::filesystem::remove(path);
        std::cout << "Wannier analytic regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
