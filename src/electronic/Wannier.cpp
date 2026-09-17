#include "electronic/Wannier.h"
#include "electronic/Volume.h"  // pi
#include "util/TaskControl.h"

#include <Eigen/Dense>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
using Complex = std::complex<double>;
using MatrixXc = Eigen::MatrixXcd;

std::string line(std::istream& in)
{
    atomforge::taskCheckpoint();
    std::string s;
    if (!std::getline(in, s)) throw std::invalid_argument("Truncated Wannier90 Hamiltonian file");
    return s;
}
std::vector<std::string> words(const std::string& s)
{
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string token; in >> token;) out.push_back(token);
    return out;
}
bool sameVector(const glm::ivec3& a, const glm::ivec3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

// H(k), or dH(k)/dk_derivative when derivative is 0, 1 or 2 (-1 for H(k) itself).
MatrixXc blochMatrix(const WannierHamiltonian& hr, glm::dvec3 kpoint, int derivative)
{
    if (!std::isfinite(kpoint.x) || !std::isfinite(kpoint.y) || !std::isfinite(kpoint.z))
        throw std::invalid_argument("kpoint must have three finite fractional components");
    MatrixXc result = MatrixXc::Zero(hr.numWann, hr.numWann);
    for (std::size_t i = 0; i < hr.vectors.size(); ++i)
    {
        atomforge::taskCheckpoint();
        const auto& R = hr.vectors[i];
        const double dot = kpoint.x * R.x + kpoint.y * R.y + kpoint.z * R.z;
        Complex weight = std::exp(Complex(0, 2 * pi * dot)) / static_cast<double>(hr.degeneracies[i]);
        if (derivative >= 0) weight *= Complex(0, 2 * pi * (derivative == 0 ? R.x : derivative == 1 ? R.y : R.z));
        const auto& h = hr.hamiltonians[i];
        for (int r = 0; r < hr.numWann; ++r)
            for (int c = 0; c < hr.numWann; ++c)
                result(r, c) += weight * h[static_cast<std::size_t>(r) * hr.numWann + c];
    }
    return result;
}

Eigen::SelfAdjointEigenSolver<MatrixXc> diagonalize(const MatrixXc& h)
{
    if ((h - h.adjoint()).cwiseAbs().maxCoeff() > 1e-8)
        throw std::invalid_argument("Interpolated Hamiltonian is not Hermitian; check R/-R hopping consistency");
    Eigen::SelfAdjointEigenSolver<MatrixXc> solver(h);
    if (solver.info() != Eigen::Success) throw std::runtime_error("Eigenvalue solver failed to converge");
    return solver;
}

void checkPlane(int a, int b)
{
    if (a == b || a < 0 || a > 2 || b < 0 || b > 2)
        throw std::invalid_argument("plane must select two distinct reciprocal directions");
}
}

WannierHamiltonian readWannierHr(const std::string& path)
{
    std::ifstream in{std::filesystem::u8path(path)};
    if (!in) throw std::invalid_argument("Cannot open Wannier90 Hamiltonian file");
    line(in);  // comment/date line, unused
    const auto header1 = words(line(in));
    const auto header2 = words(line(in));
    if (header1.size() != 1 || header2.size() != 1) throw std::invalid_argument("Invalid Wannier90 Hamiltonian header");
    const int numWann = std::stoi(header1[0]);
    const int numR = std::stoi(header2[0]);
    if (numWann < 1 || numWann > 2000 || numR < 1 || numR > 1000000)
        throw std::invalid_argument("Invalid Wannier90 Hamiltonian dimensions");

    std::vector<int> degeneracyValues;
    while (static_cast<int>(degeneracyValues.size()) < numR)
        for (const auto& token : words(line(in))) degeneracyValues.push_back(std::stoi(token));
    if (static_cast<int>(degeneracyValues.size()) != numR) throw std::invalid_argument("Invalid Wannier90 degeneracy block");
    for (int d : degeneracyValues) if (d < 1) throw std::invalid_argument("Invalid Wannier90 degeneracy block");

    WannierHamiltonian result;
    result.numWann = numWann;
    const std::size_t entriesPerR = static_cast<std::size_t>(numWann) * numWann;
    for (int degeneracy : degeneracyValues)
    {
        std::vector<Complex> matrix(entriesPerR);
        glm::ivec3 R{};
        bool first = true;
        for (std::size_t entry = 0; entry < entriesPerR; ++entry)
        {
            const auto tokens = words(line(in));
            if (tokens.size() != 7) throw std::invalid_argument("Invalid Wannier90 Hamiltonian row");
            const glm::ivec3 r{std::stoi(tokens[0]), std::stoi(tokens[1]), std::stoi(tokens[2])};
            const int m = std::stoi(tokens[3]) - 1, n = std::stoi(tokens[4]) - 1;
            if (first) { R = r; first = false; }
            else if (!sameVector(r, R)) throw std::invalid_argument("Hamiltonian block mixes R vectors");
            if (m < 0 || m >= numWann || n < 0 || n >= numWann) throw std::invalid_argument("Wannier index outside num_wann");
            matrix[static_cast<std::size_t>(m) * numWann + n] = Complex(std::stod(tokens[5]), std::stod(tokens[6]));
        }
        for (const auto& existing : result.vectors)
            if (sameVector(existing, R)) throw std::invalid_argument("Repeated R vector in Wannier90 Hamiltonian");
        result.vectors.push_back(R);
        result.hamiltonians.push_back(std::move(matrix));
        result.degeneracies.push_back(degeneracy);
    }
    return result;
}

std::vector<std::complex<double>> blochHamiltonian(const WannierHamiltonian& hr, glm::dvec3 kpoint)
{
    const auto m = blochMatrix(hr, kpoint, -1);
    std::vector<Complex> out(static_cast<std::size_t>(hr.numWann) * hr.numWann);
    for (int r = 0; r < hr.numWann; ++r)
        for (int c = 0; c < hr.numWann; ++c)
            out[static_cast<std::size_t>(r) * hr.numWann + c] = m(r, c);
    return out;
}

std::vector<double> wannierBands(const WannierHamiltonian& hr, const std::vector<glm::dvec3>& kpoints)
{
    std::vector<double> energies(kpoints.size() * static_cast<std::size_t>(hr.numWann));
    for (std::size_t i = 0; i < kpoints.size(); ++i)
    {
        atomforge::taskProgress(kpoints.empty() ? 1.0 : double(i) / kpoints.size());
        const auto solver = diagonalize(blochMatrix(hr, kpoints[i], -1));
        const auto& values = solver.eigenvalues();
        for (int b = 0; b < hr.numWann; ++b) energies[i * static_cast<std::size_t>(hr.numWann) + b] = values(b);
    }
    return energies;
}

std::vector<double> berryCurvature(const WannierHamiltonian& hr, glm::dvec3 kpoint, int planeA, int planeB)
{
    checkPlane(planeA, planeB);
    const auto solver = diagonalize(blochMatrix(hr, kpoint, -1));
    const auto values = solver.eigenvalues();
    const auto& states = solver.eigenvectors();
    const MatrixXc va = states.adjoint() * blochMatrix(hr, kpoint, planeA) * states;
    const MatrixXc vb = states.adjoint() * blochMatrix(hr, kpoint, planeB) * states;
    const int n = hr.numWann;
    std::vector<double> curvature(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i)
    {
        double sum = 0;
        for (int j = 0; j < n; ++j)
        {
            if (i == j) continue;
            const double gap = values(i) - values(j);
            const double contribution = 2 * (va(i, j) * vb(j, i)).imag() / (gap * gap);
            if (!std::isfinite(contribution))
                throw std::invalid_argument("Band " + std::to_string(i) + " is degenerate at this k point; Kubo curvature diverges");
            sum += contribution;
        }
        curvature[static_cast<std::size_t>(i)] = sum;
    }
    return curvature;
}

double chernNumber(const WannierHamiltonian& hr, int band, int planeA, int planeB, int grid, double fixed)
{
    checkPlane(planeA, planeB);
    if (grid < 2) throw std::invalid_argument("grid must be an integer >= 2");
    if (band < 0 || band >= hr.numWann) throw std::invalid_argument("band index is outside num_wann");
    std::vector<Eigen::VectorXcd> states(static_cast<std::size_t>(grid) * grid);
    for (int i = 0; i < grid; ++i)
        for (int j = 0; j < grid; ++j)
        {
            atomforge::taskProgress(double(i * grid + j) / (grid * grid));
            glm::dvec3 k(fixed, fixed, fixed);
            k[planeA] = double(i) / grid;
            k[planeB] = double(j) / grid;
            states[static_cast<std::size_t>(i) * grid + j] = diagonalize(blochMatrix(hr, k, -1)).eigenvectors().col(band);
        }
    double flux = 0.0;
    for (int i = 0; i < grid; ++i)
        for (int j = 0; j < grid; ++j)
        {
            const int i1 = (i + 1) % grid, j1 = (j + 1) % grid;
            const auto& s00 = states[static_cast<std::size_t>(i) * grid + j];
            const auto& s10 = states[static_cast<std::size_t>(i1) * grid + j];
            const auto& s11 = states[static_cast<std::size_t>(i1) * grid + j1];
            const auto& s01 = states[static_cast<std::size_t>(i) * grid + j1];
            Complex product(1, 0);
            for (Complex link : {s00.dot(s10), s10.dot(s11), s11.dot(s01), s01.dot(s00)})
            {
                const double magnitude = std::abs(link);
                if (magnitude < 1e-10)
                    throw std::invalid_argument("Band " + std::to_string(band) + " touches a degeneracy on this mesh; refine grid or shift fixed");
                product *= link / magnitude;
            }
            flux += std::atan2(product.imag(), product.real());
        }
    return flux / (2 * pi);
}
}
