#include "electronic/Volume.h"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace atomforge::electronic
{
EwaldResult ewald(const glm::dmat3& cell, const std::vector<glm::dvec3>& sites,
    const std::vector<double>& charges, double alpha, double realCutoff, double reciprocalCutoff)
{
    Grid check;
    check.cell = cell;
    check.values.resize(8);
    check.validate();
    if (sites.empty() || sites.size() != charges.size()) throw std::invalid_argument("Ewald requires one explicit charge per site");
    if (!std::isfinite(alpha) || !std::isfinite(realCutoff) || !std::isfinite(reciprocalCutoff) || alpha <= 0 || realCutoff <= 0 || reciprocalCutoff <= 0)
        throw std::invalid_argument("Ewald alpha and both cutoffs must be positive and finite");
    double total = 0, magnitude = 0;
    for (std::size_t i = 0; i < sites.size(); ++i)
    {
        if (!std::isfinite(charges[i])) throw std::invalid_argument("Non-finite charge");
        for (int a = 0; a < 3; ++a) if (!std::isfinite(sites[i][a])) throw std::invalid_argument("Non-finite site position");
        total += charges[i]; magnitude += std::abs(charges[i]);
    }
    if (std::abs(total) > 1e-10 * std::max(1.0, magnitude)) throw std::invalid_argument("Ewald requires a neutral cell; no implicit compensating background is used");
    const auto inverse = glm::inverse(cell);
    const auto reciprocal = 2 * pi * glm::transpose(inverse);
    const double volume = std::abs(glm::determinant(cell));
    glm::ivec3 realBounds, reciprocalBounds;
    double realWork = static_cast<double>(sites.size()) * sites.size(), reciprocalWork = sites.size();
    for (int a = 0; a < 3; ++a)
    {
        const double rb = std::ceil(realCutoff * glm::length(glm::dvec3(inverse[0][a], inverse[1][a], inverse[2][a]))) + 1;
        const double kb = std::ceil(reciprocalCutoff * glm::length(cell[a]) / (2 * pi));
        if (rb > 1000 || kb > 1000) throw std::invalid_argument("Excessive Ewald image bounds");
        realBounds[a] = static_cast<int>(rb); reciprocalBounds[a] = static_cast<int>(kb);
        realWork *= 2 * rb + 1; reciprocalWork *= 2 * kb + 1;
    }
    if (realWork + reciprocalWork > 2e8) throw std::invalid_argument("Ewald workload too large; reduce cutoffs or use a dedicated solver");
    std::vector<glm::dvec3> positions;
    for (auto site : sites)
    {
        auto f = inverse * site;
        positions.push_back(cell * (f - glm::floor(f)));
    }
    EwaldResult out;
    out.potentials.assign(sites.size(), 0);
    for (std::size_t i = 0; i < sites.size(); ++i)
    {
        for (std::size_t j = 0; j < sites.size(); ++j)
            for (int z = -realBounds.z; z <= realBounds.z; ++z)
                for (int y = -realBounds.y; y <= realBounds.y; ++y)
                    for (int x = -realBounds.x; x <= realBounds.x; ++x)
                    {
                        if (i == j && x == 0 && y == 0 && z == 0) continue;
                        const double r = glm::length(positions[i] - positions[j] + cell * glm::dvec3(x, y, z));
                        if (r < 1e-10) throw std::invalid_argument("Coincident periodic Ewald sites");
                        if (r <= realCutoff) out.potentials[i] += charges[j] * std::erfc(alpha * r) / r;
                    }
        out.potentials[i] -= 2 * alpha * charges[i] / std::sqrt(pi);
    }
    for (int z = -reciprocalBounds.z; z <= reciprocalBounds.z; ++z)
        for (int y = -reciprocalBounds.y; y <= reciprocalBounds.y; ++y)
            for (int x = -reciprocalBounds.x; x <= reciprocalBounds.x; ++x)
            {
                if (x == 0 && y == 0 && z == 0) continue;
                const auto k = reciprocal * glm::dvec3(x, y, z);
                const double k2 = glm::dot(k, k);
                if (k2 > reciprocalCutoff * reciprocalCutoff) continue;
                std::complex<double> rho = 0;
                for (std::size_t j = 0; j < sites.size(); ++j) rho += charges[j] * std::polar(1.0, glm::dot(k, positions[j]));
                const double factor = 4 * pi / volume * std::exp(-k2 / (4 * alpha * alpha)) / k2;
                for (std::size_t i = 0; i < sites.size(); ++i)
                    out.potentials[i] += factor * (rho * std::polar(1.0, -glm::dot(k, positions[i]))).real();
            }
    // Conducting (tin-foil) boundary: no dipole surface term. Potentials are V.
    const double coulomb = hartree * bohr;
    for (std::size_t i = 0; i < sites.size(); ++i)
    {
        out.potentials[i] *= coulomb;
        out.energy += 0.5 * charges[i] * out.potentials[i];
    }
    return out;
}
}
