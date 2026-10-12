#include "math/NeighbourSearch.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace neighbours
{
namespace
{
struct BinHash
{
    std::size_t operator()(const std::array<long long, 3>& key) const
    {
        return static_cast<std::size_t>((key[0] * 73856093LL) ^ (key[1] * 19349663LL) ^ (key[2] * 83492791LL));
    }
};
}

void visitWithin(const Structure& structure, double cutoff, bool periodic, const std::function<void(const Pair&)>& visit,
                 int maxImages)
{
    const std::size_t n = structure.atoms.size();
    if (n == 0 || !(cutoff > 0.0)) return;

    glm::dmat3 cell(1.0);  // columns: a, b, c
    glm::dmat3 inverse(1.0);
    int range[3] = {0, 0, 0};
    if (periodic) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) cell[r][c] = structure.cellVectors[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
        if (std::abs(glm::determinant(cell)) < 1e-12) throw std::runtime_error("The unit cell has no volume");
        inverse = glm::inverse(cell);
        // Planes of constant fractional coordinate k lie 1/|row k of the inverse| apart.
        for (int k = 0; k < 3; ++k) {
            const glm::dvec3 row(inverse[0][k], inverse[1][k], inverse[2][k]);
            range[k] = static_cast<int>(std::ceil(cutoff * glm::length(row)));
            if (range[k] > maxImages)
                throw std::runtime_error("The cutoff spans more than " + std::to_string(maxImages) +
                                         " periodic images of this cell; use a smaller cutoff or a supercell");
        }
    }

    std::vector<glm::dvec3> position(n);
    for (std::size_t i = 0; i < n; ++i) {
        glm::dvec3 p(structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z);
        if (periodic) {
            glm::dvec3 f = inverse * p;
            f -= glm::floor(f);
            p = cell * f;
        }
        position[i] = p;
    }

    const auto binOf = [cutoff](const glm::dvec3& p) {
        return std::array<long long, 3>{static_cast<long long>(std::floor(p.x / cutoff)),
                                        static_cast<long long>(std::floor(p.y / cutoff)),
                                        static_cast<long long>(std::floor(p.z / cutoff))};
    };
    std::unordered_map<std::array<long long, 3>, std::vector<int>, BinHash> bins;
    for (std::size_t j = 0; j < n; ++j) bins[binOf(position[j])].push_back(static_cast<int>(j));

    const double cutoff2 = cutoff * cutoff;
    for (std::size_t i = 0; i < n; ++i)
        for (int a = -range[0]; a <= range[0]; ++a)
            for (int b = -range[1]; b <= range[1]; ++b)
                for (int c = -range[2]; c <= range[2]; ++c) {
                    const glm::dvec3 shift = periodic ? cell * glm::dvec3(a, b, c) : glm::dvec3(0.0);
                    // Atom j's image at position[j] + shift is near i when position[j] is near position[i] - shift.
                    const auto centre = binOf(position[i] - shift);
                    for (long long dx = -1; dx <= 1; ++dx)
                        for (long long dy = -1; dy <= 1; ++dy)
                            for (long long dz = -1; dz <= 1; ++dz) {
                                const auto found = bins.find({centre[0] + dx, centre[1] + dy, centre[2] + dz});
                                if (found == bins.end()) continue;
                                for (int j : found->second) {
                                    const glm::dvec3 d = position[static_cast<std::size_t>(j)] + shift - position[i];
                                    const double d2 = glm::dot(d, d);
                                    if (d2 > cutoff2 || d2 < 1e-12) continue;  // also skips i itself
                                    visit(Pair{static_cast<int>(i), j, d, std::sqrt(d2)});
                                }
                            }
                }
}

std::vector<Pair> within(const Structure& structure, double cutoff, bool periodic)
{
    std::vector<Pair> pairs;
    // Stored pairs grow with the images searched, so at most 20.
    visitWithin(structure, cutoff, periodic, [&pairs](const Pair& pair) { pairs.push_back(pair); }, 20);
    return pairs;
}
}
