#include "electronic/BaderPartition.h"
#include "electronic/Volume.h"
#include "util/TaskControl.h"

#include <array>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace atomforge::electronic
{
namespace
{
std::array<glm::ivec3, 26> neighborOffsets()
{
    std::array<glm::ivec3, 26> offsets{};
    int n = 0;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
                if (dx || dy || dz) offsets[static_cast<std::size_t>(n++)] = {dx, dy, dz};
    return offsets;
}
}

BaderPartition baderOnGrid(const Grid& density)
{
    density.validate();
    if (density.unit != "e/A^3") throw std::invalid_argument("Bader partitioning requires electron density in e/A^3");
    for (double v : density.values) if (v < 0) throw std::invalid_argument("Bader partitioning requires nonnegative density");
    const int nx = density.shape[0], ny = density.shape[1], nz = density.shape[2];
    const auto offsets = neighborOffsets();
    const auto steps = density.steps();
    std::array<double, 26> neighborDistance{};
    for (std::size_t i = 0; i < 26; ++i) neighborDistance[i] = glm::length(steps * glm::dvec3(offsets[i]));

    auto valueAt = [&](int x, int y, int z) { return density.values[density.index(x, y, z)]; };
    auto wrap = [](int v, int n) { return ((v % n) + n) % n; };

    // next[i]: linear index of the steepest-ascent neighbor of voxel i, or
    // -1 if voxel i is a local density maximum among its 26 neighbors.
    const std::size_t total = static_cast<std::size_t>(nx) * ny * nz;
    std::vector<int> next(total, -1);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                atomforge::taskCheckpoint();
                const double here = valueAt(x, y, z);
                double best = 0.0;
                int bestNeighbor = -1;
                for (std::size_t i = 0; i < 26; ++i)
                {
                    int nxp = x + offsets[i].x, nyp = y + offsets[i].y, nzp = z + offsets[i].z;
                    if (density.periodic) { nxp = wrap(nxp, nx); nyp = wrap(nyp, ny); nzp = wrap(nzp, nz); }
                    else if (nxp < 0 || nxp >= nx || nyp < 0 || nyp >= ny || nzp < 0 || nzp >= nz) continue;
                    const double slope = (valueAt(nxp, nyp, nzp) - here) / neighborDistance[i];
                    if (slope > best) { best = slope; bestNeighbor = density.index(nxp, nyp, nzp); }
                }
                next[density.index(x, y, z)] = bestNeighbor;
            }

    // Follow each voxel's ascent chain to its basin root (a local maximum),
    // memoizing every visited voxel along the way (path compression) so the
    // whole pass is amortized linear in the number of voxels.
    std::vector<int> rootOf(total, -1);
    for (std::size_t start = 0; start < total; ++start)
    {
        std::vector<std::size_t> path;
        std::size_t current = start;
        while (rootOf[current] == -1 && next[current] != -1)
        {
            path.push_back(current);
            current = static_cast<std::size_t>(next[current]);
        }
        const int root = rootOf[current] != -1 ? rootOf[current] : static_cast<int>(current);
        rootOf[current] = root;
        for (auto voxel : path) rootOf[voxel] = root;
    }

    std::unordered_map<int, int> rootToBasin;
    std::vector<int> rootList;
    BaderPartition result;
    result.basin.resize(total);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                const auto i = density.index(x, y, z);
                const auto [it, inserted] = rootToBasin.try_emplace(rootOf[i], static_cast<int>(rootList.size()));
                if (inserted) { rootList.push_back(rootOf[i]); result.basins.emplace_back(); }
                const int id = it->second;
                result.basin[i] = id;
                const double w = density.weight(x, y, z);
                result.basins[static_cast<std::size_t>(id)].charge += w * valueAt(x, y, z);
                result.basins[static_cast<std::size_t>(id)].volume += w;
            }
    for (std::size_t i = 0; i < rootList.size(); ++i)
    {
        const int r = rootList[i];
        const int x = r % nx, y = (r / nx) % ny, z = r / (nx * ny);
        result.basins[i].maximum = density.position(x, y, z);
    }
    return result;
}

std::vector<double> baderPopulations(const Grid& density, const std::vector<glm::dvec3>& sites)
{
    if (sites.empty()) throw std::invalid_argument("Supply at least one site");
    const auto partition = baderOnGrid(density);
    std::vector<double> populations(sites.size(), 0.0);
    for (const auto& basin : partition.basins)
    {
        std::size_t nearest = 0;
        double bestDistance2 = std::numeric_limits<double>::infinity();
        for (std::size_t s = 0; s < sites.size(); ++s)
        {
            glm::dvec3 delta = basin.maximum - sites[s];
            if (density.periodic) delta = minimumImage(delta, density.cell);
            const double d2 = glm::dot(delta, delta);
            if (d2 < bestDistance2) { bestDistance2 = d2; nearest = s; }
        }
        populations[nearest] += basin.charge;
    }
    return populations;
}
}
