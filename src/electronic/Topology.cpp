#include "electronic/Topology.h"
#include "electronic/Volume.h"
#include "util/TaskControl.h"

#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace atomforge::electronic
{
namespace
{
struct UnionFind
{
    std::vector<int> parent;
    explicit UnionFind(std::size_t n) : parent(n) { for (std::size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i); }
    int find(int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) parent[a] = b; }
};
}

BettiNumbers bettiNumbers(const Grid& field, double threshold)
{
    field.validate();
    if (field.periodic)
        throw std::invalid_argument("Betti-number topology currently supports finite (non-periodic) grids only");
    if (!std::isfinite(threshold)) throw std::invalid_argument("threshold must be finite");
    const int nx = field.shape[0], ny = field.shape[1], nz = field.shape[2];
    auto solidAt = [&](int x, int y, int z)
    {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) return false;
        return field.values[field.index(x, y, z)] >= threshold;
    };
    auto voxelIndex = [&](int x, int y, int z) { return (z * ny + y) * nx + x; };
    const std::size_t total = static_cast<std::size_t>(nx) * ny * nz;

    // Betti0: face-connected components of the solid set (6-connectivity,
    // matching the cubical complex's actual face-gluing topology).
    // Betti2: face-connected components of the background that never touch
    // the grid boundary (Alexander duality: a bounded complementary
    // component of a compact set in R^3 is an enclosed cavity).
    UnionFind solidUF(total), bgUF(total);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                atomforge::taskCheckpoint();
                if (solidAt(x, y, z))
                {
                    if (solidAt(x + 1, y, z)) solidUF.unite(voxelIndex(x, y, z), voxelIndex(x + 1, y, z));
                    if (solidAt(x, y + 1, z)) solidUF.unite(voxelIndex(x, y, z), voxelIndex(x, y + 1, z));
                    if (solidAt(x, y, z + 1)) solidUF.unite(voxelIndex(x, y, z), voxelIndex(x, y, z + 1));
                }
                else
                {
                    if (x + 1 < nx && !solidAt(x + 1, y, z)) bgUF.unite(voxelIndex(x, y, z), voxelIndex(x + 1, y, z));
                    if (y + 1 < ny && !solidAt(x, y + 1, z)) bgUF.unite(voxelIndex(x, y, z), voxelIndex(x, y + 1, z));
                    if (z + 1 < nz && !solidAt(x, y, z + 1)) bgUF.unite(voxelIndex(x, y, z), voxelIndex(x, y, z + 1));
                }
            }
    std::unordered_set<int> solidRoots, allBackgroundRoots, unboundedBackgroundRoots;
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                if (solidAt(x, y, z)) { solidRoots.insert(solidUF.find(voxelIndex(x, y, z))); continue; }
                const int root = bgUF.find(voxelIndex(x, y, z));
                allBackgroundRoots.insert(root);
                if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1 || z == 0 || z == nz - 1)
                    unboundedBackgroundRoots.insert(root);
            }
    int betti2 = 0;
    for (int root : allBackgroundRoots) if (!unboundedBackgroundRoots.count(root)) ++betti2;

    // Euler characteristic via the doubled-coordinate encoding: every
    // sub-cell (vertex/edge/face/cube) of a solid voxel's closed unit cube is
    // a triple (2x+dx, 2y+dy, 2z+dz), dx/dy/dz in {0,1,2}; the number of odd
    // coordinates (0..3) gives its dimension. A hash set of these triples
    // deduplicates sub-cells shared between face-adjacent solid voxels.
    const long long strideX = 2LL * nx + 3, strideY = 2LL * ny + 3;
    std::unordered_set<long long> seen;
    long long counts[4] = {0, 0, 0, 0};
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                if (!solidAt(x, y, z)) continue;
                atomforge::taskCheckpoint();
                for (int dz = 0; dz < 3; ++dz)
                    for (int dy = 0; dy < 3; ++dy)
                        for (int dx = 0; dx < 3; ++dx)
                        {
                            const int X = 2 * x + dx, Y = 2 * y + dy, Z = 2 * z + dz;
                            const long long key = X + strideX * (Y + strideY * Z);
                            if (seen.insert(key).second) ++counts[(X & 1) + (Y & 1) + (Z & 1)];
                        }
            }
    // counts[3] (the all-odd cells) equals the number of solid voxels
    // exactly, since (2x+1,2y+1,2z+1) is unique to voxel (x,y,z); this is
    // asserted by the "single solid voxel" case in the native test.
    const long long euler = counts[0] - counts[1] + counts[2] - counts[3];

    BettiNumbers result;
    result.betti0 = static_cast<int>(solidRoots.size());
    result.betti2 = betti2;
    result.betti1 = result.betti0 + result.betti2 - static_cast<int>(euler);
    return result;
}

std::vector<glm::dvec4> bettiCurve(const Grid& field, const std::vector<double>& thresholds)
{
    if (thresholds.empty()) throw std::invalid_argument("Supply at least one threshold");
    std::vector<glm::dvec4> result;
    result.reserve(thresholds.size());
    for (double threshold : thresholds)
    {
        const auto numbers = bettiNumbers(field, threshold);
        result.push_back({threshold, numbers.betti0, numbers.betti1, numbers.betti2});
    }
    return result;
}
}
