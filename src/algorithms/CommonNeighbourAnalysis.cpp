#include "algorithms/CommonNeighbourAnalysis.h"
#include "util/ElementData.h"
#include "math/StructureMath.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace atomforge::analysis
{
namespace
{
constexpr float kSpatialHashCellSize = 4.0f;

struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u) ^ ((size_t)value.y * 19349663u) ^ ((size_t)value.z * 83492791u);
    }
};

glm::ivec3 getGridCell(const glm::vec3& position)
{
    return glm::ivec3((int)std::floor(position.x / kSpatialHashCellSize),
                      (int)std::floor(position.y / kSpatialHashCellSize),
                      (int)std::floor(position.z / kSpatialHashCellSize));
}

std::string classifyEnvironment(const std::map<Signature, int>& signatures, int coordination)
{
    const auto count = [&](Signature signature) {
        const auto found = signatures.find(signature);
        return found == signatures.end() ? 0 : found->second;
    };
    if (coordination == 12 && count({4, 2, 1}) == 12) return "FCC-like";
    if (coordination == 12 && count({4, 2, 1}) == 6 && count({4, 2, 2}) == 6) return "HCP-like";
    if (coordination == 14 && count({4, 4, 4}) == 6 && count({6, 6, 6}) == 8) return "BCC-like";
    if (coordination == 12 && count({5, 5, 5}) == 12) return "ICO-like";
    return "Unknown";
}

uint64_t edgeKey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

int longestChainLength(const std::vector<int>& commonNodes,
                       const std::unordered_set<uint64_t>& edgeSet)
{
    // Count bonds, not vertices. A ring includes its closing bond.
    // Bound the exponential search for highly coordinated, non-crystalline graphs.
    if (commonNodes.size() > 16) return -1;
    std::vector<std::vector<int>> adjacency(commonNodes.size());
    for (int i = 0; i < (int)commonNodes.size(); ++i)
        for (int j = i + 1; j < (int)commonNodes.size(); ++j)
            if (edgeSet.count(edgeKey(commonNodes[i], commonNodes[j]))) {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
    int best = 0;
    int visits = 0;
    bool exhausted = false;
    for (int start = 0; start < (int)commonNodes.size(); ++start) {
        const auto visit = [&](const auto& self, int node, unsigned visited, int length) -> void {
            if (++visits > 200000) { exhausted=true; return; }
            best = std::max(best, length);
            for (int next : adjacency[node]) {
                if (exhausted) return;
                if (next == start && length >= 2) best = std::max(best, length + 1);
                else if (!(visited & (1u << next)))
                    self(self, next, visited | (1u << next), length + 1);
            }
        };
        visit(visit, start, 1u << start, 0);
    }

    return exhausted ? -1 : best;
}

CnaResult runCna(const Structure& structure, float cutoffScale, bool usePbcRequest, float cutoffRadius)
{
    CnaResult result;
    result.atomCount = (int)structure.atoms.size();

    if (structure.atoms.empty())
    {
        result.valid = false;
        result.message = "No atoms available.";
        return result;
    }

    std::vector<float> radii = makeLiteratureCovalentRadii();

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    bool usePbc = false;
    if (usePbcRequest && structure.hasUnitCell)
    {
        usePbc = tryMakeCellMatrices(structure, cell, invCell);
    }
    result.pbcUsed = usePbc;

    std::vector<glm::vec3> positions(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        positions[i] = glm::vec3((float)structure.atoms[i].x,
                                 (float)structure.atoms[i].y,
                                 (float)structure.atoms[i].z);
    }

    std::vector<std::vector<int>> neighbors(structure.atoms.size());
    std::unordered_set<uint64_t> edgeSet;

    auto tryAddBond = [&](int i, int j)
    {
        if (j <= i)
            return;

        int zi = structure.atoms[i].atomicNumber;
        float ri = (zi >= 0 && zi < (int)radii.size()) ? radii[zi] : 1.0f;
        int zj = structure.atoms[j].atomicNumber;
        float rj = (zj >= 0 && zj < (int)radii.size()) ? radii[zj] : 1.0f;

        glm::vec3 delta = minimumImageDelta(positions[j] - positions[i], usePbc, cell, invCell);
        float d = glm::length(delta);
        if (d <= kMinBondDistance)
            return;

        float cutoff = cutoffRadius > 0.0f ? cutoffRadius : (ri + rj) * cutoffScale;
        if (d > cutoff)
            return;

        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        edgeSet.insert(edgeKey(i, j));
    };

    float maximumCutoff = cutoffRadius;
    if (maximumCutoff == 0.0f)
        for (const auto& atom : structure.atoms) {
            const int z = atom.atomicNumber;
            const float radius = z >= 0 && z < (int)radii.size() ? radii[z] : 1.0f;
            maximumCutoff = std::max(maximumCutoff, 2.0f * radius * cutoffScale);
        }
    const int reach = (int)std::ceil(maximumCutoff / kSpatialHashCellSize);
    if (!usePbc && reach <= 8)
    {
        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
        grid.reserve(positions.size());
        for (int i = 0; i < (int)positions.size(); ++i)
            grid[getGridCell(positions[i])].push_back(i);

        for (int i = 0; i < (int)positions.size(); ++i)
        {
            const glm::ivec3 cellCoord = getGridCell(positions[i]);
            for (int dx = -reach; dx <= reach; ++dx)
            {
                for (int dy = -reach; dy <= reach; ++dy)
                {
                    for (int dz = -reach; dz <= reach; ++dz)
                    {
                        const glm::ivec3 neighborCell(cellCoord.x + dx,
                                                      cellCoord.y + dy,
                                                      cellCoord.z + dz);
                        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash>::const_iterator it = grid.find(neighborCell);
                        if (it == grid.end())
                            continue;

                        const std::vector<int>& candidates = it->second;
                        for (int index = 0; index < (int)candidates.size(); ++index)
                            tryAddBond(i, candidates[index]);
                    }
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            for (int j = i + 1; j < (int)structure.atoms.size(); ++j)
                tryAddBond(i, j);
        }
    }

    for (int i = 0; i < (int)neighbors.size(); ++i)
        std::sort(neighbors[i].begin(), neighbors[i].end());

    std::vector<std::map<Signature, int>> atomSignatures(structure.atoms.size());

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        for (int t = 0; t < (int)neighbors[i].size(); ++t)
        {
            int j = neighbors[i][t];
            if (j <= i)
                continue;

            std::vector<int> common;
            common.reserve(std::min(neighbors[i].size(), neighbors[j].size()));
            std::set_intersection(neighbors[i].begin(), neighbors[i].end(),
                                  neighbors[j].begin(), neighbors[j].end(),
                                  std::back_inserter(common));

            int commonCount = (int)common.size();
            int bondCount = 0;
            for (int a = 0; a < (int)common.size(); ++a)
            {
                for (int b = a + 1; b < (int)common.size(); ++b)
                {
                    if (edgeSet.find(edgeKey(common[a], common[b])) != edgeSet.end())
                        ++bondCount;
                }
            }

            int chainLen = longestChainLength(common, edgeSet);

            Signature sig;
            sig.common = commonCount;
            sig.bonds = bondCount;
            sig.chain = chainLen;

            ++result.signatureCounts[sig];
            ++atomSignatures[i][sig];
            ++atomSignatures[j][sig];
            ++result.pairCount;
        }
    }

    result.atomRows.reserve(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        AtomRow row;
        row.index = i;
        row.atomicNumber = structure.atoms[i].atomicNumber;
        row.symbol = structure.atoms[i].symbol;
        row.coordination = (int)neighbors[i].size();

        Signature dominant;
        int dominantCount = 0;
        for (std::map<Signature, int>::const_iterator it = atomSignatures[i].begin(); it != atomSignatures[i].end(); ++it)
        {
            if (it->second > dominantCount)
            {
                dominant = it->first;
                dominantCount = it->second;
            }
        }

        row.dominantSignature = dominant;
        row.dominantSignatureCount = dominantCount;
        row.environment = classifyEnvironment(atomSignatures[i], row.coordination);

        ++result.environmentCounts[row.environment];
        result.atomRows.push_back(row);
    }

    result.valid = true;
    result.message = "CNA completed.";
    return result;
}

} // namespace

CnaResult computeCna(const Structure& structure, const CnaParams& params)
{
    if (!std::isfinite(params.cutoffScale) || params.cutoffScale <= 0.0f ||
        !std::isfinite(params.cutoffRadius) || params.cutoffRadius < 0.0f ||
        params.cutoffRadius > 100000.0f || params.cutoffScale > 100.0f)
    {
        CnaResult result;
        result.message = "Cutoff scale must be finite and positive.";
        return result;
    }
    for (const auto& atom : structure.atoms)
        if (!std::isfinite(atom.x) || !std::isfinite(atom.y) || !std::isfinite(atom.z) ||
            std::abs(atom.x) > 1e8 || std::abs(atom.y) > 1e8 || std::abs(atom.z) > 1e8) {
            CnaResult result;
            result.message = "Atom coordinates must be finite and within 1e8 A.";
            return result;
        }
    return runCna(structure, params.cutoffScale, params.usePbc, params.cutoffRadius);
}
} // namespace atomforge::analysis
