#include "algorithms/CommonNeighbourAnalysis.h"
#include "math/NeighbourSearch.h"
#include <stdexcept>
#include "util/ElementData.h"
#include "math/StructureMath.h"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace atomforge::analysis
{
namespace
{
// A neighbour: atom j in periodic image n (cell translations), relative to
// the central atom in image 0. Without periodicity n is always 0.
struct Neighbour
{
    int j = 0;
    glm::ivec3 n{0};
    bool operator<(const Neighbour& o) const
    {
        if (j != o.j) return j < o.j;
        if (n.x != o.n.x) return n.x < o.n.x;
        if (n.y != o.n.y) return n.y < o.n.y;
        return n.z < o.n.z;
    }
    bool operator==(const Neighbour& o) const { return j == o.j && n == o.n; }
};

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

template <typename Bonded>
int longestChainLength(int nodeCount, Bonded&& bonded)
{
    // Count bonds, not vertices. A ring includes its closing bond.
    // Bound the exponential search for highly coordinated, non-crystalline graphs.
    if (nodeCount > 16) return -1;
    std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(nodeCount));
    for (int i = 0; i < nodeCount; ++i)
        for (int j = i + 1; j < nodeCount; ++j)
            if (bonded(i, j)) {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
    int best = 0;
    int visits = 0;
    bool exhausted = false;
    for (int start = 0; start < nodeCount; ++start) {
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

    float maximumCutoff = cutoffRadius;
    if (maximumCutoff == 0.0f)
        for (const auto& atom : structure.atoms) {
            const int z = atom.atomicNumber;
            const float radius = z >= 0 && z < (int)radii.size() ? radii[z] : 1.0f;
            maximumCutoff = std::max(maximumCutoff, 2.0f * radius * cutoffScale);
        }

    // Every image of every neighbour within the cutoff, so a small periodic
    // cell gives the same environments as a supercell of it.
    std::vector<neighbours::Pair> pairs;
    try {
        pairs = neighbours::within(structure, maximumCutoff, usePbc);
    } catch (const std::runtime_error& error) {
        result.valid = false;
        result.message = error.what();
        return result;
    }
    glm::dmat3 cellD(1.0), inverseD(1.0);
    if (usePbc) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) cellD[r][c] = structure.cellVectors[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
        inverseD = glm::inverse(cellD);
    }
    std::vector<std::vector<Neighbour>> neighbors(structure.atoms.size());
    for (const neighbours::Pair& pair : pairs) {
        if (pair.r <= kMinBondDistance) continue;
        const int zi = structure.atoms[pair.i].atomicNumber;
        const int zj = structure.atoms[pair.j].atomicNumber;
        const float ri = (zi >= 0 && zi < (int)radii.size()) ? radii[zi] : 1.0f;
        const float rj = (zj >= 0 && zj < (int)radii.size()) ? radii[zj] : 1.0f;
        const double cutoff = cutoffRadius > 0.0f ? cutoffRadius : (ri + rj) * cutoffScale;
        if (pair.r > cutoff) continue;
        Neighbour neighbour;
        neighbour.j = pair.j;
        if (usePbc) {
            const auto& from = structure.atoms[pair.i];
            const auto& to = structure.atoms[pair.j];
            const glm::dvec3 image = inverseD * (pair.d - glm::dvec3(to.x - from.x, to.y - from.y, to.z - from.z));
            neighbour.n = glm::ivec3((int)std::lround(image.x), (int)std::lround(image.y), (int)std::lround(image.z));
        }
        neighbors[pair.i].push_back(neighbour);
    }
    for (auto& list : neighbors) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
    // Is atom b in image n a neighbour of atom a in image 0?
    const auto bonded = [&](int a, int b, const glm::ivec3& n) {
        const auto& list = neighbors[a];
        Neighbour key;
        key.j = b;
        key.n = n;
        return std::binary_search(list.begin(), list.end(), key);
    };

    std::vector<std::map<Signature, int>> atomSignatures(structure.atoms.size());
    std::vector<Neighbour> common;
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        for (const Neighbour& partner : neighbors[i])
        {
            // Each bonded pair once: from the lower index (or, for an atom
            // and its own image, from the image with the positive offset).
            const int j = partner.j;
            if (j < i) continue;
            if (j == i && !(Neighbour{i, glm::ivec3(0)} < Neighbour{i, partner.n})) continue;

            common.clear();
            for (const Neighbour& k : neighbors[i])
                if (bonded(j, k.j, k.n - partner.n)) common.push_back(k);

            int bondCount = 0;
            for (std::size_t x = 0; x < common.size(); ++x)
                for (std::size_t y = x + 1; y < common.size(); ++y)
                    if (bonded(common[x].j, common[y].j, common[y].n - common[x].n)) ++bondCount;

            const int chainLen = longestChainLength((int)common.size(), [&](int x, int y) {
                return bonded(common[x].j, common[y].j, common[y].n - common[x].n);
            });

            Signature sig;
            sig.common = (int)common.size();
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
