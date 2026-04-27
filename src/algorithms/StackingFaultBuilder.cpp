#include "algorithms/StackingFaultBuilder.h"

#include "graphics/StructureInstanceBuilder.h"
#include "math/StructureMath.h"
#include "util/ElementData.h"

#include <glm/glm.hpp>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{
struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u)
             ^ ((size_t)value.y * 19349663u)
             ^ ((size_t)value.z * 83492791u);
    }
};

struct Signature
{
    int common = 0;
    int bonds = 0;
    int chain = 0;

    bool operator<(const Signature& other) const
    {
        if (common != other.common) return common < other.common;
        if (bonds != other.bonds) return bonds < other.bonds;
        return chain < other.chain;
    }

    bool operator==(const Signature& other) const
    {
        return common == other.common && bonds == other.bonds && chain == other.chain;
    }
};

struct DetectionSummary
{
    StackingFaultFamily family = StackingFaultFamily::Unknown;
    int fccCount = 0;
    int hcpCount = 0;
    int bccCount = 0;
    int recognizedCount = 0;
};

glm::ivec3 getGridCell(const glm::vec3& position)
{
    constexpr float kSpatialHashCellSize = 4.0f;
    return glm::ivec3((int)std::floor(position.x / kSpatialHashCellSize),
                      (int)std::floor(position.y / kSpatialHashCellSize),
                      (int)std::floor(position.z / kSpatialHashCellSize));
}

uint64_t edgeKey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

int longestChainLength(const std::vector<int>& commonNodes,
                       const std::unordered_set<uint64_t>& edgeSet)
{
    if (commonNodes.empty())
        return 0;
    if (commonNodes.size() == 1)
        return 1;

    std::vector<std::vector<int>> adjacency(commonNodes.size());
    for (int i = 0; i < (int)commonNodes.size(); ++i)
    {
        for (int j = i + 1; j < (int)commonNodes.size(); ++j)
        {
            if (edgeSet.find(edgeKey(commonNodes[i], commonNodes[j])) != edgeSet.end())
            {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
        }
    }

    int best = 1;
    for (int src = 0; src < (int)commonNodes.size(); ++src)
    {
        std::vector<int> dist(commonNodes.size(), -1);
        std::queue<int> q;
        dist[src] = 0;
        q.push(src);

        while (!q.empty())
        {
            const int u = q.front();
            q.pop();
            for (int v : adjacency[u])
            {
                if (dist[v] >= 0)
                    continue;
                dist[v] = dist[u] + 1;
                best = std::max(best, dist[v] + 1);
                q.push(v);
            }
        }
    }

    return best;
}

DetectionSummary detectFamily(const Structure& structure, bool usePbcRequest)
{
    DetectionSummary summary;
    if (structure.atoms.size() < 4)
        return summary;

    std::vector<float> radii = makeLiteratureCovalentRadii();
    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    const bool usePbc = usePbcRequest && structure.hasUnitCell
                     && tryMakeCellMatrices(structure, cell, invCell);
    constexpr size_t kMaxPairwisePbcAtoms = 5000;
    const bool usePairwisePbc = usePbc && structure.atoms.size() <= kMaxPairwisePbcAtoms;

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

        const int zi = structure.atoms[i].atomicNumber;
        const int zj = structure.atoms[j].atomicNumber;
        const float ri = (zi >= 0 && zi < (int)radii.size()) ? radii[zi] : 1.0f;
        const float rj = (zj >= 0 && zj < (int)radii.size()) ? radii[zj] : 1.0f;

        const glm::vec3 delta = minimumImageDelta(positions[j] - positions[i], usePairwisePbc, cell, invCell);
        const float distance = glm::length(delta);
        if (distance <= kMinBondDistance)
            return;

        const float cutoff = (ri + rj) * kBondToleranceFactor;
        if (distance > cutoff)
            return;

        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        edgeSet.insert(edgeKey(i, j));
    };

    if (!usePairwisePbc)
    {
        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
        grid.reserve(positions.size());
        for (int i = 0; i < (int)positions.size(); ++i)
            grid[getGridCell(positions[i])].push_back(i);

        for (int i = 0; i < (int)positions.size(); ++i)
        {
            const glm::ivec3 cellCoord = getGridCell(positions[i]);
            for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
            {
                const glm::ivec3 neighborCell(cellCoord.x + dx, cellCoord.y + dy, cellCoord.z + dz);
                auto it = grid.find(neighborCell);
                if (it == grid.end())
                    continue;
                for (int idx : it->second)
                    tryAddBond(i, idx);
            }
        }
    }
    else
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
            for (int j = i + 1; j < (int)structure.atoms.size(); ++j)
                tryAddBond(i, j);
    }

    for (std::vector<int>& nbrs : neighbors)
        std::sort(nbrs.begin(), nbrs.end());

    const Signature fcc = {4, 2, 1};
    const Signature hcp = {4, 2, 2};
    const Signature bccA = {4, 4, 1};
    const Signature bccB = {6, 6, 1};

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        std::map<Signature, int> counts;
        for (int neighbor : neighbors[i])
        {
            if (neighbor <= i)
                continue;

            std::vector<int> common;
            std::set_intersection(neighbors[i].begin(), neighbors[i].end(),
                                  neighbors[neighbor].begin(), neighbors[neighbor].end(),
                                  std::back_inserter(common));

            Signature sig;
            sig.common = (int)common.size();
            for (int a = 0; a < (int)common.size(); ++a)
            {
                for (int b = a + 1; b < (int)common.size(); ++b)
                {
                    if (edgeSet.find(edgeKey(common[a], common[b])) != edgeSet.end())
                        sig.bonds++;
                }
            }
            sig.chain = longestChainLength(common, edgeSet);
            counts[sig]++;
        }

        Signature dominant;
        int dominantCount = 0;
        for (const auto& entry : counts)
        {
            if (entry.second > dominantCount)
            {
                dominant = entry.first;
                dominantCount = entry.second;
            }
        }

        if (dominantCount <= 0)
            continue;

        if (dominant == fcc)
            summary.fccCount++;
        else if (dominant == hcp)
            summary.hcpCount++;
        else if (dominant == bccA || dominant == bccB)
            summary.bccCount++;
    }

    summary.recognizedCount = summary.fccCount + summary.hcpCount + summary.bccCount;
    if (summary.recognizedCount <= 0)
        return summary;

    if (summary.fccCount >= summary.hcpCount && summary.fccCount >= summary.bccCount)
        summary.family = StackingFaultFamily::Fcc;
    else if (summary.hcpCount >= summary.bccCount)
        summary.family = StackingFaultFamily::Hcp;
    else
        summary.family = StackingFaultFamily::Bcc;

    return summary;
}

glm::vec3 directionFromUvw(const glm::mat3& cell, const glm::ivec3& uvw)
{
    return cell[0] * (float)uvw.x + cell[1] * (float)uvw.y + cell[2] * (float)uvw.z;
}

glm::vec3 planeNormalFromHkl(const glm::mat3& cell, const glm::ivec3& hkl)
{
    const float volume = glm::dot(cell[0], glm::cross(cell[1], cell[2]));
    if (std::abs(volume) <= 1e-8f)
        return glm::vec3(0.0f, 0.0f, 1.0f);

    const glm::vec3 astar = glm::cross(cell[1], cell[2]) / volume;
    const glm::vec3 bstar = glm::cross(cell[2], cell[0]) / volume;
    const glm::vec3 cstar = glm::cross(cell[0], cell[1]) / volume;
    return astar * (float)hkl.x + bstar * (float)hkl.y + cstar * (float)hkl.z;
}

void resolvePlaneAndDirection(StackingFaultFamily family,
                              StackingFaultPlane requestedPlane,
                              glm::ivec3& planeHkl,
                              glm::ivec3& directionUvw,
                              std::string& planeLabel)
{
    switch (family)
    {
    case StackingFaultFamily::Fcc:
        planeHkl = glm::ivec3(1, 1, 2);
        directionUvw = glm::ivec3(1, 1, 1);
        planeLabel = "(112) <111>";
        return;
    case StackingFaultFamily::Hcp:
        if (requestedPlane == StackingFaultPlane::HcpPrismatic)
        {
            planeHkl = glm::ivec3(1, 0, 0);
            directionUvw = glm::ivec3(0, 1, 0);
            planeLabel = "prismatic";
        }
        else if (requestedPlane == StackingFaultPlane::HcpPyramidal)
        {
            planeHkl = glm::ivec3(1, 0, 1);
            directionUvw = glm::ivec3(1, -1, 1);
            planeLabel = "pyramidal";
        }
        else
        {
            planeHkl = glm::ivec3(0, 0, 1);
            directionUvw = glm::ivec3(1, -1, 0);
            planeLabel = "basal";
        }
        return;
    case StackingFaultFamily::Bcc:
        if (requestedPlane == StackingFaultPlane::Bcc112)
        {
            planeHkl = glm::ivec3(1, 1, 2);
            directionUvw = glm::ivec3(1, 1, -1);
            planeLabel = "(112) <111>";
        }
        else
        {
            planeHkl = glm::ivec3(1, 1, 0);
            directionUvw = glm::ivec3(1, -1, 1);
            planeLabel = "(110) <111>";
        }
        return;
    default:
        planeHkl = glm::ivec3(0, 0, 1);
        directionUvw = glm::ivec3(1, 0, 0);
        planeLabel = "unknown";
        return;
    }
}

std::string planeLabelForSelection(StackingFaultFamily family,
                                   StackingFaultPlane requestedPlane)
{
    if (family == StackingFaultFamily::Fcc)
        return "(112) <111>";
    if (family == StackingFaultFamily::Bcc)
        return (requestedPlane == StackingFaultPlane::Bcc112) ? "(112) <111>" : "(110) <111>";
    if (family == StackingFaultFamily::Hcp)
    {
        if (requestedPlane == StackingFaultPlane::HcpPrismatic)
            return "prismatic";
        if (requestedPlane == StackingFaultPlane::HcpPyramidal)
            return "pyramidal";
        return "basal";
    }
    return "unknown";
}

float familyPartialScale(StackingFaultFamily family)
{
    switch (family)
    {
    case StackingFaultFamily::Fcc:
        return 1.0f / 6.0f;
    case StackingFaultFamily::Hcp:
        return 1.0f / 3.0f;
    case StackingFaultFamily::Bcc:
        return 0.5f;
    default:
        return 0.25f;
    }
}

std::vector<int> assignLayerIds(const std::vector<float>& projections,
                                float& layerSpacing)
{
    std::vector<float> sorted = projections;
    std::sort(sorted.begin(), sorted.end());

    float minPositiveDiff = std::numeric_limits<float>::max();
    for (size_t i = 1; i < sorted.size(); ++i)
    {
        const float diff = sorted[i] - sorted[i - 1];
        if (diff > 1e-4f)
            minPositiveDiff = std::min(minPositiveDiff, diff);
    }

    if (!std::isfinite(minPositiveDiff))
        minPositiveDiff = 1.0f;

    const float tolerance = std::max(1e-3f, minPositiveDiff * 0.35f);

    std::vector<float> layerCenters;
    layerCenters.push_back(sorted.front());
    for (size_t i = 1; i < sorted.size(); ++i)
    {
        if (std::abs(sorted[i] - layerCenters.back()) > tolerance)
            layerCenters.push_back(sorted[i]);
    }

    layerSpacing = (layerCenters.size() >= 2)
        ? std::max(1e-3f, layerCenters[1] - layerCenters[0])
        : minPositiveDiff;

    std::vector<int> layerIds(projections.size(), 0);
    for (size_t i = 0; i < projections.size(); ++i)
    {
        int bestIdx = 0;
        float bestDist = std::abs(projections[i] - layerCenters[0]);
        for (int layer = 1; layer < (int)layerCenters.size(); ++layer)
        {
            const float dist = std::abs(projections[i] - layerCenters[layer]);
            if (dist < bestDist)
            {
                bestDist = dist;
                bestIdx = layer;
            }
        }
        layerIds[i] = bestIdx;
    }

    return layerIds;
}

std::vector<glm::ivec3> planeCandidates(StackingFaultFamily family,
                                        StackingFaultPlane requestedPlane)
{
    std::vector<glm::ivec3> candidates;

    if (family == StackingFaultFamily::Fcc)
    {
        candidates = {
            {1, 1, 2}, {1, 2, 1}, {2, 1, 1},
            {-1, 1, 2}, {1, -1, 2}, {1, 1, -2},
            {-1, -1, 2}, {-1, 1, -2}, {1, -1, -2}
        };
        return candidates;
    }

    if (family == StackingFaultFamily::Bcc)
    {
        if (requestedPlane == StackingFaultPlane::Bcc112)
        {
            candidates = {
                {1, 1, 2}, {1, 2, 1}, {2, 1, 1},
                {-1, 1, 2}, {1, -1, 2}, {1, 1, -2},
                {-1, -1, 2}, {-1, 1, -2}, {1, -1, -2}
            };
        }
        else
        {
            candidates = {
                {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
                {-1, 1, 0}, {1, -1, 0},
                {-1, 0, 1}, {1, 0, -1},
                {0, -1, 1}, {0, 1, -1}
            };
        }
        return candidates;
    }

    if (family == StackingFaultFamily::Hcp)
    {
        if (requestedPlane == StackingFaultPlane::HcpPrismatic)
        {
            candidates = {{1, 0, 0}, {0, 1, 0}, {-1, 1, 0},
                          {-1, 0, 0}, {0, -1, 0}, {1, -1, 0}};
        }
        else if (requestedPlane == StackingFaultPlane::HcpPyramidal)
        {
            candidates = {{1, 0, 1}, {0, 1, 1}, {-1, 1, 1},
                          {-1, 0, -1}, {0, -1, -1}, {1, -1, -1}};
        }
        else
        {
            candidates = {{0, 0, 1}, {0, 0, -1}};
        }
        return candidates;
    }

    candidates.push_back(glm::ivec3(0, 0, 1));
    return candidates;
}

std::vector<glm::ivec3> directionCandidates(StackingFaultFamily family,
                                            StackingFaultPlane requestedPlane)
{
    if (family == StackingFaultFamily::Fcc)
    {
        return {
            {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1},
            {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1}, {-1, -1, -1}
        };
    }

    if (family == StackingFaultFamily::Bcc)
    {
        return {
            {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1},
            {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1}, {-1, -1, -1}
        };
    }

    if (family == StackingFaultFamily::Hcp)
    {
        if (requestedPlane == StackingFaultPlane::HcpPyramidal)
        {
            return {
                {1, 0, 1}, {0, 1, 1}, {-1, 1, 1},
                {-1, 0, -1}, {0, -1, -1}, {1, -1, -1}
            };
        }
        return {
            {1, 0, 0}, {0, 1, 0}, {1, -1, 0},
            {-1, 0, 0}, {0, -1, 0}, {-1, 1, 0}
        };
    }

    return {glm::ivec3(1, 0, 0)};
}

float evaluatePlaneForStructure(const Structure& structure,
                                const glm::vec3& normal,
                                int& layerCountOut)
{
    layerCountOut = 0;
    if (structure.atoms.empty() || glm::length(normal) <= 1e-6f)
        return -1e9f;

    std::vector<float> projection(structure.atoms.size(), 0.0f);
    for (size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const glm::vec3 pos((float)structure.atoms[i].x,
                            (float)structure.atoms[i].y,
                            (float)structure.atoms[i].z);
        projection[i] = glm::dot(pos, normal);
    }

    float spacing = 0.0f;
    const std::vector<int> layerIds = assignLayerIds(projection, spacing);
    const int maxLayer = *std::max_element(layerIds.begin(), layerIds.end());
    const int layers = maxLayer + 1;
    layerCountOut = layers;
    if (layers < 2)
        return -1e9f;

    const int shiftStart = layers / 2;
    int shiftedCount = 0;
    for (int id : layerIds)
    {
        if (id >= shiftStart)
            shiftedCount++;
    }

    const float half = std::max(1.0f, 0.5f * (float)structure.atoms.size());
    const float imbalance = std::abs((float)shiftedCount - half) / half;
    return (float)layers - 0.5f * imbalance;
}

float estimateNearestNeighborDistance(const Structure& structure,
                                      bool usePbc,
                                      const glm::mat3& cell,
                                      const glm::mat3& invCell)
{
    if (structure.atoms.size() < 2)
        return 0.0f;

    const int n = (int)std::min<size_t>(structure.atoms.size(), 1024);
    float best = std::numeric_limits<float>::max();
    for (int i = 0; i < n; ++i)
    {
        const glm::vec3 pi((float)structure.atoms[i].x,
                           (float)structure.atoms[i].y,
                           (float)structure.atoms[i].z);
        for (int j = i + 1; j < n; ++j)
        {
            const glm::vec3 pj((float)structure.atoms[j].x,
                               (float)structure.atoms[j].y,
                               (float)structure.atoms[j].z);
            const glm::vec3 delta = minimumImageDelta(pj - pi, usePbc, cell, invCell);
            const float d = glm::length(delta);
            if (d > kMinBondDistance && d < best)
                best = d;
        }
    }

    if (std::isfinite(best) && best > 0.0f)
        return best;

    const float la = glm::length(cell[0]);
    const float lb = glm::length(cell[1]);
    const float lc = glm::length(cell[2]);
    best = std::min(la, std::min(lb, lc));
    return (best > 0.0f && std::isfinite(best)) ? best : 1.0f;
}

AtomSite representativeAtom(const Structure& structure)
{
    if (structure.atoms.empty())
        return AtomSite{};

    std::map<int, int> counts;
    for (const AtomSite& atom : structure.atoms)
    {
        if (atom.atomicNumber > 0)
            counts[atom.atomicNumber]++;
    }

    int bestZ = structure.atoms.front().atomicNumber;
    int bestCount = -1;
    for (const auto& [z, c] : counts)
    {
        if (c > bestCount)
        {
            bestCount = c;
            bestZ = z;
        }
    }

    for (const AtomSite& atom : structure.atoms)
    {
        if (atom.atomicNumber == bestZ)
            return atom;
    }

    return structure.atoms.front();
}

Structure buildFccPrimitiveLayeredCellLikeReference(const Structure& source,
                                                    int layers,
                                                    float nearestNeighbor)
{
    Structure out;
    out.hasUnitCell = true;
    out.cellOffset = {0.0, 0.0, 0.0};

    const float a = std::max(1e-5f, nearestNeighbor);
    const float b = a;
    const float cPerLayer = a * std::sqrt(2.0f / 3.0f);
    const float c = cPerLayer * (float)layers;

    out.cellVectors = {{
        {{a, 0.0, 0.0}},
        {{0.5f * b, 0.5f * std::sqrt(3.0f) * b, 0.0}},
        {{0.0, 0.0, c}}
    }};

    const AtomSite seed = representativeAtom(source);
    out.atoms.reserve((size_t)layers);
    for (int layer = 0; layer < layers; ++layer)
    {
        AtomSite atom = seed;
        const int mod = layer % 3;
        if (mod == 0)
        {
            atom.x = 0.0;
            atom.y = 0.0;
        }
        else if (mod == 1)
        {
            atom.x = 0.5 * a;
            atom.y = a / (2.0 * std::sqrt(3.0));
        }
        else
        {
            atom.x = a;
            atom.y = a / std::sqrt(3.0);
        }
        atom.z = ((double)layer / (double)layers) * (double)c;
        out.atoms.push_back(atom);
    }

    return out;
}

void refinePlaneAndDirection(const Structure& structure,
                             const glm::mat3& cell,
                             StackingFaultFamily family,
                             StackingFaultPlane requestedPlane,
                             glm::ivec3& planeHkl,
                             glm::ivec3& directionUvw)
{
    const std::vector<glm::ivec3> planes = planeCandidates(family, requestedPlane);
    glm::ivec3 bestPlane = planeHkl;
    float bestPlaneScore = -1e9f;
    glm::vec3 bestNormal(0.0f, 0.0f, 1.0f);

    for (const glm::ivec3& hkl : planes)
    {
        glm::vec3 normal = planeNormalFromHkl(cell, hkl);
        if (glm::length(normal) <= 1e-6f)
            continue;
        normal = glm::normalize(normal);

        int layers = 0;
        const float score = evaluatePlaneForStructure(structure, normal, layers);
        if (score > bestPlaneScore)
        {
            bestPlaneScore = score;
            bestPlane = hkl;
            bestNormal = normal;
        }
    }

    planeHkl = bestPlane;

    const std::vector<glm::ivec3> directions = directionCandidates(family, requestedPlane);
    glm::ivec3 bestDirection = directionUvw;
    float bestOrthogonality = std::numeric_limits<float>::max();
    float bestLength = std::numeric_limits<float>::max();
    for (const glm::ivec3& uvw : directions)
    {
        const glm::vec3 dir = directionFromUvw(cell, uvw);
        const float len = glm::length(dir);
        if (len <= 1e-6f)
            continue;

        const float orthogonality = std::abs(glm::dot(glm::normalize(dir), bestNormal));
        if (orthogonality < bestOrthogonality - 1e-5f
            || (std::abs(orthogonality - bestOrthogonality) <= 1e-5f && len < bestLength))
        {
            bestOrthogonality = orthogonality;
            bestLength = len;
            bestDirection = uvw;
        }
    }

    directionUvw = bestDirection;
}

Structure makeOrthogonalCellStructure(const Structure& source,
                                      const glm::vec3& slipDir,
                                      const glm::vec3& planeNormal)
{
    Structure output = source;

    glm::vec3 e1 = glm::normalize(slipDir);
    glm::vec3 e3 = glm::normalize(planeNormal);
    glm::vec3 e2 = glm::normalize(glm::cross(e3, e1));
    if (glm::length(e2) <= 1e-6f)
        e2 = glm::vec3(0.0f, 1.0f, 0.0f);

    std::array<float, 3> minCoord = {
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()};
    std::array<float, 3> maxCoord = {
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max()};

    std::vector<glm::vec3> localCoords;
    localCoords.reserve(source.atoms.size());
    for (const AtomSite& atom : source.atoms)
    {
        const glm::vec3 pos((float)atom.x, (float)atom.y, (float)atom.z);
        const glm::vec3 local(glm::dot(pos, e1), glm::dot(pos, e2), glm::dot(pos, e3));
        localCoords.push_back(local);
        minCoord[0] = std::min(minCoord[0], local.x);
        minCoord[1] = std::min(minCoord[1], local.y);
        minCoord[2] = std::min(minCoord[2], local.z);
        maxCoord[0] = std::max(maxCoord[0], local.x);
        maxCoord[1] = std::max(maxCoord[1], local.y);
        maxCoord[2] = std::max(maxCoord[2], local.z);
    }

    const glm::vec3 basisX = e1 * std::max(1.0f, maxCoord[0] - minCoord[0]);
    const glm::vec3 basisY = e2 * std::max(1.0f, maxCoord[1] - minCoord[1]);
    const glm::vec3 basisZ = e3 * std::max(1.0f, maxCoord[2] - minCoord[2]);

    output.hasUnitCell = true;
    output.cellVectors = {{
        {{basisX.x, basisX.y, basisX.z}},
        {{basisY.x, basisY.y, basisY.z}},
        {{basisZ.x, basisZ.y, basisZ.z}}
    }};
    output.cellOffset = {0.0, 0.0, 0.0};

    for (size_t i = 0; i < output.atoms.size(); ++i)
    {
        const glm::vec3 shifted = localCoords[i]
            - glm::vec3(minCoord[0], minCoord[1], minCoord[2]);
        const glm::vec3 world = e1 * shifted.x + e2 * shifted.y + e3 * shifted.z;
        output.atoms[i].x = world.x;
        output.atoms[i].y = world.y;
        output.atoms[i].z = world.z;
    }

    return output;
}

std::string buildSequenceLabel(StackingFaultFamily family, float factor, float interval)
{
    const float tol = std::max(0.02f, interval * 0.5f);
    if ((family == StackingFaultFamily::Fcc || family == StackingFaultFamily::Hcp)
        && std::abs(factor - 1.0f) <= tol)
        return "ISF";
    if ((family == StackingFaultFamily::Fcc || family == StackingFaultFamily::Hcp)
        && std::abs(factor - 2.0f) <= tol)
        return "ESF";

    std::ostringstream out;
    out.setf(std::ios::fixed);
    out.precision(2);
    out << "Shift " << factor;
    return out.str();
}

#ifdef ATOMS_ENABLE_SPGLIB
bool classifyFamilyWithSpglibHomogeneous(const Structure& structure,
                                         StackingFaultFamily& family,
                                         int& spaceGroup,
                                         std::string& internationalSymbol,
                                         std::string& error)
{
    family = StackingFaultFamily::Unknown;
    spaceGroup = 0;
    internationalSymbol.clear();
    error.clear();

    if (!structure.hasUnitCell)
    {
        error = "No unit cell available for spglib detection.";
        return false;
    }
    if (structure.atoms.empty())
    {
        error = "No atoms available for spglib detection.";
        return false;
    }

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(structure, cell, invCell))
    {
        error = "Invalid unit cell for spglib detection.";
        return false;
    }

    const glm::vec3 origin((float)structure.cellOffset[0],
                           (float)structure.cellOffset[1],
                           (float)structure.cellOffset[2]);

    double latticeData[3][3] = {
        {structure.cellVectors[0][0], structure.cellVectors[0][1], structure.cellVectors[0][2]},
        {structure.cellVectors[1][0], structure.cellVectors[1][1], structure.cellVectors[1][2]},
        {structure.cellVectors[2][0], structure.cellVectors[2][1], structure.cellVectors[2][2]}
    };

    std::vector<std::array<double, 3>> positions(structure.atoms.size());
    // Force all atoms into a single species so structure type is identified
    // from lattice+geometry, independent of chemical ordering.
    std::vector<int> types(structure.atoms.size(), 1);

    for (size_t index = 0; index < structure.atoms.size(); ++index)
    {
        const AtomSite& atom = structure.atoms[index];
        glm::vec3 cart((float)atom.x, (float)atom.y, (float)atom.z);
        glm::vec3 frac = invCell * (cart - origin);

        frac.x -= std::floor(frac.x);
        frac.y -= std::floor(frac.y);
        frac.z -= std::floor(frac.z);

        positions[index][0] = frac.x;
        positions[index][1] = frac.y;
        positions[index][2] = frac.z;
    }

    const std::array<double, 5> symprecs = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3};
    const SpglibDataset* dataset = nullptr;
    for (double symprec : symprecs)
    {
        dataset = spg_get_dataset(latticeData,
                                  reinterpret_cast<double (*)[3]>(positions.data()),
                                  types.data(),
                                  (int)positions.size(),
                                  symprec);
        if (dataset)
            break;
    }

    if (!dataset)
    {
        error = "spglib could not identify symmetry for this structure.";
        return false;
    }

    spaceGroup = dataset->spacegroup_number;
    internationalSymbol = dataset->international_symbol;

    char centering = '\0';
    if (!internationalSymbol.empty())
        centering = (char)std::toupper((unsigned char)internationalSymbol[0]);

    if (spaceGroup >= 195 && spaceGroup <= 230)
    {
        if (centering == 'F')
            family = StackingFaultFamily::Fcc;
        else if (centering == 'I')
            family = StackingFaultFamily::Bcc;
    }
    else if (spaceGroup >= 168 && spaceGroup <= 194)
    {
        family = StackingFaultFamily::Hcp;
    }

    spg_free_dataset(const_cast<SpglibDataset*>(dataset));

    if (family == StackingFaultFamily::Unknown)
    {
        std::ostringstream out;
        out << "spglib identified SG " << spaceGroup << " (" << internationalSymbol
            << "), but it does not map to FCC/HCP/BCC.";
        error = out.str();
        return false;
    }

    return true;
}
#endif
}

const char* stackingFaultFamilyName(StackingFaultFamily family)
{
    switch (family)
    {
    case StackingFaultFamily::Fcc:
        return "FCC-like";
    case StackingFaultFamily::Hcp:
        return "HCP-like";
    case StackingFaultFamily::Bcc:
        return "BCC-like";
    default:
        return "Unknown";
    }
}

StackingFaultDetectionResult detectStackingFaultFamily(const Structure& structure,
                                                       bool usePbcForDetection)
{
    StackingFaultDetectionResult result;
    if (structure.atoms.empty())
    {
        result.message = "Input structure has no atoms.";
        return result;
    }

#ifdef ATOMS_ENABLE_SPGLIB
    if (structure.hasUnitCell)
    {
        StackingFaultFamily spglibFamily = StackingFaultFamily::Unknown;
        int spglibSpaceGroup = 0;
        std::string spglibSymbol;
        std::string spglibError;
        if (classifyFamilyWithSpglibHomogeneous(structure,
                                                spglibFamily,
                                                spglibSpaceGroup,
                                                spglibSymbol,
                                                spglibError))
        {
            result.success = true;
            result.family = spglibFamily;
            result.recognizedCount = (int)structure.atoms.size();

            std::ostringstream message;
            message << "Detected " << stackingFaultFamilyName(spglibFamily)
                    << " using spglib (homogeneous species), SG "
                    << spglibSpaceGroup << " (" << spglibSymbol << ").";
            result.message = message.str();
            return result;
        }
    }
#endif

    const DetectionSummary detection = detectFamily(structure, usePbcForDetection);
    result.family = detection.family;
    result.fccCount = detection.fccCount;
    result.hcpCount = detection.hcpCount;
    result.bccCount = detection.bccCount;
    result.recognizedCount = detection.recognizedCount;

    if (detection.family == StackingFaultFamily::Unknown)
    {
        result.message = "Could not identify an FCC, HCP, or BCC-like environment.";
        return result;
    }

    std::ostringstream message;
    message << "Detected " << stackingFaultFamilyName(detection.family)
            << " environment (FCC=" << detection.fccCount
            << ", HCP=" << detection.hcpCount
            << ", BCC=" << detection.bccCount << ").";
    result.message = message.str();
    result.success = true;
    return result;
}

StackingFaultResult buildStackingFaultSequence(const Structure& base,
                                               const StackingFaultParams& params)
{
    StackingFaultResult result;
    if (base.atoms.empty())
    {
        result.message = "Input structure has no atoms.";
        return result;
    }
    if (!base.hasUnitCell)
    {
        result.message = "Stacking faults builder requires a structure with a unit cell.";
        return result;
    }
    if (params.layerCount < 2)
    {
        result.message = "Layer count must be at least 2.";
        return result;
    }
    if (params.interval <= 0.0f || params.maxDisplacementFactor < 0.0f)
    {
        result.message = "Interval and maximum displacement must be positive.";
        return result;
    }

    const StackingFaultDetectionResult detection = detectStackingFaultFamily(base,
                                                                             params.usePbcForDetection);
    result.family = detection.family;
    result.detectedPhase = stackingFaultFamilyName(detection.family);
    if (detection.family == StackingFaultFamily::Unknown)
    {
        result.message = "Could not identify an FCC, HCP, or BCC-like environment from the loaded structure.";
        return result;
    }

    if (detection.family == StackingFaultFamily::Fcc
        && params.cellMode == StackingFaultCellMode::SmallestUnitCell)
    {
        glm::mat3 baseCell(1.0f);
        glm::mat3 baseInvCell(1.0f);
        if (!tryMakeCellMatrices(base, baseCell, baseInvCell))
        {
            result.message = "Input unit cell is singular.";
            return result;
        }

        const float nearestNeighbor = estimateNearestNeighborDistance(base,
                                                                      params.usePbcForDetection,
                                                                      baseCell,
                                                                      baseInvCell);

        const int layers = std::max(3, params.layerCount);
        const int startLayerBase = std::max(0, layers - 3);

        Structure reference = buildFccPrimitiveLayeredCellLikeReference(base, layers, nearestNeighbor);

        const glm::vec3 aVec((float)reference.cellVectors[0][0],
                             (float)reference.cellVectors[0][1],
                             (float)reference.cellVectors[0][2]);
        const glm::vec3 bVec((float)reference.cellVectors[1][0],
                             (float)reference.cellVectors[1][1],
                             (float)reference.cellVectors[1][2]);
        const glm::vec3 cVec((float)reference.cellVectors[2][0],
                             (float)reference.cellVectors[2][1],
                             (float)reference.cellVectors[2][2]);

        const glm::vec3 slipDirection = glm::normalize(
            glm::vec3((float)(-aVec.x + 2.0 * bVec.x),
                      (float)(-aVec.y + 2.0 * bVec.y),
                      0.0f));

        const float bTwin = nearestNeighbor / std::sqrt(3.0f);
        result.partialDisplacement = bTwin;
        result.detectedLayerCount = layers;

        const int stepCount = std::max(1, (int)std::floor(params.maxDisplacementFactor / params.interval + 0.5f));
        int shiftedCountAtOne = 0;
        const int stageAtOne = std::max(0, (int)std::floor(1.0f + 1e-6f));
        const int startAtOne = std::min(layers - 1, startLayerBase + stageAtOne);
        for (int layer = startAtOne; layer < layers; ++layer)
            shiftedCountAtOne++;
        result.shiftedAtomCount = shiftedCountAtOne;

        result.sequence.resize((size_t)(stepCount + 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int step = 0; step <= stepCount; ++step)
        {
            const float factor = std::min(params.maxDisplacementFactor, step * params.interval);
            const glm::vec3 displacement = slipDirection * (result.partialDisplacement * factor);

            Structure generated = reference;

            const int stage = std::max(0, (int)std::floor(factor + 1e-6f));
            const int activeStart = std::min(layers - 1, startLayerBase + stage);
            for (int layer = activeStart; layer < layers; ++layer)
            {
                generated.atoms[(size_t)layer].x += displacement.x;
                generated.atoms[(size_t)layer].y += displacement.y;
                generated.atoms[(size_t)layer].z += displacement.z;
            }

            generated.cellVectors[2][0] += displacement.x;
            generated.cellVectors[2][1] += displacement.y;
            generated.cellVectors[2][2] += displacement.z;

            StackingFaultSequenceItem& item = result.sequence[(size_t)step];
            item.structure = std::move(generated);
            item.displacementFactor = factor;
            item.label = buildSequenceLabel(detection.family, factor, params.interval);
        }

        result.message = std::string("Generated ") + std::to_string(result.sequence.size())
            + " FCC-like stacking-fault structures on (111) <-1-12> (reference primitive mode).";
        result.success = true;
        return result;
    }

    Structure working = base;
    const bool useFccAlignedCell = detection.family == StackingFaultFamily::Fcc;
    if (useFccAlignedCell)
    {
        // Build an FCC-oriented working cell with one lattice axis along [111]
        // and an in-plane shear axis along [11-2].
        int fccTransform[3][3] = {
            { 1,  1, -2},
            { 1, -1,  0},
            { 1,  1,  1}
        };
        working = buildSupercell(base, fccTransform);
    }

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(working, cell, invCell))
    {
        result.message = "Input unit cell is singular.";
        return result;
    }

    glm::ivec3 planeHkl(0, 0, 1);
    glm::ivec3 directionUvw(1, 0, 0);
    std::string planeLabel;
    glm::vec3 planeNormal(0.0f, 0.0f, 1.0f);
    glm::vec3 rawDirection(1.0f, 0.0f, 0.0f);
    int normalVectorIndex = 0;

    if (useFccAlignedCell)
    {
        planeNormal = glm::normalize(cell[2]);
        rawDirection = cell[0];
        planeLabel = "(111) <-1-12>";
        normalVectorIndex = 2;
    }
    else
    {
        resolvePlaneAndDirection(detection.family, params.plane, planeHkl, directionUvw, planeLabel);
        refinePlaneAndDirection(working, cell, detection.family, params.plane, planeHkl, directionUvw);
        planeLabel = planeLabelForSelection(detection.family, params.plane);

        planeNormal = planeNormalFromHkl(cell, planeHkl);
        rawDirection = directionFromUvw(cell, directionUvw);
        if (glm::length(planeNormal) <= 1e-6f || glm::length(rawDirection) <= 1e-6f)
        {
            result.message = "Could not derive stacking-fault plane vectors from the current cell.";
            return result;
        }

        planeNormal = glm::normalize(planeNormal);
        float bestAlignment = 0.0f;
        for (int axis = 0; axis < 3; ++axis)
        {
            const glm::vec3 cellVec = glm::normalize(cell[axis]);
            const float alignment = std::abs(glm::dot(cellVec, planeNormal));
            if (alignment > bestAlignment)
            {
                bestAlignment = alignment;
                normalVectorIndex = axis;
            }
        }
    }

    glm::vec3 slipDirection = rawDirection - planeNormal * glm::dot(rawDirection, planeNormal);
    if (glm::length(slipDirection) <= 1e-6f)
    {
        result.message = "Chosen slip direction is not compatible with the detected plane.";
        return result;
    }
    slipDirection = glm::normalize(slipDirection);

    std::vector<float> normalProjection(working.atoms.size(), 0.0f);
    for (size_t i = 0; i < working.atoms.size(); ++i)
    {
        const glm::vec3 pos((float)working.atoms[i].x,
                            (float)working.atoms[i].y,
                            (float)working.atoms[i].z);
        normalProjection[i] = glm::dot(pos, planeNormal);
    }

    float layerSpacing = 0.0f;
    std::vector<int> layerIds = assignLayerIds(normalProjection, layerSpacing);
    int maxLayerId = *std::max_element(layerIds.begin(), layerIds.end());
    int detectedLayers = maxLayerId + 1;

    if (detectedLayers >= 1 && detectedLayers < params.layerCount)
    {
        const int repeats = std::max(2, (params.layerCount + detectedLayers - 1) / detectedLayers);
        int transform[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        transform[normalVectorIndex][normalVectorIndex] = repeats;
        working = buildSupercell(working, transform);

        if (!tryMakeCellMatrices(working, cell, invCell))
        {
            result.message = "Could not build a valid supercell for stacking-fault generation.";
            return result;
        }

        if (useFccAlignedCell)
        {
            planeNormal = glm::normalize(cell[2]);
            rawDirection = cell[0];
            normalVectorIndex = 2;
        }
        else
        {
            planeNormal = planeNormalFromHkl(cell, planeHkl);
            rawDirection = directionFromUvw(cell, directionUvw);
            if (glm::length(planeNormal) <= 1e-6f || glm::length(rawDirection) <= 1e-6f)
            {
                result.message = "Could not derive stacking-fault vectors after supercell expansion.";
                return result;
            }
            planeNormal = glm::normalize(planeNormal);

            normalVectorIndex = 0;
            float bestAlignment = 0.0f;
            for (int axis = 0; axis < 3; ++axis)
            {
                const glm::vec3 cellVec = glm::normalize(cell[axis]);
                const float alignment = std::abs(glm::dot(cellVec, planeNormal));
                if (alignment > bestAlignment)
                {
                    bestAlignment = alignment;
                    normalVectorIndex = axis;
                }
            }
        }

        slipDirection = rawDirection - planeNormal * glm::dot(rawDirection, planeNormal);
        if (glm::length(slipDirection) <= 1e-6f)
        {
            result.message = "Slip direction became invalid after supercell expansion.";
            return result;
        }
        slipDirection = glm::normalize(slipDirection);

        normalProjection.assign(working.atoms.size(), 0.0f);
        for (size_t i = 0; i < working.atoms.size(); ++i)
        {
            const glm::vec3 pos((float)working.atoms[i].x,
                                (float)working.atoms[i].y,
                                (float)working.atoms[i].z);
            normalProjection[i] = glm::dot(pos, planeNormal);
        }
        layerIds = assignLayerIds(normalProjection, layerSpacing);
        maxLayerId = *std::max_element(layerIds.begin(), layerIds.end());
        detectedLayers = maxLayerId + 1;
    }

    result.detectedLayerCount = detectedLayers;
    if (detectedLayers < 2)
    {
        result.message = "The selected structure does not contain enough distinct layers along the fault normal.";
        return result;
    }

    const int layerCountUsed = std::min(params.layerCount, detectedLayers);
    const int startLayer = std::max(0, (detectedLayers - layerCountUsed) / 2);
    const int endLayer = startLayer + layerCountUsed;

    std::vector<int> candidateAtomIndices;
    candidateAtomIndices.reserve(working.atoms.size());
    for (int i = 0; i < (int)layerIds.size(); ++i)
    {
        if (layerIds[i] >= startLayer && layerIds[i] < endLayer)
            candidateAtomIndices.push_back(i);
    }

    if (candidateAtomIndices.empty())
    {
        result.message = "No atoms fell into the shifted layer window. Increase the number of layers or load a thicker slab.";
        return result;
    }

    int referenceShiftedCount = 0;
    if (detection.family == StackingFaultFamily::Fcc)
    {
        const int fccBaseStart = std::max(startLayer, endLayer - 4);
        for (int atomIndex : candidateAtomIndices)
        {
            if (layerIds[atomIndex] >= fccBaseStart)
                referenceShiftedCount++;
        }
    }
    else
    {
        const int shiftStartLayer = startLayer + std::max(1, layerCountUsed / 2);
        for (int atomIndex : candidateAtomIndices)
        {
            if (layerIds[atomIndex] >= shiftStartLayer)
                referenceShiftedCount++;
        }
    }

    if (referenceShiftedCount <= 0)
    {
        result.message = "No atoms available in the selected shift region.";
        return result;
    }

    result.partialDisplacement = glm::length(rawDirection) * familyPartialScale(detection.family);
    result.shiftedAtomCount = referenceShiftedCount;

    const int stepCount = std::max(1, (int)std::floor(params.maxDisplacementFactor / params.interval + 0.5f));
    result.sequence.resize((size_t)(stepCount + 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int step = 0; step <= stepCount; ++step)
    {
        const float factor = std::min(params.maxDisplacementFactor, step * params.interval);
        const glm::vec3 displacement = slipDirection * (result.partialDisplacement * factor);

        int activeShiftStartLayer = startLayer + std::max(1, layerCountUsed / 2);
        if (detection.family == StackingFaultFamily::Fcc)
        {
            // Match the reference primitive script behavior:
            // start layer advances by one at factors 1.0 and 2.0.
            const int fccBaseStart = std::max(startLayer, endLayer - 4);
            const int stage = std::max(0, (int)std::floor(factor + 1e-6f));
            activeShiftStartLayer = std::min(endLayer - 1, fccBaseStart + stage);
        }

        Structure generated = working;
        for (int atomIndex : candidateAtomIndices)
        {
            if (layerIds[atomIndex] >= activeShiftStartLayer)
            {
                generated.atoms[atomIndex].x += displacement.x;
                generated.atoms[atomIndex].y += displacement.y;
                generated.atoms[atomIndex].z += displacement.z;
            }
        }

        generated.cellVectors[normalVectorIndex][0] += displacement.x;
        generated.cellVectors[normalVectorIndex][1] += displacement.y;
        generated.cellVectors[normalVectorIndex][2] += displacement.z;

        if (params.cellMode == StackingFaultCellMode::OrthogonalCell)
            generated = makeOrthogonalCellStructure(generated, slipDirection, planeNormal);

        StackingFaultSequenceItem& item = result.sequence[(size_t)step];
        item.structure = std::move(generated);
        item.displacementFactor = factor;
        item.label = buildSequenceLabel(detection.family, factor, params.interval);
    }

    std::ostringstream message;
    message << "Generated " << result.sequence.size() << " "
            << stackingFaultFamilyName(detection.family) << " stacking-fault structures on "
            << planeLabel << ". Shifted atoms: " << result.shiftedAtomCount
            << ", detected layers: " << detectedLayers << ".";

    result.message = message.str();
    result.success = true;
    return result;
}