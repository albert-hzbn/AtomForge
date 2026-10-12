#include "algorithms/InterstitialVoidAnalysis.h"

#include "math/StructureMath.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <array>
#include <cmath>
#include <cstdint>
#include <climits>
#include <limits>
#include <set>
#include <thread>
#include <utility>
#include <vector>

namespace
{
constexpr float kPlaneTol = 2e-3f;
constexpr float kDegenerateTol = 1e-8f;
constexpr float kEmptySphereTol = 1e-4f;

struct PlaneGroup
{
    glm::vec3 normal = glm::vec3(0.0f);
    float d = 0.0f;
    std::vector<int> localVertices;
};

struct NeighborShellPoint
{
    int atomIndex = -1;
    glm::vec3 position = glm::vec3(0.0f);
};

struct PeriodicPoint
{
    int atomIndex = -1;
    glm::vec3 position = glm::vec3(0.0f);
    glm::ivec3 shift = glm::ivec3(0);
};

struct LocalNeighbor
{
    int atomIndex = -1;
    glm::vec3 position = glm::vec3(0.0f);
    float d2 = 0.0f;
};

void pushUnique(std::vector<int>& out, int value)
{
    if (std::find(out.begin(), out.end(), value) == out.end())
        out.push_back(value);
}

glm::vec3 wrapPointToCell(const glm::vec3& p,
                          bool usePbc,
                          const glm::vec3& origin,
                          const glm::mat3& cell,
                          const glm::mat3& invCell)
{
    if (!usePbc)
        return p;

    glm::vec3 frac = invCell * (p - origin);
    frac -= glm::floor(frac);
    return origin + cell * frac;
}

bool computeCircumsphere(const glm::vec3& a,
                        const glm::vec3& b,
                        const glm::vec3& c,
                        const glm::vec3& d,
                        glm::vec3& center,
                        float& radius)
{
    const glm::vec3 ba = b - a;
    const glm::vec3 ca = c - a;
    const glm::vec3 da = d - a;

    glm::mat3 m;
    m[0] = 2.0f * ba;
    m[1] = 2.0f * ca;
    m[2] = 2.0f * da;

    const glm::mat3 mt = glm::transpose(m);
    const float det = glm::determinant(mt);
    if (std::abs(det) < 1e-9f)
        return false;

    const glm::vec3 rhs(glm::dot(b, b) - glm::dot(a, a),
                        glm::dot(c, c) - glm::dot(a, a),
                        glm::dot(d, d) - glm::dot(a, a));

    center = glm::inverse(mt) * rhs;
    radius = glm::length(center - a);
    return std::isfinite(radius);
}

std::vector<NeighborShellPoint> buildFirstNeighborShell(
    const glm::vec3& center,
    const std::vector<glm::vec3>& atomPos,
    bool usePbc,
    const glm::mat3& cell,
    const glm::mat3& invCell,
    float shellFactor,
    int maxNeighborCount)
{
    struct DistAtom
    {
        float d2 = 0.0f;
        int atomIndex = -1;
        glm::vec3 pos = glm::vec3(0.0f);
    };

    std::vector<DistAtom> near;
    near.reserve(atomPos.size());
    float minD2 = std::numeric_limits<float>::max();
    for (int ai = 0; ai < (int)atomPos.size(); ++ai)
    {
        const glm::vec3 delta = minimumImageDelta(atomPos[(size_t)ai] - center, usePbc, cell, invCell);
        const float d2 = glm::dot(delta, delta);
        minD2 = std::min(minD2, d2);
        near.push_back({d2, ai, center + delta});
    }

    std::vector<NeighborShellPoint> shell;
    if (near.empty())
        return shell;

    const float factor = std::max(1.02f, shellFactor);
    const float limit2 = std::max(0.0f, minD2) * factor * factor;

    // In a cell narrower than the shell, one atom can sit in the shell more
    // than once (the octahedral site of a single FCC cell has six neighbours
    // but three atoms): add the other images within reach.
    if (usePbc)
    {
        float width = std::numeric_limits<float>::max();
        for (int k = 0; k < 3; ++k)
        {
            const glm::vec3 row(invCell[0][k], invCell[1][k], invCell[2][k]);
            const float length = glm::length(row);
            if (length > 0.0f) width = std::min(width, 1.0f / length);
        }
        if (2.0f * std::sqrt(limit2) > width)
        {
            const std::size_t count = near.size();
            for (std::size_t n = 0; n < count; ++n)
                for (int a = -1; a <= 1; ++a)
                    for (int b = -1; b <= 1; ++b)
                        for (int c = -1; c <= 1; ++c)
                        {
                            if (a == 0 && b == 0 && c == 0) continue;
                            const glm::vec3 position = near[n].pos + cell * glm::vec3((float)a, (float)b, (float)c);
                            const glm::vec3 delta = position - center;
                            const float d2 = glm::dot(delta, delta);
                            if (d2 <= limit2) near.push_back({d2, near[n].atomIndex, position});
                        }
        }
    }

    std::vector<DistAtom> candidates;
    candidates.reserve(std::min((int)atomPos.size(), maxNeighborCount * 4));
    for (const DistAtom& da : near)
    {
        if (da.d2 <= limit2)
            candidates.push_back(da);
    }

    if ((int)candidates.size() > maxNeighborCount)
    {
        std::nth_element(candidates.begin(), candidates.begin() + maxNeighborCount, candidates.end(),
                         [](const DistAtom& a, const DistAtom& b) { return a.d2 < b.d2; });
        candidates.resize((size_t)maxNeighborCount);
        std::sort(candidates.begin(), candidates.end(),
                  [](const DistAtom& a, const DistAtom& b) { return a.d2 < b.d2; });
    }

    shell.reserve(candidates.size());
    for (const DistAtom& da : candidates)
        shell.push_back({da.atomIndex, da.pos});

    // Ensure nearest atom is always present in pathological tolerance cases.
    if (shell.empty())
    {
        auto it = std::min_element(near.begin(), near.end(),
                                   [](const DistAtom& a, const DistAtom& b) { return a.d2 < b.d2; });
        if (it != near.end())
            shell.push_back({it->atomIndex, it->pos});
    }

    return shell;
}

// Whether the origin lies in the convex hull of the unit vectors: the
// minimum-norm point of the hull (Frank-Wolfe iterations) is at the origin.
bool originInHull(const std::vector<glm::vec3>& directions)
{
    if (directions.size() < 2)
        return false;
    glm::vec3 x = directions[0];
    for (int iteration = 0; iteration < 500; ++iteration)
    {
        if (glm::dot(x, x) < 1e-6f)
            return true;
        std::size_t best = 0;
        for (std::size_t i = 1; i < directions.size(); ++i)
            if (glm::dot(x, directions[i]) < glm::dot(x, directions[best])) best = i;
        const glm::vec3 d = directions[best] - x;
        const float gap = -glm::dot(x, d);
        if (gap < 1e-7f)
            break;  // x is the minimum-norm point
        const float step = std::min(1.0f, gap / std::max(1e-12f, glm::dot(d, d)));
        x += step * d;
    }
    return glm::dot(x, x) < 1e-4f;
}

bool buildPolyhedralFacesFromFirstNeighbors(const glm::vec3& center,
                                            const std::vector<NeighborShellPoint>& shellPoints,
                                            std::vector<std::array<int, 3>>& outFaces,
                                            std::vector<std::array<glm::vec3, 3>>& outFaceVertices,
                                            float& outVolume)
{
    outFaces.clear();
    outFaceVertices.clear();
    outVolume = 0.0f;

    const int n = (int)shellPoints.size();
    if (n < 4)
        return false;

    std::vector<glm::vec3> shell;
    shell.reserve((size_t)n);
    for (const NeighborShellPoint& p : shellPoints)
        shell.push_back(p.position);

    std::vector<PlaneGroup> planes;
    planes.reserve((size_t)n * 2);

    for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
    for (int k = j + 1; k < n; ++k)
    {
        glm::vec3 a = shell[(size_t)i];
        glm::vec3 b = shell[(size_t)j];
        glm::vec3 c = shell[(size_t)k];

        glm::vec3 normal = glm::cross(b - a, c - a);
        const float len2 = glm::dot(normal, normal);
        if (len2 < kDegenerateTol)
            continue;

        if (glm::dot(normal, center - a) > 0.0f)
        {
            std::swap(b, c);
            normal = -normal;
        }

        bool supportsHull = true;
        for (int t = 0; t < n; ++t)
        {
            if (t == i || t == j || t == k)
                continue;
            const float side = glm::dot(normal, shell[(size_t)t] - a);
            if (side > kPlaneTol)
            {
                supportsHull = false;
                break;
            }
        }
        if (!supportsHull)
            continue;

        const glm::vec3 nrm = normal / std::sqrt(len2);
        const float d = glm::dot(nrm, a);

        int groupIdx = -1;
        for (int p = 0; p < (int)planes.size(); ++p)
        {
            const float parallel = glm::dot(planes[(size_t)p].normal, nrm);
            if (std::abs(parallel - 1.0f) > 1e-4f)
                continue;
            if (std::abs(planes[(size_t)p].d - d) > 1e-3f)
                continue;
            groupIdx = p;
            break;
        }

        if (groupIdx < 0)
        {
            PlaneGroup pg;
            pg.normal = nrm;
            pg.d = d;
            pg.localVertices = {i, j, k};
            planes.push_back(std::move(pg));
        }
        else
        {
            pushUnique(planes[(size_t)groupIdx].localVertices, i);
            pushUnique(planes[(size_t)groupIdx].localVertices, j);
            pushUnique(planes[(size_t)groupIdx].localVertices, k);
        }
    }

    if (planes.empty())
        return false;

    std::set<std::array<int, 3>> uniqueFaces;
    for (const PlaneGroup& plane : planes)
    {
        std::vector<int> verts = plane.localVertices;
        if (verts.size() < 3)
            continue;

        glm::vec3 centroid(0.0f);
        for (int li : verts)
            centroid += shell[(size_t)li];
        centroid /= (float)verts.size();

        glm::vec3 u = shell[(size_t)verts[0]] - centroid;
        if (glm::dot(u, u) < kDegenerateTol)
            u = glm::cross(plane.normal, glm::vec3(1.0f, 0.0f, 0.0f));
        if (glm::dot(u, u) < kDegenerateTol)
            u = glm::cross(plane.normal, glm::vec3(0.0f, 1.0f, 0.0f));
        if (glm::dot(u, u) < kDegenerateTol)
            continue;
        u = glm::normalize(u);

        glm::vec3 v = glm::normalize(glm::cross(plane.normal, u));

        std::sort(verts.begin(), verts.end(), [&](int a, int b) {
            const glm::vec3 ra = shell[(size_t)a] - centroid;
            const glm::vec3 rb = shell[(size_t)b] - centroid;
            const float aa = std::atan2(glm::dot(ra, v), glm::dot(ra, u));
            const float ab = std::atan2(glm::dot(rb, v), glm::dot(rb, u));
            return aa < ab;
        });

        const int first = verts[0];
        for (size_t t = 1; t + 1 < verts.size(); ++t)
        {
            int li0 = first;
            int li1 = verts[t];
            int li2 = verts[t + 1];

            const glm::vec3& a = shell[(size_t)li0];
            const glm::vec3& b = shell[(size_t)li1];
            const glm::vec3& c = shell[(size_t)li2];
            if (glm::dot(glm::cross(b - a, c - a), plane.normal) < 0.0f)
                std::swap(li1, li2);

            const int gi0 = shellPoints[(size_t)li0].atomIndex;
            const int gi1 = shellPoints[(size_t)li1].atomIndex;
            const int gi2 = shellPoints[(size_t)li2].atomIndex;

            std::array<int, 3> key = {gi0, gi1, gi2};
            // Unique by shell vertex, not atom: in a small cell two vertices
            // can be images of one atom.
            std::array<int, 3> sortedKey = {li0, li1, li2};
            std::sort(sortedKey.begin(), sortedKey.end());
            if (uniqueFaces.insert(sortedKey).second)
            {
                outFaces.push_back(key);
                outFaceVertices.push_back({a, b, c});
            }
        }
    }

    if (outFaces.empty())
        return false;

    // From the face vertices themselves: an atom can appear in the shell
    // more than once (periodic images), so its index does not name a vertex.
    double volume = 0.0;
    for (const auto& tri : outFaceVertices)
    {
        const glm::vec3 a = tri[0] - center;
        const glm::vec3 b = tri[1] - center;
        const glm::vec3 c = tri[2] - center;
        // The centre is inside the convex shell: every face adds a positive cone.
        volume += std::abs((double)glm::dot(a, glm::cross(b, c))) / 6.0;
    }
    outVolume = (float)std::abs(volume);
    return outVolume > 1e-6f;
}

void assignVolumeBins(std::vector<InterstitialVoidRegion>& regions)
{
    if (regions.empty())
        return;

    std::vector<int> order(regions.size());
    for (int i = 0; i < (int)regions.size(); ++i)
        order[(size_t)i] = i;

    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return regions[(size_t)a].volume < regions[(size_t)b].volume;
    });

    float representative = -1.0f;
    int currentBin = -1;
    for (int idx : order)
    {
        const float v = regions[(size_t)idx].volume;
        if (currentBin < 0)
        {
            currentBin = 0;
            representative = v;
        }
        else
        {
            const float tol = std::max(1e-4f, representative * 0.05f);
            if (std::abs(v - representative) > tol)
            {
                ++currentBin;
                representative = v;
            }
            else
            {
                representative = 0.5f * (representative + v);
            }
        }

        regions[(size_t)idx].volumeBin = currentBin;
    }
}

} // namespace

InterstitialVoidDetectionResult detectInterstitialVoidRegions(
    const Structure& structure,
    const InterstitialVoidDetectionParams& params)
{
    InterstitialVoidDetectionResult result;

    if (structure.atoms.empty())
    {
        result.message = "Source structure has no atoms.";
        return result;
    }

    const int gridN = std::max(6, std::min(48, params.gridResolution));
    const bool limitVoids = params.maxVoids > 0;
    const int maxVoids = limitVoids ? std::max(1, params.maxVoids) : INT_MAX;
    const float minClear = std::max(0.05f, params.minClearance);
    const float minSep = std::max(0.05f, params.minSeparation);
    const float neighborShell = std::max(1.02f, params.firstNeighborShellFactor);
    const int maxNeighbor = std::max(4, std::min(24, params.maxNeighborCount));

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    const bool usePbc = tryMakeCellMatrices(structure, cell, invCell);
    const glm::vec3 origin((float)structure.cellOffset[0],
                           (float)structure.cellOffset[1],
                           (float)structure.cellOffset[2]);

    std::vector<glm::vec3> atomPos;
    atomPos.reserve(structure.atoms.size());
    for (const AtomSite& a : structure.atoms)
        atomPos.push_back(glm::vec3((float)a.x, (float)a.y, (float)a.z));

    const int neighborPool = std::max(16, std::min(40, gridN + 10));

    // A cell narrower than about two pool radii (some 2.2 atomic spacings
    // each) holds neighbours whose other images are as near as the minimum
    // image: list the 26 adjacent images of every atom too.
    std::vector<glm::vec3> imageShifts{glm::vec3(0.0f)};
    if (usePbc)
    {
        float width = std::numeric_limits<float>::max();
        for (int k = 0; k < 3; ++k)
        {
            const glm::vec3 row(invCell[0][k], invCell[1][k], invCell[2][k]);
            if (glm::length(row) > 0.0f) width = std::min(width, 1.0f / glm::length(row));
        }
        const float spacing = std::cbrt(std::abs(glm::determinant(cell)) / (float)std::max<std::size_t>(1, atomPos.size()));
        if (width < 4.4f * spacing)
            for (int a = -1; a <= 1; ++a)
                for (int b = -1; b <= 1; ++b)
                    for (int c = -1; c <= 1; ++c)
                        if (a != 0 || b != 0 || c != 0) imageShifts.push_back(cell * glm::vec3((float)a, (float)b, (float)c));
    }

    std::vector<InterstitialVoidRegion> regions;
    if (limitVoids)
        regions.reserve((size_t)maxVoids * 2);
    else
        regions.reserve(std::max((size_t)4096, atomPos.size() * (size_t)2));

    struct CandidateCenter
    {
        glm::vec3 center = glm::vec3(0.0f);
        float radius = 0.0f;
        int a = -1;
        int b = -1;
        int c = -1;
        int d = -1;
    };

    std::vector<CandidateCenter> centers;
    centers.reserve(atomPos.size() * 6);

    const int atomCount = (int)atomPos.size();
    const unsigned int hw = std::thread::hardware_concurrency();
    const int workerCount = std::max(1, std::min(atomCount, (int)((hw == 0u) ? 8u : hw)));
    std::vector<std::vector<CandidateCenter>> centersPerWorker((size_t)workerCount);
    std::vector<std::thread> workers;
    workers.reserve((size_t)workerCount);

    for (int wi = 0; wi < workerCount; ++wi)
    {
        workers.emplace_back([&, wi]()
        {
            std::vector<CandidateCenter>& localCenters = centersPerWorker[(size_t)wi];
            localCenters.reserve((size_t)std::max(64, atomCount / workerCount * 4));

            for (int ai = wi; ai < atomCount; ai += workerCount)
            {
                struct Near
                {
                    float first = 0.0f;  // squared distance
                    int second = -1;     // atom
                    glm::vec3 delta = glm::vec3(0.0f);
                };
                std::vector<Near> near;
                near.reserve(atomPos.size() * imageShifts.size());
                for (int aj = 0; aj < atomCount; ++aj)
                {
                    const glm::vec3 base = minimumImageDelta(atomPos[(size_t)aj] - atomPos[(size_t)ai], usePbc, cell, invCell);
                    for (const glm::vec3& shift : imageShifts)
                    {
                        const glm::vec3 d = base + shift;
                        const float d2 = glm::dot(d, d);
                        if (d2 < 1e-8f)
                            continue;  // the atom itself
                        near.push_back({d2, aj, d});
                    }
                }
                if ((int)near.size() < 4)
                    continue;

                const int poolN = std::min(neighborPool, (int)near.size());
                if (poolN < 4)
                    continue;

                if (poolN < (int)near.size())
                {
                    std::nth_element(near.begin(), near.begin() + poolN, near.end(),
                                     [](const auto& a, const auto& b) { return a.first < b.first; });
                    near.resize((size_t)poolN);
                }
                std::sort(near.begin(), near.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

                std::vector<LocalNeighbor> pool;
                pool.reserve((size_t)poolN);
                for (int i = 0; i < poolN; ++i)
                    pool.push_back({near[(size_t)i].second, atomPos[(size_t)ai] + near[(size_t)i].delta, near[(size_t)i].first});

                const float d1 = std::sqrt(std::max(0.0f, near[0].first));
                const float candidateSpan2 = (d1 * 1.9f) * (d1 * 1.9f);
                std::vector<int> tetraPts;
                tetraPts.reserve(pool.size());
                for (int pi = 0; pi < (int)pool.size(); ++pi)
                {
                    if (pool[(size_t)pi].d2 <= candidateSpan2)
                        tetraPts.push_back(pi);
                }
                if ((int)tetraPts.size() < 3)
                    continue;

                for (int ib = 0; ib < (int)tetraPts.size(); ++ib)
                for (int ic = ib + 1; ic < (int)tetraPts.size(); ++ic)
                for (int id = ic + 1; id < (int)tetraPts.size(); ++id)
                {
                    const LocalNeighbor& pb = pool[(size_t)tetraPts[(size_t)ib]];
                    const LocalNeighbor& pc = pool[(size_t)tetraPts[(size_t)ic]];
                    const LocalNeighbor& pd = pool[(size_t)tetraPts[(size_t)id]];

                    glm::vec3 cc(0.0f);
                    float rr = 0.0f;
                    if (!computeCircumsphere(atomPos[(size_t)ai], pb.position, pc.position, pd.position, cc, rr))
                        continue;
                    if (!(rr > 0.0f) || !std::isfinite(rr))
                        continue;

                    bool empty = true;
                    const float rrTol = rr - kEmptySphereTol;
                    const float rrTol2 = rrTol * rrTol;
                    for (const LocalNeighbor& p : pool)
                    {
                        if (&p == &pb || &p == &pc || &p == &pd)
                            continue;

                        const glm::vec3 d = p.position - cc;
                        if (glm::dot(d, d) < rrTol2)
                        {
                            empty = false;
                            break;
                        }
                    }
                    if (!empty)
                        continue;

                    glm::vec3 centerWrapped = wrapPointToCell(cc, usePbc, origin, cell, invCell);
                    localCenters.push_back({centerWrapped, rr, ai, pb.atomIndex, pc.atomIndex, pd.atomIndex});
                }

                // Midpoints of nearby atom pairs recover edge/bridge-centered void candidates
                // (important for BCC octahedral families under periodic boundaries).
                const int pairN = std::min(poolN, 14);
                for (int pi = 0; pi < pairN; ++pi)
                {
                    const LocalNeighbor& pb = pool[(size_t)pi];
                    const glm::vec3 mid = atomPos[(size_t)ai] + 0.5f * (pb.position - atomPos[(size_t)ai]);
                    const glm::vec3 wrappedMid = wrapPointToCell(mid, usePbc, origin, cell, invCell);

                    // Cheap radius estimate here; exact clearance is recomputed later.
                    const float rr = std::sqrt(std::max(0.0f, pb.d2)) * 0.5f;
                    if (rr >= minClear)
                        localCenters.push_back({wrappedMid, rr, ai, pb.atomIndex, -1, -1});
                }
            }
        });
    }

    for (std::thread& t : workers)
        t.join();

    for (std::vector<CandidateCenter>& local : centersPerWorker)
        centers.insert(centers.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));

    // Supplemental periodic clearance peaks improve octahedral-site recovery
    // in lattices where pure tetra circumsphere candidates under-sample maxima.
    glm::vec3 bbMin( std::numeric_limits<float>::max());
    glm::vec3 bbMax(-std::numeric_limits<float>::max());
    for (const glm::vec3& p : atomPos)
    {
        bbMin = glm::min(bbMin, p);
        bbMax = glm::max(bbMax, p);
    }

    struct GridSample
    {
        glm::vec3 p = glm::vec3(0.0f);
        float clearance = 0.0f;
    };

    std::vector<GridSample> samples((size_t)gridN * (size_t)gridN * (size_t)gridN);
    auto sampleIdx = [gridN](int i, int j, int k) {
        return (size_t)i * (size_t)gridN * (size_t)gridN + (size_t)j * (size_t)gridN + (size_t)k;
    };
    auto wrapGrid = [gridN](int v) {
        int r = v % gridN;
        if (r < 0)
            r += gridN;
        return r;
    };

    const int totalSamples = gridN * gridN * gridN;
    std::vector<std::thread> sampleWorkers;
    sampleWorkers.reserve((size_t)workerCount);
    for (int wi = 0; wi < workerCount; ++wi)
    {
        sampleWorkers.emplace_back([&, wi]()
        {
            for (int idx = wi; idx < totalSamples; idx += workerCount)
            {
                const int i = idx / (gridN * gridN);
                const int rem = idx % (gridN * gridN);
                const int j = rem / gridN;
                const int k = rem % gridN;

                const glm::vec3 frac(((float)i + 0.5f) / (float)gridN,
                                     ((float)j + 0.5f) / (float)gridN,
                                     ((float)k + 0.5f) / (float)gridN);

                glm::vec3 p(0.0f);
                if (usePbc)
                    p = origin + cell * frac;
                else
                    p = bbMin + frac * (bbMax - bbMin);

                float nearest2 = std::numeric_limits<float>::max();
                for (const glm::vec3& a : atomPos)
                {
                    const glm::vec3 d = minimumImageDelta(p - a, usePbc, cell, invCell);
                    nearest2 = std::min(nearest2, glm::dot(d, d));
                }

                samples[sampleIdx(i, j, k)] = {p, std::sqrt(std::max(0.0f, nearest2))};
            }
        });
    }
    for (std::thread& t : sampleWorkers)
        t.join();

    const int iBeg = usePbc ? 0 : 1;
    const int jBeg = usePbc ? 0 : 1;
    const int kBeg = usePbc ? 0 : 1;
    const int iEnd = usePbc ? gridN : (gridN - 1);
    const int jEnd = usePbc ? gridN : (gridN - 1);
    const int kEnd = usePbc ? gridN : (gridN - 1);

    for (int i = iBeg; i < iEnd; ++i)
    for (int j = jBeg; j < jEnd; ++j)
    for (int k = kBeg; k < kEnd; ++k)
    {
        const GridSample& s = samples[sampleIdx(i, j, k)];
        if (s.clearance < minClear)
            continue;

        bool isPeak = true;
        for (int di = -1; di <= 1 && isPeak; ++di)
        for (int dj = -1; dj <= 1 && isPeak; ++dj)
        for (int dk = -1; dk <= 1 && isPeak; ++dk)
        {
            if (di == 0 && dj == 0 && dk == 0)
                continue;

            const int ni = usePbc ? wrapGrid(i + di) : (i + di);
            const int nj = usePbc ? wrapGrid(j + dj) : (j + dj);
            const int nk = usePbc ? wrapGrid(k + dk) : (k + dk);
            if (!usePbc)
            {
                if (ni <= 0 || nj <= 0 || nk <= 0 || ni >= gridN - 1 || nj >= gridN - 1 || nk >= gridN - 1)
                    continue;
            }

            if (samples[sampleIdx(ni, nj, nk)].clearance > s.clearance)
                isPeak = false;
        }
        if (!isPeak)
            continue;

        centers.push_back({wrapPointToCell(s.p, usePbc, origin, cell, invCell), s.clearance, -1, -1, -1, -1});
    }

    const size_t candidateBudget = std::max((size_t)20000,
                                            std::max((size_t)maxVoids * 500,
                                                     (size_t)atomCount * 120));
    if (centers.size() > candidateBudget)
    {
        std::nth_element(centers.begin(), centers.begin() + (int)candidateBudget, centers.end(),
                         [](const CandidateCenter& a, const CandidateCenter& b) {
                             return a.radius > b.radius;
                         });
        centers.resize(candidateBudget);
    }

    std::sort(centers.begin(), centers.end(), [](const CandidateCenter& a, const CandidateCenter& b) {
        return a.radius > b.radius;
    });

    const auto clearanceAt = [&](const glm::vec3& point) {
        float nearest2 = std::numeric_limits<float>::max();
        for (const glm::vec3& a : atomPos)
        {
            const glm::vec3 d = minimumImageDelta(point - a, usePbc, cell, invCell);
            nearest2 = std::min(nearest2, glm::dot(d, d));
        }
        return std::sqrt(std::max(0.0f, nearest2));
    };
    const auto shellAt = [&](const glm::vec3& point) {
        std::vector<NeighborShellPoint> shell = buildFirstNeighborShell(
            point, atomPos, usePbc, cell, invCell, neighborShell, maxNeighbor);
        if ((int)shell.size() < 4)
            shell = buildFirstNeighborShell(
                point, atomPos, usePbc, cell, invCell,
                std::max(1.5f, neighborShell), std::max(8, maxNeighbor));
        return shell;
    };

    // Accepted sites by bin of size minSep (fractional bins in a periodic
    // cell), so "is a site already within minSep" looks at 27 bins rather than
    // every site: most candidates repeat a site already found.
    int binCounts[3] = {1, 1, 1};
    if (usePbc)
        for (int k = 0; k < 3; ++k)
        {
            const glm::vec3 row(invCell[0][k], invCell[1][k], invCell[2][k]);
            const float width = glm::length(row) > 0.0f ? 1.0f / glm::length(row) : minSep;
            binCounts[k] = std::max(1, std::min(1024, (int)std::floor(width / minSep)));
        }
    const auto binOf = [&](const glm::vec3& p) {
        glm::ivec3 bin;
        if (usePbc)
        {
            const glm::vec3 f = invCell * (p - origin);
            for (int k = 0; k < 3; ++k)
            {
                int b = (int)std::floor((f[k] - std::floor(f[k])) * (float)binCounts[k]);
                bin[k] = std::min(std::max(b, 0), binCounts[k] - 1);
            }
        }
        else
            bin = glm::ivec3(glm::floor(p / minSep));
        return bin;
    };
    const auto binKey = [](const glm::ivec3& b) {
        return ((std::uint64_t)(std::uint32_t)(b.x + 1048576) << 42) ^ ((std::uint64_t)(std::uint32_t)(b.y + 1048576) << 21) ^
               (std::uint64_t)(std::uint32_t)(b.z + 1048576);
    };
    std::unordered_map<std::uint64_t, std::vector<int>> siteBins;
    const auto nearSite = [&](const glm::vec3& p) {
        const glm::ivec3 centre = binOf(p);
        // In a periodic direction with fewer than three bins every bin is adjacent.
        int reach[3];
        for (int k = 0; k < 3; ++k) reach[k] = usePbc && binCounts[k] < 3 ? 0 : 1;
        for (int dx = -reach[0]; dx <= reach[0]; ++dx)
            for (int dy = -reach[1]; dy <= reach[1]; ++dy)
                for (int dz = -reach[2]; dz <= reach[2]; ++dz)
                {
                    glm::ivec3 b = centre + glm::ivec3(dx, dy, dz);
                    if (usePbc)
                        for (int k = 0; k < 3; ++k) b[k] = ((b[k] % binCounts[k]) + binCounts[k]) % binCounts[k];
                    const auto found = siteBins.find(binKey(b));
                    if (found == siteBins.end()) continue;
                    for (int index : found->second)
                    {
                        const glm::vec3 d = minimumImageDelta(p - regions[(size_t)index].position, usePbc, cell, invCell);
                        if (glm::dot(d, d) < minSep * minSep) return true;
                    }
                }
        if (usePbc && (binCounts[0] < 3 || binCounts[1] < 3 || binCounts[2] < 3))
        {
            // Small cells: bins wrap onto themselves; check every site.
            for (const InterstitialVoidRegion& c : regions)
            {
                const glm::vec3 d = minimumImageDelta(p - c.position, usePbc, cell, invCell);
                if (glm::dot(d, d) < minSep * minSep) return true;
            }
        }
        return false;
    };

    for (const CandidateCenter& cand : centers)
    {
        glm::vec3 center = cand.center;
        if (nearSite(center))
            continue;
        float clearance = clearanceAt(center);
        if (clearance < minClear)
            continue;

        std::vector<NeighborShellPoint> shell = shellAt(center);
        if ((int)shell.size() < 4)
            continue;

        // Grid and bond-midpoint candidates lie near, not at, the site: move to
        // the centre of the surrounding atoms when that is roomier (the
        // octahedral and tetrahedral sites of a crystal are such centres).
        // Repeat while it helps: the new place can have a fuller shell.
        for (int step = 0; step < 8; ++step)
        {
            glm::vec3 centroid(0.0f);
            for (const NeighborShellPoint& sp : shell) centroid += sp.position;
            centroid = wrapPointToCell(centroid / (float)shell.size(), usePbc, origin, cell, invCell);
            const float refined = clearanceAt(centroid);
            if (refined > clearance + 1e-4f)
            {
                std::vector<NeighborShellPoint> refinedShell = shellAt(centroid);
                if ((int)refinedShell.size() < 4)
                    break;
                clearance = refined;
                shell = std::move(refinedShell);
                center = centroid;
            }
            else break;
        }
        // A clearance peak on the sampling grid approximates a largest empty
        // sphere, whose centre is equidistant from its four nearest atoms.
        if (cand.a < 0 && shell.size() >= 4)
        {
            std::vector<NeighborShellPoint> nearest = shell;
            std::sort(nearest.begin(), nearest.end(), [&](const NeighborShellPoint& x, const NeighborShellPoint& y) {
                return glm::dot(x.position - center, x.position - center) < glm::dot(y.position - center, y.position - center);
            });
            glm::vec3 cc(0.0f);
            float rr = 0.0f;
            if (computeCircumsphere(nearest[0].position, nearest[1].position, nearest[2].position, nearest[3].position, cc, rr) &&
                glm::length(cc - center) < clearance)
            {
                const glm::vec3 wrapped = wrapPointToCell(cc, usePbc, origin, cell, invCell);
                const float refined = clearanceAt(wrapped);
                if (refined > clearance + 1e-4f && std::abs(refined - rr) < 1e-3f)
                {
                    std::vector<NeighborShellPoint> refinedShell = shellAt(wrapped);
                    if ((int)refinedShell.size() >= 4)
                    {
                        clearance = refined;
                        shell = std::move(refinedShell);
                        center = wrapped;
                    }
                }
            }
        }

        // A site is a critical point of the clearance: the atoms it touches
        // surround it (a largest empty sphere, or a bridge between two
        // opposite atoms such as the BCC octahedral site). Elsewhere the
        // sphere could still grow by moving, so the point is no site.
        {
            std::vector<glm::vec3> touching;
            for (const NeighborShellPoint& sp : shell)
            {
                const glm::vec3 d = sp.position - center;
                const float r = glm::length(d);
                if (r > 0.0f && r <= clearance * 1.002f + 2e-3f) touching.push_back(d / r);
            }
            if (!originInHull(touching))
                continue;
        }

        std::vector<std::array<int, 3>> faces;
        std::vector<std::array<glm::vec3, 3>> faceVertices;
        float volume = 0.0f;
        if (!buildPolyhedralFacesFromFirstNeighbors(center, shell, faces, faceVertices, volume))
            continue;

        if (volume <= 1e-6f)
            continue;

        // Reject obvious oversized false cavities: compare against nearest-shell characteristic scale.
        float r1 = std::numeric_limits<float>::max();
        for (const NeighborShellPoint& sp : shell)
        {
            const float d = glm::length(sp.position - center);
            r1 = std::min(r1, d);
        }
        const float maxReasonableVolume = (4.0f / 3.0f) * 3.14159265f * std::pow(std::max(0.0f, r1 * 1.5f), 3.0f);
        if (volume > maxReasonableVolume)
            continue;

        InterstitialVoidRegion region;
        region.position = center;
        region.clearance = clearance;
        region.coordination = (int)shell.size();
        region.volume = volume;
        for (const NeighborShellPoint& sp : shell)
            region.neighborAtomIndices.push_back(sp.atomIndex);
        region.faces = std::move(faces);
        region.faceVertices = std::move(faceVertices);

        if (region.coordination <= 4)
            region.kind = InterstitialVoidKind::Tetrahedral;
        else if (region.coordination <= 6)
            region.kind = InterstitialVoidKind::Octahedral;
        else
            region.kind = InterstitialVoidKind::Irregular;

        if (!nearSite(region.position))
        {
            siteBins[binKey(binOf(region.position))].push_back((int)regions.size());
            regions.push_back(std::move(region));
        }
    }

    std::sort(regions.begin(), regions.end(), [](const InterstitialVoidRegion& a, const InterstitialVoidRegion& b) {
        return a.clearance > b.clearance;
    });

    if (limitVoids && (int)regions.size() > maxVoids)
        regions.resize((size_t)maxVoids);

    assignVolumeBins(regions);

    result.success = !regions.empty();
    if (result.success)
        result.message = "Detected " + std::to_string(regions.size()) + " void regions.";
    else
        result.message = "No void regions detected.";
    result.regions = std::move(regions);
    return result;
}
