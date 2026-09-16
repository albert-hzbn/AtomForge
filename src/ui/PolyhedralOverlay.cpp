#include "ui/PolyhedralOverlay.h"

#include "math/StructureMath.h"
#include "util/ElementData.h"

#include "imgui.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <queue>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace
{
constexpr float kNearClip = 0.03f;
constexpr int kMinNeighborCount = 4;
constexpr float kFacetEpsilon = 1e-4f;

bool projectToScreen(const glm::vec3& p,
                     const glm::mat4& projection,
                     const glm::mat4& view,
                     int w,
                     int h,
                     ImVec2& out)
{
    const glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
    if (clip.w <= 1e-5f)
        return false;

    const float invW = 1.0f / clip.w;
    const float nx = clip.x * invW;
    const float ny = clip.y * invW;

    out.x = (nx * 0.5f + 0.5f) * (float)w;
    out.y = (1.0f - (ny * 0.5f + 0.5f)) * (float)h;
    return true;
}

std::vector<glm::vec3> clipPolygonAgainstNearPlane(const std::vector<glm::vec3>& worldPoints,
                                                   const glm::mat4& view,
                                                   float nearDepth)
{
    if (worldPoints.size() < 3)
        return {};

    std::vector<glm::vec3> clipped;
    clipped.reserve(worldPoints.size() + 2);

    for (size_t i = 0; i < worldPoints.size(); ++i)
    {
        const glm::vec3& currentWorld = worldPoints[i];
        const glm::vec3& nextWorld = worldPoints[(i + 1) % worldPoints.size()];

        const glm::vec4 currentView = view * glm::vec4(currentWorld, 1.0f);
        const glm::vec4 nextView = view * glm::vec4(nextWorld, 1.0f);
        const bool currentInside = currentView.z <= -nearDepth;
        const bool nextInside = nextView.z <= -nearDepth;

        if (currentInside)
            clipped.push_back(currentWorld);

        if (currentInside != nextInside)
        {
            const float denom = nextView.z - currentView.z;
            if (std::abs(denom) > 1e-6f)
            {
                float t = (-nearDepth - currentView.z) / denom;
                t = std::max(0.0f, std::min(1.0f, t));
                clipped.push_back(currentWorld + (nextWorld - currentWorld) * t);
            }
        }
    }

    return (clipped.size() >= 3) ? clipped : std::vector<glm::vec3>();
}

bool isOutsideUnitCell(const Structure& structure, const glm::vec3& position)
{
    glm::vec3 fractional(0.0f);
    if (!tryCartesianToFractional(structure, position, fractional))
        return false;

    constexpr float kFracTol = 1e-4f;
    return (fractional.x < -kFracTol || fractional.x > 1.0f + kFracTol
         || fractional.y < -kFracTol || fractional.y > 1.0f + kFracTol
         || fractional.z < -kFracTol || fractional.z > 1.0f + kFracTol);
}

std::vector<int> gatherCenters(const Structure& structure,
                               const std::vector<int>& selectedInstanceIndices,
                               const std::vector<int>& instanceToAtomIndex,
                               const PolyhedralOverlaySettings& settings)
{
    (void)selectedInstanceIndices;
    (void)instanceToAtomIndex;

    // Elements that are almost always ligands / anions, never coordination centres.
    // When no explicit element filter is set we exclude these by default so that
    // the overlay "just works" for perovskites, oxides, halides, etc.
    static const int kDefaultExclude[] = {
        1,  // H
        2,  // He
        7,  // N
        8,  // O
        9,  // F
        10, // Ne
        16, // S
        17, // Cl
        18, // Ar
        34, // Se
        35, // Br
        36, // Kr
        52, // Te
        53, // I
        54, // Xe
        85, // At
        86, // Rn
    };
    auto isDefaultExcluded = [](int z) {
        for (int ex : kDefaultExclude)
            if (ex == z) return true;
        return false;
    };

    std::vector<int> centers;
    const int centerCap = std::max(1, settings.maxDisplayedCenters);

    auto centerAllowed = [&](int atomIdx) -> bool {
        if (atomIdx < 0 || atomIdx >= (int)structure.atoms.size())
            return false;
        const int z = structure.atoms[(size_t)atomIdx].atomicNumber;
        if (z < 1 || z > 118)
            return false;
        if (settings.centerElementFilterEnabled)
            return settings.centerElementMask[(size_t)z];
        // Default: exclude typical anion / ligand species.
        return !isDefaultExcluded(z);
    };

    if (settings.centerAtomIndexFilterEnabled && !settings.centerAtomIndices.empty())
    {
        centers.reserve(settings.centerAtomIndices.size());
        for (const int atomIdx : settings.centerAtomIndices)
        {
            if (centerAllowed(atomIdx))
                centers.push_back(atomIdx);
        }
        return centers;
    }

    centers.reserve((size_t)centerCap);
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        if (!centerAllowed(i))
            continue;
        centers.push_back(i);
        if ((int)centers.size() >= centerCap)
            break;
    }

    return centers;
}

// Returned by buildNeighborCloud.
struct NeighborCloud
{
    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> colors;      // per-position source-atom display color (r,g,b)
    int ligandAtomicNumber = 0;
};

NeighborCloud buildNeighborCloud(const Structure& structure,
                                 const glm::vec3& center,
                                 int centerAtomicNumber,
                                 bool usePbc,
                                 const glm::mat3& cell,
                                 const glm::mat3& invCell,
                                 const PolyhedralOverlaySettings& settings,
                                 const std::vector<glm::vec3>& elementColors)
{
    (void)invCell;

    struct Candidate
    {
        float dist = 0.0f;
        glm::vec3 pos = glm::vec3(0.0f);
        int atomicNumber = 0;
        glm::vec3 color = glm::vec3(0.8f); // source atom r,g,b
    };

    const int maxNeighborCandidates = std::max(kMinNeighborCount, settings.maxNeighborCandidates);

    const auto nearer=[](const Candidate& a,const Candidate& b) { return a.dist<b.dist; };
    using Nearest=std::priority_queue<Candidate,std::vector<Candidate>,decltype(nearer)>;
    std::map<int,Nearest> byElement;
    const auto retain=[&](const Candidate& candidate) {
        auto found=byElement.find(candidate.atomicNumber);
        if (found==byElement.end()) found=byElement.emplace(candidate.atomicNumber,Nearest(nearer)).first;
        auto& nearest=found->second;
        if (nearest.size()<static_cast<std::size_t>(maxNeighborCandidates)) nearest.push(candidate);
        else if (candidate.dist<nearest.top().dist) { nearest.pop(); nearest.push(candidate); }
    };
    std::vector<Candidate> candidates;

    const glm::vec3 cellA = usePbc ? cell[0] : glm::vec3(0.0f);
    const glm::vec3 cellB = usePbc ? cell[1] : glm::vec3(0.0f);
    const glm::vec3 cellC = usePbc ? cell[2] : glm::vec3(0.0f);

    for (int j = 0; j < (int)structure.atoms.size(); ++j)
    {
        const auto& atom = structure.atoms[(size_t)j];
        if (atom.atomicNumber == centerAtomicNumber)
            continue;

        if (settings.ligandElementFilterEnabled)
        {
            const int z = atom.atomicNumber;
            if (z < 1 || z > 118 || !settings.ligandElementMask[(size_t)z])
                continue;
        }

        const glm::vec3 atomPos((float)atom.x, (float)atom.y, (float)atom.z);

        if (!usePbc)
        {
            const glm::vec3 delta = atomPos - center;
            const float dist = glm::length(delta);
            if (dist <= 1e-5f)
                continue;

            Candidate candidate;
            candidate.dist = dist;
            candidate.pos = atomPos;
            candidate.atomicNumber = atom.atomicNumber;
            candidate.color = (atom.atomicNumber >= 1 && atom.atomicNumber < (int)elementColors.size())
                ? elementColors[(size_t)atom.atomicNumber]
                : glm::vec3(atom.r, atom.g, atom.b);
            retain(candidate);
            continue;
        }

        // Include neighboring periodic images so centers near boundaries can
        // still form complete ligand cages (e.g., O/N octahedra in perovskites).
        for (int sx = -1; sx <= 1; ++sx)
        {
            for (int sy = -1; sy <= 1; ++sy)
            {
                for (int sz = -1; sz <= 1; ++sz)
                {
                    const glm::vec3 shiftedPos = atomPos
                        + (float)sx * cellA
                        + (float)sy * cellB
                        + (float)sz * cellC;

                    const glm::vec3 delta = shiftedPos - center;
                    const float dist = glm::length(delta);
                    if (dist <= 1e-5f)
                        continue;

                    Candidate candidate;
                    candidate.dist = dist;
                    candidate.pos = shiftedPos;
                    candidate.atomicNumber = atom.atomicNumber;
                    candidate.color = (atom.atomicNumber >= 1 && atom.atomicNumber < (int)elementColors.size())
                        ? elementColors[(size_t)atom.atomicNumber]
                        : glm::vec3(atom.r, atom.g, atom.b);
                    retain(candidate);
                }
            }
        }
    }

    for (auto& entry:byElement) while (!entry.second.empty()) {
        candidates.push_back(entry.second.top()); entry.second.pop();
    }
    if (candidates.size() < (size_t)kMinNeighborCount)
        return {};

    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& lhs, const Candidate& rhs) { return lhs.dist < rhs.dist; });

    // Use the dominant ligand species (most common among nearest neighbours).
    // In a perovskite the 6 nearest ligands to a B-site are all O, so this
    // naturally picks O even if a second ligand type appears further out.
    const int ligandType = candidates.front().atomicNumber;

    // First-shell cutoff: include only ligands within kFirstShellFactor of the
    // nearest ligand distance.  This eliminates second-shell atoms that would
    // otherwise produce long diagonal faces across the polyhedron.
    constexpr float kFirstShellFactor = 1.20f;
    const float nearestDist = candidates.front().dist;
    const float shellCutoff = nearestDist * kFirstShellFactor;

    NeighborCloud cloud;
    cloud.ligandAtomicNumber = ligandType;
    cloud.positions.reserve((size_t)maxNeighborCandidates);
    cloud.colors.reserve((size_t)maxNeighborCandidates);

    constexpr float kUniquePointEpsilonSq = 1e-8f;
    for (const auto& candidate : candidates)
    {
        if (candidate.atomicNumber != ligandType)
            continue;
        if (candidate.dist > shellCutoff)
            break;
        if ((int)cloud.positions.size() >= maxNeighborCandidates)
            break;

        bool duplicate = false;
        for (const auto& existing : cloud.positions)
        {
            const glm::vec3 diff = existing - candidate.pos;
            if (glm::dot(diff, diff) <= kUniquePointEpsilonSq)
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate)
        {
            cloud.positions.push_back(candidate.pos);
            cloud.colors.push_back(candidate.color);
        }
    }

    return cloud;
}

std::vector<std::array<int, 3>> buildConvexHullTriangles(const glm::vec3& center,
                                                         const std::vector<glm::vec3>& neighbors)
{
    std::vector<std::array<int, 3>> triangles;
    const int n = (int)neighbors.size();
    if (n < 4)
        return triangles;

    std::set<std::tuple<int, int, int>> uniqueFacets;

    for (int a = 0; a < n - 2; ++a)
    {
        for (int b = a + 1; b < n - 1; ++b)
        {
            for (int c = b + 1; c < n; ++c)
            {
                const glm::vec3 pa = neighbors[(size_t)a] - center;
                const glm::vec3 pb = neighbors[(size_t)b] - center;
                const glm::vec3 pc = neighbors[(size_t)c] - center;

                glm::vec3 normal = glm::cross(pb - pa, pc - pa);
                const float area2 = glm::length(normal);
                if (area2 < 1e-6f)
                    continue;

                normal /= area2;
                float posCount = 0.0f;
                float negCount = 0.0f;
                for (int d = 0; d < n; ++d)
                {
                    if (d == a || d == b || d == c)
                        continue;

                    const float side = glm::dot(normal, neighbors[(size_t)d] - neighbors[(size_t)a]);
                    if (side > kFacetEpsilon)
                        posCount += 1.0f;
                    else if (side < -kFacetEpsilon)
                        negCount += 1.0f;

                    if (posCount > 0.0f && negCount > 0.0f)
                        break;
                }

                if (posCount > 0.0f && negCount > 0.0f)
                    continue;

                int ia = a;
                int ib = b;
                int ic = c;
                if (glm::dot(normal, (pa + pb + pc) * (1.0f / 3.0f)) < 0.0f)
                    std::swap(ib, ic);

                std::array<int, 3> sortedIdx = {ia, ib, ic};
                std::sort(sortedIdx.begin(), sortedIdx.end());
                const auto key = std::make_tuple(sortedIdx[0], sortedIdx[1], sortedIdx[2]);
                if (uniqueFacets.insert(key).second)
                    triangles.push_back({ia, ib, ic});
            }
        }
    }

    return triangles;
}

// buildMergedHullFaces removed — raw triangles from buildConvexHullTriangles
// are used directly so each face connects nearest-neighbor ligand atoms.

struct ProjectedFace
{
    ImVector<ImVec2> points;
    float depth   = 0.0f;
    ImU32 fillCol = IM_COL32_WHITE;
    ImU32 edgeCol = IM_COL32_WHITE;
};
}

void drawPolyhedralOverlay(ImDrawList* drawList,
                           const glm::mat4& projection,
                           const glm::mat4& view,
                           int framebufferWidth,
                           int framebufferHeight,
                           const Structure& structure,
                           const std::vector<int>& selectedInstanceIndices,
                           const std::vector<int>& instanceToAtomIndex,
                           const PolyhedralOverlaySettings& settings,
                           const std::vector<glm::vec3>& elementColors,
                           bool enabled)
{
    if (!enabled || drawList == nullptr || structure.atoms.size() < 4)
        return;

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    const bool usePbc = tryMakeCellMatrices(structure, cell, invCell);

    const std::vector<int> centers = gatherCenters(structure,
                                                   selectedInstanceIndices,
                                                   instanceToAtomIndex,
                                                   settings);
    if (centers.empty())
        return;

    const glm::vec3 cellOrigin((float)structure.cellOffset[0],
                               (float)structure.cellOffset[1],
                               (float)structure.cellOffset[2]);

    const float faceOpacity = std::max(0.0f, std::min(1.0f, settings.faceOpacity));
    const float edgeOpacity = std::max(0.0f, std::min(1.0f, settings.edgeOpacity));

    // Build the list of centre positions to draw polyhedra around.
    struct CenterEntry
    {
        glm::vec3 pos;
        int       atomicNumber = 0;
        glm::vec3 color = glm::vec3(0.8f);
    };
    std::vector<CenterEntry> centerEntries;
    centerEntries.reserve(centers.size() * 4);

    for (const int idx : centers)
    {
        if (idx < 0 || idx >= (int)structure.atoms.size())
            continue;
        const auto& atom = structure.atoms[(size_t)idx];
        const glm::vec3 cartPos((float)atom.x, (float)atom.y, (float)atom.z);
        // Use elementColors (live table from EditMenuDialogs) so the face color
        // updates immediately when the user edits an element color.
        const glm::vec3 col = (atom.atomicNumber >= 1 && atom.atomicNumber < (int)elementColors.size())
            ? elementColors[(size_t)atom.atomicNumber]
            : glm::vec3(atom.r, atom.g, atom.b);
        centerEntries.push_back({cartPos, atom.atomicNumber, col});

        if (!usePbc)
            continue;

        const glm::vec3 frac = invCell * (cartPos - cellOrigin);

        constexpr float kBndTol = 0.002f;
        for (int sx = -1; sx <= 1; ++sx)
        {
            for (int sy = -1; sy <= 1; ++sy)
            {
                for (int sz = -1; sz <= 1; ++sz)
                {
                    if (sx == 0 && sy == 0 && sz == 0)
                        continue;
                    const glm::vec3 imageFrac = frac + glm::vec3((float)sx, (float)sy, (float)sz);
                    if (imageFrac.x < -kBndTol || imageFrac.x > 1.0f + kBndTol) continue;
                    if (imageFrac.y < -kBndTol || imageFrac.y > 1.0f + kBndTol) continue;
                    if (imageFrac.z < -kBndTol || imageFrac.z > 1.0f + kBndTol) continue;
                    centerEntries.push_back({cellOrigin + cell * imageFrac, atom.atomicNumber, col});
                }
            }
        }
    }

    // Out-of-box ligand atoms (periodic images outside [0,1]^3) need to be
    // drawn as element-coloured spheres so they look identical to in-box atoms.
    struct OutsideLigand
    {
        glm::vec3 pos;
        glm::vec3 color;      // source atom display color (r,g,b)
        int atomicNumber = 0; // kept for reference
    };
    static std::vector<OutsideLigand> outsideLigands;
    outsideLigands.reserve(256);

    auto appendOutsideLigand = [&](const glm::vec3& p, const glm::vec3& col, int z)
    {
        if (!isOutsideUnitCell(structure, p))
            return;
        constexpr float kMergeSq = 0.001f * 0.001f;
        for (const auto& ex : outsideLigands)
        {
            const glm::vec3 d = ex.pos - p;
            if (glm::dot(d, d) <= kMergeSq)
                return;
        }
        outsideLigands.push_back({p, col, z});
    };

    // --- first pass: collect hull geometry + out-of-box ligands ---
    struct PolyhedronData
    {
        std::vector<std::array<int, 3>> facets;
        std::vector<glm::vec3>          neighbors;
        ImU32 fillCol = IM_COL32_WHITE;
        ImU32 edgeCol = IM_COL32_WHITE;
    };
    // World-space geometry survives camera moves. A complete content hash also
    // invalidates it after edits, tab switches, recoloring or selection changes.
    static std::vector<PolyhedronData> polyhedra;
    static std::size_t previousKey=0, nextCenter=0;
    std::size_t key=structure.atoms.size();
    const auto hash=[&](double value) {
        key ^= std::hash<double>{}(value)+0x9e3779b9+(key<<6)+(key>>2);
    };
    hash(usePbc); hash(faceOpacity); hash(edgeOpacity); hash(settings.maxNeighborCandidates);
    hash(settings.ligandElementFilterEnabled);
    for (bool value:settings.ligandElementMask) hash(value);
    for (const auto& v:structure.cellVectors) for (double value:v) hash(value);
    for (double value:structure.cellOffset) hash(value);
    for (const auto& a:structure.atoms) {
        hash(a.x); hash(a.y); hash(a.z); hash(a.atomicNumber); hash(a.r); hash(a.g); hash(a.b);
    }
    for (const auto& color:elementColors) { hash(color.r); hash(color.g); hash(color.b); }
    for (const auto& entry:centerEntries) {
        hash(entry.pos.x); hash(entry.pos.y); hash(entry.pos.z); hash(entry.atomicNumber);
    }
    if (key!=previousKey) {
        previousKey=key; nextCenter=0; polyhedra.clear(); outsideLigands.clear();
    }
    polyhedra.reserve(centerEntries.size());
    // Large structures fill progressively; orbit and pan never rebuild hulls.
    const auto end=std::min(centerEntries.size(),nextCenter+4);
    for (;nextCenter<end;++nextCenter)
    {
        const auto& entry=centerEntries[nextCenter];
        const NeighborCloud cloud = buildNeighborCloud(structure,
                                                       entry.pos,
                                                       entry.atomicNumber,
                                                       usePbc,
                                                       cell,
                                                       invCell,
                                                       settings,
                                                       elementColors);
        if ((int)cloud.positions.size() < kMinNeighborCount)
            continue;

        for (size_t k = 0; k < cloud.positions.size(); ++k)
        {
            const glm::vec3& col = cloud.colors.size() > k ? cloud.colors[k] : glm::vec3(0.8f);
            appendOutsideLigand(cloud.positions[k], col, cloud.ligandAtomicNumber);
        }

        const std::vector<std::array<int, 3>> facets = buildConvexHullTriangles(entry.pos, cloud.positions);
        if (facets.empty())
            continue;

        // Face color = center atom CPK color, darkened slightly for edge.
        const glm::vec3& c = entry.color;
        const glm::vec3 dark = c * 0.55f;
        const ImU32 fillCol = ImColor(c.r, c.g, c.b, faceOpacity);
        const ImU32 edgeCol = ImColor(dark.r, dark.g, dark.b, edgeOpacity);

        polyhedra.push_back({facets, cloud.positions, fillCol, edgeCol});
    }

    // --- second pass: draw out-of-box ligand atom markers (before polyhedra) ---
    // Uses the same source-atom display color as in-box atoms, plus a simulated
    // sphere highlight so the shading matches the OpenGL-lit in-box spheres.
    constexpr float kAtomWorldRadius = 0.65f;
    const float projYScale = projection[1][1]; // = 1/tan(fovY/2)

    for (const auto& lig : outsideLigands)
    {
        const glm::vec4 viewPos4 = view * glm::vec4(lig.pos, 1.0f);
        const float viewZ = viewPos4.z;
        if (viewZ >= -kNearClip)
            continue;

        const float depth = -viewZ;
        const float screenRadius = kAtomWorldRadius * projYScale * ((float)framebufferHeight * 0.5f) / depth;
        if (screenRadius < 1.0f)
            continue;

        ImVec2 screen;
        if (!projectToScreen(lig.pos, projection, view, framebufferWidth, framebufferHeight, screen))
            continue;

        const glm::vec3& col  = lig.color;
        const glm::vec3  dark = col * 0.45f;      // rim / shadow
        const glm::vec3  mid  = col;               // base
        const glm::vec3  lite = col * 0.55f + glm::vec3(0.45f); // highlight tint

        // 1. Dark rim (1 px border)
        drawList->AddCircleFilled(screen, screenRadius + 1.0f,
                                  ImColor(dark.r, dark.g, dark.b, 1.0f), 24);
        // 2. Base color fill
        drawList->AddCircleFilled(screen, screenRadius,
                                  ImColor(mid.r, mid.g, mid.b, 1.0f), 24);
        // 3. Specular highlight: small bright circle offset upper-left
        const float hlR = screenRadius * 0.38f;
        const ImVec2 hlCenter(screen.x - screenRadius * 0.28f,
                              screen.y - screenRadius * 0.28f);
        drawList->AddCircleFilled(hlCenter, hlR,
                                  ImColor(lite.r, lite.g, lite.b, 0.70f), 16);
    }

    // --- third pass: draw polyhedra ---
    // ALL triangles from ALL polyhedra are collected into one flat list and
    // sorted globally back-to-front.  This is the painter's-algorithm equivalent
    // of depth testing: a triangle from polyhedron B can correctly appear in
    // front of one from polyhedron A without any per-object ordering bias.
    std::vector<ProjectedFace> allFaces;
    allFaces.reserve(polyhedra.size() * 8);

    for (const auto& poly : polyhedra)
    {
        const std::vector<std::array<int, 3>>& facets    = poly.facets;
        const std::vector<glm::vec3>&           neighbors = poly.neighbors;

        for (const auto& tri : facets)
        {
            const glm::vec3 p0 = neighbors[(size_t)tri[0]];
            const glm::vec3 p1 = neighbors[(size_t)tri[1]];
            const glm::vec3 p2 = neighbors[(size_t)tri[2]];

            const std::vector<glm::vec3> polygon = {p0, p1, p2};
            const auto clipped = clipPolygonAgainstNearPlane(polygon, view, kNearClip);
            if (clipped.size() < 3)
                continue;

            ProjectedFace projectedFace;
            projectedFace.fillCol = poly.fillCol;
            projectedFace.edgeCol = poly.edgeCol;
            projectedFace.points.reserve((int)clipped.size());
            bool projected = true;
            float depthAccum = 0.0f;
            for (const auto& vertex : clipped)
            {
                ImVec2 screen;
                if (!projectToScreen(vertex,
                                     projection,
                                     view,
                                     framebufferWidth,
                                     framebufferHeight,
                                     screen))
                {
                    projected = false;
                    break;
                }
                projectedFace.points.push_back(screen);
                depthAccum += (view * glm::vec4(vertex, 1.0f)).z;
            }

            if (!projected || projectedFace.points.Size < 3)
                continue;

            projectedFace.depth = depthAccum / (float)clipped.size();
            allFaces.push_back(std::move(projectedFace));
        }
    }

    // Global back-to-front sort (view-space z < 0 in front of camera).
    std::sort(allFaces.begin(), allFaces.end(),
              [](const ProjectedFace& lhs, const ProjectedFace& rhs)
              { return lhs.depth < rhs.depth; });

    for (const auto& face : allFaces)
    {
        drawList->AddConvexPolyFilled(face.points.Data, face.points.Size, face.fillCol);
        if (settings.showEdges)
            drawList->AddPolyline(face.points.Data, face.points.Size, face.edgeCol, ImDrawFlags_Closed, 1.2f);
    }
}
