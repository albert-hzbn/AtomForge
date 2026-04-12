#include "algorithms/NanoCrystalBuilder.h"

#include "ElementData.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <unordered_map>

namespace
{
bool rayIntersectsTriangle(const glm::vec3& orig,
                           const glm::vec3& dir,
                           const glm::vec3& v0,
                           const glm::vec3& v1,
                           const glm::vec3& v2)
{
    const float eps = 1e-8f;
    const glm::vec3 e1 = v1 - v0;
    const glm::vec3 e2 = v2 - v0;
    const glm::vec3 p = glm::cross(dir, e2);
    const float det = glm::dot(e1, p);
    if (std::abs(det) < eps)
        return false;

    const float invDet = 1.0f / det;
    const glm::vec3 t = orig - v0;
    const float u = glm::dot(t, p) * invDet;
    if (u < 0.0f || u > 1.0f)
        return false;

    const glm::vec3 q = glm::cross(t, e1);
    const float v = glm::dot(dir, q) * invDet;
    if (v < 0.0f || (u + v) > 1.0f)
        return false;

    const float hitT = glm::dot(e2, q) * invDet;
    return hitT > eps;
}

bool isInsideMeshModel(const glm::vec3& p,
                       const NanoParams& params,
                       const std::vector<glm::vec3>& modelVertices,
                       const std::vector<unsigned int>& modelIndices)
{
    if (modelVertices.empty() || modelIndices.size() < 3)
        return false;

    const float invScale = (params.modelScale > 1e-8f) ? (1.0f / params.modelScale) : 1.0f;
    const glm::vec3 pm = p * invScale;

    // Use 3 ray directions and take majority vote for robustness
    // against meshes with edge/vertex coincidence along a single axis.
    static const glm::vec3 rayDirs[3] = {
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 1.0f),
    };

    int insideVotes = 0;
    for (int r = 0; r < 3; ++r)
    {
        int hitCount = 0;
        for (size_t i = 0; i + 2 < modelIndices.size(); i += 3)
        {
            const unsigned int i0 = modelIndices[i + 0];
            const unsigned int i1 = modelIndices[i + 1];
            const unsigned int i2 = modelIndices[i + 2];
            if (i0 >= modelVertices.size() || i1 >= modelVertices.size() || i2 >= modelVertices.size())
                continue;

            if (rayIntersectsTriangle(pm, rayDirs[r], modelVertices[i0], modelVertices[i1], modelVertices[i2]))
                ++hitCount;
        }
        if ((hitCount % 2) == 1)
            ++insideVotes;
    }

    return insideVotes >= 2;
}
}

// -- Shape label -------------------------------------------------------------

const char* shapeLabel(NanoShape s)
{
    switch (s) {
        case NanoShape::Sphere:              return "Sphere";
        case NanoShape::Ellipsoid:           return "Ellipsoid";
        case NanoShape::Box:                 return "Box";
        case NanoShape::Cylinder:            return "Cylinder";
        case NanoShape::Octahedron:          return "Octahedron";
        case NanoShape::TruncatedOctahedron: return "Truncated Octahedron";
        case NanoShape::Cuboctahedron:       return "Cuboctahedron";
        case NanoShape::MeshModel:           return "3D Model Fill";
    }
    return "Unknown";
}

// -- Geometry helpers --------------------------------------------------------

float computeBoundingRadius(const NanoParams& p)
{
    switch (p.shape) {
        case NanoShape::Sphere:
            return p.sphereRadius;
        case NanoShape::Ellipsoid:
            return std::max({p.ellipRx, p.ellipRy, p.ellipRz});
        case NanoShape::Box:
            return std::sqrt(p.boxHx*p.boxHx + p.boxHy*p.boxHy + p.boxHz*p.boxHz);
        case NanoShape::Cylinder:
            return std::sqrt(p.cylRadius*p.cylRadius
                             + (p.cylHeight*0.5f)*(p.cylHeight*0.5f));
        case NanoShape::Octahedron:
            return p.octRadius;
        case NanoShape::TruncatedOctahedron:
            return std::max(p.truncOctRadius, p.truncOctTrunc * std::sqrt(3.0f));
        case NanoShape::Cuboctahedron:
            return p.cuboRadius;
        case NanoShape::MeshModel:
            return p.modelScale * std::sqrt(p.modelHx*p.modelHx + p.modelHy*p.modelHy + p.modelHz*p.modelHz);
    }
    return 10.0f;
}

HalfExtents computeShapeHalfExtents(const NanoParams& p)
{
    switch (p.shape) {
        case NanoShape::Sphere:
            return {p.sphereRadius, p.sphereRadius, p.sphereRadius};
        case NanoShape::Ellipsoid:
            return {p.ellipRx, p.ellipRy, p.ellipRz};
        case NanoShape::Box:
            return {p.boxHx, p.boxHy, p.boxHz};
        case NanoShape::Cylinder:
            if (p.cylAxis == 0) return {p.cylHeight*0.5f, p.cylRadius, p.cylRadius};
            if (p.cylAxis == 1) return {p.cylRadius, p.cylHeight*0.5f, p.cylRadius};
            return {p.cylRadius, p.cylRadius, p.cylHeight*0.5f};
        case NanoShape::Octahedron:
            return {p.octRadius, p.octRadius, p.octRadius};
        case NanoShape::TruncatedOctahedron: {
            float r = std::min(p.truncOctTrunc, p.truncOctRadius);
            return {r, r, r};
        }
        case NanoShape::Cuboctahedron:
            return {p.cuboRadius, p.cuboRadius, p.cuboRadius};
        case NanoShape::MeshModel:
            return {p.modelHx * p.modelScale,
                    p.modelHy * p.modelScale,
                    p.modelHz * p.modelScale};
    }
    return {10.0f, 10.0f, 10.0f};
}

bool isInsideShape(const glm::vec3& p,
                   const NanoParams& params,
                   const std::vector<glm::vec3>& modelVertices,
                   const std::vector<unsigned int>& modelIndices)
{
    switch (params.shape) {
        case NanoShape::Sphere:
            return glm::dot(p, p) <= params.sphereRadius * params.sphereRadius;
        case NanoShape::Ellipsoid: {
            float fx = p.x / params.ellipRx;
            float fy = p.y / params.ellipRy;
            float fz = p.z / params.ellipRz;
            return fx*fx + fy*fy + fz*fz <= 1.0f;
        }
        case NanoShape::Box:
            return std::abs(p.x) <= params.boxHx &&
                   std::abs(p.y) <= params.boxHy &&
                   std::abs(p.z) <= params.boxHz;
        case NanoShape::Cylinder: {
            float r2, ax;
            if      (params.cylAxis == 0) { r2 = p.y*p.y + p.z*p.z; ax = p.x; }
            else if (params.cylAxis == 1) { r2 = p.x*p.x + p.z*p.z; ax = p.y; }
            else                          { r2 = p.x*p.x + p.y*p.y; ax = p.z; }
            return r2 <= params.cylRadius * params.cylRadius &&
                   std::abs(ax) <= params.cylHeight * 0.5f;
        }
        case NanoShape::Octahedron:
            return std::abs(p.x) + std::abs(p.y) + std::abs(p.z) <= params.octRadius;
        case NanoShape::TruncatedOctahedron:
            return (std::abs(p.x) + std::abs(p.y) + std::abs(p.z) <= params.truncOctRadius) &&
                   std::abs(p.x) <= params.truncOctTrunc &&
                   std::abs(p.y) <= params.truncOctTrunc &&
                   std::abs(p.z) <= params.truncOctTrunc;
        case NanoShape::Cuboctahedron:
            return std::abs(p.x) + std::abs(p.y) <= params.cuboRadius &&
                   std::abs(p.y) + std::abs(p.z) <= params.cuboRadius &&
                   std::abs(p.x) + std::abs(p.z) <= params.cuboRadius;
        case NanoShape::MeshModel:
            return isInsideMeshModel(p, params, modelVertices, modelIndices);
    }
    return false;
}

glm::vec3 computeAtomCentroid(const std::vector<AtomSite>& atoms)
{
    if (atoms.empty()) return glm::vec3(0.0f);
    glm::vec3 sum(0.0f);
    for (const AtomSite& a : atoms)
        sum += glm::vec3((float)a.x, (float)a.y, (float)a.z);
    return sum / (float)atoms.size();
}

static float safeLen3(const glm::vec3& v)
{
    return std::sqrt(glm::dot(v, v));
}

// -- Builder -----------------------------------------------------------------

NanoBuildResult buildNanocrystal(Structure& structure,
                                 const Structure& reference,
                                 const NanoParams& params,
                                 const std::vector<glm::vec3>& elementColors,
                                 const std::vector<glm::vec3>& modelVertices,
                                 const std::vector<unsigned int>& modelIndices)
{
    NanoBuildResult result;
    result.shape      = params.shape;
    result.inputAtoms = (int)reference.atoms.size();

    if (params.shape == NanoShape::MeshModel)
    {
        if (modelVertices.empty() || modelIndices.size() < 3)
        {
            result.message = "No 3D model loaded. Drop an OBJ/STL file in the builder first.";
            return result;
        }
        if (params.modelScale <= 0.0f)
        {
            result.message = "Model scale must be positive.";
            return result;
        }
    }

    if (reference.atoms.empty()) {
        result.message = "Reference structure has no atoms.";
        return result;
    }

    glm::vec3 center;
    if (params.autoCenterFromAtoms)
        center = computeAtomCentroid(reference.atoms);
    else
        center = glm::vec3(params.cx, params.cy, params.cz);

    const float maxR = computeBoundingRadius(params);
    std::vector<AtomSite> generatedAtoms;

    if (reference.hasUnitCell) {
        const auto& cv = reference.cellVectors;
        const glm::vec3 a((float)cv[0][0], (float)cv[0][1], (float)cv[0][2]);
        const glm::vec3 b((float)cv[1][0], (float)cv[1][1], (float)cv[1][2]);
        const glm::vec3 c((float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
        const float la = safeLen3(a), lb = safeLen3(b), lc = safeLen3(c);

        if (la < 1e-8f || lb < 1e-8f || lc < 1e-8f) {
            result.message = "Reference structure has degenerate lattice vectors.";
            return result;
        }

        const int kMaxReps = 40;
        bool clamped = false;
        auto safeRep = [&](float L) -> int {
            int n = (int)std::ceil(maxR / L) + 2;
            if (n > kMaxReps) { clamped = true; n = kMaxReps; }
            return n;
        };

        int nA, nB, nC;
        if (params.autoReplicate) {
            nA = safeRep(la); nB = safeRep(lb); nC = safeRep(lc);
        } else {
            nA = std::max(1, params.repA);
            nB = std::max(1, params.repB);
            nC = std::max(1, params.repC);
            clamped = false;
        }

        result.repA = nA; result.repB = nB; result.repC = nC;
        result.tilingUsed  = true;
        result.repClamped  = clamped;

        const long long total =
            (long long)(2*nA+1) * (long long)(2*nB+1) * (long long)(2*nC+1)
            * (long long)reference.atoms.size();

        if (total > 8000000LL) {
            std::ostringstream msg;
            msg << "Tiling would test " << total
                << " atoms (limit 8M). Reduce shape size or use manual replication.";
            result.message = msg.str();
            return result;
        }

        generatedAtoms.reserve((size_t)std::min(total, (long long)2000000));

        for (int ia = -nA; ia <= nA; ++ia)
        for (int ib = -nB; ib <= nB; ++ib)
        for (int ic = -nC; ic <= nC; ++ic) {
            const glm::vec3 offset = (float)ia*a + (float)ib*b + (float)ic*c;
            for (const AtomSite& atom : reference.atoms) {
                const glm::vec3 pos(
                    (float)atom.x + offset.x,
                    (float)atom.y + offset.y,
                    (float)atom.z + offset.z);
                if (!isInsideShape(pos - center, params, modelVertices, modelIndices)) continue;

                AtomSite out = atom;
                out.x = (double)pos.x;
                out.y = (double)pos.y;
                out.z = (double)pos.z;
                int z = out.atomicNumber;
                if (z >= 0 && z < (int)elementColors.size()) {
                    out.r = elementColors[z].r;
                    out.g = elementColors[z].g;
                    out.b = elementColors[z].b;
                } else {
                    getDefaultElementColor(z, out.r, out.g, out.b);
                }
                generatedAtoms.push_back(out);
            }
        }
    } else {
        result.tilingUsed = false;
        generatedAtoms.reserve(reference.atoms.size());
        for (const AtomSite& atom : reference.atoms) {
            const glm::vec3 pos((float)atom.x, (float)atom.y, (float)atom.z);
            if (!isInsideShape(pos - center, params, modelVertices, modelIndices)) continue;
            AtomSite out = atom;
            int z = out.atomicNumber;
            if (z >= 0 && z < (int)elementColors.size()) {
                out.r = elementColors[z].r;
                out.g = elementColors[z].g;
                out.b = elementColors[z].b;
            } else {
                getDefaultElementColor(z, out.r, out.g, out.b);
            }
            generatedAtoms.push_back(out);
        }
    }

    if (generatedAtoms.empty()) {
        result.message =
            "No atoms within the specified shape. "
            "Try increasing the size parameter(s) or adjusting Model Scale.";
        return result;
    }

    // Remove duplicate atoms at unit-cell boundary overlaps using spatial hash.
    {
        const float dedupTol = 0.01f; // 0.01 A tolerance
        const float cellSize = 0.5f;  // hash cell size in A
        const float invCell  = 1.0f / cellSize;

        struct Vec3Hash {
            float invC;
            Vec3Hash(float ic) : invC(ic) {}
            size_t operator()(const glm::ivec3& v) const {
                size_t h = (size_t)(v.x * 73856093u);
                h ^= (size_t)(v.y * 19349663u);
                h ^= (size_t)(v.z * 83492791u);
                return h;
            }
        };
        struct Vec3Eq {
            bool operator()(const glm::ivec3& a, const glm::ivec3& b) const {
                return a.x == b.x && a.y == b.y && a.z == b.z;
            }
        };
        std::unordered_map<glm::ivec3, std::vector<size_t>, Vec3Hash, Vec3Eq> grid(
            generatedAtoms.size(), Vec3Hash(invCell));

        std::vector<AtomSite> unique;
        unique.reserve(generatedAtoms.size());
        const float dedupTol2 = dedupTol * dedupTol;

        for (size_t i = 0; i < generatedAtoms.size(); ++i)
        {
            const glm::vec3 pi((float)generatedAtoms[i].x,
                               (float)generatedAtoms[i].y,
                               (float)generatedAtoms[i].z);
            const glm::ivec3 cell((int)std::floor(pi.x * invCell),
                                  (int)std::floor(pi.y * invCell),
                                  (int)std::floor(pi.z * invCell));
            bool isDup = false;
            // Check 3x3x3 neighborhood of hash cells
            for (int dx = -1; dx <= 1 && !isDup; ++dx)
            for (int dy = -1; dy <= 1 && !isDup; ++dy)
            for (int dz = -1; dz <= 1 && !isDup; ++dz)
            {
                const glm::ivec3 nc(cell.x + dx, cell.y + dy, cell.z + dz);
                auto it = grid.find(nc);
                if (it == grid.end()) continue;
                for (size_t idx : it->second)
                {
                    const glm::vec3 pj((float)unique[idx].x,
                                       (float)unique[idx].y,
                                       (float)unique[idx].z);
                    const glm::vec3 d = pi - pj;
                    if (glm::dot(d, d) < dedupTol2)
                    {
                        isDup = true;
                        break;
                    }
                }
            }
            if (!isDup)
            {
                const size_t uid = unique.size();
                unique.push_back(generatedAtoms[i]);
                grid[cell].push_back(uid);
            }
        }
        generatedAtoms.swap(unique);
    }

    result.estimatedDiameter = 2.0f * maxR;
    result.outputAtoms = (int)generatedAtoms.size();
    structure.atoms.swap(generatedAtoms);

    // Clear stale per-atom metadata from any previous structure.
    structure.grainColors.clear();
    structure.grainRegionIds.clear();
    structure.pbcBoundaryTol = 0.0f;
    structure.ipfLoadStatus.clear();

    if (params.setOutputCell) {
        const HalfExtents he = computeShapeHalfExtents(params);
        const float pad = params.vacuumPadding;
        structure.hasUnitCell = true;
        structure.cellOffset  = {{
            (double)(center.x - he.hx - pad),
            (double)(center.y - he.hy - pad),
            (double)(center.z - he.hz - pad) }};
        structure.cellVectors = {{
            {{ 2.0*(he.hx+pad), 0.0, 0.0 }},
            {{ 0.0, 2.0*(he.hy+pad), 0.0 }},
            {{ 0.0, 0.0, 2.0*(he.hz+pad) }} }};
    } else {
        structure.hasUnitCell = false;
    }

    result.success = true;
    if (params.shape == NanoShape::MeshModel)
    {
        std::ostringstream msg;
        msg << "Custom structure built: " << result.outputAtoms << " atoms";
        if (result.tilingUsed)
            msg << " (reps=" << result.repA << "x" << result.repB << "x" << result.repC << ")";
        result.message = msg.str();
    }
    else
    {
        result.message = "Nanocrystal built successfully.";
    }
    return result;
}
