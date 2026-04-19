#include "algorithms/NanoCrystalBuilder.h"

// The Wulff-construction path in this builder follows the equilibrium-shape
// formulation reviewed in: Georgios D. Barmparis, Zbigniew Lodziana, Nuria
// Lopez, and Ioannis N. Remediakis, "Nanoparticle shapes by using Wulff
// constructions and first-principles calculations," Beilstein J.
// Nanotechnol. 6, 361-368 (2015).

#include "ElementData.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
#include <limits>
#include <numeric>
#include <sstream>
#include <unordered_map>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

namespace
{
struct SymmetryRotationSet
{
    bool success = false;
    std::string error;
    std::vector<std::array<std::array<int, 3>, 3>> rotations;
};

struct MillerKey
{
    int h = 0;
    int k = 0;
    int l = 0;

    bool operator==(const MillerKey& other) const
    {
        return h == other.h && k == other.k && l == other.l;
    }
};

struct MillerKeyHash
{
    std::size_t operator()(const MillerKey& key) const
    {
        std::size_t value = (std::size_t)(key.h * 73856093);
        value ^= (std::size_t)(key.k * 19349663);
        value ^= (std::size_t)(key.l * 83492791);
        return value;
    }
};

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

static glm::mat3 rotationAligningDirectionToZ(const glm::vec3& fromDir);
void buildReciprocalBasis(const glm::vec3& a,
                         const glm::vec3& b,
                         const glm::vec3& c,
                         glm::vec3& aStar,
                         glm::vec3& bStar,
                         glm::vec3& cStar);

bool invertIntegerMatrix(const std::array<std::array<int, 3>, 3>& matrix,
                        std::array<std::array<int, 3>, 3>& inverse)
{
    const int a = matrix[0][0], b = matrix[0][1], c = matrix[0][2];
    const int d = matrix[1][0], e = matrix[1][1], f = matrix[1][2];
    const int g = matrix[2][0], h = matrix[2][1], i = matrix[2][2];

    const int c00 =  (e * i - f * h);
    const int c01 = -(d * i - f * g);
    const int c02 =  (d * h - e * g);
    const int c10 = -(b * i - c * h);
    const int c11 =  (a * i - c * g);
    const int c12 = -(a * h - b * g);
    const int c20 =  (b * f - c * e);
    const int c21 = -(a * f - c * d);
    const int c22 =  (a * e - b * d);

    const int determinant = a * c00 + b * c01 + c * c02;
    if (determinant != 1 && determinant != -1)
        return false;

    inverse = {{
        {{ c00 / determinant, c10 / determinant, c20 / determinant }},
        {{ c01 / determinant, c11 / determinant, c21 / determinant }},
        {{ c02 / determinant, c12 / determinant, c22 / determinant }}
    }};
    return true;
}

glm::vec3 buildMillerPlaneNormal(const Structure& reference,
                                 int h,
                                 int k,
                                 int l,
                                 std::string* errorMessage = nullptr)
{
    if (!reference.hasUnitCell)
    {
        if (errorMessage)
            *errorMessage = "Miller orientation requires a reference structure with a unit cell.";
        return glm::vec3(0.0f);
    }

    const glm::vec3 a((float)reference.cellVectors[0][0], (float)reference.cellVectors[0][1], (float)reference.cellVectors[0][2]);
    const glm::vec3 b((float)reference.cellVectors[1][0], (float)reference.cellVectors[1][1], (float)reference.cellVectors[1][2]);
    const glm::vec3 c((float)reference.cellVectors[2][0], (float)reference.cellVectors[2][1], (float)reference.cellVectors[2][2]);
    const float volume = glm::dot(a, glm::cross(b, c));
    if (std::abs(volume) <= 1e-8f)
    {
        if (errorMessage)
            *errorMessage = "Reference structure has a degenerate lattice.";
        return glm::vec3(0.0f);
    }

    glm::vec3 aStar(0.0f), bStar(0.0f), cStar(0.0f);
    buildReciprocalBasis(a, b, c, aStar, bStar, cStar);
    const glm::vec3 normal = (float)h * aStar + (float)k * bStar + (float)l * cStar;
    const float normalLength = glm::length(normal);
    if (normalLength <= 1e-8f)
    {
        if (errorMessage)
            *errorMessage = "Miller indices do not define a valid plane normal for the current lattice.";
        return glm::vec3(0.0f);
    }

    return normal / normalLength;
}

glm::mat3 buildOrientationMatrix(const Structure& reference,
                                 const NanoParams& params,
                                 std::string* errorMessage = nullptr)
{
    glm::mat3 orientation(1.0f);
    if (!params.applyCrystalOrientation)
        return orientation;

    if (params.useMillerOrientation)
    {
        if (params.millerH == 0 && params.millerK == 0 && params.millerL == 0)
        {
            if (errorMessage)
                *errorMessage = "Miller indices [h k l] cannot all be zero.";
            return glm::mat3(0.0f);
        }

        std::string normalError;
        const glm::vec3 planeNormal = buildMillerPlaneNormal(reference,
                                                             params.millerH,
                                                             params.millerK,
                                                             params.millerL,
                                                             &normalError);
        if (!normalError.empty())
        {
            if (errorMessage)
                *errorMessage = normalError;
            return glm::mat3(0.0f);
        }

        return rotationAligningDirectionToZ(planeNormal);
    }

    const float rx = glm::radians(params.orientXDeg);
    const float ry = glm::radians(params.orientYDeg);
    const float rz = glm::radians(params.orientZDeg);
    const glm::mat4 rot4 = glm::rotate(glm::mat4(1.0f), rz, glm::vec3(0.0f, 0.0f, 1.0f))
                         * glm::rotate(glm::mat4(1.0f), ry, glm::vec3(0.0f, 1.0f, 0.0f))
                         * glm::rotate(glm::mat4(1.0f), rx, glm::vec3(1.0f, 0.0f, 0.0f));
    return glm::mat3(rot4);
}

glm::vec3 buildClipCenter(const Structure& reference, const NanoParams& params)
{
    if (params.autoCenterFromAtoms)
        return computeAtomCentroid(reference.atoms);
    return glm::vec3(params.cx, params.cy, params.cz);
}

glm::ivec3 reduceMillerIndices(const glm::ivec3& hkl)
{
    const int divisor = std::gcd(std::abs(hkl.x), std::gcd(std::abs(hkl.y), std::abs(hkl.z)));
    if (divisor <= 1)
        return hkl;
    return glm::ivec3(hkl.x / divisor, hkl.y / divisor, hkl.z / divisor);
}

bool isZeroMiller(const glm::ivec3& hkl)
{
    return hkl.x == 0 && hkl.y == 0 && hkl.z == 0;
}

glm::ivec3 transformMillerIndices(const std::array<std::array<int, 3>, 3>& rotation,
                                  const glm::ivec3& hkl)
{
    std::array<std::array<int, 3>, 3> inverse = {{
        {{1, 0, 0}},
        {{0, 1, 0}},
        {{0, 0, 1}}
    }};
    if (!invertIntegerMatrix(rotation, inverse))
        return hkl;

    return glm::ivec3(
        inverse[0][0] * hkl.x + inverse[1][0] * hkl.y + inverse[2][0] * hkl.z,
        inverse[0][1] * hkl.x + inverse[1][1] * hkl.y + inverse[2][1] * hkl.z,
        inverse[0][2] * hkl.x + inverse[1][2] * hkl.y + inverse[2][2] * hkl.z);
}

void buildReciprocalBasis(const glm::vec3& a,
                         const glm::vec3& b,
                         const glm::vec3& c,
                         glm::vec3& aStar,
                         glm::vec3& bStar,
                         glm::vec3& cStar)
{
    const float volume = glm::dot(a, glm::cross(b, c));
    aStar = glm::cross(b, c) / volume;
    bStar = glm::cross(c, a) / volume;
    cStar = glm::cross(a, b) / volume;
}

SymmetryRotationSet collectSymmetryRotations(const Structure& reference)
{
    SymmetryRotationSet result;
    result.success = false;

#ifndef ATOMS_ENABLE_SPGLIB
    result.error = "Wulff construction requires spglib support at build time.";
    return result;
#else
    if (!reference.hasUnitCell)
    {
        result.error = "Wulff construction requires a reference structure with a unit cell.";
        return result;
    }
    if (reference.atoms.empty())
    {
        result.error = "Reference structure has no atoms.";
        return result;
    }

    double lattice[3][3] = {
        {reference.cellVectors[0][0], reference.cellVectors[0][1], reference.cellVectors[0][2]},
        {reference.cellVectors[1][0], reference.cellVectors[1][1], reference.cellVectors[1][2]},
        {reference.cellVectors[2][0], reference.cellVectors[2][1], reference.cellVectors[2][2]}
    };

    const glm::vec3 a((float)reference.cellVectors[0][0], (float)reference.cellVectors[0][1], (float)reference.cellVectors[0][2]);
    const glm::vec3 b((float)reference.cellVectors[1][0], (float)reference.cellVectors[1][1], (float)reference.cellVectors[1][2]);
    const glm::vec3 c((float)reference.cellVectors[2][0], (float)reference.cellVectors[2][1], (float)reference.cellVectors[2][2]);
    const glm::mat3 latticeMat(a, b, c);
    const float det = glm::determinant(latticeMat);
    if (std::abs(det) <= 1e-10f)
    {
        result.error = "Reference structure has a degenerate lattice.";
        return result;
    }

    const glm::mat3 invLattice = glm::inverse(latticeMat);
    const glm::vec3 origin((float)reference.cellOffset[0],
                           (float)reference.cellOffset[1],
                           (float)reference.cellOffset[2]);

    std::vector<std::array<double, 3>> positions(reference.atoms.size());
    std::vector<int> types(reference.atoms.size(), 0);
    for (size_t index = 0; index < reference.atoms.size(); ++index)
    {
        const AtomSite& atom = reference.atoms[index];
        glm::vec3 frac = invLattice * (glm::vec3((float)atom.x, (float)atom.y, (float)atom.z) - origin);
        frac.x -= std::floor(frac.x);
        frac.y -= std::floor(frac.y);
        frac.z -= std::floor(frac.z);
        positions[index] = {{frac.x, frac.y, frac.z}};
        types[index] = (atom.atomicNumber > 0) ? atom.atomicNumber : (int)index + 1000;
    }

    const SpglibDataset* dataset = spg_get_dataset(
        lattice,
        reinterpret_cast<double (*)[3]>(positions.data()),
        types.data(),
        (int)positions.size(),
        1e-5);

    if (!dataset)
    {
        result.error = "spglib could not determine symmetry for the reference structure.";
        return result;
    }

    result.rotations.reserve((size_t)dataset->n_operations);
    std::unordered_map<std::string, bool> seen;
    for (int operationIndex = 0; operationIndex < dataset->n_operations; ++operationIndex)
    {
        std::array<std::array<int, 3>, 3> rotation = {{
            {{dataset->rotations[operationIndex][0][0], dataset->rotations[operationIndex][0][1], dataset->rotations[operationIndex][0][2]}},
            {{dataset->rotations[operationIndex][1][0], dataset->rotations[operationIndex][1][1], dataset->rotations[operationIndex][1][2]}},
            {{dataset->rotations[operationIndex][2][0], dataset->rotations[operationIndex][2][1], dataset->rotations[operationIndex][2][2]}}
        }};

        std::ostringstream key;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                key << rotation[row][col] << ',';
        if (seen.emplace(key.str(), true).second)
            result.rotations.push_back(rotation);
    }

    spg_free_dataset(const_cast<SpglibDataset*>(dataset));
    result.success = !result.rotations.empty();
    if (!result.success)
        result.error = "No symmetry rotations were returned by spglib.";
    return result;
#endif
}

bool appendUniquePoint(std::vector<glm::vec3>& points,
                       const glm::vec3& candidate,
                       float tolerance)
{
    const float tolerance2 = tolerance * tolerance;
    for (const glm::vec3& point : points)
    {
        const glm::vec3 delta = point - candidate;
        if (glm::dot(delta, delta) <= tolerance2)
            return false;
    }
    points.push_back(candidate);
    return true;
}

void sortFaceVertices(std::vector<glm::vec3>& vertices, const glm::vec3& normal)
{
    if (vertices.size() < 3)
        return;

    glm::vec3 centroid(0.0f);
    for (const glm::vec3& vertex : vertices)
        centroid += vertex;
    centroid /= (float)vertices.size();

    glm::vec3 referenceAxis = (std::abs(normal.z) < 0.9f)
        ? glm::vec3(0.0f, 0.0f, 1.0f)
        : glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 basisU = glm::normalize(glm::cross(referenceAxis, normal));
    glm::vec3 basisV = glm::normalize(glm::cross(normal, basisU));

    std::sort(vertices.begin(), vertices.end(), [&](const glm::vec3& left, const glm::vec3& right) {
        const glm::vec3 dl = left - centroid;
        const glm::vec3 dr = right - centroid;
        const float angleL = std::atan2(glm::dot(dl, basisV), glm::dot(dl, basisU));
        const float angleR = std::atan2(glm::dot(dr, basisV), glm::dot(dr, basisU));
        return angleL < angleR;
    });
}

bool isInsideWulffPolyhedron(const glm::vec3& localPos, const WulffPreview& preview, float tolerance)
{
    for (const WulffPlane& plane : preview.planes)
    {
        if (glm::dot(plane.normal, localPos) > plane.distance + tolerance)
            return false;
    }
    return true;
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
    if (p.generationMode == NanoGenerationMode::WulffConstruction)
        return std::max(p.wulffMaxRadius, 0.1f);

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
    if (p.generationMode == NanoGenerationMode::WulffConstruction)
    {
        const float radius = std::max(p.wulffMaxRadius, 0.1f);
        return {radius, radius, radius};
    }

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

WulffPreview computeWulffPreview(const Structure& reference, const NanoParams& params)
{
    WulffPreview preview;
    preview.center = buildClipCenter(reference, params);

    if (params.wulffPlanes.empty())
    {
        preview.message = "Add at least one Miller-index and surface-energy entry.";
        return preview;
    }
    if (!reference.hasUnitCell)
    {
        preview.message = "Wulff construction requires a periodic reference structure with a unit cell.";
        return preview;
    }

    std::string orientationError;
    const glm::mat3 orientation = buildOrientationMatrix(reference, params, &orientationError);
    if (!orientationError.empty())
    {
        preview.message = orientationError;
        return preview;
    }

    const glm::vec3 a = orientation * glm::vec3((float)reference.cellVectors[0][0], (float)reference.cellVectors[0][1], (float)reference.cellVectors[0][2]);
    const glm::vec3 b = orientation * glm::vec3((float)reference.cellVectors[1][0], (float)reference.cellVectors[1][1], (float)reference.cellVectors[1][2]);
    const glm::vec3 c = orientation * glm::vec3((float)reference.cellVectors[2][0], (float)reference.cellVectors[2][1], (float)reference.cellVectors[2][2]);
    const float volume = glm::dot(a, glm::cross(b, c));
    if (std::abs(volume) <= 1e-8f)
    {
        preview.message = "Reference structure has a degenerate lattice.";
        return preview;
    }

    float maxEnergy = 0.0f;
    for (const WulffPlaneInput& input : params.wulffPlanes)
    {
        if (input.surfaceEnergy <= 0.0f)
        {
            preview.message = "Surface energies must be positive.";
            return preview;
        }
        if (input.h == 0 && input.k == 0 && input.l == 0)
        {
            preview.message = "Miller indices [h k l] cannot all be zero.";
            return preview;
        }
        maxEnergy = std::max(maxEnergy, input.surfaceEnergy);
    }

    if (params.wulffMaxRadius <= 0.0f)
    {
        preview.message = "Maximum Wulff radius must be positive.";
        return preview;
    }

    const SymmetryRotationSet symmetry = collectSymmetryRotations(reference);
    if (!symmetry.success)
    {
        preview.message = symmetry.error;
        return preview;
    }

    glm::vec3 aStar(0.0f), bStar(0.0f), cStar(0.0f);
    buildReciprocalBasis(a, b, c, aStar, bStar, cStar);

    std::unordered_map<MillerKey, size_t, MillerKeyHash> planeIndexByHkl;
    for (size_t familyIndex = 0; familyIndex < params.wulffPlanes.size(); ++familyIndex)
    {
        const WulffPlaneInput& input = params.wulffPlanes[familyIndex];
        const glm::ivec3 base = reduceMillerIndices(glm::ivec3(input.h, input.k, input.l));
        const float distance = params.wulffMaxRadius * (input.surfaceEnergy / maxEnergy);

        for (const auto& rotation : symmetry.rotations)
        {
            const glm::ivec3 transformed = reduceMillerIndices(transformMillerIndices(rotation, base));
            if (isZeroMiller(transformed))
                continue;

            const glm::vec3 reciprocalNormal = (float)transformed.x * aStar
                                             + (float)transformed.y * bStar
                                             + (float)transformed.z * cStar;
            const float normalLength = glm::length(reciprocalNormal);
            if (normalLength <= 1e-8f)
                continue;

            const MillerKey key{transformed.x, transformed.y, transformed.z};
            const auto found = planeIndexByHkl.find(key);
            if (found == planeIndexByHkl.end())
            {
                WulffPlane plane;
                plane.familyIndex = (int)familyIndex;
                plane.h = transformed.x;
                plane.k = transformed.y;
                plane.l = transformed.z;
                plane.distance = distance;
                plane.normal = reciprocalNormal / normalLength;
                planeIndexByHkl.emplace(key, preview.planes.size());
                preview.planes.push_back(plane);
            }
            else if (distance < preview.planes[found->second].distance)
            {
                WulffPlane& plane = preview.planes[found->second];
                plane.familyIndex = (int)familyIndex;
                plane.distance = distance;
                plane.normal = reciprocalNormal / normalLength;
            }
        }
    }

    if (preview.planes.size() < 4)
    {
        preview.message = "Wulff construction produced too few symmetry-equivalent planes.";
        return preview;
    }

    const float intersectionTolerance = std::max(1e-3f, params.wulffMaxRadius * 1e-4f);
    const float faceTolerance = std::max(2e-3f, params.wulffMaxRadius * 2e-4f);

    std::vector<glm::vec3> vertices;
    for (size_t i = 0; i < preview.planes.size(); ++i)
    {
        for (size_t j = i + 1; j < preview.planes.size(); ++j)
        {
            for (size_t k = j + 1; k < preview.planes.size(); ++k)
            {
                glm::mat3 system(1.0f);
                system[0][0] = preview.planes[i].normal.x;
                system[1][0] = preview.planes[i].normal.y;
                system[2][0] = preview.planes[i].normal.z;
                system[0][1] = preview.planes[j].normal.x;
                system[1][1] = preview.planes[j].normal.y;
                system[2][1] = preview.planes[j].normal.z;
                system[0][2] = preview.planes[k].normal.x;
                system[1][2] = preview.planes[k].normal.y;
                system[2][2] = preview.planes[k].normal.z;

                const float determinant = glm::determinant(system);
                if (std::abs(determinant) <= 1e-6f)
                    continue;

                const glm::vec3 distances(preview.planes[i].distance,
                                          preview.planes[j].distance,
                                          preview.planes[k].distance);
                const glm::vec3 localPoint = glm::inverse(system) * distances;
                if (!isInsideWulffPolyhedron(localPoint, preview, intersectionTolerance))
                    continue;

                appendUniquePoint(vertices, preview.center + localPoint, intersectionTolerance);
            }
        }
    }

    if (vertices.size() < 4)
    {
        preview.message = "No valid Wulff polyhedron could be formed from the supplied planes and energies.";
        return preview;
    }

    preview.minPoint = glm::vec3(std::numeric_limits<float>::max());
    preview.maxPoint = glm::vec3(-std::numeric_limits<float>::max());
    preview.boundingRadius = 0.0f;
    preview.maxPlaneDistance = 0.0f;
    for (const glm::vec3& vertex : vertices)
    {
        preview.minPoint = glm::min(preview.minPoint, vertex);
        preview.maxPoint = glm::max(preview.maxPoint, vertex);
        preview.boundingRadius = std::max(preview.boundingRadius, glm::length(vertex - preview.center));
    }
    for (const WulffPlane& plane : preview.planes)
        preview.maxPlaneDistance = std::max(preview.maxPlaneDistance, plane.distance);

    for (const WulffPlane& plane : preview.planes)
    {
        std::vector<glm::vec3> faceVertices;
        for (const glm::vec3& vertex : vertices)
        {
            const float signedDistance = glm::dot(plane.normal, vertex - preview.center);
            if (std::abs(signedDistance - plane.distance) <= faceTolerance)
                faceVertices.push_back(vertex);
        }

        if (faceVertices.size() < 3)
            continue;

        sortFaceVertices(faceVertices, plane.normal);

        WulffFace face;
        face.familyIndex = plane.familyIndex;
        face.h = plane.h;
        face.k = plane.k;
        face.l = plane.l;
        face.distance = plane.distance;
        face.normal = plane.normal;
        face.vertices.swap(faceVertices);
        preview.faces.push_back(face);
    }

    if (preview.faces.empty())
    {
        preview.message = "Wulff construction generated vertices but no valid faces.";
        return preview;
    }

    preview.success = true;
    std::ostringstream message;
    message << "Wulff preview: " << preview.faces.size() << " faces from "
            << preview.planes.size() << " symmetry-expanded planes.";
    preview.message = message.str();
    return preview;
}

static float safeLen3(const glm::vec3& v)
{
    return std::sqrt(glm::dot(v, v));
}

static glm::mat3 rotationAligningDirectionToZ(const glm::vec3& fromDir)
{
    const glm::vec3 source = glm::normalize(fromDir);
    const glm::vec3 target(0.0f, 0.0f, 1.0f);
    const float dotST = glm::dot(source, target);

    if (dotST > 0.999999f)
        return glm::mat3(1.0f);

    if (dotST < -0.999999f)
    {
        glm::vec3 axis = glm::cross(source, glm::vec3(1.0f, 0.0f, 0.0f));
        if (glm::length(axis) < 1e-6f)
            axis = glm::cross(source, glm::vec3(0.0f, 1.0f, 0.0f));
        axis = glm::normalize(axis);
        return glm::mat3(glm::rotate(glm::mat4(1.0f), 3.14159265358979323846f, axis));
    }

    glm::vec3 axis = glm::cross(source, target);
    const float axisLen = glm::length(axis);
    if (axisLen < 1e-6f)
        return glm::mat3(1.0f);

    axis /= axisLen;
    const float clampedDot = std::max(-1.0f, std::min(1.0f, dotST));
    const float angle = std::acos(clampedDot);
    return glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
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
    result.mode       = params.generationMode;
    result.shape      = params.shape;
    result.inputAtoms = (int)reference.atoms.size();

    if (params.generationMode == NanoGenerationMode::WulffConstruction)
    {
        if (reference.atoms.empty())
        {
            result.message = "Reference structure has no atoms.";
            return result;
        }

        const WulffPreview preview = computeWulffPreview(reference, params);
        if (!preview.success)
        {
            result.message = preview.message;
            return result;
        }

        std::string orientationError;
        const glm::mat3 orientation = buildOrientationMatrix(reference, params, &orientationError);
        if (!orientationError.empty())
        {
            result.message = orientationError;
            return result;
        }

        const glm::vec3 center = preview.center;
        const float maxR = std::max(preview.boundingRadius, 0.1f);
        std::vector<AtomSite> generatedAtoms;

        if (reference.hasUnitCell)
        {
            const glm::vec3 a = orientation * glm::vec3((float)reference.cellVectors[0][0], (float)reference.cellVectors[0][1], (float)reference.cellVectors[0][2]);
            const glm::vec3 b = orientation * glm::vec3((float)reference.cellVectors[1][0], (float)reference.cellVectors[1][1], (float)reference.cellVectors[1][2]);
            const glm::vec3 c = orientation * glm::vec3((float)reference.cellVectors[2][0], (float)reference.cellVectors[2][1], (float)reference.cellVectors[2][2]);
            const float la = safeLen3(a);
            const float lb = safeLen3(b);
            const float lc = safeLen3(c);
            if (la < 1e-8f || lb < 1e-8f || lc < 1e-8f)
            {
                result.message = "Reference structure has degenerate lattice vectors.";
                return result;
            }

            const int kMaxReps = 40;
            bool clamped = false;
            auto safeRep = [&](float length) -> int {
                int count = (int)std::ceil(maxR / length) + 2;
                if (count > kMaxReps)
                {
                    clamped = true;
                    count = kMaxReps;
                }
                return count;
            };

            const int nA = params.autoReplicate ? safeRep(la) : std::max(1, params.repA);
            const int nB = params.autoReplicate ? safeRep(lb) : std::max(1, params.repB);
            const int nC = params.autoReplicate ? safeRep(lc) : std::max(1, params.repC);

            result.repA = nA;
            result.repB = nB;
            result.repC = nC;
            result.tilingUsed = true;
            result.repClamped = clamped;

            const long long total = (long long)(2 * nA + 1)
                                  * (long long)(2 * nB + 1)
                                  * (long long)(2 * nC + 1)
                                  * (long long)reference.atoms.size();
            if (total > 8000000LL)
            {
                std::ostringstream message;
                message << "Tiling would test " << total
                        << " atoms (limit 8M). Reduce the Wulff radius or use manual replication.";
                result.message = message.str();
                return result;
            }

            generatedAtoms.reserve((size_t)std::min(total, (long long)2000000));
            const float insideTolerance = std::max(1e-3f, preview.maxPlaneDistance * 1e-4f);
            for (int ia = -nA; ia <= nA; ++ia)
            for (int ib = -nB; ib <= nB; ++ib)
            for (int ic = -nC; ic <= nC; ++ic)
            {
                const glm::vec3 offset = (float)ia * a + (float)ib * b + (float)ic * c;
                for (const AtomSite& atom : reference.atoms)
                {
                    const glm::vec3 localPos((float)atom.x - center.x,
                                             (float)atom.y - center.y,
                                             (float)atom.z - center.z);
                    const glm::vec3 pos = center + (orientation * localPos) + offset;
                    if (!isInsideWulffPolyhedron(pos - center, preview, insideTolerance))
                        continue;

                    AtomSite out = atom;
                    out.x = (double)pos.x;
                    out.y = (double)pos.y;
                    out.z = (double)pos.z;
                    const int z = out.atomicNumber;
                    if (z >= 0 && z < (int)elementColors.size())
                    {
                        out.r = elementColors[z].r;
                        out.g = elementColors[z].g;
                        out.b = elementColors[z].b;
                    }
                    else
                    {
                        getDefaultElementColor(z, out.r, out.g, out.b);
                    }
                    generatedAtoms.push_back(out);
                }
            }
        }
        else
        {
            result.tilingUsed = false;
            generatedAtoms.reserve(reference.atoms.size());
            const float insideTolerance = std::max(1e-3f, preview.maxPlaneDistance * 1e-4f);
            for (const AtomSite& atom : reference.atoms)
            {
                const glm::vec3 localPos((float)atom.x - center.x,
                                         (float)atom.y - center.y,
                                         (float)atom.z - center.z);
                const glm::vec3 pos = center + (orientation * localPos);
                if (!isInsideWulffPolyhedron(pos - center, preview, insideTolerance))
                    continue;

                AtomSite out = atom;
                out.x = (double)pos.x;
                out.y = (double)pos.y;
                out.z = (double)pos.z;
                const int z = out.atomicNumber;
                if (z >= 0 && z < (int)elementColors.size())
                {
                    out.r = elementColors[z].r;
                    out.g = elementColors[z].g;
                    out.b = elementColors[z].b;
                }
                else
                {
                    getDefaultElementColor(z, out.r, out.g, out.b);
                }
                generatedAtoms.push_back(out);
            }
        }

        if (generatedAtoms.empty())
        {
            result.message = "No atoms fell inside the generated Wulff polyhedron.";
            return result;
        }

        {
            const float dedupTol = 0.01f;
            const float cellSize = 0.5f;
            const float invCell = 1.0f / cellSize;

            struct Vec3Hash {
                float invC;
                explicit Vec3Hash(float ic) : invC(ic) {}
                size_t operator()(const glm::ivec3& v) const {
                    size_t h = (size_t)(v.x * 73856093u);
                    h ^= (size_t)(v.y * 19349663u);
                    h ^= (size_t)(v.z * 83492791u);
                    return h;
                }
            };
            struct Vec3Eq {
                bool operator()(const glm::ivec3& left, const glm::ivec3& right) const {
                    return left.x == right.x && left.y == right.y && left.z == right.z;
                }
            };

            std::unordered_map<glm::ivec3, std::vector<size_t>, Vec3Hash, Vec3Eq> grid(
                generatedAtoms.size(), Vec3Hash(invCell));
            std::vector<AtomSite> unique;
            unique.reserve(generatedAtoms.size());
            const float dedupTol2 = dedupTol * dedupTol;

            for (size_t index = 0; index < generatedAtoms.size(); ++index)
            {
                const glm::vec3 position((float)generatedAtoms[index].x,
                                         (float)generatedAtoms[index].y,
                                         (float)generatedAtoms[index].z);
                const glm::ivec3 cell((int)std::floor(position.x * invCell),
                                      (int)std::floor(position.y * invCell),
                                      (int)std::floor(position.z * invCell));
                bool duplicate = false;
                for (int dx = -1; dx <= 1 && !duplicate; ++dx)
                for (int dy = -1; dy <= 1 && !duplicate; ++dy)
                for (int dz = -1; dz <= 1 && !duplicate; ++dz)
                {
                    const glm::ivec3 neighbor(cell.x + dx, cell.y + dy, cell.z + dz);
                    if (auto it = grid.find(neighbor); it != grid.end())
                    for (size_t existingIndex : it->second)
                    {
                        const glm::vec3 existing((float)unique[existingIndex].x,
                                                 (float)unique[existingIndex].y,
                                                 (float)unique[existingIndex].z);
                        const glm::vec3 delta = position - existing;
                        if (glm::dot(delta, delta) < dedupTol2)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                }

                if (!duplicate)
                {
                    const size_t uniqueIndex = unique.size();
                    unique.push_back(generatedAtoms[index]);
                    grid[cell].push_back(uniqueIndex);
                }
            }

            generatedAtoms.swap(unique);
        }

        result.estimatedDiameter = 2.0f * preview.maxPlaneDistance;
        result.outputAtoms = (int)generatedAtoms.size();
        result.wulffFamilyCount = (int)params.wulffPlanes.size();
        result.wulffFaceCount = (int)preview.faces.size();
        structure.atoms.swap(generatedAtoms);

        structure.grainColors.clear();
        structure.grainRegionIds.clear();
        structure.pbcBoundaryTol = 0.0f;
        structure.ipfLoadStatus.clear();

        if (params.setOutputCell)
        {
            const float pad = params.vacuumPadding;
            structure.hasUnitCell = true;
            structure.cellOffset = {{
                (double)(preview.minPoint.x - pad),
                (double)(preview.minPoint.y - pad),
                (double)(preview.minPoint.z - pad)
            }};
            structure.cellVectors = {{
                {{(double)(preview.maxPoint.x - preview.minPoint.x + 2.0f * pad), 0.0, 0.0}},
                {{0.0, (double)(preview.maxPoint.y - preview.minPoint.y + 2.0f * pad), 0.0}},
                {{0.0, 0.0, (double)(preview.maxPoint.z - preview.minPoint.z + 2.0f * pad)}}
            }};
        }
        else
        {
            structure.hasUnitCell = false;
        }

        result.success = true;
        std::ostringstream message;
        message << "Wulff nanocrystal built: " << result.outputAtoms << " atoms, "
                << result.wulffFaceCount << " faces.";
        result.message = message.str();
        return result;
    }

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

    const glm::vec3 center = buildClipCenter(reference, params);

    const float maxR = computeBoundingRadius(params);
    std::vector<AtomSite> generatedAtoms;

    std::string orientationError;
    const glm::mat3 orientation = buildOrientationMatrix(reference, params, &orientationError);
    if (!orientationError.empty())
    {
        result.message = orientationError;
        return result;
    }

    if (reference.hasUnitCell) {
        const auto& cv = reference.cellVectors;
        const glm::vec3 a = orientation * glm::vec3((float)cv[0][0], (float)cv[0][1], (float)cv[0][2]);
        const glm::vec3 b = orientation * glm::vec3((float)cv[1][0], (float)cv[1][1], (float)cv[1][2]);
        const glm::vec3 c = orientation * glm::vec3((float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
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
                const glm::vec3 localPos((float)atom.x - center.x,
                                         (float)atom.y - center.y,
                                         (float)atom.z - center.z);
                const glm::vec3 pos = center + (orientation * localPos) + offset;
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
            const glm::vec3 localPos((float)atom.x - center.x,
                                     (float)atom.y - center.y,
                                     (float)atom.z - center.z);
            const glm::vec3 pos = center + (orientation * localPos);
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
                if (auto it = grid.find(nc); it != grid.end())
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
