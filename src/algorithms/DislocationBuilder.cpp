#include "algorithms/DislocationBuilder.h"

#include "algorithms/StackingFaultBuilder.h"
#include "math/StructureMath.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

// Dislocation generation overview
// 1) Detect lattice family (FCC/HCP/BCC) from the input structure.
// 2) Build crystallographic directions (auto presets or user-provided hkl/uvw).
// 3) Construct a local orthonormal frame around the dislocation line.
// 4) Compute isotropic-elastic displacement fields (edge/screw/mixed) in that frame.
// 5) Apply optional shape/cutoff masks to limit where displacements are active.
// 6) Apply atom displacements (parallel on large systems), then recenter COM for edge/mixed.
// 7) Run validation checks and return a rich result for UI/CLI reporting.

namespace
{
constexpr float kPi = 3.14159265358979323846f;

struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u)
             ^ ((size_t)value.y * 19349663u)
             ^ ((size_t)value.z * 83492791u);
    }
};

DislocationLatticeFamily toDislocationFamily(StackingFaultFamily family)
{
    // Reuse stacking-fault local-environment classification for dislocation defaults.
    switch (family)
    {
    case StackingFaultFamily::Fcc: return DislocationLatticeFamily::Fcc;
    case StackingFaultFamily::Hcp: return DislocationLatticeFamily::Hcp;
    case StackingFaultFamily::Bcc: return DislocationLatticeFamily::Bcc;
    default: return DislocationLatticeFamily::Unknown;
    }
}

float nearestDistanceApprox(const Structure& structure);

float safeNorm(const glm::vec3& v)
{
    return glm::length(v);
}

glm::vec3 safeNormalize(const glm::vec3& v, const glm::vec3& fallback)
{
    const float n = safeNorm(v);
    if (n <= 1e-6f)
        return fallback;
    return v / n;
}

glm::vec3 directionFromUvw(const glm::mat3& cell, const glm::ivec3& uvw)
{
    // Convert crystallographic direction [u v w] to Cartesian using direct lattice basis.
    const glm::vec3 d = cell[0] * (float)uvw.x
                      + cell[1] * (float)uvw.y
                      + cell[2] * (float)uvw.z;
    return safeNormalize(d, glm::vec3(1.0f, 0.0f, 0.0f));
}

glm::vec3 planeNormalFromHkl(const glm::mat3& cell, const glm::ivec3& hkl)
{
    // Convert plane normal (h k l) using reciprocal basis vectors a*, b*, c*.
    const float volume = glm::dot(cell[0], glm::cross(cell[1], cell[2]));
    if (std::abs(volume) <= 1e-8f)
        return glm::vec3(0.0f, 0.0f, 1.0f);

    const glm::vec3 astar = glm::cross(cell[1], cell[2]) / volume;
    const glm::vec3 bstar = glm::cross(cell[2], cell[0]) / volume;
    const glm::vec3 cstar = glm::cross(cell[0], cell[1]) / volume;
    return safeNormalize(astar * (float)hkl.x + bstar * (float)hkl.y + cstar * (float)hkl.z,
                         glm::vec3(0.0f, 0.0f, 1.0f));
}

void getAutoLatticeVectors(DislocationLatticeFamily family,
                           glm::ivec3& planeHkl,
                           glm::ivec3& burgersUvw,
                           glm::ivec3& lineUvw)
{
    // Typical textbook defaults used when the user enables automatic directions.
    switch (family)
    {
    case DislocationLatticeFamily::Fcc:
        planeHkl = glm::ivec3(1, 1, 1);
        burgersUvw = glm::ivec3(1, -1, 0);
        lineUvw = glm::ivec3(1, 1, -2);
        break;
    case DislocationLatticeFamily::Hcp:
        planeHkl = glm::ivec3(0, 0, 1);
        burgersUvw = glm::ivec3(1, -1, 0);
        lineUvw = glm::ivec3(1, 1, 0);
        break;
    case DislocationLatticeFamily::Bcc:
        planeHkl = glm::ivec3(1, 1, 0);
        burgersUvw = glm::ivec3(1, 1, 1);
        lineUvw = glm::ivec3(1, -1, 0);
        break;
    default:
        planeHkl = glm::ivec3(0, 0, 1);
        burgersUvw = glm::ivec3(1, 0, 0);
        lineUvw = glm::ivec3(0, 1, 0);
        break;
    }
}

float estimateLatticeParameter(const Structure& structure, DislocationLatticeFamily family)
{
    // Prefer nearest-neighbor distance because supercell vectors can be very large
    // and are not suitable for recovering primitive lattice parameter a.
    const float dnn = nearestDistanceApprox(structure);
    if (dnn > 1e-4f)
    {
        switch (family)
        {
        case DislocationLatticeFamily::Fcc:
            return std::max(1e-3f, dnn * std::sqrt(2.0f));
        case DislocationLatticeFamily::Bcc:
            return std::max(1e-3f, dnn * (2.0f / std::sqrt(3.0f)));
        case DislocationLatticeFamily::Hcp:
            return std::max(1e-3f, dnn);
        default:
            return std::max(1e-3f, dnn);
        }
    }

    if (!structure.hasUnitCell)
        return 1.0f;

    const glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
    const glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
    const glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);
    return std::max(1e-3f, (glm::length(a) + glm::length(b) + glm::length(c)) / 3.0f);
}

float baseBurgersMagnitude(DislocationLatticeFamily family, float latticeA)
{
    // Canonical full-dislocation magnitudes used as defaults per lattice type.
    switch (family)
    {
    case DislocationLatticeFamily::Fcc:
        return latticeA / std::sqrt(2.0f);
    case DislocationLatticeFamily::Bcc:
        return latticeA * std::sqrt(3.0f) * 0.5f;
    case DislocationLatticeFamily::Hcp:
        return latticeA;
    default:
        return latticeA;
    }
}

bool pointInPolygon(const glm::vec2& point, const std::vector<glm::vec2>& polygon)
{
    // Ray-casting 2D point-in-polygon test for freeform region masking.
    if (polygon.size() < 3)
        return false;

    bool inside = false;
    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
    {
        const glm::vec2& pi = polygon[i];
        const glm::vec2& pj = polygon[j];

        const bool intersect = ((pi.y > point.y) != (pj.y > point.y))
                            && (point.x < (pj.x - pi.x) * (point.y - pi.y)
                                            / std::max(1e-8f, (pj.y - pi.y)) + pi.x);
        if (intersect)
            inside = !inside;
    }
    return inside;
}

float shapeMask(const DislocationParams& params, float x, float y, float z)
{
    // Geometry gate in local dislocation coordinates; returns 1 inside, 0 outside.
    switch (params.shape)
    {
    case DislocationShape::HalfPlane:
        return (y >= 0.0f) ? 1.0f : 0.0f;
    case DislocationShape::Cylinder:
        return ((x * x + y * y) <= params.cylinderRadius * params.cylinderRadius) ? 1.0f : 0.0f;
    case DislocationShape::Sphere:
        return ((x * x + y * y + z * z) <= params.sphereRadius * params.sphereRadius) ? 1.0f : 0.0f;
    case DislocationShape::Ellipsoid:
    {
        const float rx = std::max(1e-4f, params.ellipsoidRadii.x);
        const float ry = std::max(1e-4f, params.ellipsoidRadii.y);
        const float rz = std::max(1e-4f, params.ellipsoidRadii.z);
        const float q = (x * x) / (rx * rx) + (y * y) / (ry * ry) + (z * z) / (rz * rz);
        return (q <= 1.0f) ? 1.0f : 0.0f;
    }
    case DislocationShape::Freeform2D:
        return pointInPolygon(glm::vec2(x, y), params.freeformPoints) ? 1.0f : 0.0f;
    default:
        return 0.0f;
    }
}

float attenuationMask(const DislocationParams& params, float radialDistance)
{
    // Smooth radial attenuation to avoid hard discontinuities at cutoff boundary.
    if (params.cutoffRadius <= 0.0f)
        return 1.0f;

    const float r = radialDistance / std::max(1e-4f, params.cutoffRadius);
    const float t = r * r;
    return std::exp(-(t * t));
}

float nearestDistanceApprox(const Structure& structure)
{
    // Uniform-grid neighbor search used as a fast approximate nearest-neighbor estimate.
    if (structure.atoms.size() < 2)
        return 0.0f;

    constexpr float kCell = 2.5f;
    std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
    grid.reserve(structure.atoms.size());

    float best = std::numeric_limits<float>::max();

    auto gridKey = [](const glm::vec3& p) -> glm::ivec3 {
        return glm::ivec3((int)std::floor(p.x / kCell),
                          (int)std::floor(p.y / kCell),
                          (int)std::floor(p.z / kCell));
    };

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const glm::vec3 p((float)structure.atoms[i].x,
                          (float)structure.atoms[i].y,
                          (float)structure.atoms[i].z);
        const glm::ivec3 k = gridKey(p);

        for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
        {
            const glm::ivec3 nk(k.x + dx, k.y + dy, k.z + dz);
            auto it = grid.find(nk);
            if (it == grid.end())
                continue;

            for (int j : it->second)
            {
                const glm::vec3 q((float)structure.atoms[j].x,
                                  (float)structure.atoms[j].y,
                                  (float)structure.atoms[j].z);
                best = std::min(best, glm::length(p - q));
            }
        }

        grid[k].push_back(i);
    }

    return std::isfinite(best) ? best : 0.0f;
}

glm::vec3 structureCenter(const Structure& structure)
{
    // Arithmetic center used as default dislocation-line anchor when no fractional point is requested.
    if (structure.atoms.empty())
        return glm::vec3(0.0f);

    glm::vec3 center(0.0f);
    for (const AtomSite& atom : structure.atoms)
        center += glm::vec3((float)atom.x, (float)atom.y, (float)atom.z);
    return center / (float)structure.atoms.size();
}

void getLinePoint(const Structure& structure,
                  const DislocationParams& params,
                  const glm::mat3& cell,
                  glm::vec3& linePoint)
{
    // Resolve dislocation line point either in fractional cell coordinates or from scene center.
    if (params.useFractionalLinePoint && structure.hasUnitCell)
    {
        const glm::vec3 origin((float)structure.cellOffset[0],
                               (float)structure.cellOffset[1],
                               (float)structure.cellOffset[2]);
        linePoint = origin + cell * params.linePointFractional + params.linePointCartesianOffset;
        return;
    }

    linePoint = structureCenter(structure) + params.linePointCartesianOffset;
}

DislocationValidationReport buildValidationReport(const Structure& before,
                                                  const Structure& after,
                                                  DislocationLatticeFamily familyBefore)
{
    // Post-process sanity checks: finite coords, nearest distance, family persistence, displacement stats.
    DislocationValidationReport report;
    report.familyBefore = familyBefore;

    report.finiteCoordinates = true;
    for (const AtomSite& atom : after.atoms)
    {
        if (!std::isfinite(atom.x) || !std::isfinite(atom.y) || !std::isfinite(atom.z))
        {
            report.finiteCoordinates = false;
            break;
        }
    }

    report.minInteratomicDistance = nearestDistanceApprox(after);
    report.minDistanceSafe = report.minInteratomicDistance >= 0.20f;

    const DislocationDetectionResult afterDetection = detectDislocationLattice(after, true);
    report.familyAfter = afterDetection.family;
    report.latticeFamilyPreserved = (familyBefore == DislocationLatticeFamily::Unknown)
        ? (report.familyAfter != DislocationLatticeFamily::Unknown)
        : (familyBefore == report.familyAfter);

    float sumSq = 0.0f;
    float maxDisp = 0.0f;
    const size_t n = std::min(before.atoms.size(), after.atoms.size());
    for (size_t i = 0; i < n; ++i)
    {
        const glm::vec3 p0((float)before.atoms[i].x, (float)before.atoms[i].y, (float)before.atoms[i].z);
        const glm::vec3 p1((float)after.atoms[i].x, (float)after.atoms[i].y, (float)after.atoms[i].z);
        const float d = glm::length(p1 - p0);
        sumSq += d * d;
        maxDisp = std::max(maxDisp, d);
    }

    report.rmsDisplacement = (n > 0) ? std::sqrt(sumSq / (float)n) : 0.0f;
    report.maxDisplacement = maxDisp;

    report.passed = report.finiteCoordinates
                 && report.minDistanceSafe
                 && report.latticeFamilyPreserved;

    std::ostringstream oss;
    if (report.passed)
    {
        oss << "Validation passed. Family preserved, coordinates finite, and min distance = "
            << report.minInteratomicDistance << " A.";
    }
    else
    {
        oss << "Validation warning: ";
        bool emitted = false;
        if (!report.finiteCoordinates)
        {
            oss << "non-finite coordinates detected";
            emitted = true;
        }
        if (!report.minDistanceSafe)
        {
            if (emitted) oss << "; ";
            oss << "very short interatomic distance (" << report.minInteratomicDistance << " A)";
            emitted = true;
        }
        if (!report.latticeFamilyPreserved)
        {
            if (emitted) oss << "; ";
            oss << "lattice family changed from "
                << dislocationLatticeFamilyName(report.familyBefore)
                << " to "
                << dislocationLatticeFamilyName(report.familyAfter);
        }
        oss << ".";
    }

    report.message = oss.str();
    return report;
}

struct DislocationKernelContext
{
    glm::vec3 linePoint = glm::vec3(0.0f);
    glm::vec3 lineDir = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 e1 = glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 e2 = glm::vec3(0.0f, 1.0f, 0.0f);
    float prefEdge = 0.0f;
    float prefScrew = 0.0f;
    float poisson = 0.33f;
    float core = 1.0f;
    float lineHalfLength = 1.0e6f;
};

struct WorkerAccum
{
    int shifted = 0;
    double sumDx = 0.0;
    double sumDy = 0.0;
    double sumDz = 0.0;
};

glm::vec3 computeDisplacement(const DislocationParams& params,
                              const DislocationKernelContext& ctx,
                              const glm::vec3& position)
{
    // Evaluate isotropic elasticity displacement at one atom in local coordinates.
    const glm::vec3 rel = position - ctx.linePoint;
    const float z = glm::dot(rel, ctx.lineDir);
    if (std::abs(z) > ctx.lineHalfLength)
        return glm::vec3(0.0f);

    const float x = glm::dot(rel, ctx.e1);
    const float y = glm::dot(rel, ctx.e2);
    const float r2 = x * x + y * y + ctx.core * ctx.core;

    const float sh = shapeMask(params, x, y, z);
    if (sh <= 0.0f)
        return glm::vec3(0.0f);

    const float attenuation = attenuationMask(params, std::sqrt(x * x + y * y));
    const float weight = sh * attenuation;
    if (weight <= 1e-6f)
        return glm::vec3(0.0f);

    float ux = 0.0f;
    float uy = 0.0f;
    float uz = 0.0f;

    if (std::abs(ctx.prefEdge) > 0.0f)
    {
        // Isotropic edge displacement in local (e1,e2) coordinates.
        const float atanTerm = std::atan2(y, x);
        const float inv2OneMinusNuR2 = 1.0f / std::max(1e-7f, 2.0f * (1.0f - ctx.poisson) * r2);
        ux = ctx.prefEdge * (atanTerm + (x * y) * inv2OneMinusNuR2);

        const float logTerm = std::log(std::max(1e-7f, r2));
        const float q = std::max(1e-7f, 4.0f * (1.0f - ctx.poisson));
        uy = -ctx.prefEdge * (((1.0f - 2.0f * ctx.poisson) / q) * logTerm
                              + (x * x - y * y) / (q * std::max(1e-7f, r2)));
    }

    if (std::abs(ctx.prefScrew) > 0.0f)
    {
        // Isotropic screw displacement along line direction.
        uz = ctx.prefScrew * std::atan2(y, x);
    }

    return weight * (ux * ctx.e1 + uy * ctx.e2 + uz * ctx.lineDir);
}
}

const char* dislocationLatticeFamilyName(DislocationLatticeFamily family)
{
    // User-facing label helper for dialogs/log messages.
    switch (family)
    {
    case DislocationLatticeFamily::Fcc: return "FCC-like";
    case DislocationLatticeFamily::Hcp: return "HCP-like";
    case DislocationLatticeFamily::Bcc: return "BCC-like";
    default: return "Unknown";
    }
}

DislocationDetectionResult detectDislocationLattice(const Structure& structure,
                                                    bool usePbcForDetection)
{
    // Wrapper around existing local-order detection pipeline.
    const StackingFaultDetectionResult sf = detectStackingFaultFamily(structure, usePbcForDetection);

    DislocationDetectionResult out;
    out.success = sf.success;
    out.message = sf.message;
    out.family = toDislocationFamily(sf.family);
    out.fccCount = sf.fccCount;
    out.hcpCount = sf.hcpCount;
    out.bccCount = sf.bccCount;
    out.recognizedCount = sf.recognizedCount;
    return out;
}

DislocationResult buildDislocation(const Structure& base,
                                   const DislocationParams& params)
{
    // Main pipeline: detect, configure local frame, apply displacement field, validate result.
    DislocationResult result;

    if (base.atoms.empty())
    {
        result.message = "Input structure has no atoms.";
        return result;
    }

    if (!base.hasUnitCell)
    {
        result.message = "Dislocation insertion requires a structure with a unit cell.";
        return result;
    }

    const DislocationDetectionResult detection = detectDislocationLattice(base, true);
    result.family = detection.family;

    if (params.autoDirections && detection.family == DislocationLatticeFamily::Unknown)
    {
        result.message = "Could not identify lattice family (FCC/HCP/BCC). Disable auto vectors and set hkl/uvw manually.";
        return result;
    }

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(base, cell, invCell))
    {
        result.message = "Could not build cell matrices from the structure.";
        return result;
    }

    glm::ivec3 planeHkl = params.planeHkl;
    glm::ivec3 burgersUvw = params.burgersUvw;
    glm::ivec3 lineUvw = params.lineUvw;
    if (params.autoDirections)
        getAutoLatticeVectors(detection.family, planeHkl, burgersUvw, lineUvw);

    const glm::vec3 planeNormal = planeNormalFromHkl(cell, planeHkl);
    const glm::vec3 burgersDir = directionFromUvw(cell, burgersUvw);

    glm::vec3 lineDir;
    if (params.character == DislocationCharacter::Screw)
    {
        lineDir = burgersDir;
    }
    else
    {
        if (params.autoDirections)
            lineDir = safeNormalize(glm::cross(planeNormal, burgersDir), directionFromUvw(cell, lineUvw));
        else
            lineDir = directionFromUvw(cell, lineUvw);
    }

    const glm::vec3 e2 = safeNormalize(glm::cross(lineDir, burgersDir), planeNormal);

    if (params.character != DislocationCharacter::Screw && glm::length(glm::cross(lineDir, burgersDir)) <= 1e-6f)
    {
        // Edge and mixed modes require distinct line and Burgers directions.
        result.message = "For edge/mixed dislocations, line direction must not be parallel to Burgers direction.";
        return result;
    }

    const glm::vec3 e1 = safeNormalize(glm::cross(e2, lineDir), burgersDir);

    const float latticeA = estimateLatticeParameter(base, detection.family);
    float bmag;
    if (params.burgersOverrideMagnitude > 0.0f)
        bmag = std::max(1e-5f, params.burgersScale * params.burgersOverrideMagnitude);
    else
        bmag = std::max(1e-5f, params.burgersScale * baseBurgersMagnitude(detection.family, latticeA));

    const float poisson = std::clamp(params.poissonRatio, 0.02f, 0.49f);
    const float core = std::max(1e-4f, params.coreRadius);

    float edgeFraction = 0.0f;
    float screwFraction = 0.0f;
    if (params.character == DislocationCharacter::Edge)
    {
        edgeFraction = 1.0f;
    }
    else if (params.character == DislocationCharacter::Screw)
    {
        screwFraction = 1.0f;
    }
    else
    {
        const float angle = params.mixedCharacterAngleDeg * (kPi / 180.0f);
        edgeFraction = std::cos(angle);
        screwFraction = std::sin(angle);
    }

    glm::vec3 linePoint(0.0f);
    getLinePoint(base, params, cell, linePoint);

    // Build an orthonormal local frame around the dislocation line:
    // e1: edge displacement axis, e2: in-plane normal to e1 and lineDir,
    // lineDir: dislocation line direction (z-like axis in local formulas).
    DislocationKernelContext kernel;
    kernel.linePoint = linePoint;
    kernel.lineDir = lineDir;
    kernel.e1 = e1;
    kernel.e2 = e2;
    kernel.poisson = poisson;
    kernel.core = core;
    kernel.lineHalfLength = params.lineHalfLength;
    kernel.prefEdge = (edgeFraction * bmag) / (2.0f * kPi);
    kernel.prefScrew = (screwFraction * bmag) / (2.0f * kPi);

    Structure output = base;
    const size_t atomCount = output.atoms.size();
    const unsigned hwThreads = std::max(1u, std::thread::hardware_concurrency());
    // Use multithreading only when there is enough work to amortize thread overhead.
    const unsigned workerCount = (atomCount >= 50000)
        ? std::min(8u, hwThreads)
        : 1u;

    std::vector<WorkerAccum> accum(workerCount);
    auto applyRange = [&](size_t begin, size_t end, WorkerAccum& local) {
        for (size_t i = begin; i < end; ++i)
        {
            AtomSite& atom = output.atoms[i];
            const glm::vec3 p((float)atom.x, (float)atom.y, (float)atom.z);
            const glm::vec3 disp = computeDisplacement(params, kernel, p);
            const float disp2 = glm::dot(disp, disp);
            if (disp2 <= 1e-14f)
                continue;

            atom.x += (double)disp.x;
            atom.y += (double)disp.y;
            atom.z += (double)disp.z;
            local.sumDx += (double)disp.x;
            local.sumDy += (double)disp.y;
            local.sumDz += (double)disp.z;
            local.shifted++;
        }
    };

    if (workerCount == 1)
    {
        applyRange(0, atomCount, accum[0]);
    }
    else
    {
        std::vector<std::thread> workers;
        workers.reserve(workerCount);
        const size_t chunk = (atomCount + workerCount - 1) / workerCount;
        for (unsigned w = 0; w < workerCount; ++w)
        {
            const size_t begin = (size_t)w * chunk;
            const size_t end = std::min(atomCount, begin + chunk);
            if (begin >= end)
                break;
            workers.emplace_back([&, begin, end, w]() {
                applyRange(begin, end, accum[w]);
            });
        }
        for (std::thread& worker : workers)
            worker.join();
    }

    int shifted = 0;
    double sumDx = 0.0;
    double sumDy = 0.0;
    double sumDz = 0.0;
    for (const WorkerAccum& local : accum)
    {
        shifted += local.shifted;
        sumDx += local.sumDx;
        sumDy += local.sumDy;
        sumDz += local.sumDz;
    }

    // Center-of-mass conservation for edge and mixed dislocations: the
    // logarithmic term in u_y introduces a non-zero mean displacement that
    // would translate the whole crystal.  Subtract the mean so the result is
    // consistent with codes that enforce this (e.g. Atomsk).
    // Not applied to pure screw (atan2 is antisymmetric, mean ≈ 0 by symmetry).
    if (shifted > 0 && output.atoms.size() > 0
        && params.character != DislocationCharacter::Screw)
    {
        const size_t n = output.atoms.size();
        const double invN = 1.0 / (double)n;
        const double meanDx = sumDx * invN;
        const double meanDy = sumDy * invN;
        const double meanDz = sumDz * invN;

        auto recenterRange = [&](size_t begin, size_t end) {
            for (size_t i = begin; i < end; ++i)
            {
                AtomSite& atom = output.atoms[i];
                atom.x -= meanDx;
                atom.y -= meanDy;
                atom.z -= meanDz;
            }
        };

        if (workerCount == 1)
        {
            recenterRange(0, n);
        }
        else
        {
            std::vector<std::thread> workers;
            workers.reserve(workerCount);
            const size_t chunk = (n + workerCount - 1) / workerCount;
            for (unsigned w = 0; w < workerCount; ++w)
            {
                const size_t begin = (size_t)w * chunk;
                const size_t end = std::min(n, begin + chunk);
                if (begin >= end)
                    break;
                workers.emplace_back([&, begin, end]() {
                    recenterRange(begin, end);
                });
            }
            for (std::thread& worker : workers)
                worker.join();
        }
    }

    result.success = shifted > 0;
    result.shiftedAtomCount = shifted;
    result.burgersMagnitude = bmag;
    result.linePoint = linePoint;
    result.lineDirection = lineDir;
    result.burgersDirection = burgersDir;
    result.e1 = e1;
    result.e2 = e2;

    // Structure-derived dislocation boundary overlay.
    //
    // Algorithm (works for any crystal, any shape mask, any dislocation character):
    //   1. Compare base atom positions with displaced output atom positions.
    //   2. Collect atoms where |Δr| exceeds a small threshold (fraction of |b|).
    //   3. Project those atom positions onto the slip plane (e1–lineDir axes).
    //   4. Compute the 2D convex hull of the projected points (Jarvis march).
    //   5. Map hull vertices back to 3D world space.
    //
    // This gives the actual boundary of the displaced region as seen in TEM,
    // independent of any geometric approximation or shape-specific parameter.
    {
        result.loopPoints.clear();

        if (params.shape == DislocationShape::Freeform2D)
        {
            // User explicitly specified the polygon boundary — use it verbatim.
            for (const auto& pt : params.freeformPoints)
                result.loopPoints.push_back(linePoint + pt.x * e1 + pt.y * e2);
        }
        else
        {
            // Threshold: flag an atom as displaced if |Δr| > 5 % of |b| (min 0.05 Å).
            const float threshold = std::max(0.05f, 0.05f * bmag);

            // Collect 2D slip-plane projections of all significantly displaced atoms.
            // Slip-plane axes: e1 (Burgers direction) and lineDir (line direction).
            const size_t natoms = std::min(base.atoms.size(), output.atoms.size());
            std::vector<glm::vec2> pts2d;
            pts2d.reserve(natoms / 4 + 1);

            for (size_t i = 0; i < natoms; ++i)
            {
                const glm::vec3 bpos(static_cast<float>(base.atoms[i].x),
                                     static_cast<float>(base.atoms[i].y),
                                     static_cast<float>(base.atoms[i].z));
                const glm::vec3 opos(static_cast<float>(output.atoms[i].x),
                                     static_cast<float>(output.atoms[i].y),
                                     static_cast<float>(output.atoms[i].z));
                if (glm::length(opos - bpos) < threshold)
                    continue;
                const glm::vec3 rel = bpos - linePoint;
                pts2d.push_back({ glm::dot(rel, e1), glm::dot(rel, lineDir) });
            }

            if (pts2d.size() >= 3)
            {
                // 2D convex hull via gift wrapping (Jarvis march).
                // Produces hull in counter-clockwise order; GL_LINE_LOOP closes it.
                const int n = static_cast<int>(pts2d.size());

                // Start from the leftmost point (tie-break: lowest y).
                int start = 0;
                for (int i = 1; i < n; ++i)
                    if (pts2d[i].x < pts2d[start].x ||
                        (pts2d[i].x == pts2d[start].x && pts2d[i].y < pts2d[start].y))
                        start = i;

                std::vector<int> hullIdx;
                hullIdx.reserve(32);
                int cur = start;
                do {
                    hullIdx.push_back(cur);
                    // Initial candidate: first index that isn't cur.
                    int nxt = (cur == 0) ? 1 : 0;
                    for (int i = 0; i < n; ++i)
                    {
                        if (i == nxt) continue;
                        const glm::vec2 a = pts2d[nxt] - pts2d[cur];
                        const glm::vec2 b = pts2d[i]   - pts2d[cur];
                        const float cross = a.x * b.y - a.y * b.x;
                        // cross > 0: i is more counter-clockwise than current nxt.
                        // cross == 0: collinear — keep the farther point.
                        if (cross > 0.0f ||
                            (std::abs(cross) < 1e-9f &&
                             glm::dot(b, b) > glm::dot(a, a)))
                            nxt = i;
                    }
                    cur = nxt;
                } while (cur != start && static_cast<int>(hullIdx.size()) <= n);

                for (int idx : hullIdx)
                    result.loopPoints.push_back(
                        linePoint + pts2d[idx].x * e1 + pts2d[idx].y * lineDir);
            }
            else
            {
                // Degenerate: fewer than 3 displaced atoms — just emit their positions.
                for (const auto& p : pts2d)
                    result.loopPoints.push_back(linePoint + p.x * e1 + p.y * lineDir);
            }
        }
    }

    result.output = std::move(output);

    if (!result.success)
    {
        result.message = "No atoms were shifted. Try increasing the shape radius/cutoff or moving the line position.";
        return result;
    }

    result.validation = buildValidationReport(base, result.output, detection.family);

    std::ostringstream oss;
    oss << "Inserted " << ((params.character == DislocationCharacter::Edge) ? "edge" :
                             (params.character == DislocationCharacter::Screw) ? "screw" : "mixed")
        << " dislocation on " << dislocationLatticeFamilyName(detection.family)
        << ". Shifted atoms: " << shifted
        << ", |b| = " << bmag << " A. "
        << result.validation.message;

    result.message = oss.str();
    return result;
}
