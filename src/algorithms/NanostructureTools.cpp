#include "algorithms/NanostructureTools.h"
#include "util/ElementData.h"

#include <cstdio>
#include <cmath>
#include <glm/glm.hpp>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace atomforge
{

namespace
{
int resolveAtomicNumber(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    return 0;
}
}

NanowireResult buildNanowire(const Structure& source, const NanowireParams& params)
{
    NanowireResult result;
    if (!source.hasUnitCell || source.atoms.empty())
    {
        result.message = "Structure needs a unit cell and at least one atom.";
        return result;
    }
    if (params.axis < 0 || params.axis > 2)
    {
        result.message = "axis must be 0 (a), 1 (b), or 2 (c).";
        return result;
    }
    if (params.radius <= 0.0)
    {
        result.message = "radius must be positive.";
        return result;
    }
    if (params.axisRepeats < 1)
    {
        result.message = "axisRepeats must be at least 1.";
        return result;
    }

    const int axis = params.axis;
    const int inPlane[2] = {(axis + 1) % 3, (axis + 2) % 3};

    const glm::dvec3 cellVecs[3] = {
        glm::dvec3(source.cellVectors[0][0], source.cellVectors[0][1], source.cellVectors[0][2]),
        glm::dvec3(source.cellVectors[1][0], source.cellVectors[1][1], source.cellVectors[1][2]),
        glm::dvec3(source.cellVectors[2][0], source.cellVectors[2][1], source.cellVectors[2][2])};
    const double axisLen = glm::length(cellVecs[axis]);
    if (axisLen < 1e-9)
    {
        result.message = "The wire-axis cell vector is degenerate.";
        return result;
    }
    const glm::dvec3 axisDir = cellVecs[axis] / axisLen;

    // Orthonormal frame spanning the plane perpendicular to the wire axis,
    // used for both the radial cut and the replacement (orthogonal, vacuum
    // padded) cell -- Cartesian, so it is correct even for oblique cells.
    glm::dvec3 e1 = std::abs(glm::dot(axisDir, glm::dvec3(1, 0, 0))) < 0.9
        ? glm::cross(axisDir, glm::dvec3(1, 0, 0))
        : glm::cross(axisDir, glm::dvec3(0, 1, 0));
    e1 = glm::normalize(e1);
    const glm::dvec3 e2 = glm::normalize(glm::cross(axisDir, e1));

    // Bounding radius large enough to guarantee the requested cross-section
    // is fully covered by repeated copies of the source cell.
    const double coverage = params.radius + params.vacuum + axisLen;
    const int reps1 = (int)std::ceil(coverage / std::max(1e-6, glm::length(cellVecs[inPlane[0]]))) + 1;
    const int reps2 = (int)std::ceil(coverage / std::max(1e-6, glm::length(cellVecs[inPlane[1]]))) + 1;

    // Centroid of the source cell's own atoms, used as the wire's central
    // axis line origin so the cut is centered on the material, not the cell
    // origin (which may sit on a corner with no atoms nearby).
    glm::dvec3 centroid(0.0);
    for (const auto& atom : source.atoms) centroid += glm::dvec3(atom.x, atom.y, atom.z);
    centroid /= (double)source.atoms.size();

    const double faceAngleStep = params.sides >= 3 ? (2.0 * M_PI / params.sides) : 0.0;

    std::vector<AtomSite> kept;
    for (int i1 = -reps1; i1 <= reps1; ++i1)
    for (int i2 = -reps2; i2 <= reps2; ++i2)
    for (int ia = 0; ia < params.axisRepeats; ++ia)
    {
        const glm::dvec3 shift = (double)i1 * cellVecs[inPlane[0]] + (double)i2 * cellVecs[inPlane[1]]
            + (double)ia * cellVecs[axis];
        for (const auto& atom : source.atoms)
        {
            const glm::dvec3 pos = glm::dvec3(atom.x, atom.y, atom.z) + shift;
            const glm::dvec3 rel = pos - centroid;
            const double alongAxis = glm::dot(rel, axisDir);
            // Keep exactly one axis-length window so axisRepeats controls
            // wire length without also re-cutting the cross-section window.
            if (alongAxis < -1e-6 || alongAxis >= axisLen * params.axisRepeats - 1e-6) continue;

            const double c1 = glm::dot(rel, e1);
            const double c2 = glm::dot(rel, e2);

            bool inside;
            if (params.sides >= 3)
            {
                inside = true;
                for (int k = 0; k < params.sides && inside; ++k)
                {
                    const double faceAngle = k * faceAngleStep;
                    const double proj = c1 * std::cos(faceAngle) + c2 * std::sin(faceAngle);
                    if (proj > params.radius) inside = false;
                }
            }
            else
            {
                inside = std::sqrt(c1 * c1 + c2 * c2) <= params.radius;
            }
            if (!inside) continue;

            // Re-express the kept atom in the new (e1, e2, axisDir) frame,
            // shifting the two non-periodic directions by (radius+vacuum)
            // so the cross-section is centered inside the new box instead
            // of sitting on the original, now-irrelevant lattice origin.
            AtomSite placed = atom;
            const glm::dvec3 finalPos = alongAxis * axisDir
                + (c1 + params.radius + params.vacuum) * e1
                + (c2 + params.radius + params.vacuum) * e2;
            placed.x = finalPos.x; placed.y = finalPos.y; placed.z = finalPos.z;
            kept.push_back(placed);
        }
    }

    if (kept.empty())
    {
        result.message = "No atoms fell inside the requested cross-section; increase radius.";
        return result;
    }

    Structure out;
    out.hasUnitCell = true;
    out.atoms = std::move(kept);
    const double boxSize = 2.0 * (params.radius + params.vacuum);
    glm::dvec3 newCell[3];
    newCell[inPlane[0]] = e1 * boxSize;
    newCell[inPlane[1]] = e2 * boxSize;
    newCell[axis] = axisDir * axisLen * (double)params.axisRepeats;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out.cellVectors[i][j] = newCell[i][j];

    result.structure = out;
    result.atomCount = (int)out.atoms.size();
    result.success = true;
    result.message = "Built a " + std::string(params.sides >= 3 ? std::to_string(params.sides) + "-sided" : "circular")
        + " nanowire with " + std::to_string(result.atomCount) + " atoms, periodic along axis "
        + std::to_string(axis) + ".";
    return result;
}

CoreShellResult applyCoreShell(const Structure& source, const CoreShellParams& params)
{
    CoreShellResult result;
    if (source.atoms.empty())
    {
        result.message = "Structure has no atoms.";
        return result;
    }
    if (params.coreElement.empty() || params.shellElement.empty())
    {
        result.message = "Both coreElement and shellElement must be specified.";
        return result;
    }
    const int coreZ = resolveAtomicNumber(params.coreElement);
    const int shellZ = resolveAtomicNumber(params.shellElement);
    if (coreZ == 0 || shellZ == 0)
    {
        result.message = "Unknown element symbol ('" + params.coreElement + "' or '" + params.shellElement + "').";
        return result;
    }

    glm::dvec3 center(params.center[0], params.center[1], params.center[2]);
    if (params.useCentroid)
    {
        center = glm::dvec3(0.0);
        for (const auto& atom : source.atoms) center += glm::dvec3(atom.x, atom.y, atom.z);
        center /= (double)source.atoms.size();
    }

    Structure out = source;
    for (auto& atom : out.atoms)
    {
        const double dist = glm::length(glm::dvec3(atom.x, atom.y, atom.z) - center);
        const bool isCore = dist < params.coreRadius;
        const int z = isCore ? coreZ : shellZ;
        atom.atomicNumber = z;
        atom.symbol = elementSymbol(z);
        getDefaultElementColor(z, atom.r, atom.g, atom.b);
        if (isCore) ++result.coreCount; else ++result.shellCount;
    }

    result.structure = out;
    result.success = true;
    result.message = "Assigned " + std::to_string(result.coreCount) + " core (" + params.coreElement
        + ") and " + std::to_string(result.shellCount) + " shell (" + params.shellElement + ") atoms.";
    return result;
}

}
