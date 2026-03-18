#pragma once

#include "io/StructureLoader.h"

#include <cmath>

#include <glm/glm.hpp>

constexpr float kBondToleranceFactor = 1.18f;
constexpr float kMinBondDistance = 0.10f;

inline glm::mat3 makeCellMatrix(const Structure& structure)
{
    return glm::mat3(
        glm::vec3((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]),
        glm::vec3((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]),
        glm::vec3((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]));
}

inline bool tryMakeCellMatrices(const Structure& structure, glm::mat3& cell, glm::mat3& invCell)
{
    if (!structure.hasUnitCell)
        return false;

    cell = makeCellMatrix(structure);
    const float det = glm::determinant(cell);
    if (std::abs(det) <= 1e-8f)
        return false;

    invCell = glm::inverse(cell);
    return true;
}

inline glm::vec3 minimumImageDelta(const glm::vec3& deltaCartesian,
                                   bool usePbc,
                                   const glm::mat3& cell,
                                   const glm::mat3& invCell)
{
    if (!usePbc)
        return deltaCartesian;

    glm::vec3 frac = invCell * deltaCartesian;
    frac -= glm::round(frac);
    return cell * frac;
}

inline bool tryCartesianToFractional(const Structure& structure,
                                     const glm::vec3& cart,
                                     glm::vec3& frac)
{
    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(structure, cell, invCell))
        return false;

    const glm::vec3 origin((float)structure.cellOffset[0],
                           (float)structure.cellOffset[1],
                           (float)structure.cellOffset[2]);
    frac = invCell * (cart - origin);
    return true;
}

inline bool tryFractionalToCartesian(const Structure& structure,
                                     const glm::vec3& frac,
                                     glm::vec3& cart)
{
    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(structure, cell, invCell))
        return false;

    const glm::vec3 origin((float)structure.cellOffset[0],
                           (float)structure.cellOffset[1],
                           (float)structure.cellOffset[2]);
    cart = origin + cell * frac;
    return true;
}
