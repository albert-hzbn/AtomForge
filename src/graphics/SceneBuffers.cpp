#include "SceneBuffers.h"
#include "math/StructureMath.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <iostream>

namespace
{
constexpr float kMinBondRadius = 0.06f;
constexpr float kMaxBondRadius = 0.16f;
constexpr float kBondInsetFactor = 0.55f;
constexpr size_t kMaxBondCount = 1000000;  // Prevent memory exhaustion
constexpr float kSpatialHashCellSize = 4.0f;  // Cell size for spatial hashing

float clampValue(float value, float minValue, float maxValue)
{
    return std::max(minValue, std::min(maxValue, value));
}

// Hash function for spatial grid cell coordinates
struct Vec3iHash
{
    size_t operator()(const glm::ivec3& v) const
    {
        return ((size_t)v.x * 73856093) ^ ((size_t)v.y * 19349663) ^ ((size_t)v.z * 83492791);
    }
};

// Get grid cell coordinates for a position
glm::ivec3 getGridCell(const glm::vec3& pos, float cellSize)
{
    return glm::ivec3(
        (int)std::floor(pos.x / cellSize),
        (int)std::floor(pos.y / cellSize),
        (int)std::floor(pos.z / cellSize)
    );
}

bool isLikelyCovalentElement(int atomicNumber)
{
    switch (atomicNumber)
    {
        case 1:   // H
        case 5:   // B
        case 6:   // C
        case 7:   // N
        case 8:   // O
        case 9:   // F
        case 14:  // Si
        case 15:  // P
        case 16:  // S
        case 17:  // Cl
        case 32:  // Ge
        case 33:  // As
        case 34:  // Se
        case 35:  // Br
        case 52:  // Te
        case 53:  // I
            return true;
        default:
            return false;
    }
}
}

void SceneBuffers::init(GLuint sphereVAO, GLuint cylinderVAO)
{
    glGenBuffers(1, &instanceVBO);
    glGenBuffers(1, &colorVBO);
    glGenBuffers(1, &scaleVBO);
    glGenBuffers(1, &shininessVBO);
    glGenBuffers(1, &bondStartVBO);
    glGenBuffers(1, &bondEndVBO);
    glGenBuffers(1, &bondColorAVBO);
    glGenBuffers(1, &bondColorBVBO);
    glGenBuffers(1, &bondRadiusVBO);
    glGenBuffers(1, &bondShininessAVBO);
    glGenBuffers(1, &bondShininessBVBO);
    glGenVertexArrays(1, &lineVAO);
    glGenBuffers(1, &lineVBO);

    // Wire instance buffers into the sphere VAO so the atom draw call
    // picks them up as instanced attributes.
    glBindVertexArray(sphereVAO);

    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(1, 1);

    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(2, 1);

    glBindBuffer(GL_ARRAY_BUFFER, scaleVBO);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glVertexAttribDivisor(3, 1);

    glBindBuffer(GL_ARRAY_BUFFER, shininessVBO);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glVertexAttribDivisor(4, 1);

    glBindVertexArray(0);

    glBindVertexArray(cylinderVAO);

    glBindBuffer(GL_ARRAY_BUFFER, bondStartVBO);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(1, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondEndVBO);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(2, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondColorAVBO);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(3, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondColorBVBO);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glVertexAttribDivisor(4, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondRadiusVBO);
    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glVertexAttribDivisor(5, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondShininessAVBO);
    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glVertexAttribDivisor(6, 1);

    glBindBuffer(GL_ARRAY_BUFFER, bondShininessBVBO);
    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glVertexAttribDivisor(7, 1);

    glBindVertexArray(0);

    // Wire the line VBO into the dedicated line VAO.
    glBindVertexArray(lineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);
}

void SceneBuffers::upload(const StructureInstanceData& data,
                          bool bondElementFilterEnabled,
                          const std::array<bool, 119>& bondElementFilterMask)
{
    atomCount     = data.positions.size();
    bondCount     = 0;
    orbitCenter   = atomCount > 0 ? data.orbitCenter : glm::vec3(0.0f);
    boxLines      = data.boxLines;

    // Determine if this is a large structure
    constexpr size_t kCacheCutoff = 100000;
    cpuCachesDisabled = (atomCount > kCacheCutoff);

    if (!cpuCachesDisabled)
    {
        // Small structure: keep CPU caches for picking and undo/redo
        atomPositions = data.positions;
        atomColors    = data.colors;
        atomRadii     = data.scales;
        atomShininess = data.shininess;
        atomIndices   = data.atomIndices;
    }
    else
    {
        // Large structure: disable CPU caches to save memory (~80% reduction)
        // GPU buffers are still uploaded (used for rendering)
        std::cerr << "Large structure (" << atomCount << " atoms): CPU caches disabled. "
                  << "Atom picking and per-atom color editing disabled." << std::endl;
        atomPositions.clear();
        atomColors.clear();
        atomRadii.clear();
        atomShininess.clear();
        atomIndices.clear();
        atomPositions.shrink_to_fit();
        atomColors.shrink_to_fit();
        atomRadii.shrink_to_fit();
        atomShininess.shrink_to_fit();
        atomIndices.shrink_to_fit();
    }

    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 data.positions.size() * sizeof(glm::vec3),
                 data.positions.empty() ? nullptr : data.positions.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 data.colors.size() * sizeof(glm::vec3),
                 data.colors.empty() ? nullptr : data.colors.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, scaleVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 data.scales.size() * sizeof(float),
                 data.scales.empty() ? nullptr : data.scales.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, shininessVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 data.shininess.size() * sizeof(float),
                 data.shininess.empty() ? nullptr : data.shininess.data(),
                 GL_STATIC_DRAW);

    std::vector<glm::vec3> bondStarts;
    std::vector<glm::vec3> bondEnds;
    std::vector<glm::vec3> bondColorA;
    std::vector<glm::vec3> bondColorB;
    std::vector<float> bondRadii;
    std::vector<float> bondShininessA;
    std::vector<float> bondShininessB;
    const std::vector<float>& bondSourceRadii =
        (data.bondRadii.size() == data.positions.size()) ? data.bondRadii : data.scales;
    const std::vector<int>& bondSourceAtomicNumbers = data.atomicNumbers;

    // Disable bond detection for very large structures (>500k atoms)
    constexpr size_t kBondDisableCutoff = 500000;
    if (atomPositions.size() >= 2 && atomCount < kBondDisableCutoff)
    {
        // Pre-allocate for expected bonds (not squared!)
        size_t expectedBonds = std::min(atomPositions.size() * 6, kMaxBondCount);
        bondStarts.reserve(std::min(expectedBonds, kMaxBondCount));
        bondEnds.reserve(std::min(expectedBonds, kMaxBondCount));
        bondColorA.reserve(std::min(expectedBonds, kMaxBondCount));
        bondColorB.reserve(std::min(expectedBonds, kMaxBondCount));
        bondRadii.reserve(std::min(expectedBonds, kMaxBondCount));
        bondShininessA.reserve(std::min(expectedBonds, kMaxBondCount));
        bondShininessB.reserve(std::min(expectedBonds, kMaxBondCount));

        // Build spatial grid for efficient neighbor lookup
        std::unordered_map<glm::ivec3, std::vector<size_t>, Vec3iHash> grid;
        for (size_t i = 0; i < atomPositions.size(); ++i)
        {
            glm::ivec3 cell = getGridCell(atomPositions[i], kSpatialHashCellSize);
            grid[cell].push_back(i);
        }

        // Check bonds only within neighboring cells
        const int neighborRange = 2;  // Check 3x3x3 neighborhood (includes diagonals)
        for (size_t i = 0; i < atomPositions.size() && bondStarts.size() < kMaxBondCount; ++i)
        {
            glm::ivec3 centerCell = getGridCell(atomPositions[i], kSpatialHashCellSize);
            float radiusA = (i < bondSourceRadii.size()) ? bondSourceRadii[i] : 1.0f;

            // Check all neighboring cells
            for (int dx = -neighborRange; dx <= neighborRange; ++dx)
            {
                for (int dy = -neighborRange; dy <= neighborRange; ++dy)
                {
                    for (int dz = -neighborRange; dz <= neighborRange; ++dz)
                    {
                        glm::ivec3 neighborCell = centerCell + glm::ivec3(dx, dy, dz);
                        auto it = grid.find(neighborCell);
                        if (it == grid.end())
                            continue;

                        for (size_t j : it->second)
                        {
                            if (j <= i)
                                continue;  // Only check each pair once

                            if (bondElementFilterEnabled)
                            {
                                const int zA = (i < bondSourceAtomicNumbers.size()) ? bondSourceAtomicNumbers[i] : 0;
                                const int zB = (j < bondSourceAtomicNumbers.size()) ? bondSourceAtomicNumbers[j] : 0;

                                const bool allowedA = (zA >= 1 && zA <= 118) ? bondElementFilterMask[(size_t)zA] : false;
                                const bool allowedB = (zB >= 1 && zB <= 118) ? bondElementFilterMask[(size_t)zB] : false;
                                // Filter semantics: include bonds touching selected elements
                                // (e.g. "C,F" keeps C-*, F-* and C-F bonds).
                                if (!allowedA && !allowedB)
                                    continue;
                            }
                            else
                            {
                                const int zA = (i < bondSourceAtomicNumbers.size()) ? bondSourceAtomicNumbers[i] : 0;
                                const int zB = (j < bondSourceAtomicNumbers.size()) ? bondSourceAtomicNumbers[j] : 0;

                                // Default mode is covalent-focused: suppress metal-metal bonds.
                                if (zA >= 1 && zA <= 118 && zB >= 1 && zB <= 118 &&
                                    !isLikelyCovalentElement(zA) && !isLikelyCovalentElement(zB))
                                    continue;
                            }

                            glm::vec3 delta = atomPositions[j] - atomPositions[i];
                            float distance = glm::length(delta);
                            if (distance <= kMinBondDistance)
                                continue;

                            float radiusB = (j < bondSourceRadii.size()) ? bondSourceRadii[j] : 1.0f;
                            float maxBondDistance = (radiusA + radiusB) * kBondToleranceFactor;
                            if (distance > maxBondDistance)
                                continue;

                            glm::vec3 direction = delta / distance;
                            float insetA = std::min(radiusA * kBondInsetFactor, distance * 0.30f);
                            float insetB = std::min(radiusB * kBondInsetFactor, distance * 0.30f);
                            glm::vec3 start = atomPositions[i] + direction * insetA;
                            glm::vec3 end = atomPositions[j] - direction * insetB;
                            float visibleLength = glm::length(end - start);
                            if (visibleLength <= kMinBondDistance)
                                continue;

                            float bondRadius = clampValue(0.12f * (radiusA + radiusB) * 0.5f,
                                                          kMinBondRadius,
                                                          kMaxBondRadius);

                            bondStarts.push_back(start);
                            bondEnds.push_back(end);
                            bondColorA.push_back((i < atomColors.size()) ? atomColors[i] : glm::vec3(0.8f));
                            bondColorB.push_back((j < atomColors.size()) ? atomColors[j] : glm::vec3(0.8f));
                            bondRadii.push_back(bondRadius);
                            bondShininessA.push_back((i < atomShininess.size()) ? atomShininess[i] : 32.0f);
                            bondShininessB.push_back((j < atomShininess.size()) ? atomShininess[j] : 32.0f);

                            if (bondStarts.size() >= kMaxBondCount)
                            {
                                std::cerr << "Warning: Bond count capped at " << kMaxBondCount
                                          << " to prevent memory exhaustion. Structure has "
                                          << atomPositions.size() << " atoms." << std::endl;
                                goto bonds_done;
                            }
                        }
                    }
                }
            }
        }
        bonds_done:;
    }
    else if (atomCount >= kBondDisableCutoff)
    {
        std::cerr << "Bond detection disabled for structure with " << atomCount 
                  << " atoms (>500k cutoff). Bonds will not be rendered." << std::endl;
    }

    bondCount = bondStarts.size();

    glBindBuffer(GL_ARRAY_BUFFER, bondStartVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondStarts.size() * sizeof(glm::vec3),
                 bondStarts.empty() ? nullptr : bondStarts.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondEndVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondEnds.size() * sizeof(glm::vec3),
                 bondEnds.empty() ? nullptr : bondEnds.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondColorAVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondColorA.size() * sizeof(glm::vec3),
                 bondColorA.empty() ? nullptr : bondColorA.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondColorBVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondColorB.size() * sizeof(glm::vec3),
                 bondColorB.empty() ? nullptr : bondColorB.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondRadiusVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondRadii.size() * sizeof(float),
                 bondRadii.empty() ? nullptr : bondRadii.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondShininessAVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondShininessA.size() * sizeof(float),
                 bondShininessA.empty() ? nullptr : bondShininessA.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, bondShininessBVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 bondShininessB.size() * sizeof(float),
                 bondShininessB.empty() ? nullptr : bondShininessB.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 data.boxLines.size() * sizeof(glm::vec3),
                 data.boxLines.empty() ? nullptr : data.boxLines.data(),
                 GL_STATIC_DRAW);
}

void SceneBuffers::highlightAtom(int idx, glm::vec3 color)
{
    if (idx < 0 || (size_t)idx >= atomCount)
        return;
    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferSubData(GL_ARRAY_BUFFER,
                    (GLintptr)(idx * (GLintptr)sizeof(glm::vec3)),
                    sizeof(glm::vec3),
                    &color);
}

void SceneBuffers::restoreAtomColor(int idx)
{
    if (idx < 0 || (size_t)idx >= atomColors.size())
        return;
    glm::vec3 orig = atomColors[idx];
    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferSubData(GL_ARRAY_BUFFER,
                    (GLintptr)(idx * (GLintptr)sizeof(glm::vec3)),
                    sizeof(glm::vec3),
                    &orig);
}
