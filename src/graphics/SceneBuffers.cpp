#include "SceneBuffers.h"

#include <algorithm>
#include <cmath>

namespace
{
constexpr float kBondToleranceFactor = 1.18f;
constexpr float kMinBondDistance = 0.10f;
constexpr float kMinBondRadius = 0.10f;
constexpr float kMaxBondRadius = 0.24f;
constexpr float kBondInsetFactor = 0.55f;

float clampValue(float value, float minValue, float maxValue)
{
    return std::max(minValue, std::min(maxValue, value));
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

    glBindVertexArray(0);

    // Wire the line VBO into the dedicated line VAO.
    glBindVertexArray(lineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);
}

void SceneBuffers::upload(const StructureInstanceData& data)
{
    atomCount     = data.positions.size();
    bondCount     = 0;
    orbitCenter   = atomCount > 0 ? data.orbitCenter : glm::vec3(0.0f);
    boxLines      = data.boxLines;
    atomPositions = data.positions;
    atomColors    = data.colors;
    atomRadii     = data.scales;
    atomShininess = data.shininess;
    atomIndices   = data.atomIndices;

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

    if (atomPositions.size() >= 2)
    {
        bondStarts.reserve(atomPositions.size());
        bondEnds.reserve(atomPositions.size());
        bondColorA.reserve(atomPositions.size());
        bondColorB.reserve(atomPositions.size());
        bondRadii.reserve(atomPositions.size());

        for (size_t i = 0; i < atomPositions.size(); ++i)
        {
            for (size_t j = i + 1; j < atomPositions.size(); ++j)
            {
                glm::vec3 delta = atomPositions[j] - atomPositions[i];
                float distance = glm::length(delta);
                if (distance <= kMinBondDistance)
                    continue;

                float radiusA = (i < atomRadii.size()) ? atomRadii[i] : 1.0f;
                float radiusB = (j < atomRadii.size()) ? atomRadii[j] : 1.0f;
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

                float bondRadius = clampValue(0.18f * (radiusA + radiusB) * 0.5f,
                                              kMinBondRadius,
                                              kMaxBondRadius);

                bondStarts.push_back(start);
                bondEnds.push_back(end);
                bondColorA.push_back((i < atomColors.size()) ? atomColors[i] : glm::vec3(0.8f));
                bondColorB.push_back((j < atomColors.size()) ? atomColors[j] : glm::vec3(0.8f));
                bondRadii.push_back(bondRadius);
            }
        }
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
