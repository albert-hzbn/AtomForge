#pragma once

#include "StructureLoader.h"

#include <glm/glm.hpp>
#include <vector>

struct StructureInstanceData
{
    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> colors;
    std::vector<float>     scales;
    std::vector<float>     shininess;
    std::vector<glm::vec3> boxLines;
    std::vector<int>       atomIndices; // maps each instance to its index in Structure::atoms
    glm::vec3 orbitCenter = glm::vec3(0.0f, 0.0f, 0.0f);
};

StructureInstanceData buildStructureInstanceData(
    const Structure& structure,
    bool useTransformMatrix,
    const int (&transformMatrix)[3][3],
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess
);

    // Expand `structure` according to `transformMatrix` and return a new Structure
    // whose atoms list contains every atom in the resulting supercell with Cartesian
    // coordinates, and whose cellVectors reflect the supercell unit cell.
    // If the transform is identity / singular or the structure has no unit cell,
    // the original structure is returned unchanged.
    Structure buildSupercell(
        const Structure& structure,
        const int (&transformMatrix)[3][3]
    );
