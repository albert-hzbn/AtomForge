#pragma once

#include "io/StructureLoader.h"

#include <glm/glm.hpp>

#include <vector>

struct VoronoiFace
{
    std::vector<glm::vec3> vertices; // ordered convex polygon vertices in world space
};

struct VoronoiCell
{
    std::vector<VoronoiFace> faces;
};

struct VoronoiDiagram
{
    std::vector<VoronoiCell> cells;
};

// Compute Voronoi tessellation for all atoms in the structure.
// Uses periodic boundary conditions when the structure has a unit cell.
// Returns cells clipped to the unit cell bounding box.
VoronoiDiagram computeVoronoi(const Structure& structure);
