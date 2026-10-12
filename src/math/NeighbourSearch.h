#pragma once
// Neighbour pairs within a cutoff, found with a cell list instead of all N^2
// atom pairs.

#include "model/Structure.h"

#include <glm/glm.hpp>
#include <vector>

namespace neighbours
{
struct Pair
{
    int i = 0;
    int j = 0;
    glm::dvec3 d{0.0};  // from atom i to (the image of) atom j
    double r = 0.0;     // |d|
};

// Every ordered pair (i, j) with 0 < |d| <= cutoff, grouped by i in
// increasing order. Periodic: over all images of the cell, so a cell smaller
// than the cutoff lists each image of a neighbour (as a larger supercell
// would); positions are wrapped into the cell first. Not periodic: the plain
// positions. Throws std::runtime_error for a singular cell or a cutoff that
// spans more than 20 images of the cell.
std::vector<Pair> within(const Structure& structure, double cutoff, bool periodic);
}
