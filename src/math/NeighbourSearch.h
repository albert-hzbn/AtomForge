#pragma once
// Neighbour pairs within a cutoff, found with a cell list instead of all N^2
// atom pairs.

#include "model/Structure.h"

#include <glm/glm.hpp>
#include <algorithm>
#include <functional>
#include <stdexcept>
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

// The same pairs, handed to visit one at a time instead of stored (for
// large structures and cutoffs). Nothing is stored, so the cutoff may span
// up to maxImages periodic images of the cell (a long RDF of a small cell).
void visitWithin(const Structure& structure, double cutoff, bool periodic, const std::function<void(const Pair&)>& visit,
                 int maxImages = 200);

// Calls visit(i, j) once for each unordered pair i < j that may lie within
// `cutoff`, as a candidate: the caller applies its own distance test (for
// example a minimum-image bond length). Uses within(); for a cell too small
// to search, every pair is a candidate.
template <typename Visit>
void forEachCandidatePair(const Structure& structure, double cutoff, bool periodic, Visit&& visit)
{
    std::vector<Pair> pairs;
    try {
        pairs = within(structure, cutoff * 1.001 + 1e-3, periodic);
    } catch (const std::runtime_error&) {
        const int n = static_cast<int>(structure.atoms.size());
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) visit(i, j);
        return;
    }
    std::vector<int> partners;
    for (std::size_t k = 0; k < pairs.size();) {
        const int i = pairs[k].i;
        partners.clear();
        for (; k < pairs.size() && pairs[k].i == i; ++k)
            if (pairs[k].j > i) partners.push_back(pairs[k].j);
        std::sort(partners.begin(), partners.end());
        partners.erase(std::unique(partners.begin(), partners.end()), partners.end());
        for (int j : partners) visit(i, j);
    }
}
}
