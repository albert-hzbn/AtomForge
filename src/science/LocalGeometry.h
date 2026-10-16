#pragma once
// Building blocks shared by the local-environment analyses (template matching,
// Voronoi, CNA, Chill+, rings, ...): reading the atoms of a request, neighbour
// lists, spherical-harmonic factors, union-find and classification results.

#include "science/ScienceData.h"

#include <string>
#include <vector>

namespace atomforge::science
{
// The atoms of a per-atom request: "positions", "cell" (optional) and "pbc"
// (periodic axes, which need the cell).
struct AtomSet
{
    std::vector<Vec3> positions;
    Mat3 cell{};
    Pbc pbc{};
    bool hasCell = false;

    bool periodic() const { return pbc[0] || pbc[1] || pbc[2]; }
    bool fullyPeriodic() const { return pbc[0] && pbc[1] && pbc[2]; }
};

AtomSet readAtomSet(const Parameters& p, const char* positions = "positions");

// (V/N)^(1/3): from the cell, or without one from the bounding box padded by
// 1 Angstrom on each side.
double meanAtomicSpacing(const AtomSet& atoms);

// For each atom, the indices into `neighbors` of the pairs centred on it.
std::vector<std::vector<std::size_t>> neighborsByCenter(const std::vector<Neighbor>& neighbors, std::size_t atoms);

// Each atom's neighbour vectors within the cutoff, nearest first. Coincident
// atoms make local environments undefined and are an error.
std::vector<std::vector<Vec3>> sortedNeighborVectors(const AtomSet& atoms, double cutoff);

// True for one of the two directions of each bond of a neighbour list (an
// atom bonded to its own periodic image keeps one of the two shifts).
inline bool firstOfPair(const Neighbor& n) { return n.i < n.j || (n.i == n.j && n.shift > std::array<int, 3>{0, 0, 0}); }

// Orthonormal associated Legendre factor of Y_lm(theta, 0) for x = cos theta
// (Condon-Shortley phase).
double normalizedLegendre(int l, int m, double x);

// Root of i in a union-find forest (path halving).
std::size_t unionRoot(std::vector<std::size_t>& parent, std::size_t i);

// Adds a per-atom classification to a result: result[key] (one type index per
// atom) and, by type name, the counts and fractions. Without a prefix the
// tables are "type_names", "counts" and "fractions"; with one (e.g. "ordering_")
// "ordering_names", "ordering_counts" and "ordering_fractions".
void addClassification(Json& result, const std::string& key, const std::vector<int>& types,
                       const std::vector<std::string>& names, const std::string& prefix = "");
}
