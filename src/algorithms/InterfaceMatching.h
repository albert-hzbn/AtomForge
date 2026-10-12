#pragma once
// Interfaces between two crystals by geometric lattice matching, after
// Stradi et al., J. Phys.: Condens. Matter 29, 185901 (2017), the method of
// QuantumATK's InterfaceBuilder: each crystal is cut along a plane, the
// surface lattice of B is rotated through a range of angles, supercells of
// the two surface lattices that coincide within a strain limit are listed,
// and a chosen match is assembled into one cell with B strained onto A (or A
// onto B, or both onto their average).

#include "model/Structure.h"

#include <cstddef>
#include <string>
#include <vector>

namespace atomforge::interfaces
{
// One side of the interface: a crystal cut along (h k l), or a structure that
// is already a slab (its a and b vectors span the interface plane).
struct Side
{
    Structure crystal;
    bool cleave = true;
    int h = 0, k = 0, l = 1;
    int layers = 4;  // surface layers (repeats of the stacking translation) in the slab
};

struct SearchParams
{
    int maxRepeats = 6;          // |i|, |j| of B's surface vectors i b1 + j b2
    double maxStrain = 0.03;     // mean absolute strain, (|e11| + |e22| + |e12|) / 3
    double maxArea = 150.0;      // interface cell area (A^2)
    double minCellAngle = 20.0;  // angle between the interface cell vectors (deg)
    double angleMin = 0.0;       // rotation of B about the interface normal (deg)
    double angleMax = 180.0;
    double angleStep = 1.0;
    std::size_t maxMatches = 2000;
};

struct Match
{
    int a[2][2] = {};       // rows: the interface cell in A's surface vectors
    int b[2][2] = {};       // rows: the interface cell in B's surface vectors
    double thetaDeg = 0;    // rotation of B about the normal
    double strain[3] = {};  // e11, e22, e12 taking B's cell onto A's
    double meanStrain = 0;  // (|e11| + |e22| + |e12|) / 3
    double area = 0;        // A^2
    double cellAngle = 0;   // deg
    int cellsA = 0, cellsB = 0;  // surface cells of each side in the interface cell
    int atoms = 0;          // atoms of the interface with the requested layers
    bool pareto = false;    // no other match has as few atoms and less strain
};

enum class StrainOn { B, A, Both };

struct Assembly
{
    double gap = 2.0;          // A between the facing atomic layers of A and B
    double vacuum = 0.0;       // 0: periodic along the normal (two interfaces, each with the gap)
    double shift[2] = {0, 0};  // in-plane shift of B, fractions of the interface cell
    StrainOn strainOn = StrainOn::B;
};

// The slab of one side with its layers, in its surface cell: a along x, b in
// the xy plane, the stacking direction along +z. Layers are stacked along the
// crystal's own translation, so fcc (111) keeps its ABC stacking. False with
// a message when it cannot be made.
bool prepareSlab(const Side& side, Structure& slab, std::string& error);

// Every match of the two surface cells within the limits, sorted by atoms
// and then strain, with the Pareto-optimal ones marked.
std::vector<Match> findMatches(const Structure& slabA, const Structure& slabB, const SearchParams& params);

// A suggested match: of the Pareto-optimal ones, the one with the lowest
// mean strain (%) plus area / 25 A^2 (1% strain weighs like 25 A^2). -1 when none.
int suggestedMatch(const std::vector<Match>& matches);

// The interface of a match: slab A below, slab B above, B rotated by the
// match's angle and strained as asked.
Structure buildInterface(const Structure& slabA, const Structure& slabB, const Match& match, const Assembly& assembly);

// The slab repeated in its plane by the integer rows of m (in-plane
// translations only, so the stacking vector is never used).
Structure inPlaneSupercell(const Structure& slab, const int m[2][2]);
}
