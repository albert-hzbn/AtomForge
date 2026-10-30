#pragma once
// Building blocks shared by the builder and edit steps: covalent bonds,
// hydrogen saturation and appending atoms with their per-atom metadata.

#include "model/Structure.h"
#include "science/ScienceCore.h"

#include <set>
#include <string>
#include <vector>

namespace atomforge::pipeline
{
// A covalent bond i-j: vector from i to the bonded image of j.
struct Bond
{
    std::size_t i = 0, j = 0;
    science::Vec3 vector{};
};

// Covalent radius of an element (Angstrom; Cordero et al. 2008).
double covalentRadius(int atomicNumber);

// Pairs closer than scale x (sum of covalent radii), across periodic cell
// vectors; each pair once (i < j, or i == j for a bond to an own image).
std::vector<Bond> covalentBonds(const Structure& s, double scale = 1.2);
// Neighbour lists from those bonds: for each atom the bonded atoms and the
// vectors to them.
struct Neighbours
{
    std::vector<std::vector<std::size_t>> atoms;
    std::vector<std::vector<science::Vec3>> vectors;
};
Neighbours bondNeighbours(const Structure& s, double scale = 1.2);

// Rotation matrices act on column vectors: v' = R v.
// The rotation taking direction `from` to direction `to` (shortest arc).
science::Mat3 rotationBetween(const science::Vec3& from, const science::Vec3& to);
// The rotation by `angle` (radians, right-handed) about `axis`.
science::Mat3 rotationAbout(const science::Vec3& axis, double angle);
science::Vec3 rotate(const science::Mat3& r, const science::Vec3& v);

// A new atom (default colour) at the end, extending the aligned per-atom
// metadata and the selection.
void appendAtom(Structure& s, std::vector<char>& selected, const std::string& symbol, const science::Vec3& position, bool select);

// Adds hydrogen to atoms of the given elements that have fewer bonds than
// their valence (C, Si, Ge 4; N, P, B 3; O, S 2), along the directions of the
// missing bonds: tetrahedral (sp3), or in the plane of two bonds at about
// 120 degrees (sp2, for planar carbon, boron and nitrogen, unless sp3 is
// forced). Returns the number of hydrogen atoms added.
std::size_t saturateWithHydrogen(Structure& s, std::vector<char>& selected, const std::set<std::string>& elements,
                                 bool forceSp3 = false, double bondScale = 1.2, bool selectNew = true);
}
