#pragma once
// Small molecules with experimental gas-phase geometries (NIST CCCBDB bond
// lengths and angles), for building molecular structures, packings and
// adsorbates. Positions are in Angstrom, centred on the centre of mass.

#include "science/ScienceCore.h"

#include <array>
#include <string>
#include <vector>

namespace atomforge::pipeline
{
struct MoleculeAtom
{
    std::string symbol;
    science::Vec3 position;
};

// The names of the library's molecules (formulas or common names).
std::vector<std::string> moleculeNames();
// A molecule by name (case-insensitive), a structure file (XYZ, POSCAR, ...)
// or a single element symbol (one atom). Throws for an unknown name, listing
// the known ones.
std::vector<MoleculeAtom> molecule(const std::string& name);
// The mass of a molecule in atomic mass units.
double moleculeMass(const std::vector<MoleculeAtom>& atoms);
// Principal axes of inertia (unit vectors, smallest moment first, right-handed):
// for a linear or long molecule the first is its long axis; for a planar one
// the last is its plane normal.
std::array<science::Vec3, 3> principalAxes(const std::vector<MoleculeAtom>& atoms);
}
