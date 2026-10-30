#pragma once
// Small molecules with experimental gas-phase geometries (NIST CCCBDB bond
// lengths and angles), for building molecular structures, packings and
// adsorbates. Positions are in Angstrom, centred on the centre of mass.

#include "science/ScienceCore.h"

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
// A molecule by name (case-insensitive); a single element symbol gives one
// atom. Throws for an unknown name, listing the known ones.
std::vector<MoleculeAtom> molecule(const std::string& name);
// The mass of a molecule in atomic mass units.
double moleculeMass(const std::vector<MoleculeAtom>& atoms);
}
