#pragma once
#include "model/Structure.h"
#include <istream>
#include <string>
#include <vector>

// XYZ/extXYZ frames, preserving a complete Lattice matrix when supplied.
std::vector<Structure> loadXyzTrajectory(const std::string& path);

// LAMMPS custom-style dump frames ("ITEM: TIMESTEP" ... "ITEM: ATOMS ...").
// Requires an "element" column (from `dump_modify ... element <symbols>`) to
// map LAMMPS's arbitrary per-simulation type IDs to chemical symbols; a
// LAMMPS "type" integer is not reliably an atomic number, so it is never
// inferred from "type" alone. Accepts wrapped (x y z) or unwrapped
// (xu yu zu) coordinates, and orthogonal or triclinic box bounds.
std::vector<Structure> loadLammpsTrajectory(const std::string& path);

// The elements of LAMMPS atom types 1, 2, ... ("Cu Ni" or "Cu,Ni"), used for a
// dump with a type column but no element column (the default "dump atom").
// Process-wide; empty restores the requirement for an element column.
// Throws for an unknown element symbol.
void setLammpsTypeElements(const std::string& symbols);
std::string lammpsTypeElements();
// Up to maxFrames dump frames from the current stream position.
std::vector<Structure> loadLammpsFrames(std::istream& input, std::size_t maxFrames);

// Auto-detects XYZ/extXYZ vs. LAMMPS dump from the first non-blank line
// ("ITEM:" marks a LAMMPS dump) and dispatches to the matching loader.
std::vector<Structure> loadTrajectory(const std::string& path);
