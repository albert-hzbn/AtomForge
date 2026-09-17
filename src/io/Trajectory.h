#pragma once
#include "model/Structure.h"
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

// Auto-detects XYZ/extXYZ vs. LAMMPS dump from the first non-blank line
// ("ITEM:" marks a LAMMPS dump) and dispatches to the matching loader.
std::vector<Structure> loadTrajectory(const std::string& path);
