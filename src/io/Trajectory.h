#pragma once
#include "model/Structure.h"
#include <string>
#include <vector>

// XYZ/extXYZ frames, preserving a complete Lattice matrix when supplied.
std::vector<Structure> loadXyzTrajectory(const std::string& path);
