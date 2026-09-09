#pragma once

#include "electronic/Volume.h"

namespace atomforge::electronic
{
Grid densityDifference(const Grid& a, const Grid& b, double weight = 1);
Grid thresholdMask(const Grid& grid, double low, double high);
Grid booleanMask(const Grid& a, const Grid& b, const std::string& operation);
Grid applyMask(const Grid& grid, const Grid& mask);
std::array<Grid,2> splitDensity(const Grid& grid);
// Positive accumulation, positive depletion magnitude, net electron change.
std::array<double,3> chargeSummary(const Grid& grid);
// Normal distance and cumulative electron count from the cell's lower face.
std::vector<glm::dvec2> cumulativeCharge(const Grid& grid, int axis);
}
