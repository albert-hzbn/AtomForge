#pragma once

#include "science/Json.h"

#include <array>
#include <string>
#include <vector>

namespace atomforge::science
{
// A per-atom scalar from a tool result, aligned with the analysed atoms.
// Non-finite values mark invalid environments.
struct AtomProperty
{
    std::string name;
    std::vector<double> values;
};

std::vector<AtomProperty> perAtomProperties(const std::string& tool, const Json& result);

// Names of properties that callers pick from a tool's several properties.
namespace propertyName
{
inline constexpr const char* coordination = "Coordination";
inline constexpr const char* voronoiCoordination = "Voronoi coordination";
inline constexpr const char* orientationAngle = "Orientation angle (degrees)";
inline constexpr const char* ptmOrdering =
    "PTM ordering (0 none, 1 pure, 2 L1_0, 3 L1_2 minority, 4 L1_2 majority, 5 B2, 6 zincblende/wurtzite, 7 other)";
}
// Tool ids that report per-atom properties (for registry consistency checks).
std::vector<std::string> toolsWithAtomProperties();

// Perceptually uniform viridis colour map (polynomial fit), t clamped to [0, 1].
std::array<float, 3> viridis(double t);

struct PropertyRange
{
    double low = 0.0;
    double high = 1.0;
};

// Finite minimum and maximum (0..1 when no value is finite).
PropertyRange finiteRange(const std::vector<double>& values);

struct PropertyDisplay
{
    bool autoRange = true;
    PropertyRange range;
    bool hideOutside = false;  // hide atoms outside [low, high]
    bool hideInvalid = false;  // hide atoms without a finite value
};

// Colour for each atom (grey for invalid) and whether it is shown.
void colourByProperty(const std::vector<double>& values, const PropertyDisplay& display,
                      std::vector<std::array<float, 3>>& colours, std::vector<bool>& visible);
}
