#pragma once

// Element colour schemes for the Display Settings (Element Colours tab).

#include <glm/glm.hpp>

#include <vector>

namespace atomforge::display
{
struct ColourScheme
{
    const char* name;
    const char* description;
};
// AtomForge (the defaults), Jmol, RasMol CPK and Greyscale, in that order.
const std::vector<ColourScheme>& colourSchemes();
// The colour of every element (indexed by atomic number, 0-118) in a scheme.
std::vector<glm::vec3> schemeColours(int scheme);
// The scheme a colour table is, or -1 after individual changes.
int matchingScheme(const std::vector<glm::vec3>& colours);
}  // namespace atomforge::display
