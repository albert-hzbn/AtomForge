#pragma once

#include <vector>
#include <glm/glm.hpp>

// Return the chemical symbol for atomic number z (1–118), or "?" otherwise.
const char* elementSymbol(int z);

// Return the full element name for atomic number z (1–118), or "Unknown" otherwise.
const char* elementName(int z);

// Covalent radii (Angstrom) from Cordero et al., Dalton Trans. 2008.
// Returns a vector of size 119 indexed by atomic number (index 0 unused).
std::vector<float> makeLiteratureCovalentRadii();

// CPK-style default colours for each element.
// Returns a vector of size 119 indexed by atomic number (index 0 unused).
std::vector<glm::vec3> makeDefaultElementColors();
