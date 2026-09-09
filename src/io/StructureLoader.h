#pragma once

#include "model/Structure.h"
#include "util/ElementData.h"

// Returns true if the filename has a supported input extension.
bool isSupportedStructureFile(const std::string& filename);

// Load structure from file and return detailed error text on failure.
bool loadStructureFromFile(const std::string& filename, Structure& structure, std::string& errorMessage);

Structure loadStructure(const std::string& filename);

// Save structure to file. format is an OpenBabel format string (e.g. "xyz", "cif", "vasp").
// Returns true on success.
bool saveStructure(const Structure& structure, const std::string& filename, const std::string& format);
