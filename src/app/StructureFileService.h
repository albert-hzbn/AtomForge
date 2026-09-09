#pragma once

#include "model/Structure.h"

#include <cstddef>
#include <string>

bool loadStructureFromPath(
    const std::string& filePath,
    Structure& structure,
    std::string& errorMessage);

bool saveStructureWithOptionalSupercell(
    const Structure& structure,
    bool useTransformMatrix,
    const int (&transformMatrix)[3][3],
    const std::string& outputPath,
    const std::string& format,
    std::size_t& savedAtomCount);
