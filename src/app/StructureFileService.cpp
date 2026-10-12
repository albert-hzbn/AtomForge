#include "app/StructureFileService.h"
#include "io/StructureLoader.h"

#include "graphics/StructureInstanceBuilder.h"

#include <algorithm>
#include <cctype>

bool loadStructureFromPath(
    const std::string& filePath,
    Structure& structure,
    std::string& errorMessage)
{
    return loadStructureFromFile(filePath, structure, errorMessage);
}

bool saveStructureWithOptionalSupercell(
    const Structure& structure,
    bool useTransformMatrix,
    const int (&transformMatrix)[3][3],
    const std::string& outputPath,
    const std::string& format,
    std::size_t& savedAtomCount)
{
    const bool shouldExpandToSupercell = useTransformMatrix && structure.hasUnitCell;
    const Structure structureToSave = shouldExpandToSupercell
        ? buildSupercell(structure, transformMatrix)
        : structure;

    // XYZ keeps a periodic cell as extended XYZ, so its atoms are written as
    // they are (not the periodic images the view adds at the cell faces).
    savedAtomCount = structureToSave.atoms.size();
    return saveStructure(structureToSave, outputPath, format);
}
