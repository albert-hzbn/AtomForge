#include "app/StructureFileService.h"

#include "graphics/StructureInstanceBuilder.h"

#include <algorithm>
#include <cctype>

namespace
{
bool isXyzFormat(const std::string& format)
{
    std::string lower = format;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char ch) { return (char)std::tolower(ch); });
    return lower == "xyz";
}

Structure buildPeriodicXyzExportStructure(const Structure& structure)
{
    // Reuse the renderer's periodic-boundary image logic so XYZ export
    // matches what users see in the viewport for periodic cells.
    static const int kIdentityTransform[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    std::vector<float> elementRadii(119, 1.0f);
    std::vector<float> elementShininess(119, 32.0f);
    const StructureInstanceData instances = buildStructureInstanceData(
        structure,
        false,
        kIdentityTransform,
        elementRadii,
        elementShininess);

    Structure exportStructure = structure;
    exportStructure.atoms.clear();
    exportStructure.atoms.reserve(instances.positions.size());

    for (std::size_t i = 0; i < instances.positions.size(); ++i)
    {
        AtomSite atom;
        atom.atomicNumber = instances.atomicNumbers[i];
        atom.x = instances.positions[i].x;
        atom.y = instances.positions[i].y;
        atom.z = instances.positions[i].z;
        atom.r = 1.0f;
        atom.g = 1.0f;
        atom.b = 1.0f;
        exportStructure.atoms.push_back(atom);
    }

    return exportStructure;
}
}

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

    const Structure exportStructure = isXyzFormat(format)
        ? buildPeriodicXyzExportStructure(structureToSave)
        : structureToSave;

    savedAtomCount = exportStructure.atoms.size();
    return saveStructure(exportStructure, outputPath, format);
}
