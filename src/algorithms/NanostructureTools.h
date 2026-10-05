#pragma once

#include "model/Structure.h"

#include <string>

// Nanostructure tools that complement the existing Wulff-construction
// Nanocrystal Builder (NanoCrystalBuilder.h), which produces finite,
// compositionally uniform particles:
//  - buildNanowire: a 1D-periodic rod (periodic along one lattice axis,
//    bounded by a circular or regular-polygon cross-section elsewhere).
//  - applyCoreShell: recolors/relabels an existing finite structure (e.g. a
//    Nanocrystal Builder output) into two concentric compositions by radius.
namespace atomforge
{

struct NanowireParams
{
    // Which of the source structure's lattice vectors (0=a, 1=b, 2=c) is
    // kept periodic as the wire axis; the other two directions are cut to a
    // finite cross-section.
    int axis = 2;

    // Cross-section size, in Angstrom: the apothem (center-to-face
    // distance) when sides >= 3, or the plain radius for a circular
    // cross-section (sides < 3).
    double radius = 10.0;

    // Regular-polygon facet count for the cross-section; < 3 means a
    // circular cross-section instead.
    int sides = 0;

    // Extra bounding-box padding (Angstrom) beyond the cross-section, so the
    // wire does not self-interact across the non-periodic directions'
    // implicit cell boundary.
    double vacuum = 10.0;

    // Number of unit cells stacked along the wire axis before cutting.
    int axisRepeats = 1;
};

struct NanowireResult
{
    bool success = false;
    std::string message;
    Structure structure;
    int atomCount = 0;
};

NanowireResult buildNanowire(const Structure& source, const NanowireParams& params);

struct CoreShellParams
{
    double coreRadius = 5.0;
    std::string coreElement;
    std::string shellElement;

    // When true, the core/shell boundary is centered on the atoms' own
    // centroid; otherwise `center` is used directly.
    bool useCentroid = true;
    double center[3] = {0.0, 0.0, 0.0};
};

struct CoreShellResult
{
    bool success = false;
    std::string message;
    Structure structure;
    int coreCount = 0;
    int shellCount = 0;
};

CoreShellResult applyCoreShell(const Structure& source, const CoreShellParams& params);

}
