#pragma once
#include <array>
#include <vector>
#include <string>

struct AtomSite
{
    std::string symbol;
    int atomicNumber;

    double x;
    double y;
    double z;

    float r;
    float g;
    float b;
};

struct Structure
{
    std::vector<AtomSite> atoms;

    // Optional unit cell information; when present we can render
    // an accurate lattice-registered bounding box and show all atoms
    // in the full unit cell.
    bool hasUnitCell = false;
    std::array<std::array<double, 3>, 3> cellVectors;
    std::array<double, 3> cellOffset = {0.0, 0.0, 0.0};

    // When > 0, appendPbcBoundaryImages uses this instead of the
    // default tight tolerance so that atoms slightly off a cell face
    // still get periodic-boundary copies.  Set by CSL grain-boundary
    // builder to half the minimum layer spacing.
    float pbcBoundaryTol = 0.0f;

    // Per-atom grain orientation colors (IPF-Z). Populated by the polycrystal
    // builder. When non-empty (same size as atoms), "Crystal Orientation"
    // coloring mode is available.
    std::vector<std::array<float, 3>> grainColors;

    // User-visible note set on load to indicate how IPF was obtained
    // (metadata sidecar, geometry fallback, or unavailable).
    std::string ipfLoadStatus;
};

void getDefaultElementColor(int atomicNumber, float& r, float& g, float& b);

// Returns true if the filename has a supported input extension.
bool isSupportedStructureFile(const std::string& filename);

// Load structure from file and return detailed error text on failure.
bool loadStructureFromFile(const std::string& filename, Structure& structure, std::string& errorMessage);

Structure loadStructure(const std::string& filename);

// Save structure to file. format is an OpenBabel format string (e.g. "xyz", "cif", "vasp").
// Returns true on success.
bool saveStructure(const Structure& structure, const std::string& filename, const std::string& format);