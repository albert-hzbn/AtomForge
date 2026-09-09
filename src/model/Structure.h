#pragma once
#include <array>
#include <vector>
#include <string>
#include <glm/glm.hpp>

struct AtomSite
{
    std::string symbol;
    int atomicNumber = 0;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    float r = 1.0f;
    float g = 1.0f;
    float b = 1.0f;
};

struct Structure
{
    std::vector<AtomSite> atoms;

    // Optional unit cell information; when present we can render
    // an accurate lattice-registered bounding box and show all atoms
    // in the full unit cell.
    bool hasUnitCell = false;
    std::array<std::array<double, 3>, 3> cellVectors{};
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

    // Optional per-atom grain/region ids for builder-generated bicrystals and
    // polycrystals. When present, Grain Boundary coloring uses bonded
    // cross-region neighbors instead of inferring boundaries from local order.
    std::vector<int> grainRegionIds;

    // User-visible note set on load to indicate how IPF was obtained
    // (metadata sidecar, geometry fallback, or unavailable).
    std::string ipfLoadStatus;

    // 3D convex-hull boundary of the displaced atom region, populated by the
    // dislocation builder.  Projected onto the slip plane (e1–lineDir), so it
    // forms a closed polygon that can be rendered with GL_LINE_LOOP.
    // Empty unless a dislocation was built and applied to this structure.
    std::vector<glm::vec3> dislocationLoopPoints;

    // Set to true after the auto-detection algorithm in refreshSceneBuffers has
    // run once for this structure (regardless of whether it found anything).
    // Prevents the O(N) detection from re-running on every updateBuffers call.
    bool dislocationDetectionDone = false;

    // Keep per-atom metadata aligned when an editor removes a site.
    void eraseAtom(std::size_t index)
    {
        if (index >= atoms.size()) return;
        if (grainColors.size() == atoms.size())
            grainColors.erase(grainColors.begin() + index);
        else
            grainColors.clear();
        if (grainRegionIds.size() == atoms.size())
            grainRegionIds.erase(grainRegionIds.begin() + index);
        else
            grainRegionIds.clear();
        atoms.erase(atoms.begin() + index);
        dislocationLoopPoints.clear();
        dislocationDetectionDone = false;
    }
};
