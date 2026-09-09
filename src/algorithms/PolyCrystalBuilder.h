#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <string>
#include <vector>

// Orientation mode for each grain in the polycrystal.
enum class GrainOrientationMode
{
    AllRandom = 0,     // Every grain gets a random orientation
    AllSpecified,      // Every grain has user-specified Euler angles
    PartialSpecified,  // Some grains are specified, rest are random
};

// Per-grain orientation override (used when mode != AllRandom).
struct GrainOrientation
{
    int   grainIndex = 0;
    float phi1       = 0.0f;   // Bunge Euler angle (degrees)
    float Phi        = 0.0f;
    float phi2       = 0.0f;
};

struct PolyParams
{
    // Box dimensions in Angstroms
    float sizeX = 50.0f;
    float sizeY = 50.0f;
    float sizeZ = 50.0f;

    // Number of grains (Voronoi seeds)
    int numGrains = 8;

    // Random seed for reproducibility
    int seed = 42;

    // Orientation assignment
    GrainOrientationMode orientationMode = GrainOrientationMode::AllRandom;
    std::vector<GrainOrientation> specifiedOrientations;
};

struct PolyBuildResult
{
    bool        success      = false;
    std::string message;
    int         inputAtoms   = 0;
    int         outputAtoms  = 0;
    int         numGrains    = 0;
};

// Build a polycrystalline microstructure by Voronoi tessellation.
// The reference structure must have a unit cell.
[[nodiscard]] PolyBuildResult buildPolycrystal(Structure& structure,
                                 const Structure& reference,
                                 const PolyParams& params,
                                 const std::vector<glm::vec3>& elementColors);
