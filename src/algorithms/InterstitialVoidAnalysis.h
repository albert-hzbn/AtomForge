#pragma once

#include "io/StructureLoader.h"

#include <array>
#include <string>
#include <vector>

#include <glm/glm.hpp>

enum class InterstitialVoidKind
{
    Tetrahedral = 0,
    Octahedral,
    Irregular,
};

struct InterstitialVoidRegion
{
    glm::vec3 position = glm::vec3(0.0f);
    float clearance = 0.0f;
    int coordination = 0;
    float volume = 0.0f;
    int volumeBin = 0;
    InterstitialVoidKind kind = InterstitialVoidKind::Irregular;
    std::vector<int> neighborAtomIndices;
    std::vector<std::array<int, 3>> faces;
    std::vector<std::array<glm::vec3, 3>> faceVertices;
};

struct InterstitialVoidDetectionParams
{
    int gridResolution = 14;
    // 0 or negative means uncapped (detect all unique void regions).
    int maxVoids = 0;
    float minClearance = 0.6f;
    float minSeparation = 0.8f;
    float firstNeighborShellFactor = 1.22f;
    int maxNeighborCount = 12;
};

struct InterstitialVoidDetectionResult
{
    bool success = false;
    std::string message;
    std::vector<InterstitialVoidRegion> regions;
};

InterstitialVoidDetectionResult detectInterstitialVoidRegions(
    const Structure& structure,
    const InterstitialVoidDetectionParams& params);
