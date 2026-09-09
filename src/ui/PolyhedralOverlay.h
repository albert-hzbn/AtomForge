#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <array>
#include <vector>

struct ImDrawList;

struct PolyhedralOverlaySettings
{
    bool centerAtomIndexFilterEnabled = false;
    std::vector<int> centerAtomIndices;
    int maxDisplayedCenters = 350;
    int maxNeighborCandidates = 24;
    bool showEdges = true;
    float faceOpacity = 0.08f;
    float edgeOpacity = 0.82f;
    bool centerElementFilterEnabled = false;
    bool ligandElementFilterEnabled = false;
    std::array<bool, 119> centerElementMask = {};
    std::array<bool, 119> ligandElementMask = {};
};

void drawPolyhedralOverlay(ImDrawList* drawList,
                           const glm::mat4& projection,
                           const glm::mat4& view,
                           int framebufferWidth,
                           int framebufferHeight,
                           const Structure& structure,
                           const std::vector<int>& selectedInstanceIndices,
                           const std::vector<int>& instanceToAtomIndex,
                           const PolyhedralOverlaySettings& settings,
                           const std::vector<glm::vec3>& elementColors,
                           bool enabled);
