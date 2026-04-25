#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "io/StructureLoader.h"
#include "graphics/SceneBuffers.h"

struct ImDrawList;

struct MeasurementOverlayState
{
    bool showDistancePopup = false;
    char distanceMessage[256] = {0};
    bool showDistanceLine = false;
    int distanceLineIdx0 = -1;
    int distanceLineIdx1 = -1;

    bool showAnglePopup = false;
    char angleMessage[256] = {0};
    bool showAngleLines = false;
    int angleLineIdx0 = -1; // first atom
    int angleLineIdx1 = -1; // vertex atom
    int angleLineIdx2 = -1; // third atom

    bool showAtomInfoPopup = false;
    char atomInfoMessage[512] = {0};

    void clearVisuals();
};

void processMeasurementRequests(MeasurementOverlayState& state,
                                bool requestMeasureDistance,
                                bool requestMeasureAngle,
                                bool requestAtomInfo,
                                const std::vector<int>& selectedInstanceIndices,
                                const SceneBuffers& sceneBuffers,
                                const Structure& structure,
                                const std::vector<float>& elementRadii);

void drawMeasurementPopups(MeasurementOverlayState& state);

void drawMeasurementOverlays(const MeasurementOverlayState& state,
                             ImDrawList* drawList,
                             const glm::mat4& projection,
                             const glm::mat4& view,
                             int viewportWidth,
                             int viewportHeight,
                             const SceneBuffers& sceneBuffers);

void drawElementLabelsOverlay(ImDrawList* drawList,
                              const glm::mat4& projection,
                              const glm::mat4& view,
                              int viewportWidth,
                              int viewportHeight,
                              const SceneBuffers& sceneBuffers,
                              const Structure& structure);

void drawOrientationAxesOverlay(ImDrawList* drawList,
                                const glm::mat4& view,
                                int viewportWidth,
                                int viewportHeight,
                                const Structure& structure);

// Draw the IPF (Inverse Pole Figure) standard triangle legend for cubic
// crystal orientation coloring.  Placed in the bottom-right corner.
void drawIPFTriangleLegend(ImDrawList* drawList,
                           int viewportWidth,
                           int viewportHeight);
