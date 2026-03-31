#pragma once

#include "algorithms/VoronoiComputation.h"
#include "io/StructureLoader.h"

#include <glm/glm.hpp>

#include <vector>

struct ImDrawList;

void drawVoronoiOverlay(ImDrawList* drawList,
                        const glm::mat4& projection,
                        const glm::mat4& view,
                        int framebufferWidth,
                        int framebufferHeight,
                        const VoronoiDiagram& diagram,
                        const std::vector<int>& selectedIndices,
                        bool enabled);
