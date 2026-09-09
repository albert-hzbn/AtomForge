#pragma once

#include "model/Structure.h"
#include "ui/FileBrowser.h"

#include <glm/glm.hpp>
#include <vector>

struct ImDrawList;

void drawLatticePlanesOverlay(ImDrawList* drawList,
                              const glm::mat4& projection,
                              const glm::mat4& view,
                              int framebufferWidth,
                              int framebufferHeight,
                              const Structure& structure,
                              const std::vector<LatticePlane>& planes,
                              bool enabled);
