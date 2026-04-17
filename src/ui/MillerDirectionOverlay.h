#pragma once

#include "io/StructureLoader.h"
#include "ui/FileBrowser.h"

#include <glm/glm.hpp>
#include <vector>

struct ImDrawList;

void drawMillerDirectionsOverlay(ImDrawList* drawList,
                                 const glm::mat4& projection,
                                 const glm::mat4& view,
                                 int framebufferWidth,
                                 int framebufferHeight,
                                 const Structure& structure,
                                 const std::vector<MillerDirection>& directions,
                                 bool enabled);
