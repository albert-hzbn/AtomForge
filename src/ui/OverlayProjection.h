#pragma once
// Projection of world points to the screen for the ImDrawList overlays drawn
// over the structure view (surface meshes, trajectory lines).

#include "imgui.h"

#include <glm/glm.hpp>

// The screen point of `world` (false behind the camera), in the coordinates
// of a view `width` x `height` pixels, as the other overlays use.
inline bool projectToView(const glm::mat4& viewProjection, const glm::vec3& world, int width, int height, ImVec2& screen)
{
    const glm::vec4 clip = viewProjection * glm::vec4(world, 1.0f);
    if (clip.w <= 1e-6f) return false;
    const glm::vec3 ndc = glm::vec3(clip) / clip.w;
    screen = ImVec2((ndc.x * 0.5f + 0.5f) * static_cast<float>(width), (1.0f - (ndc.y * 0.5f + 0.5f)) * static_cast<float>(height));
    return true;
}
