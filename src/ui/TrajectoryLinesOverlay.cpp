#include "ui/TrajectoryLinesOverlay.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>

TrajectoryLinesOverlay& trajectoryLinesOverlay()
{
    static TrajectoryLinesOverlay overlay;
    return overlay;
}

void drawTrajectoryLinesOverlay(ImDrawList* drawList, const glm::mat4& projection, const glm::mat4& view,
                                int framebufferWidth, int framebufferHeight)
{
    const TrajectoryLinesOverlay& overlay = trajectoryLinesOverlay();
    if (!overlay.visible || overlay.paths.empty() || !drawList) return;
    // Points up to the frame shown.
    std::size_t count = overlay.frames.size();
    if (!overlay.wholePath)
        count = static_cast<std::size_t>(std::upper_bound(overlay.frames.begin(), overlay.frames.end(), overlay.currentFrame) - overlay.frames.begin());
    if (count < 2) return;
    const glm::mat4 viewProjection = projection * view;
    std::vector<ImVec2> screen;
    for (std::size_t p = 0; p < overlay.paths.size(); ++p) {
        const auto& path = overlay.paths[p];
        // A colour per path, spread around the hue circle.
        const float hue = std::fmod(static_cast<float>(p) * 0.61803399f, 1.0f);
        const ImU32 colour = ImColor::HSV(hue, 0.75f, 0.85f, 0.9f);
        screen.clear();
        for (std::size_t k = 0; k < std::min(count, path.size()); ++k) {
            const glm::vec4 clip = viewProjection * glm::vec4(path[k], 1.0f);
            if (clip.w <= 1e-6f) {
                if (screen.size() > 1) drawList->AddPolyline(screen.data(), static_cast<int>(screen.size()), colour, ImDrawFlags_None, overlay.width);
                screen.clear();
                continue;
            }
            const glm::vec3 ndc = glm::vec3(clip) / clip.w;
            screen.push_back(ImVec2((ndc.x * 0.5f + 0.5f) * static_cast<float>(framebufferWidth),
                                    (1.0f - (ndc.y * 0.5f + 0.5f)) * static_cast<float>(framebufferHeight)));
        }
        if (screen.size() > 1) drawList->AddPolyline(screen.data(), static_cast<int>(screen.size()), colour, ImDrawFlags_None, overlay.width);
    }
}
