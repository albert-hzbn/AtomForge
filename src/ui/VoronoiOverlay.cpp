#include "ui/VoronoiOverlay.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>

namespace
{
bool projectToScreen(const glm::vec3& p,
                     const glm::mat4& projection,
                     const glm::mat4& view,
                     int w,
                     int h,
                     ImVec2& out)
{
    glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
    if (clip.w <= 1e-5f)
        return false;

    const float invW = 1.0f / clip.w;
    const float nx = clip.x * invW;
    const float ny = clip.y * invW;

    out.x = (nx * 0.5f + 0.5f) * (float)w;
    out.y = (1.0f - (ny * 0.5f + 0.5f)) * (float)h;
    return true;
}

std::vector<glm::vec3> clipPolygonAgainstNearPlane(const std::vector<glm::vec3>& worldPoints,
                                                   const glm::mat4& view,
                                                   float nearDepth)
{
    if (worldPoints.size() < 3)
        return {};

    std::vector<glm::vec3> clipped;
    clipped.reserve(worldPoints.size() + 2);

    for (size_t i = 0; i < worldPoints.size(); ++i)
    {
        const glm::vec3& currentWorld = worldPoints[i];
        const glm::vec3& nextWorld = worldPoints[(i + 1) % worldPoints.size()];

        glm::vec4 currentView = view * glm::vec4(currentWorld, 1.0f);
        glm::vec4 nextView = view * glm::vec4(nextWorld, 1.0f);
        bool currentInside = currentView.z <= -nearDepth;
        bool nextInside = nextView.z <= -nearDepth;

        if (currentInside)
            clipped.push_back(currentWorld);

        if (currentInside != nextInside)
        {
            float denom = nextView.z - currentView.z;
            if (std::abs(denom) > 1e-6f)
            {
                float t = (-nearDepth - currentView.z) / denom;
                t = std::max(0.0f, std::min(1.0f, t));
                clipped.push_back(currentWorld + (nextWorld - currentWorld) * t);
            }
        }
    }

    return (clipped.size() >= 3) ? clipped : std::vector<glm::vec3>();
}
}

void drawVoronoiOverlay(ImDrawList* drawList,
                        const glm::mat4& projection,
                        const glm::mat4& view,
                        int framebufferWidth,
                        int framebufferHeight,
                        const VoronoiDiagram& diagram,
                        const std::vector<int>& selectedIndices,
                        bool enabled)
{
    if (!enabled || drawList == nullptr || diagram.cells.empty())
        return;

    const ImU32 defaultFill = ImColor(1.0f, 1.0f, 1.0f, 0.06f);
    const ImU32 defaultEdge = ImColor(1.0f, 1.0f, 1.0f, 0.7f);
    const ImU32 selectedFill = ImColor(1.0f, 0.2f, 0.2f, 0.15f);
    const ImU32 selectedEdge = ImColor(1.0f, 0.2f, 0.2f, 0.85f);

    // Build a fast lookup for selected cell indices.
    std::vector<bool> isSelected(diagram.cells.size(), false);
    for (int idx : selectedIndices)
    {
        if (idx >= 0 && idx < (int)isSelected.size())
            isSelected[idx] = true;
    }

    for (size_t ci = 0; ci < diagram.cells.size(); ++ci)
    {
        const auto& cell = diagram.cells[ci];
        const ImU32 fillColor = isSelected[ci] ? selectedFill : defaultFill;
        const ImU32 edgeColor = isSelected[ci] ? selectedEdge : defaultEdge;
        for (const auto& face : cell.faces)
        {
            if (face.vertices.size() < 3)
                continue;

            const auto clippedWorld =
                clipPolygonAgainstNearPlane(face.vertices, view, 0.03f);
            if (clippedWorld.size() < 3)
                continue;

            ImVector<ImVec2> screenPoints;
            screenPoints.reserve((int)clippedWorld.size());
            bool allProjected = true;
            for (const auto& v : clippedWorld)
            {
                ImVec2 screen;
                if (!projectToScreen(v, projection, view,
                                     framebufferWidth, framebufferHeight, screen))
                {
                    allProjected = false;
                    break;
                }
                screenPoints.push_back(screen);
            }

            if (!allProjected || screenPoints.Size < 3)
                continue;

            drawList->AddConvexPolyFilled(screenPoints.Data, screenPoints.Size, fillColor);
            drawList->AddPolyline(screenPoints.Data, screenPoints.Size,
                                  edgeColor, ImDrawFlags_Closed, 1.5f);
        }
    }
}
