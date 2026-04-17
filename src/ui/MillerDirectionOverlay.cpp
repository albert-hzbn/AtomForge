#include "ui/MillerDirectionOverlay.h"

#include "imgui.h"
#include "math/StructureMath.h"

#include <cmath>
#include <cstdio>

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
    const float nx   = clip.x * invW;
    const float ny   = clip.y * invW;

    out.x = (nx * 0.5f + 0.5f) * (float)w;
    out.y = (1.0f - (ny * 0.5f + 0.5f)) * (float)h;
    return true;
}

// Clip a world-space segment to the camera near plane (returns false if fully clipped).
bool clipSegmentNear(glm::vec3& p0, glm::vec3& p1, const glm::mat4& view, float nearDepth)
{
    glm::vec4 v0 = view * glm::vec4(p0, 1.0f);
    glm::vec4 v1 = view * glm::vec4(p1, 1.0f);

    const bool in0 = v0.z <= -nearDepth;
    const bool in1 = v1.z <= -nearDepth;

    if (!in0 && !in1)
        return false;

    if (in0 != in1)
    {
        const float denom = v0.z - v1.z;
        if (std::abs(denom) > 1e-6f)
        {
            const float t = (v0.z + nearDepth) / denom;
            const glm::vec3 mid = p0 + (p1 - p0) * glm::clamp(t, 0.0f, 1.0f);
            if (!in0) p0 = mid;
            else      p1 = mid;
        }
    }
    return true;
}

void drawArrowHead(ImDrawList* drawList,
                   ImVec2 tip,
                   ImVec2 base,
                   float headSize,
                   ImU32 color)
{
    // Direction from base to tip
    const float dx = tip.x - base.x;
    const float dy = tip.y - base.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-3f)
        return;

    const float nx = dx / len;
    const float ny = dy / len;
    // Perpendicular
    const float px = -ny;
    const float py =  nx;

    const ImVec2 left (tip.x - nx * headSize + px * headSize * 0.45f,
                       tip.y - ny * headSize + py * headSize * 0.45f);
    const ImVec2 right(tip.x - nx * headSize - px * headSize * 0.45f,
                       tip.y - ny * headSize - py * headSize * 0.45f);

    drawList->AddTriangleFilled(tip, left, right, color);
}
} // anonymous namespace

void drawMillerDirectionsOverlay(ImDrawList* drawList,
                                 const glm::mat4& projection,
                                 const glm::mat4& view,
                                 int framebufferWidth,
                                 int framebufferHeight,
                                 const Structure& structure,
                                 const std::vector<MillerDirection>& directions,
                                 bool enabled)
{
    if (!enabled || !structure.hasUnitCell || drawList == nullptr || directions.empty())
        return;

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    if (!tryMakeCellMatrices(structure, cell, invCell))
        return;
    (void)invCell;

    // Arrow origin: cell corner (cellOffset)
    const glm::vec3 origin((float)structure.cellOffset[0],
                           (float)structure.cellOffset[1],
                           (float)structure.cellOffset[2]);

    for (size_t i = 0; i < directions.size(); ++i)
    {
        const MillerDirection& dir = directions[i];
        if (!dir.visible)
            continue;
        if (dir.u == 0 && dir.v == 0 && dir.w == 0)
            continue;

        // Cartesian direction vector
        const glm::vec3 uvw((float)dir.u, (float)dir.v, (float)dir.w);
        const glm::vec3 cartDir = cell * uvw;
        const float cartLen = glm::length(cartDir);
        if (cartLen < 1e-6f)
            continue;
        const glm::vec3 cartNorm = cartDir / cartLen;

        glm::vec3 p0 = origin;
        glm::vec3 p1 = origin + cartNorm * dir.length;

        if (!clipSegmentNear(p0, p1, view, 0.03f))
            continue;

        ImVec2 screen0, screen1;
        if (!projectToScreen(p0, projection, view, framebufferWidth, framebufferHeight, screen0))
            continue;
        if (!projectToScreen(p1, projection, view, framebufferWidth, framebufferHeight, screen1))
            continue;

        const ImU32 lineColor = ImColor(dir.color[0], dir.color[1], dir.color[2], 1.0f);

        // Shaft
        drawList->AddLine(screen0, screen1, lineColor, 2.5f);

        // Arrowhead at tip
        drawArrowHead(drawList, screen1, screen0, 10.0f, lineColor);

        // Direction label: "[uvw]"
        char label[32];
        std::snprintf(label, sizeof(label), "[%d%d%d]", dir.u, dir.v, dir.w);

        const float offsetX = 6.0f;
        const float offsetY = -4.0f;
        const ImVec2 labelPos(screen1.x + offsetX, screen1.y + offsetY);
        // Shadow for readability
        drawList->AddText(ImVec2(labelPos.x + 1.0f, labelPos.y + 1.0f),
                          IM_COL32(0, 0, 0, 180), label);
        drawList->AddText(labelPos, lineColor, label);
    }
}
