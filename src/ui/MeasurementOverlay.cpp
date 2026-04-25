#include "ui/MeasurementOverlay.h"

#include "ElementData.h"
#include "math/StructureMath.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <limits>

namespace
{
constexpr float kPi = 3.14159265358979323846f;

struct BondStats
{
    int coordinationNumber = 0;
    float averageDistance = 0.0f;
    float minDistance = 0.0f;
    float maxDistance = 0.0f;
};

BondStats computeBondStatsForAtom(int centerInstanceIdx,
                                  int centerBaseIdx,
                                  const SceneBuffers& sceneBuffers,
                                  const Structure& structure,
                                  const std::vector<float>& elementRadii)
{
    BondStats stats;

    if (centerInstanceIdx < 0 || centerInstanceIdx >= (int)sceneBuffers.atomPositions.size() ||
        centerInstanceIdx >= (int)sceneBuffers.atomRadii.size())
    {
        return stats;
    }

    bool canUseBaseAtoms = centerBaseIdx >= 0 && centerBaseIdx < (int)structure.atoms.size();
    const glm::vec3 centerPos = canUseBaseAtoms
        ? glm::vec3((float)structure.atoms[centerBaseIdx].x,
                    (float)structure.atoms[centerBaseIdx].y,
                    (float)structure.atoms[centerBaseIdx].z)
        : sceneBuffers.atomPositions[centerInstanceIdx];

    auto lookupAtomRadius = [&](int baseIdx, int instanceIdx) -> float
    {
        if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
        {
            const int atomicNumber = structure.atoms[baseIdx].atomicNumber;
            if (atomicNumber >= 0 && atomicNumber < (int)elementRadii.size() && elementRadii[atomicNumber] > 0.0f)
                return elementRadii[atomicNumber];
        }

        if (instanceIdx >= 0 && instanceIdx < (int)sceneBuffers.atomRadii.size())
            return sceneBuffers.atomRadii[instanceIdx];

        return 1.0f;
    };

    const float radiusCenter = lookupAtomRadius(centerBaseIdx, centerInstanceIdx);

    bool usePbc = false;
    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    usePbc = tryMakeCellMatrices(structure, cell, invCell);

    float minDistance = std::numeric_limits<float>::max();
    float maxDistance = 0.0f;
    float sumDistance = 0.0f;
    int count = 0;

    if (canUseBaseAtoms)
    {
        for (int j = 0; j < (int)structure.atoms.size(); ++j)
        {
            if (j == centerBaseIdx)
                continue;

            const glm::vec3 otherPos((float)structure.atoms[j].x,
                                     (float)structure.atoms[j].y,
                                     (float)structure.atoms[j].z);
            const glm::vec3 delta = minimumImageDelta(otherPos - centerPos, usePbc, cell, invCell);
            const float distance = glm::length(delta);
            if (distance <= kMinBondDistance)
                continue;

            const float radiusOther = lookupAtomRadius(j, -1);
            const float maxBondDistance = (radiusCenter + radiusOther) * kBondToleranceFactor;
            if (distance > maxBondDistance)
                continue;

            minDistance = std::min(minDistance, distance);
            maxDistance = std::max(maxDistance, distance);
            sumDistance += distance;
            ++count;
        }
    }
    else
    {
        for (int i = 0; i < (int)sceneBuffers.atomPositions.size(); ++i)
        {
            if (i == centerInstanceIdx || i >= (int)sceneBuffers.atomRadii.size())
                continue;

            const glm::vec3 delta = sceneBuffers.atomPositions[i] - centerPos;
            const float distance = glm::length(delta);
            if (distance <= kMinBondDistance)
                continue;

            const int otherBaseIdx = (i < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[i] : -1;
            const float radiusOther = lookupAtomRadius(otherBaseIdx, i);
            const float maxBondDistance = (radiusCenter + radiusOther) * kBondToleranceFactor;
            if (distance > maxBondDistance)
                continue;

            minDistance = std::min(minDistance, distance);
            maxDistance = std::max(maxDistance, distance);
            sumDistance += distance;
            ++count;
        }
    }

    stats.coordinationNumber = count;
    if (count > 0)
    {
        stats.averageDistance = sumDistance / (float)count;
        stats.minDistance = minDistance;
        stats.maxDistance = maxDistance;
    }

    return stats;
}

bool projectToScreen(const glm::vec3& p,
                     const glm::mat4& projection,
                     const glm::mat4& view,
                     int w,
                     int h,
                     float& sx,
                     float& sy)
{
    glm::vec4 clip = projection * view * glm::vec4(p, 1.0f);
    if (clip.w <= 0.0f)
        return false;

    float invW = 1.0f / clip.w;
    float nx = clip.x * invW;
    float ny = clip.y * invW;
    float nz = clip.z * invW;
    if (nx < -1.0f || nx > 1.0f || ny < -1.0f || ny > 1.0f || nz < -1.0f || nz > 1.0f)
        return false;

    sx = (nx * 0.5f + 0.5f) * (float)w;
    sy = (1.0f - (ny * 0.5f + 0.5f)) * (float)h;
    return true;
}

void drawDashedLine(ImDrawList* drawList, ImVec2 a, ImVec2 b, ImU32 col)
{
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1.0f)
        return;

    float ux = dx / len;
    float uy = dy / len;
    for (float t = 0.0f; t < len; t += 12.0f)
    {
        float t2 = std::min(t + 7.0f, len);
        drawList->AddLine(ImVec2(a.x + ux * t, a.y + uy * t),
                          ImVec2(a.x + ux * t2, a.y + uy * t2),
                          col, 2.0f);
    }
}

struct AxisOverlayEntry
{
    glm::vec3 dir;
    const char* label;
    ImU32 color;
};

struct GizmoEdge
{
    int a;
    int b;
};

struct ProjectedPoint
{
    ImVec2 screen = ImVec2(0.0f, 0.0f);
    float depth = 0.0f;
};

ProjectedPoint projectGizmoPoint(const glm::vec3& point,
                                 ImVec2 origin,
                                 float projectionScale)
{
    ProjectedPoint projected;
    projected.screen = ImVec2(
        origin.x + point.x * projectionScale,
        origin.y - point.y * projectionScale);
    projected.depth = point.z;
    return projected;
}

float toGizmoDepth(const glm::vec3& cameraAxisDirection)
{
    // In view space, negative Z is toward viewer. Flip for local gizmo camera.
    return -cameraAxisDirection.z;
}

void draw3DAxis(ImDrawList* drawList,
                ImVec2 origin,
                const AxisOverlayEntry& axis,
                float axisLength,
                float projectionScale)
{
    const glm::vec3 endPoint(axis.dir.x * axisLength,
                             axis.dir.y * axisLength,
                             toGizmoDepth(axis.dir) * axisLength);
    const ProjectedPoint p0 = projectGizmoPoint(glm::vec3(0.0f), origin, projectionScale);
    const ProjectedPoint p1 = projectGizmoPoint(endPoint, origin, projectionScale);

    const float depth01 = glm::clamp((p1.depth + 1.0f) * 0.5f, 0.0f, 1.0f);
    const float thickness = 1.8f + 1.2f * depth01;
    drawList->AddLine(p0.screen, p1.screen, axis.color, thickness);

    const float endpointRadius = 2.5f + 2.0f * depth01;
    drawList->AddCircleFilled(p1.screen, endpointRadius, axis.color, 16);
    drawList->AddCircle(p1.screen, endpointRadius + 1.0f,
                        isLightTheme() ? IM_COL32(30, 30, 40, 200) : IM_COL32(245, 245, 250, 200), 16, 1.0f);

    drawList->AddText(ImVec2(p1.screen.x + 5.5f, p1.screen.y + 4.0f), axis.color, axis.label);
}

void draw3DWireCube(ImDrawList* drawList,
                    ImVec2 origin,
                    const glm::mat3& viewRotation,
                    float halfExtent,
                    float projectionScale)
{
    const std::array<glm::vec3, 8> cubeLocal = {{
        {-halfExtent, -halfExtent, -halfExtent},
        { halfExtent, -halfExtent, -halfExtent},
        { halfExtent,  halfExtent, -halfExtent},
        {-halfExtent,  halfExtent, -halfExtent},
        {-halfExtent, -halfExtent,  halfExtent},
        { halfExtent, -halfExtent,  halfExtent},
        { halfExtent,  halfExtent,  halfExtent},
        {-halfExtent,  halfExtent,  halfExtent}
    }};

    const std::array<GizmoEdge, 12> edges = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    }};

    std::array<ProjectedPoint, 8> projectedPoints;
    for (std::size_t i = 0; i < cubeLocal.size(); ++i)
    {
        const glm::vec3 rotated = viewRotation * cubeLocal[i];
        const glm::vec3 gizmoPoint(rotated.x, rotated.y, toGizmoDepth(rotated));
        projectedPoints[i] = projectGizmoPoint(gizmoPoint, origin, projectionScale);
    }

    std::array<int, 12> edgeOrder = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
    std::sort(edgeOrder.begin(), edgeOrder.end(), [&](int lhs, int rhs) {
        const GizmoEdge& e0 = edges[(std::size_t)lhs];
        const GizmoEdge& e1 = edges[(std::size_t)rhs];
        const float z0 = 0.5f * (projectedPoints[(std::size_t)e0.a].depth + projectedPoints[(std::size_t)e0.b].depth);
        const float z1 = 0.5f * (projectedPoints[(std::size_t)e1.a].depth + projectedPoints[(std::size_t)e1.b].depth);
        return z0 < z1;
    });

    for (int orderIndex : edgeOrder)
    {
        const GizmoEdge& edge = edges[(std::size_t)orderIndex];
        const ProjectedPoint& a = projectedPoints[(std::size_t)edge.a];
        const ProjectedPoint& b = projectedPoints[(std::size_t)edge.b];
        const float zMid = 0.5f * (a.depth + b.depth);
        const float alphaFactor = glm::clamp(0.45f + 0.35f * ((zMid + 1.0f) * 0.5f), 0.25f, 0.9f);
        const int alpha = (int)(alphaFactor * 255.0f);
        drawList->AddLine(a.screen, b.screen,
                          isLightTheme() ? IM_COL32(60, 50, 40, alpha) : IM_COL32(205, 215, 230, alpha), 1.0f);
    }
}

void draw3DOrientationGizmo(ImDrawList* drawList,
                            ImVec2 origin,
                            const std::array<AxisOverlayEntry, 3>& axes)
{
    const float projectionScale = 53.0f;
    const float axisLength = 0.77f;

    std::array<int, 3> axisOrder = {{0, 1, 2}};
    std::sort(axisOrder.begin(), axisOrder.end(), [&](int lhs, int rhs) {
        return axes[(std::size_t)lhs].dir.z > axes[(std::size_t)rhs].dir.z;
    });

    for (int idx : axisOrder)
        draw3DAxis(drawList, origin, axes[(std::size_t)idx], axisLength, projectionScale);
}

void drawGizmoBackdrop(ImDrawList* drawList, ImVec2 origin)
{
    const float backgroundRadius = 53.0f;
    const bool light = isLightTheme();
    drawList->AddCircleFilled(origin, backgroundRadius,
                              light ? IM_COL32(235, 238, 242, 200) : IM_COL32(18, 22, 28, 190), 48);
    drawList->AddCircle(origin, backgroundRadius,
                        light ? IM_COL32(0, 0, 0, 55) : IM_COL32(255, 255, 255, 55), 48, 1.0f);
    drawList->AddCircle(origin, backgroundRadius - 7.0f,
                        light ? IM_COL32(0, 0, 0, 20) : IM_COL32(255, 255, 255, 20), 48, 1.0f);
}

void drawGizmoCenter(ImDrawList* drawList, ImVec2 origin)
{
    drawList->AddCircleFilled(origin, 3.2f,
                              isLightTheme() ? IM_COL32(20, 20, 25, 240) : IM_COL32(240, 240, 245, 240), 16);
}

void drawOrientationAxesGizmo(ImDrawList* drawList,
                              const glm::mat4& view,
                              int viewportHeight,
                              const Structure& structure)
{
    const ImVec2 origin(73.0f, (float)viewportHeight - 73.0f);
    drawGizmoBackdrop(drawList, origin);

    const bool light = isLightTheme();
    const glm::mat3 viewRotation(view);

    glm::vec3 worldAxes[3];
    const char* labels[3];
    if (structure.hasUnitCell)
    {
        worldAxes[0] = glm::normalize(glm::vec3((float)structure.cellVectors[0][0],
                                                 (float)structure.cellVectors[0][1],
                                                 (float)structure.cellVectors[0][2]));
        worldAxes[1] = glm::normalize(glm::vec3((float)structure.cellVectors[1][0],
                                                 (float)structure.cellVectors[1][1],
                                                 (float)structure.cellVectors[1][2]));
        worldAxes[2] = glm::normalize(glm::vec3((float)structure.cellVectors[2][0],
                                                 (float)structure.cellVectors[2][1],
                                                 (float)structure.cellVectors[2][2]));
        labels[0] = "a"; labels[1] = "b"; labels[2] = "c";
    }
    else
    {
        worldAxes[0] = glm::vec3(1.0f, 0.0f, 0.0f);
        worldAxes[1] = glm::vec3(0.0f, 1.0f, 0.0f);
        worldAxes[2] = glm::vec3(0.0f, 0.0f, 1.0f);
        labels[0] = "X"; labels[1] = "Y"; labels[2] = "Z";
    }

    std::array<AxisOverlayEntry, 3> axes = {{
        { glm::normalize(viewRotation * worldAxes[0]), labels[0],
          light ? IM_COL32(200, 50, 50, 255) : IM_COL32(235, 92, 92, 255) },
        { glm::normalize(viewRotation * worldAxes[1]), labels[1],
          light ? IM_COL32(40, 160, 50, 255) : IM_COL32(110, 220, 120, 255) },
        { glm::normalize(viewRotation * worldAxes[2]), labels[2],
          light ? IM_COL32(40, 110, 220, 255) : IM_COL32(110, 175, 255, 255) }
    }};

    draw3DOrientationGizmo(drawList, origin, axes);
    drawGizmoCenter(drawList, origin);
}
} // namespace

void MeasurementOverlayState::clearVisuals()
{
    showDistanceLine = false;
    showAngleLines = false;
}

void processMeasurementRequests(MeasurementOverlayState& state,
                                bool requestMeasureDistance,
                                bool requestMeasureAngle,
                                bool requestAtomInfo,
                                const std::vector<int>& selectedInstanceIndices,
                                const SceneBuffers& sceneBuffers,
                                const Structure& structure,
                                const std::vector<float>& elementRadii)
{
    // Measurements disabled for large structures (CPU caches not available)
    if (sceneBuffers.cpuCachesDisabled)
    {
        std::snprintf(state.distanceMessage, sizeof(state.distanceMessage),
                      "Measurements disabled for large structures (>100k atoms)");
        std::snprintf(state.angleMessage, sizeof(state.angleMessage),
                      "Measurements disabled for large structures (>100k atoms)");
        state.clearVisuals();
        return;
    }

    if (requestMeasureDistance)
    {
        if (selectedInstanceIndices.size() != 2)
        {
            std::snprintf(state.distanceMessage, sizeof(state.distanceMessage),
                          "Select exactly 2 atoms to measure distance.");
        }
        else
        {
            int idxA = selectedInstanceIndices[0];
            int idxB = selectedInstanceIndices[1];

            bool validA = idxA >= 0 && idxA < (int)sceneBuffers.atomPositions.size();
            bool validB = idxB >= 0 && idxB < (int)sceneBuffers.atomPositions.size();

            if (!validA || !validB)
            {
                std::snprintf(state.distanceMessage, sizeof(state.distanceMessage),
                              "Unable to measure distance for current selection.");
                state.showDistanceLine = false;
            }
            else
            {
                glm::vec3 pA = sceneBuffers.atomPositions[idxA];
                glm::vec3 pB = sceneBuffers.atomPositions[idxB];
                float distance = glm::length(pA - pB);

                int baseA = (idxA < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idxA] : -1;
                int baseB = (idxB < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idxB] : -1;
                const char* symA = (baseA >= 0 && baseA < (int)structure.atoms.size()) ? structure.atoms[baseA].symbol.c_str() : "A";
                const char* symB = (baseB >= 0 && baseB < (int)structure.atoms.size()) ? structure.atoms[baseB].symbol.c_str() : "B";

                std::snprintf(state.distanceMessage, sizeof(state.distanceMessage),
                              "Distance %s-%s: %.4f", symA, symB, distance);
                state.showDistanceLine = true;
                state.distanceLineIdx0 = idxA;
                state.distanceLineIdx1 = idxB;
            }
        }
        state.showDistancePopup = true;
    }

    if (requestMeasureAngle)
    {
        if (selectedInstanceIndices.size() != 3)
        {
            std::snprintf(state.angleMessage, sizeof(state.angleMessage),
                          "Select exactly 3 atoms to measure angle.");
            state.showAngleLines = false;
        }
        else
        {
            int idx0 = selectedInstanceIndices[0];
            int idx1 = selectedInstanceIndices[1];
            int idx2 = selectedInstanceIndices[2];

            bool ok = idx0 >= 0 && idx0 < (int)sceneBuffers.atomPositions.size() &&
                      idx1 >= 0 && idx1 < (int)sceneBuffers.atomPositions.size() &&
                      idx2 >= 0 && idx2 < (int)sceneBuffers.atomPositions.size();
            if (!ok)
            {
                std::snprintf(state.angleMessage, sizeof(state.angleMessage),
                              "Unable to measure angle for current selection.");
                state.showAngleLines = false;
            }
            else
            {
                glm::vec3 p0 = sceneBuffers.atomPositions[idx0];
                glm::vec3 p1 = sceneBuffers.atomPositions[idx1];
                glm::vec3 p2 = sceneBuffers.atomPositions[idx2];

                glm::vec3 v0 = glm::normalize(p0 - p1);
                glm::vec3 v2 = glm::normalize(p2 - p1);
                float cosA = glm::clamp(glm::dot(v0, v2), -1.0f, 1.0f);
                float angleDeg = glm::degrees(std::acos(cosA));

                int b0 = (idx0 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx0] : -1;
                int b1 = (idx1 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx1] : -1;
                int b2 = (idx2 < (int)sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[idx2] : -1;
                const char* s0 = (b0 >= 0 && b0 < (int)structure.atoms.size()) ? structure.atoms[b0].symbol.c_str() : "A";
                const char* s1 = (b1 >= 0 && b1 < (int)structure.atoms.size()) ? structure.atoms[b1].symbol.c_str() : "B";
                const char* s2 = (b2 >= 0 && b2 < (int)structure.atoms.size()) ? structure.atoms[b2].symbol.c_str() : "C";

                std::snprintf(state.angleMessage, sizeof(state.angleMessage),
                              "Angle %s-%s-%s: %.2f deg", s0, s1, s2, angleDeg);
                state.showAngleLines = true;
                state.angleLineIdx0 = idx0;
                state.angleLineIdx1 = idx1;
                state.angleLineIdx2 = idx2;
            }
        }
        state.showAnglePopup = true;
    }

    if (requestAtomInfo)
    {
        if (selectedInstanceIndices.size() != 1)
        {
            std::snprintf(state.atomInfoMessage, sizeof(state.atomInfoMessage),
                          "Select exactly 1 atom to view info.");
        }
        else
        {
            int idx = selectedInstanceIndices[0];
            bool valid = idx >= 0 && idx < (int)sceneBuffers.atomPositions.size() &&
                         idx < (int)sceneBuffers.atomIndices.size();
            if (!valid)
            {
                std::snprintf(state.atomInfoMessage, sizeof(state.atomInfoMessage),
                              "Unable to retrieve atom info.");
            }
            else
            {
                int baseIdx = sceneBuffers.atomIndices[idx];
                glm::vec3 pos = sceneBuffers.atomPositions[idx];
                int len = 0;

                if (baseIdx >= 0 && baseIdx < (int)structure.atoms.size())
                {
                    const AtomSite& atom = structure.atoms[baseIdx];
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Element:  %s (%s)\n",
                                         elementName(atom.atomicNumber), atom.symbol.c_str());
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Atomic number:  %d\n", atom.atomicNumber);
                }

                len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                     "Cartesian:  (%.6f, %.6f, %.6f) A\n",
                                     pos.x, pos.y, pos.z);

                if (structure.hasUnitCell)
                {
                    glm::mat3 cellMat(1.0f);
                    glm::mat3 invCellMat(1.0f);
                    if (tryMakeCellMatrices(structure, cellMat, invCellMat))
                    {
                        glm::vec3 origin((float)structure.cellOffset[0],
                                         (float)structure.cellOffset[1],
                                         (float)structure.cellOffset[2]);
                        glm::vec3 frac = invCellMat * (pos - origin);
                        std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                      "Direct:  (%.6f, %.6f, %.6f)",
                                      frac.x, frac.y, frac.z);
                        len = (int)std::strlen(state.atomInfoMessage);
                    }
                }

                BondStats stats = computeBondStatsForAtom(idx, baseIdx, sceneBuffers, structure, elementRadii);
                len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                     "\nCoordination number:  %d\n",
                                     stats.coordinationNumber);

                if (stats.coordinationNumber > 0)
                {
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Avg bond length:  %.4f A\n",
                                         stats.averageDistance);
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Min bond length:  %.4f A\n",
                                         stats.minDistance);
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Max bond length:  %.4f A",
                                         stats.maxDistance);
                }
                else
                {
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Avg bond length:  N/A\n");
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Min bond length:  N/A\n");
                    len += std::snprintf(state.atomInfoMessage + len, sizeof(state.atomInfoMessage) - len,
                                         "Max bond length:  N/A");
                }
            }
        }
        state.showAtomInfoPopup = true;
    }
}

void drawMeasurementPopups(MeasurementOverlayState& state)
{
    if (state.showAnglePopup)
    {
        ImGui::OpenPopup("Measure Angle");
        state.showAnglePopup = false;
    }

    ImGui::SetNextWindowSize(ImVec2(540.0f, 0.0f), ImGuiCond_Appearing);
    bool measureAngleOpen = true;
    if (ImGui::BeginPopupModal("Measure Angle", &measureAngleOpen))
    {
        ImGui::TextWrapped("%s", state.angleMessage);
        ImGui::EndPopup();
    }
    if (!measureAngleOpen)
    {
        state.showAngleLines = false;
        ImGui::CloseCurrentPopup();
    }

    if (state.showDistancePopup)
    {
        ImGui::OpenPopup("Measure Distance");
        state.showDistancePopup = false;
    }

    ImGui::SetNextWindowSize(ImVec2(540.0f, 0.0f), ImGuiCond_Appearing);
    bool measureDistanceOpen = true;
    if (ImGui::BeginPopupModal("Measure Distance", &measureDistanceOpen))
    {
        ImGui::TextWrapped("%s", state.distanceMessage);
        ImGui::EndPopup();
    }
    if (!measureDistanceOpen)
    {
        state.showDistanceLine = false;
        ImGui::CloseCurrentPopup();
    }

    if (state.showAtomInfoPopup)
    {
        ImGui::OpenPopup("Atom Info");
        state.showAtomInfoPopup = false;
    }

    ImGui::SetNextWindowSize(ImVec2(480.0f, 0.0f), ImGuiCond_Appearing);
    bool atomInfoOpen = true;
    if (ImGui::BeginPopupModal("Atom Info", &atomInfoOpen))
    {
        ImGui::TextWrapped("%s", state.atomInfoMessage);
        ImGui::EndPopup();
    }
    if (!atomInfoOpen)
        ImGui::CloseCurrentPopup();
}

void drawMeasurementOverlays(const MeasurementOverlayState& state,
                             ImDrawList* drawList,
                             const glm::mat4& projection,
                             const glm::mat4& view,
                             int viewportWidth,
                             int viewportHeight,
                             const SceneBuffers& sceneBuffers)
{
    if (state.showDistanceLine)
    {
        int idxA = state.distanceLineIdx0;
        int idxB = state.distanceLineIdx1;
        bool validA = idxA >= 0 && idxA < (int)sceneBuffers.atomPositions.size();
        bool validB = idxB >= 0 && idxB < (int)sceneBuffers.atomPositions.size();

        if (validA && validB)
        {
            float ax = 0.0f, ay = 0.0f, bx = 0.0f, by = 0.0f;
            if (projectToScreen(sceneBuffers.atomPositions[idxA], projection, view, viewportWidth, viewportHeight, ax, ay) &&
                projectToScreen(sceneBuffers.atomPositions[idxB], projection, view, viewportWidth, viewportHeight, bx, by))
            {
                float dx = bx - ax;
                float dy = by - ay;
                float len = std::sqrt(dx * dx + dy * dy);
                if (len > 1.0f)
                {
                    const float dashLen = 7.0f;
                    const float gapLen = 5.0f;
                    float ux = dx / len;
                    float uy = dy / len;

                    for (float t = 0.0f; t < len; t += (dashLen + gapLen))
                    {
                        float t2 = std::min(t + dashLen, len);
                        ImVec2 p0(ax + ux * t, ay + uy * t);
                        ImVec2 p1(ax + ux * t2, ay + uy * t2);
                        drawList->AddLine(p0, p1, IM_COL32(255, 255, 80, 230), 2.0f);
                    }
                }
            }
        }
    }

    if (state.showAngleLines)
    {
        bool ok = state.angleLineIdx0 >= 0 && state.angleLineIdx0 < (int)sceneBuffers.atomPositions.size() &&
                  state.angleLineIdx1 >= 0 && state.angleLineIdx1 < (int)sceneBuffers.atomPositions.size() &&
                  state.angleLineIdx2 >= 0 && state.angleLineIdx2 < (int)sceneBuffers.atomPositions.size();

        if (ok)
        {
            float x0, y0, x1, y1, x2, y2;
            bool v0ok = projectToScreen(sceneBuffers.atomPositions[state.angleLineIdx0], projection, view, viewportWidth, viewportHeight, x0, y0);
            bool v1ok = projectToScreen(sceneBuffers.atomPositions[state.angleLineIdx1], projection, view, viewportWidth, viewportHeight, x1, y1);
            bool v2ok = projectToScreen(sceneBuffers.atomPositions[state.angleLineIdx2], projection, view, viewportWidth, viewportHeight, x2, y2);

            ImU32 col = IM_COL32(80, 220, 255, 230);
            if (v0ok && v1ok) drawDashedLine(drawList, ImVec2(x0, y0), ImVec2(x1, y1), col);
            if (v2ok && v1ok) drawDashedLine(drawList, ImVec2(x2, y2), ImVec2(x1, y1), col);

            if (v0ok && v1ok && v2ok)
            {
                float d0x = x0 - x1;
                float d0y = y0 - y1;
                float d2x = x2 - x1;
                float d2y = y2 - y1;
                float len0 = std::sqrt(d0x * d0x + d0y * d0y);
                float len2 = std::sqrt(d2x * d2x + d2y * d2y);
                float r = std::min({len0, len2, 30.0f}) * 0.35f;
                if (r > 4.0f)
                {
                    float a0 = std::atan2(d0y, d0x);
                    float a2 = std::atan2(d2y, d2x);
                    float da = a2 - a0;
                    while (da > kPi) da -= 2.0f * kPi;
                    while (da < -kPi) da += 2.0f * kPi;
                    const int segs = 20;
                    for (int s = 0; s < segs; ++s)
                    {
                        float t0 = a0 + da * ((float)s / segs);
                        float t1 = a0 + da * ((float)(s + 1) / segs);
                        drawList->AddLine(ImVec2(x1 + r * std::cos(t0), y1 + r * std::sin(t0)),
                                          ImVec2(x1 + r * std::cos(t1), y1 + r * std::sin(t1)),
                                          col, 1.5f);
                    }
                }
            }
        }
    }
}

void drawElementLabelsOverlay(ImDrawList* drawList,
                              const glm::mat4& projection,
                              const glm::mat4& view,
                              int viewportWidth,
                              int viewportHeight,
                              const SceneBuffers& sceneBuffers,
                              const Structure& structure)
{
    if (!drawList || structure.atoms.empty() || sceneBuffers.atomPositions.empty())
        return;

    struct LabelCandidate
    {
        int baseIdx = -1;
        float sx = 0.0f;
        float sy = 0.0f;
        float depth = 0.0f;
    };

    auto intersects = [](const ImVec4& a, const ImVec4& b) {
        return a.x < b.z && a.z > b.x && a.y < b.w && a.w > b.y;
    };

    const size_t kMaxCandidates = 8000;
    const size_t kMaxLabelsPerFrame = 1200;
    const float kLabelPadding = 2.0f;

    std::vector<LabelCandidate> candidates;
    candidates.reserve(std::min(sceneBuffers.atomPositions.size(), kMaxCandidates));

    for (size_t i = 0; i < sceneBuffers.atomPositions.size() && candidates.size() < kMaxCandidates; ++i)
    {
        const int baseIdx = (i < sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[i] : -1;
        if (baseIdx < 0 || baseIdx >= (int)structure.atoms.size())
            continue;

        const std::string& label = structure.atoms[baseIdx].symbol;
        if (label.empty())
            continue;

        const glm::vec4 viewPos = view * glm::vec4(sceneBuffers.atomPositions[i], 1.0f);
        if (viewPos.z >= 0.0f)
            continue;

        float sx = 0.0f;
        float sy = 0.0f;
        if (!projectToScreen(sceneBuffers.atomPositions[i], projection, view, viewportWidth, viewportHeight, sx, sy))
            continue;

        LabelCandidate c;
        c.baseIdx = baseIdx;
        c.sx = sx;
        c.sy = sy;
        c.depth = -viewPos.z;
        candidates.push_back(c);
    }

    std::sort(candidates.begin(), candidates.end(), [](const LabelCandidate& a, const LabelCandidate& b) {
        return a.depth < b.depth;
    });

    std::vector<ImVec4> occupiedRects;
    occupiedRects.reserve(std::min(candidates.size(), kMaxLabelsPerFrame));

    size_t labelsPlaced = 0;
    for (const LabelCandidate& c : candidates)
    {
        const std::string& label = structure.atoms[c.baseIdx].symbol;
        const ImVec2 textSize = ImGui::CalcTextSize(label.c_str());

        ImVec4 rect(
            c.sx - textSize.x * 0.5f - kLabelPadding,
            c.sy - textSize.y * 0.5f - kLabelPadding,
            c.sx + textSize.x * 0.5f + kLabelPadding,
            c.sy + textSize.y * 0.5f + kLabelPadding);

        if (rect.z < 0.0f || rect.x > (float)viewportWidth ||
            rect.w < 0.0f || rect.y > (float)viewportHeight)
        {
            continue;
        }

        bool overlaps = false;
        for (const ImVec4& placed : occupiedRects)
        {
            if (intersects(rect, placed))
            {
                overlaps = true;
                break;
            }
        }
        if (overlaps)
            continue;

        drawList->AddText(ImVec2(rect.x + kLabelPadding, rect.y + kLabelPadding),
                          isLightTheme() ? IM_COL32(15, 15, 15, 255) : IM_COL32(255, 255, 255, 255),
                          label.c_str());
        occupiedRects.push_back(rect);

        ++labelsPlaced;
        if (labelsPlaced >= kMaxLabelsPerFrame)
            break;
    }
}

void drawOrientationAxesOverlay(ImDrawList* drawList,
                                const glm::mat4& view,
                                int viewportWidth,
                                int viewportHeight,
                                const Structure& structure)
{
    if (!drawList || viewportWidth < 120 || viewportHeight < 120)
        return;

    drawOrientationAxesGizmo(drawList, view, viewportHeight, structure);
}

// ---------------------------------------------------------------------------
// IPF triangle legend (cubic symmetry, Z-direction)
// ---------------------------------------------------------------------------

namespace
{
// Compute IPF color from a crystal direction in the standard triangle.
// h >= k >= l >= 0, normalized.
static ImU32 ipfDirToColor(float h, float k, float l)
{
    float r = h - k;
    float g = k - l;
    float b = l * 1.7320508f;
    float maxC = std::max({r, g, b, 1e-6f});
    r /= maxC; g /= maxC; b /= maxC;
    return IM_COL32((int)(r * 255), (int)(g * 255), (int)(b * 255), 255);
}
} // namespace

void drawIPFTriangleLegend(ImDrawList* drawList,
                           int viewportWidth,
                           int viewportHeight)
{
    if (!drawList || viewportWidth < 200 || viewportHeight < 200)
        return;

    // Triangle size and position (bottom-right corner)
    const float triSize = 140.0f;
    const float margin  = 20.0f;
    const float labelGap = 6.0f;
    const float boxPadX = 3.0f;
    const float boxPadY = 1.0f;

    ImFont* font = ImGui::GetFont();
    const float fontSize = ImGui::GetFontSize();

    const char* title = "IPF-Z (cubic)";
    const char* lbl001 = "[001]";
    const char* lbl011 = "[011]";
    const char* lbl111 = "[111]";

    ImVec2 titleSize = font->CalcTextSizeA(fontSize, FLT_MAX, 0.0f, title);
    ImVec2 sz001 = font->CalcTextSizeA(fontSize, FLT_MAX, 0.0f, lbl001);
    ImVec2 sz011 = font->CalcTextSizeA(fontSize, FLT_MAX, 0.0f, lbl011);
    ImVec2 sz111 = font->CalcTextSizeA(fontSize, FLT_MAX, 0.0f, lbl111);

    // Reserve enough room for the right-side label and centered title so nothing clips.
    const float titleOverhang = std::max(0.0f, (titleSize.x - triSize) * 0.5f);
    const float rightReserve = std::max(sz011.x + boxPadX * 2.0f, titleOverhang + 4.0f);

    // The standard triangle vertices in a 2D layout:
    //   [001] at bottom-left   -> Red
    //   [011] at bottom-right  -> Green
    //   [111] at top           -> Blue
    //
    // Stereographic projection positions (approximate):
    //   [001] = (0, 0)
    //   [011] = (1, 0)    (on equator, 45° from [001])
    //   [111] = (0.5, ~0.61)  (above midpoint)
    //
    // We'll use a simplified right triangle layout for clarity:
    //   [001] bottom-left, [011] bottom-right, [111] top-right area

    const float ox = (float)viewportWidth  - margin - rightReserve - triSize;
    const float oy = (float)viewportHeight - margin;

    // Triangle corners in screen space:
    // V0 = [001] bottom-left
    // V1 = [011] bottom-right
    // V2 = [111] top (we offset it right and up for the standard triangle shape)
    const ImVec2 v0(ox,              oy);                    // [001]
    const ImVec2 v1(ox + triSize,    oy);                    // [011]
    const ImVec2 v2(ox + triSize * 0.58f, oy - triSize * 0.82f); // [111]

    // Draw filled triangle with colored sub-triangles for smooth gradient
    // We sample the triangle on a grid and draw small colored triangles.
    const int kSteps = 30;
    for (int i = 0; i < kSteps; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            // Barycentric coordinates for grid point
            float u0 = 1.0f - (float)i / kSteps;
            float u1 = (float)j / kSteps;
            float u2 = 1.0f - u0 - u1;
            if (u2 < -0.001f) continue;

            // Next row points
            float u0b = 1.0f - (float)(i + 1) / kSteps;
            float u1b = (float)j / kSteps;
            float u2b = 1.0f - u0b - u1b;

            float u1c = (float)(j + 1) / kSteps;
            float u2c = 1.0f - u0b - u1c;

            // Map barycentric to crystal direction:
            // [001] = (0,0,1), [011] = (0,1,1)/sqrt2, [111] = (1,1,1)/sqrt3
            auto baryToDir = [](float b0, float b1, float b2, float& h, float& k, float& l) {
                // Interpolate in crystal direction space
                h = b2 * (1.0f / 1.7320508f);                      // [111] contribution
                k = b1 * (1.0f / 1.4142136f) + b2 * (1.0f / 1.7320508f); // [011] + [111]
                l = b0 * 1.0f + b1 * (1.0f / 1.4142136f) + b2 * (1.0f / 1.7320508f);
                // Sort: l >= k >= h to match convention (largest first)
                if (l < k) std::swap(l, k);
                if (l < h) std::swap(l, h);
                if (k < h) std::swap(k, h);
                float len = std::sqrt(h * h + k * k + l * l);
                if (len > 1e-10f) { h /= len; k /= len; l /= len; }
            };

            auto baryToScreen = [&](float b0, float b1, float b2) -> ImVec2 {
                return ImVec2(b0 * v0.x + b1 * v1.x + b2 * v2.x,
                              b0 * v0.y + b1 * v1.y + b2 * v2.y);
            };

            // Triangle A: (i,j) -> (i+1,j) -> (i+1,j+1)
            if (u2b >= -0.001f && u2c >= -0.001f)
            {
                ImVec2 pA = baryToScreen(u0, u1, std::max(u2, 0.0f));
                ImVec2 pB = baryToScreen(u0b, u1b, std::max(u2b, 0.0f));
                ImVec2 pC = baryToScreen(u0b, u1c, std::max(u2c, 0.0f));

                float h, k, l;
                float bc0 = (u0 + u0b + u0b) / 3.0f;
                float bc1 = (u1 + u1b + u1c) / 3.0f;
                float bc2 = 1.0f - bc0 - bc1;
                baryToDir(bc0, bc1, std::max(bc2, 0.0f), h, k, l);
                ImU32 col = ipfDirToColor(l, k, h);

                drawList->AddTriangleFilled(pA, pB, pC, col);
            }

            // Triangle B (upward): (i,j) -> (i,j+1) -> (i+1,j+1)  [only if valid]
            if (j < i)
            {
                float u1d = (float)(j + 1) / kSteps;
                float u2d = 1.0f - u0 - u1d;
                if (u2d >= -0.001f && u2c >= -0.001f)
                {
                    ImVec2 pA = baryToScreen(u0, u1, std::max(u2, 0.0f));
                    ImVec2 pB = baryToScreen(u0, u1d, std::max(u2d, 0.0f));
                    ImVec2 pC = baryToScreen(u0b, u1c, std::max(u2c, 0.0f));

                    float h, k, l;
                    float bc0 = (u0 + u0 + u0b) / 3.0f;
                    float bc1 = (u1 + u1d + u1c) / 3.0f;
                    float bc2 = 1.0f - bc0 - bc1;
                    baryToDir(bc0, bc1, std::max(bc2, 0.0f), h, k, l);
                    ImU32 col = ipfDirToColor(l, k, h);

                    drawList->AddTriangleFilled(pA, pB, pC, col);
                }
            }
        }
    }

    // Draw border
    const bool light = isLightTheme();
    drawList->AddTriangle(v0, v1, v2, light ? IM_COL32(60, 60, 60, 255) : IM_COL32(200, 200, 200, 255), 2.0f);

    // Draw vertex labels
    ImVec2 pos111(v2.x - sz111.x * 0.5f, v2.y - sz111.y - labelGap);
    ImVec2 titlePos(ox + triSize * 0.5f - titleSize.x * 0.5f,
                    pos111.y - titleSize.y - labelGap);
    const ImU32 labelBg = light ? IM_COL32(240, 240, 240, 210) : IM_COL32(30, 30, 30, 200);
    const ImU32 labelFg = light ? IM_COL32(30, 30, 30, 255) : IM_COL32(220, 220, 220, 255);
    drawList->AddRectFilled(
        ImVec2(titlePos.x - 4, titlePos.y - 2),
        ImVec2(titlePos.x + titleSize.x + 4, titlePos.y + titleSize.y + 2),
        labelBg, 3.0f);
    drawList->AddText(titlePos, labelFg, title);

    // Keep the bottom labels above the base edge so they stay inside the legend area.
    ImVec2 pos001(v0.x, v0.y - sz001.y - labelGap);
    drawList->AddRectFilled(
        ImVec2(pos001.x - boxPadX, pos001.y - boxPadY),
        ImVec2(pos001.x + sz001.x + boxPadX, pos001.y + sz001.y + boxPadY),
        labelBg, 2.0f);
    drawList->AddText(pos001, IM_COL32(255, 80, 80, 255), lbl001);

    ImVec2 pos011(v1.x - sz011.x, v1.y - sz011.y - labelGap);
    drawList->AddRectFilled(
        ImVec2(pos011.x - boxPadX, pos011.y - boxPadY),
        ImVec2(pos011.x + sz011.x + boxPadX, pos011.y + sz011.y + boxPadY),
        labelBg, 2.0f);
    drawList->AddText(pos011, IM_COL32(80, 255, 80, 255), lbl011);

    drawList->AddRectFilled(
        ImVec2(pos111.x - boxPadX, pos111.y - boxPadY),
        ImVec2(pos111.x + sz111.x + boxPadX, pos111.y + sz111.y + boxPadY),
        labelBg, 2.0f);
    drawList->AddText(pos111, IM_COL32(100, 100, 255, 255), lbl111);
}
