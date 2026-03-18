#include "ui/MeasurementOverlay.h"

#include "ElementData.h"
#include "math/StructureMath.h"
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
                                  const Structure& structure)
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

    const float radiusCenter = sceneBuffers.atomRadii[centerInstanceIdx];

    std::vector<float> baseRadii;
    baseRadii.assign(structure.atoms.size(), -1.0f);
    const int pairCount = std::min((int)sceneBuffers.atomIndices.size(), (int)sceneBuffers.atomRadii.size());
    for (int i = 0; i < pairCount; ++i)
    {
        int b = sceneBuffers.atomIndices[i];
        if (b >= 0 && b < (int)baseRadii.size() && baseRadii[b] < 0.0f)
            baseRadii[b] = sceneBuffers.atomRadii[i];
    }

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

            const float radiusOther = (j >= 0 && j < (int)baseRadii.size() && baseRadii[j] > 0.0f)
                ? baseRadii[j]
                : 1.0f;
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

            const float radiusOther = sceneBuffers.atomRadii[i];
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

void drawAxisArrow(ImDrawList* drawList,
                   ImVec2 origin,
                   const AxisOverlayEntry& axis,
                   float baseLength)
{
    glm::vec2 screenDir(axis.dir.x, -axis.dir.y);
    float screenLen = glm::length(screenDir);
    if (screenLen < 1e-4f)
        screenDir = glm::vec2(0.0f, -1.0f);
    else
        screenDir /= screenLen;

    const float towardViewer = std::max(0.0f, -axis.dir.z);
    const float axisLength = baseLength + 7.0f * towardViewer;
    const float headLength = 8.0f;
    const float headWidth = 4.0f;

    ImVec2 end(origin.x + screenDir.x * axisLength,
               origin.y + screenDir.y * axisLength);
    ImVec2 headBase(end.x - screenDir.x * headLength,
                    end.y - screenDir.y * headLength);
    ImVec2 perp(-screenDir.y, screenDir.x);

    drawList->AddLine(origin, headBase, axis.color, 2.5f);
    drawList->AddTriangleFilled(
        end,
        ImVec2(headBase.x + perp.x * headWidth, headBase.y + perp.y * headWidth),
        ImVec2(headBase.x - perp.x * headWidth, headBase.y - perp.y * headWidth),
        axis.color);

    ImVec2 labelOffset(screenDir.x * 8.0f + perp.x * 3.0f,
                       screenDir.y * 8.0f + perp.y * 3.0f);
    drawList->AddText(ImVec2(end.x + labelOffset.x, end.y + labelOffset.y), axis.color, axis.label);
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
                                const Structure& structure)
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

                BondStats stats = computeBondStatsForAtom(idx, baseIdx, sceneBuffers, structure);
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

    // Guard against excessive text allocations when many periodic instances
    // are visible.
    const size_t kMaxLabelsPerFrame = 5000;

    size_t labelCount = 0;
    for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
    {
        int baseIdx = (i < sceneBuffers.atomIndices.size()) ? sceneBuffers.atomIndices[i] : -1;
        if (baseIdx < 0 || baseIdx >= (int)structure.atoms.size())
            continue;

        float sx, sy;
        if (!projectToScreen(sceneBuffers.atomPositions[i], projection, view, viewportWidth, viewportHeight, sx, sy))
            continue;

        const std::string& label = structure.atoms[baseIdx].symbol;
        ImVec2 textSize = ImGui::CalcTextSize(label.c_str());
        drawList->AddText(ImVec2(sx - textSize.x * 0.5f, sy - textSize.y * 0.5f),
                          IM_COL32(255, 255, 255, 255),
                          label.c_str());

        ++labelCount;
        if (labelCount >= kMaxLabelsPerFrame)
            break;
    }
}

void drawOrientationAxesOverlay(ImDrawList* drawList,
                                const glm::mat4& view,
                                int viewportWidth,
                                int viewportHeight)
{
    if (!drawList || viewportWidth < 120 || viewportHeight < 120)
        return;

    const ImVec2 origin(58.0f, (float)viewportHeight - 58.0f);
    const float backgroundRadius = 34.0f;
    const float baseAxisLength = 24.0f;

    drawList->AddCircleFilled(origin, backgroundRadius, IM_COL32(20, 24, 30, 170), 32);
    drawList->AddCircle(origin, backgroundRadius, IM_COL32(255, 255, 255, 45), 32, 1.0f);

    const glm::mat3 viewRotation(view);
    std::array<AxisOverlayEntry, 3> axes = {{
        { glm::normalize(viewRotation * glm::vec3(1.0f, 0.0f, 0.0f)), "X", IM_COL32(235, 92, 92, 255) },
        { glm::normalize(viewRotation * glm::vec3(0.0f, 1.0f, 0.0f)), "Y", IM_COL32(110, 220, 120, 255) },
        { glm::normalize(viewRotation * glm::vec3(0.0f, 0.0f, 1.0f)), "Z", IM_COL32(110, 175, 255, 255) }
    }};

    std::sort(axes.begin(), axes.end(), [](const AxisOverlayEntry& a, const AxisOverlayEntry& b) {
        return a.dir.z > b.dir.z;
    });

    for (std::size_t i = 0; i < axes.size(); ++i)
        drawAxisArrow(drawList, origin, axes[i], baseAxisLength);

    drawList->AddCircleFilled(origin, 3.5f, IM_COL32(240, 240, 245, 240), 16);
}
