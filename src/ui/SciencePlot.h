#pragma once

#include "science/ResultPlots.h"
#include "ui/ResponsiveLayout.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// Numeric x/y plot of a scientific result: true x axis, ticks, legend,
// optional log y axis, reference lines and a nearest-point hover readout.
namespace uiPlot
{
inline std::vector<double> niceTicks(double low, double high, int target)
{
    std::vector<double> ticks;
    if (!(high > low)) return ticks;
    const double raw = (high - low) / std::max(1, target);
    const double magnitude = std::pow(10.0, std::floor(std::log10(raw)));
    double step = magnitude;
    for (double m : {1.0, 2.0, 5.0, 10.0})
        if (raw <= m * magnitude) { step = m * magnitude; break; }
    for (double t = std::ceil(low / step) * step; t <= high + step * 1e-9; t += step) ticks.push_back(std::abs(t) < step * 1e-9 ? 0.0 : t);
    return ticks;
}

inline void drawSciencePlot(const char* id, const atomforge::science::PlotSpec& spec, float height = 280.0f)
{
    responsive::beginChild(id, ImVec2(-1, responsive::dp(height)), true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 80.0f); avail.y = std::max(avail.y, 80.0f);
    ImGui::InvisibleButton("##canvas", avail);
    const bool hovered = ImGui::IsItemHovered();
    const bool light = isLightTheme();
    const ImU32 background = light ? IM_COL32(248, 249, 251, 255) : IM_COL32(20, 24, 31, 255);
    const ImU32 frame = light ? IM_COL32(150, 155, 170, 255) : IM_COL32(80, 90, 110, 255);
    const ImU32 grid = light ? IM_COL32(215, 218, 226, 255) : IM_COL32(45, 54, 68, 255);
    const ImU32 text = light ? IM_COL32(45, 45, 60, 255) : IM_COL32(205, 205, 210, 255);
    draw->AddRectFilled(origin, ImVec2(origin.x + avail.x, origin.y + avail.y), background);
    const float lineHeight = ImGui::GetTextLineHeight();
    auto transform = [&](double y) { return spec.logY ? (y > 0 ? std::log10(y) : NAN) : y; };
    double xMin = HUGE_VAL, xMax = -HUGE_VAL, yMin = HUGE_VAL, yMax = -HUGE_VAL;
    for (const auto& series : spec.series)
        for (std::size_t i = 0; i < std::min(series.x.size(), series.y.size()); ++i) {
            const double y = transform(series.y[i]);
            if (!std::isfinite(series.x[i]) || !std::isfinite(y)) continue;
            xMin = std::min(xMin, series.x[i]); xMax = std::max(xMax, series.x[i]);
            yMin = std::min(yMin, y); yMax = std::max(yMax, y);
        }
    for (double h : spec.horizontal) { const double y = transform(h); if (std::isfinite(y)) { yMin = std::min(yMin, y); yMax = std::max(yMax, y); } }
    if (!(xMax >= xMin) || !(yMax >= yMin)) {
        draw->AddText(ImVec2(origin.x + 10, origin.y + 10), text, "No data to plot.");
        ImGui::EndChild();
        return;
    }
    if (xMax == xMin) { xMin -= 0.5; xMax += 0.5; }
    if (yMax == yMin) { const double pad = std::max(1e-12, std::abs(yMin) * 0.05); yMin -= pad; yMax += pad; }
    const double yPad = (yMax - yMin) * 0.06;
    yMin -= yPad; yMax += yPad;
    const ImVec2 gMin(origin.x + responsive::dp(64), origin.y + lineHeight * 1.6f);
    const ImVec2 gMax(origin.x + avail.x - responsive::dp(14), origin.y + avail.y - lineHeight * 2.4f);
    auto px = [&](double x) { return gMin.x + static_cast<float>((x - xMin) / (xMax - xMin)) * (gMax.x - gMin.x); };
    auto py = [&](double y) { return gMax.y - static_cast<float>((y - yMin) / (yMax - yMin)) * (gMax.y - gMin.y); };
    char label[48];
    for (double t : niceTicks(xMin, xMax, 6)) {
        draw->AddLine(ImVec2(px(t), gMin.y), ImVec2(px(t), gMax.y), grid);
        std::snprintf(label, sizeof(label), "%.4g", t);
        draw->AddText(ImVec2(px(t) - ImGui::CalcTextSize(label).x / 2, gMax.y + 3), text, label);
    }
    for (double t : niceTicks(yMin, yMax, 5)) {
        draw->AddLine(ImVec2(gMin.x, py(t)), ImVec2(gMax.x, py(t)), grid);
        std::snprintf(label, sizeof(label), "%.4g", spec.logY ? std::pow(10.0, t) : t);
        draw->AddText(ImVec2(gMin.x - ImGui::CalcTextSize(label).x - 4, py(t) - lineHeight / 2), text, label);
    }
    draw->AddRect(gMin, gMax, frame);
    draw->AddText(ImVec2(gMin.x, origin.y + 2), text, (spec.title + "   [" + spec.yLabel + "]").c_str());
    draw->AddText(ImVec2((gMin.x + gMax.x) / 2 - ImGui::CalcTextSize(spec.xLabel.c_str()).x / 2, gMax.y + lineHeight + 4), text, spec.xLabel.c_str());
    draw->PushClipRect(gMin, gMax, true);
    for (double h : spec.horizontal) {
        const double y = transform(h);
        if (std::isfinite(y)) draw->AddLine(ImVec2(gMin.x, py(y)), ImVec2(gMax.x, py(y)), frame, 1.0f);
    }
    for (const auto& marker : spec.markers) {
        draw->AddLine(ImVec2(px(marker.x), gMin.y), ImVec2(px(marker.x), gMax.y), frame, 1.0f);
        draw->AddText(ImVec2(px(marker.x) + 2, gMin.y + 2), text, marker.label.c_str());
    }
    double bestDistance = 12.0 * 12.0, bestX = 0, bestY = 0;
    std::string bestName;
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    for (std::size_t s = 0; s < spec.series.size(); ++s) {
        const auto& series = spec.series[s];
        const ImU32 colour = ImColor::HSV(std::fmod(0.58f + s * 0.61803f, 1.0f), 0.7f, light ? 0.6f : 0.95f);
        std::vector<ImVec2> points;
        for (std::size_t i = 0; i < std::min(series.x.size(), series.y.size()); ++i) {
            const double y = transform(series.y[i]);
            if (!std::isfinite(series.x[i]) || !std::isfinite(y)) {
                if (points.size() > 1 && !series.points) draw->AddPolyline(points.data(), static_cast<int>(points.size()), colour, 0, 1.6f);
                points.clear();
                continue;
            }
            const ImVec2 p(px(series.x[i]), py(y));
            points.push_back(p);
            if (series.points) draw->AddCircleFilled(p, 3.0f, colour);
            const double dx = p.x - mouse.x, dy = p.y - mouse.y;
            if (hovered && dx * dx + dy * dy < bestDistance) { bestDistance = dx * dx + dy * dy; bestX = series.x[i]; bestY = series.y[i]; bestName = series.name; }
        }
        if (!series.points && points.size() > 1) draw->AddPolyline(points.data(), static_cast<int>(points.size()), colour, 0, 1.6f);
    }
    draw->PopClipRect();
    if (spec.series.size() > 1) {
        float y = gMin.y + 4;
        for (std::size_t s = 0; s < spec.series.size(); ++s) {
            const ImU32 colour = ImColor::HSV(std::fmod(0.58f + s * 0.61803f, 1.0f), 0.7f, light ? 0.6f : 0.95f);
            const float width = ImGui::CalcTextSize(spec.series[s].name.c_str()).x;
            draw->AddRectFilled(ImVec2(gMax.x - width - 22, y + 4), ImVec2(gMax.x - width - 12, y + lineHeight - 4), colour);
            draw->AddText(ImVec2(gMax.x - width - 6, y), text, spec.series[s].name.c_str());
            y += lineHeight;
        }
    }
    if (!bestName.empty()) ImGui::SetTooltip("%s\n%s = %.6g\n%s = %.6g", bestName.c_str(), spec.xLabel.c_str(), bestX, spec.yLabel.c_str(), bestY);
    ImGui::EndChild();
}
}
