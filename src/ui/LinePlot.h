#pragma once

#include "ui/ResponsiveLayout.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

// Shared hand-drawn multi-series line plot (ImDrawList polylines with grid
// lines and axis labels), used by analysis dialogs whose data has no single
// canonical x-axis unit ImGui::PlotLines can express (multiple overlaid
// series, or a shared but non-uniform x range).
namespace uiPlot
{
inline void drawColoredLines(const char* childId, int numSeries, int numSamples,
    const std::function<double(int series, int sample)>& valueAt, float valueMin, float valueMax, const char* xLabel)
{
    responsive::beginChild(childId, ImVec2(-1, 260), true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 avail = ImGui::GetContentRegionAvail();
    avail.x = std::max(avail.x, 50.0f); avail.y = std::max(avail.y, 50.0f);
    ImGui::InvisibleButton("##canvas", avail);
    const bool light = isLightTheme();
    const ImVec2 pMin = origin, pMax(origin.x + avail.x, origin.y + avail.y);
    draw->AddRectFilled(pMin, pMax, light ? IM_COL32(245,246,248,255) : IM_COL32(20,24,31,255));
    draw->AddRect(pMin, pMax, light ? IM_COL32(140,145,160,255) : IM_COL32(110,125,150,255));
    if (numSamples < 2 || numSeries < 1)
    {
        draw->AddText(ImVec2(pMin.x+10,pMin.y+10), light ? IM_COL32(80,80,100,255) : IM_COL32(220,220,220,255), "No data.");
        ImGui::EndChild();
        return;
    }
    const float leftPad=52, rightPad=16, topPad=16, bottomPad=28;
    const ImVec2 gMin(pMin.x+leftPad,pMin.y+topPad), gMax(pMax.x-rightPad,pMax.y-bottomPad);
    draw->AddRect(gMin, gMax, light ? IM_COL32(150,155,170,255) : IM_COL32(80,90,110,255));
    const ImU32 gridCol = light ? IM_COL32(170,175,190,180) : IM_COL32(50,60,75,255);
    const ImU32 labelCol = light ? IM_COL32(50,50,65,255) : IM_COL32(200,200,200,255);
    for (int i=0;i<=4;++i)
    {
        const float t=i/4.0f, y=gMax.y+(gMin.y-gMax.y)*t;
        draw->AddLine(ImVec2(gMin.x,y),ImVec2(gMax.x,y),gridCol);
        char label[32]; std::snprintf(label,sizeof(label),"%.3g",valueMin+(valueMax-valueMin)*t);
        draw->AddText(ImVec2(pMin.x+4,y-7),labelCol,label);
    }
    draw->AddText(ImVec2(gMin.x,gMax.y+6),labelCol,xLabel);
    const float range = std::max(valueMax-valueMin, 1e-12f);
    for (int series=0; series<numSeries; ++series)
    {
        std::vector<ImVec2> points(static_cast<std::size_t>(numSamples));
        for (int s=0;s<numSamples;++s)
        {
            const float x = gMin.x + (gMax.x-gMin.x) * (static_cast<float>(s)/(numSamples-1));
            const float y = gMax.y - (gMax.y-gMin.y) * ((static_cast<float>(valueAt(series,s))-valueMin)/range);
            points[static_cast<std::size_t>(s)] = ImVec2(x,y);
        }
        const ImColor c = ImColor::HSV(std::fmod(series*0.61803f,1.0f),0.65f,light?0.55f:0.95f);
        draw->AddPolyline(points.data(),numSamples,static_cast<ImU32>(c),0,1.5f);
    }
    ImGui::EndChild();
}
}
