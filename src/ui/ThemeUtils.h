#pragma once

#include "imgui.h"

inline bool isLightTheme()
{
    const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    return (bg.x + bg.y + bg.z) / 3.0f > 0.5f;
}

// Theme-aware status colors — visible on both dark and light backgrounds.

// Green: success / OK / loaded
inline ImVec4 themeStatusGood()
{
    return isLightTheme() ? ImVec4(0.06f, 0.50f, 0.12f, 1.0f)
                          : ImVec4(0.22f, 0.90f, 0.35f, 1.0f);
}

// Amber/orange: warning / in-progress hint
inline ImVec4 themeStatusWarn()
{
    return isLightTheme() ? ImVec4(0.68f, 0.38f, 0.00f, 1.0f)
                          : ImVec4(1.00f, 0.62f, 0.10f, 1.0f);
}

// Dark amber: computing / running (replaces yellow which is invisible on white)
inline ImVec4 themeStatusComputing()
{
    return isLightTheme() ? ImVec4(0.52f, 0.33f, 0.00f, 1.0f)
                          : ImVec4(1.00f, 0.95f, 0.10f, 1.0f);
}

// Blue accent: section headers, panel titles
inline ImVec4 themeAccentColor()
{
    return isLightTheme() ? ImVec4(0.05f, 0.38f, 0.72f, 1.0f)
                          : ImVec4(0.70f, 0.85f, 1.00f, 1.0f);
}

// Subtle blue accent: secondary section headers / sub-labels
inline ImVec4 themeAccentSubtle()
{
    return isLightTheme() ? ImVec4(0.10f, 0.44f, 0.65f, 1.0f)
                          : ImVec4(0.60f, 0.75f, 0.90f, 1.0f);
}

// Teal/cyan: informational text (replaces light-cyan which is invisible on white)
inline ImVec4 themeInfoColor()
{
    return isLightTheme() ? ImVec4(0.00f, 0.44f, 0.62f, 1.0f)
                          : ImVec4(0.50f, 0.90f, 1.00f, 1.0f);
}
