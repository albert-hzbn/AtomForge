#pragma once

#include "ui/ThemeUtils.h"
#include "ui/ResponsiveLayout.h"

namespace dialogLayout
{
inline void section(const char* title)
{
    ImGui::Spacing();
    ImGui::TextColored(themeAccentColor(),"%s",title);
    ImGui::Separator();
    ImGui::Spacing();
}

inline ImVec2 actionSize()
{
    return ImVec2(ImGui::GetFontSize()*6.5f,0);
}

inline bool primaryButton(const char* label, ImVec2 size = ImVec2(0,0))
{
    ImGui::PushStyleColor(ImGuiCol_Button,ImVec4(.10f,.39f,.48f,1));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered,ImVec4(.12f,.48f,.58f,1));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive,ImVec4(.08f,.31f,.39f,1));
    ImGui::PushStyleColor(ImGuiCol_Text,ImVec4(1,1,1,1));
    const bool pressed=responsive::button(label,size);
    ImGui::PopStyleColor(4);
    return pressed;
}
}
