#pragma once

#include "imgui.h"

// Dimensions passed to dp/size are design pixels at the 18 px base font.
// Dimensions obtained from ImGui are already scaled and must not be scaled again.
namespace responsive
{
float scale();
float dp(float value);
ImVec2 size(float width, float height);
void windowSize(ImVec2 designSize, ImGuiCond condition = ImGuiCond_None);
void windowConstraints(ImVec2 minimum, ImVec2 maximum);
bool begin(const char* name, bool* open = nullptr, ImGuiWindowFlags flags = 0);
bool beginModal(const char* name, bool* open = nullptr, ImGuiWindowFlags flags = 0);
bool beginChild(const char* name, ImVec2 size = ImVec2(0, 0),
                ImGuiChildFlags childFlags = 0, ImGuiWindowFlags flags = 0);
bool button(const char* label, ImVec2 size = ImVec2(0, 0));
// Keep the controls legible; previews give up width first, then stack above them.
bool stacked(float minimumWidth = 860.0f);
float previewWidth(float preferred, float controlsWidth = 400.0f);
float panelHeight(float preferred, float footerRows = 2.0f);
void nextPanel(bool stack);
}
