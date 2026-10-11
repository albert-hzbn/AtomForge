#pragma once
// Shared building blocks for dialogs, so that every window looks and reads
// the same: section headings, callouts (info, success, warning, error), help
// markers, primary/secondary buttons and right-aligned action rows.

#include "ui/Fonts.h"
#include "ui/ThemeUtils.h"
#include "ui/ResponsiveLayout.h"

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <initializer_list>

namespace dialogLayout
{
// The theme's accent (the check-mark colour), used for primary actions.
inline ImVec4 accent() { return ImGui::GetStyle().Colors[ImGuiCol_CheckMark]; }

inline ImVec4 mix(ImVec4 a, ImVec4 b, float t)
{
    return ImVec4(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t);
}

// A section heading: semibold title over a thin rule, with breathing room above.
inline void section(const char* title)
{
    ImGui::Dummy(ImVec2(0, responsive::dp(4)));
    uiFonts::pushHeading();
    ImGui::TextUnformatted(title);
    uiFonts::popHeading();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float width = ImGui::GetContentRegionAvail().x;
    ImGui::GetWindowDrawList()->AddLine(ImVec2(p.x, p.y - ImGui::GetStyle().ItemSpacing.y * 0.5f),
                                        ImVec2(p.x + width, p.y - ImGui::GetStyle().ItemSpacing.y * 0.5f),
                                        ImGui::GetColorU32(ImGuiCol_Separator), 1.0f);
    ImGui::Dummy(ImVec2(0, responsive::dp(2)));
}

// A quieter label for a group of fields inside a section.
inline void caption(const char* text)
{
    ImGui::TextDisabled("%s", text);
}

enum class Note { Info, Success, Warning, Error };

inline ImVec4 noteColour(Note kind)
{
    switch (kind) {
    case Note::Success: return themeStatusGood();
    case Note::Warning: return themeStatusWarn();
    case Note::Error: return themeStatusError();
    case Note::Info: default: return accent();
    }
}

// A callout: tinted, rounded box with a coloured edge, an icon and wrapped text.
inline void note(Note kind, const char* text)
{
    const ImVec4 colour = noteColour(kind);
    const float width = ImGui::GetContentRegionAvail().x;
    const float pad = responsive::dp(9);
    const float iconWidth = responsive::dp(22);
    const ImVec2 textSize = ImGui::CalcTextSize(text, nullptr, false, width - 2 * pad - iconWidth);
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float height = textSize.y + 2 * pad;
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const float rounding = ImGui::GetStyle().FrameRounding;
    draw->AddRectFilled(p, ImVec2(p.x + width, p.y + height), ImGui::GetColorU32(ImVec4(colour.x, colour.y, colour.z, isLightTheme() ? 0.09f : 0.16f)), rounding);
    draw->AddRectFilled(p, ImVec2(p.x + responsive::dp(3), p.y + height), ImGui::GetColorU32(colour), rounding, ImDrawFlags_RoundCornersLeft);
    // The icon: a glyph from the icon font, or a small drawn badge.
    const ImVec2 iconCentre(p.x + pad + iconWidth * 0.40f, p.y + pad + ImGui::GetTextLineHeight() * 0.5f);
    if (uiFonts::icons) {
        const char* glyph = kind == Note::Warning ? Icon::Warning : kind == Note::Error ? Icon::Error : kind == Note::Success ? Icon::Check : Icon::Info;
        const ImVec2 g = ImGui::CalcTextSize(glyph);
        draw->AddText(ImVec2(iconCentre.x - g.x * 0.5f, iconCentre.y - g.y * 0.5f), ImGui::GetColorU32(colour), glyph);
    } else {
        draw->AddCircleFilled(iconCentre, responsive::dp(6.5f), ImGui::GetColorU32(colour));
        const char* mark = kind == Note::Info ? "i" : kind == Note::Success ? "+" : "!";
        const ImVec2 g = ImGui::CalcTextSize(mark);
        draw->AddText(ImVec2(iconCentre.x - g.x * 0.5f, iconCentre.y - g.y * 0.5f), IM_COL32(255, 255, 255, 255), mark);
    }
    ImGui::SetCursorScreenPos(ImVec2(p.x + pad + iconWidth, p.y + pad));
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + width - 2 * pad - iconWidth);
    ImGui::TextUnformatted(text);
    ImGui::PopTextWrapPos();
    ImGui::SetCursorScreenPos(ImVec2(p.x, p.y + height));
    ImGui::Dummy(ImVec2(width, 0));
}

inline void notef(Note kind, const char* format, ...)
{
    char text[1024];
    va_list args;
    va_start(args, format);
    std::vsnprintf(text, sizeof(text), format, args);
    va_end(args);
    note(kind, text);
}

// The banner of a dialog opened on a structure-pipeline step.
inline void stepBanner(const char* step)
{
    notef(Note::Info, "Editing the pipeline step %s. Update step writes these settings into the pipeline; the structure changes through the pipeline.", step);
}

// A "(?)" marker after the previous item with a wrapped tooltip.
inline void helpMarker(const char* text)
{
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

inline ImVec2 actionSize()
{
    return ImVec2(ImGui::GetFontSize() * 6.5f, 0);
}

// The main action of a dialog, in the accent colour.
inline bool primaryButton(const char* label, ImVec2 size = ImVec2(0, 0))
{
    const ImVec4 a = accent();
    const ImVec4 white(1, 1, 1, 1);
    ImGui::PushStyleColor(ImGuiCol_Button, a);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, mix(a, white, 0.14f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, mix(a, ImVec4(0, 0, 0, 1), 0.15f));
    ImGui::PushStyleColor(ImGuiCol_Border, mix(a, ImVec4(0, 0, 0, 1), 0.10f));
    ImGui::PushStyleColor(ImGuiCol_Text, white);
    const bool pressed = responsive::button(label, size);
    ImGui::PopStyleColor(5);
    return pressed;
}

// Secondary actions keep the neutral button style.
inline bool secondaryButton(const char* label, ImVec2 size = ImVec2(0, 0))
{
    return responsive::button(label, size);
}

// A destructive action (delete, discard), in red.
inline bool dangerButton(const char* label, ImVec2 size = ImVec2(0, 0))
{
    const ImVec4 red = noteColour(Note::Error);
    ImGui::PushStyleColor(ImGuiCol_Text, red);
    const bool pressed = responsive::button(label, size);
    ImGui::PopStyleColor();
    return pressed;
}

// Width of a row of buttons with these labels (at least `minimum` wide each).
inline float buttonRowWidth(std::initializer_list<const char*> labels, float minimum = 0)
{
    const ImGuiStyle& style = ImGui::GetStyle();
    float width = 0;
    int count = 0;
    for (const char* label : labels) {
        width += std::max(minimum, ImGui::CalcTextSize(label, nullptr, true).x + 2 * style.FramePadding.x);
        ++count;
    }
    return width + style.ItemSpacing.x * static_cast<float>(count > 0 ? count - 1 : 0);
}

// Starts a footer: a rule, then the cursor placed so that a button row of
// `width` ends at the right edge.
inline void beginFooter(float width)
{
    ImGui::Dummy(ImVec2(0, responsive::dp(2)));
    ImGui::Separator();
    const float available = ImGui::GetContentRegionAvail().x;
    if (width < available) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + available - width);
}

// The common footer: [secondary] [primary], right-aligned. Returns 1 when the
// primary action is pressed, 2 for the secondary one, else 0.
inline int footer(const char* primary, const char* secondary = "Close", bool primaryEnabled = true)
{
    const float minimum = ImGui::GetFontSize() * 5.5f;
    beginFooter(buttonRowWidth({secondary, primary}, minimum));
    int result = 0;
    if (secondary && secondaryButton(secondary, ImVec2(std::max(minimum, ImGui::CalcTextSize(secondary, nullptr, true).x + 2 * ImGui::GetStyle().FramePadding.x), 0))) result = 2;
    ImGui::SameLine();
    ImGui::BeginDisabled(!primaryEnabled);
    if (primaryButton(primary, ImVec2(std::max(minimum, ImGui::CalcTextSize(primary, nullptr, true).x + 2 * ImGui::GetStyle().FramePadding.x), 0))) result = 1;
    ImGui::EndDisabled();
    return result;
}

// A placeholder centred in the remaining space of a results area.
inline void emptyState(const char* text)
{
    const ImVec2 available = ImGui::GetContentRegionAvail();
    const ImVec2 size = ImGui::CalcTextSize(text, nullptr, false, std::max(available.x * 0.8f, 1.0f));
    const ImVec2 p = ImGui::GetCursorPos();
    ImGui::SetCursorPos(ImVec2(p.x + std::max(0.0f, (available.x - size.x) * 0.5f), p.y + std::max(0.0f, (available.y - size.y) * 0.45f)));
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + available.x * 0.8f);
    ImGui::TextDisabled("%s", text);
    ImGui::PopTextWrapPos();
}

// A label in a left column with its field to the right, aligned across rows.
inline void fieldLabel(const char* label, float column = 0)
{
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    ImGui::SameLine(column > 0 ? column : ImGui::GetFontSize() * 9.0f);
}
}
