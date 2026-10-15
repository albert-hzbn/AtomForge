#pragma once
// Fonts beyond the default text font, loaded by ImGuiSetup: a semibold face
// for headings when the platform has one, and the icon glyphs compiled into
// the application, merged into the text font on every platform. Everything
// degrades to plain text when they are missing, so callers never depend on them.
#include "imgui.h"

#include <cstddef>
#include <cstdio>

namespace uiFonts
{
inline ImFont* heading = nullptr;   // semibold text font, or nullptr
inline bool icons = false;          // icon glyphs available in the text font

// Pushes the heading font (paired with popHeading()).
inline void pushHeading(float scale = 1.0f)
{
    ImGui::PushFont(heading, heading || scale != 1.0f ? ImGui::GetStyle().FontSizeBase * scale : 0.0f);
}
inline void popHeading() { ImGui::PopFont(); }

// "<icon> <label>" when icons are available, else the label alone.
// Icon codes are UTF-8 strings (see Icon below).
inline const char* withIcon(const char* icon, const char* label, char* buffer, std::size_t size)
{
    if (!icons || !icon) return label;
    std::snprintf(buffer, size, "%s  %s", icon, label);
    return buffer;
}
}

// Icon code points (UTF-8) of the AtomForge Icons font compiled into the
// application (src/third_party/icons), the same on every platform.
namespace Icon
{
inline constexpr const char* Open =       "\xe2\x9d\x90";  // U+2750
inline constexpr const char* Save =       "\xe2\xac\x87";  // U+2B07
inline constexpr const char* Undo =       "\xe2\x86\xb6";  // U+21B6
inline constexpr const char* Redo =       "\xe2\x86\xb7";  // U+21B7
inline constexpr const char* Home =       "\xe2\x8c\x82";  // U+2302 (reset view)
inline constexpr const char* Ruler =      "\xe2\x9f\xb7";  // U+27F7
inline constexpr const char* Rotate =     "\xe2\x9f\xb3";  // U+27F3
inline constexpr const char* Info =       "\xe2\x84\xb9";  // U+2139
inline constexpr const char* Warning =    "\xe2\x9a\xa0";  // U+26A0
inline constexpr const char* Error =      "\xe2\x9c\x96";  // U+2716
inline constexpr const char* Check =      "\xe2\x9c\x94";  // U+2714 (completed)
inline constexpr const char* Add =        "\xe2\x9c\x9a";  // U+271A
inline constexpr const char* Delete =     "\xe2\x9c\x98";  // U+2718
inline constexpr const char* Settings =   "\xe2\x9a\x99";  // U+2699
inline constexpr const char* Build =      "\xe2\x9c\x8e";  // U+270E (edit/pencil)
inline constexpr const char* Cube =       "\xe2\xac\xa1";  // U+2B21 (3D)
inline constexpr const char* Help =       "\x3f";  // U+003F (the text font's question mark)
inline constexpr const char* Camera =     "\xe2\x97\x89";  // U+25C9
inline constexpr const char* Folder =     "\xe2\x96\xa4";  // U+25A4
inline constexpr const char* Document =   "\xe2\x96\xa1";  // U+25A1
inline constexpr const char* Play =       "\xe2\x96\xb6";  // U+25B6
inline constexpr const char* Pipeline =   "\xe2\x87\x89";  // U+21C9
inline constexpr const char* Search =     "\xe2\x9a\xb2";  // U+26B2
inline constexpr const char* Pause =      "\xe2\x80\x96";  // U+2016
inline constexpr const char* Previous =   "\xe2\x87\xa4";  // U+21E4
inline constexpr const char* Next =       "\xe2\x87\xa5";  // U+21E5
}
