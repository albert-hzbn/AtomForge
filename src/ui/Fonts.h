#pragma once
// Fonts beyond the default text font, loaded by ImGuiSetup when the platform
// has them: a semibold face for headings, and icon glyphs (Segoe Fluent/MDL2
// icons on Windows) merged into the text font. Everything degrades to plain
// text when they are missing, so callers never depend on them.
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
// Icon codes are UTF-8 strings of Private Use Area code points (see Icon below).
inline const char* withIcon(const char* icon, const char* label, char* buffer, std::size_t size)
{
    if (!icons || !icon) return label;
    std::snprintf(buffer, size, "%s  %s", icon, label);
    return buffer;
}
}

// Segoe Fluent / MDL2 icon code points used by the interface (UTF-8).
namespace Icon
{
inline constexpr const char* Open = "\xee\xa3\xa5";       // U+E8E5
inline constexpr const char* Save = "\xee\x9d\x8e";       // U+E74E
inline constexpr const char* Undo = "\xee\x9e\xa7";       // U+E7A7
inline constexpr const char* Redo = "\xee\x9e\xa6";       // U+E7A6
inline constexpr const char* Home = "\xee\xa0\x8f";       // U+E80F (reset view)
inline constexpr const char* Ruler = "\xee\xb5\x9e";      // U+ED5E
inline constexpr const char* Rotate = "\xee\x9e\xad";     // U+E7AD
inline constexpr const char* Info = "\xee\xa5\x86";       // U+E946
inline constexpr const char* Warning = "\xee\x9e\xba";    // U+E7BA
inline constexpr const char* Error = "\xee\xa8\xb9";      // U+EA39
inline constexpr const char* Check = "\xee\xa4\xb0";      // U+E930 (completed)
inline constexpr const char* Add = "\xee\x9c\x90";        // U+E710
inline constexpr const char* Delete = "\xee\x9d\x8d";     // U+E74D
inline constexpr const char* Settings = "\xee\x9c\x93";   // U+E713
inline constexpr const char* Build = "\xee\x9c\x8f";      // U+E70F (edit/pencil)
inline constexpr const char* Cube = "\xef\x85\x98";       // U+F158 (3D)
inline constexpr const char* Help = "\xee\xa2\x97";       // U+E897
inline constexpr const char* Camera = "\xee\x9c\xa2";     // U+E722
inline constexpr const char* Folder = "\xee\xa2\xb7";     // U+E8B7
inline constexpr const char* Document = "\xee\xa2\xa5";   // U+E8A5
inline constexpr const char* Play = "\xee\x9d\xa8";       // U+E768
inline constexpr const char* Pipeline = "\xee\xa5\xa4";   // U+E964
}
