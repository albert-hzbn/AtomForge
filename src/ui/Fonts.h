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

// Icon strings (namespace Icon): the Codicons-based font compiled into the
// application (src/third_party/icons), the same on every platform.
#include "third_party/icons/AtomForgeIconCodes.h"
