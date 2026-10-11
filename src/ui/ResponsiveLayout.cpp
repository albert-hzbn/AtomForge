#include "ui/ResponsiveLayout.h"
#include "imgui_internal.h"

#include <algorithm>
#include <cfloat>

namespace responsive
{
float scale() { return ImGui::GetFontSize() / 18.0f; }
float dp(float value) { return value * scale(); }
ImVec2 size(float width, float height) { return ImVec2(dp(width), dp(height)); }

namespace
{
ImVec2 availableSize()
{
    const ImVec2 work = ImGui::GetMainViewport()->WorkSize;
    const float margin = dp(8);
    return ImVec2(std::max(1.0f, work.x - 2 * margin),
                  std::max(1.0f, work.y - 2 * margin));
}

void prepare(const char* name)
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImVec2 available = availableSize();
    ImGuiNextWindowData& next = ImGui::GetCurrentContext()->NextWindowData;
    ImVec2 minimum = size(240, 120);
    ImVec2 maximum = available;
    if (next.HasFlags & ImGuiNextWindowDataFlags_HasSizeConstraint)
    {
        minimum = next.SizeConstraintRect.Min;
        maximum = next.SizeConstraintRect.Max;
    }
    minimum.x = std::clamp(minimum.x, 0.0f, available.x);
    minimum.y = std::clamp(minimum.y, 0.0f, available.y);
    maximum.x = std::clamp(maximum.x, minimum.x, available.x);
    maximum.y = std::clamp(maximum.y, minimum.y, available.y);
    ImGui::SetNextWindowSizeConstraints(minimum, maximum);

    // Existing dialogs must also recover after the host window is made smaller.
    if (ImGuiWindow* window = ImGui::FindWindowByName(name))
    {
        ImVec2 extent(std::min(window->SizeFull.x, maximum.x),
                      std::min(window->SizeFull.y, maximum.y));
        if (window->SizeFull.x > maximum.x || window->SizeFull.y > maximum.y)
            ImGui::SetNextWindowSize(extent);
        const ImVec2 origin(viewport->WorkPos.x + dp(8), viewport->WorkPos.y + dp(8));
        const ImVec2 position(std::clamp(window->Pos.x, origin.x, origin.x + available.x - extent.x),
                              std::clamp(window->Pos.y, origin.y, origin.y + available.y - extent.y));
        if (position.x != window->Pos.x || position.y != window->Pos.y)
            ImGui::SetNextWindowPos(position);
    }
    else if (!(next.HasFlags & ImGuiNextWindowDataFlags_HasPos))
        ImGui::SetNextWindowPos(viewport->GetWorkCenter(), ImGuiCond_Appearing, ImVec2(.5f, .5f));
}
}

void windowSize(ImVec2 designSize, ImGuiCond condition)
{
    ImGui::SetNextWindowSize(size(designSize.x, designSize.y), condition);
}

void windowConstraints(ImVec2 minimum, ImVec2 maximum)
{
    const auto scaledLimit = [](float v) { return v >= FLT_MAX / 4 ? FLT_MAX : dp(v); };
    ImGui::SetNextWindowSizeConstraints(size(minimum.x, minimum.y),
        ImVec2(scaledLimit(maximum.x), scaledLimit(maximum.y)));
}

bool begin(const char* name, bool* open, ImGuiWindowFlags flags)
{
    // Toolbars and overlays deliberately occupy application-owned screen regions.
    if (!(flags & ImGuiWindowFlags_NoTitleBar))
    {
        prepare(name);
        flags = (flags & ~ImGuiWindowFlags_NoScrollbar) | ImGuiWindowFlags_HorizontalScrollbar;
    }
    const bool visible = ImGui::Begin(name, open, flags);
    if (visible && open && *open && !(flags & ImGuiWindowFlags_NoTitleBar) && escapeCloses()) *open = false;
    return visible;
}

bool beginModal(const char* name, bool* open, ImGuiWindowFlags flags)
{
    prepare(name);
    flags = (flags & ~ImGuiWindowFlags_NoScrollbar) | ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoCollapse;
    const bool visible = ImGui::BeginPopupModal(name, open, flags);
    if (visible && open && *open && escapeCloses()) {
        *open = false;
        ImGui::CloseCurrentPopup();
    }
    return visible;
}

// Esc closes the focused dialog, like its close button, unless a field is
// being edited (Esc then cancels the edit), a popup of it (a combo list) has
// the focus, or something is being dragged (Esc then cancels the drag).
bool escapeCloses()
{
    if (!ImGui::IsKeyPressed(ImGuiKey_Escape, false) || ImGui::IsAnyItemActive() || ImGui::GetDragDropPayload()) return false;
    for (int button = 0; button < ImGuiMouseButton_COUNT; ++button)
        if (ImGui::IsMouseDown(button)) return false;
    return ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows);
}

bool beginChild(const char* name, ImVec2 extent, ImGuiChildFlags childFlags, ImGuiWindowFlags flags)
{
    // Canvas children retain their intentional no-scrollbar policy.
    if (!(flags & ImGuiWindowFlags_NoScrollbar)) flags |= ImGuiWindowFlags_HorizontalScrollbar;
    return ImGui::BeginChild(name, extent, childFlags, flags);
}

bool button(const char* label, ImVec2 extent)
{
    const ImVec2 text = ImGui::CalcTextSize(label, nullptr, true);
    const ImVec2 padding = ImGui::GetStyle().FramePadding;
    if (extent.x > 0) extent.x = std::max(extent.x, text.x + 2 * padding.x);
    if (extent.y > 0) extent.y = std::max(extent.y, text.y + 2 * padding.y);
    return ImGui::Button(label, extent);
}

bool stacked(float minimumWidth) { return ImGui::GetContentRegionAvail().x < dp(minimumWidth); }
float previewWidth(float preferred, float controlsWidth)
{
    const float available = ImGui::GetContentRegionAvail().x;
    if (stacked()) return available;
    return std::max(dp(160), std::min(dp(preferred), available - dp(controlsWidth) - ImGui::GetStyle().ItemSpacing.x));
}
float panelHeight(float preferred, float footerRows)
{
    return std::max(dp(180), std::min(dp(preferred),
        ImGui::GetContentRegionAvail().y - footerRows * ImGui::GetFrameHeightWithSpacing()));
}
void nextPanel(bool stack) { if (!stack) ImGui::SameLine(); }
}
