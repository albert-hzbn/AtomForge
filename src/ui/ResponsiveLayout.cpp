#include "ui/ResponsiveLayout.h"
#include "imgui_internal.h"

#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>

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
    // A docked window (or one being docked) is placed and sized by its dock
    // node: a position or size set for it would undock it.
    if (const ImGuiWindow* window = ImGui::FindWindowByName(name); window && (window->DockId != 0 || window->DockIsActive)) {
        next.HasFlags &= ~(ImGuiNextWindowDataFlags_HasPos | ImGuiNextWindowDataFlags_HasSize | ImGuiNextWindowDataFlags_HasSizeConstraint);
        return;
    }
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
        flags = (flags & ~(ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse)) | ImGuiWindowFlags_HorizontalScrollbar;
    }
    // Docked windows are square, flush with the view and with each other.
    const ImGuiWindow* existing = ImGui::FindWindowByName(name);
    const bool docked = existing && existing->DockIsActive;
    if (docked) ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    const bool visible = ImGui::Begin(name, open, flags);
    if (docked) ImGui::PopStyleVar();
    if (visible && open && *open && !(flags & ImGuiWindowFlags_NoTitleBar) && escapeCloses()) *open = false;
    return visible;
}

namespace
{
// Dialogs worth hiding to look at (or work in) the view while they keep their settings.
bool minimizable(const std::string& name)
{
    static const char* const names[] = {
        "Build Bulk Crystal", "Amorphous Structure Builder", "Build Nanocrystal", "Build Polycrystal",
        "CSL Grain Boundary Builder", "Interface Builder", "Stacking Faults", "Custom Structure",
        "Substitutional Solid Solution", "Insert Dislocation", "Add Interstitial Atoms", "Merge Structures",
        "Cell Sculptor", "Radial Distribution Function", "Angular Distribution Function",
        "Common Neighbour Analysis", "Short-Range Order Analysis"};
    const std::string title = name.substr(0, name.find("##"));
    for (const char* candidate : names)
        if (title == candidate) return true;
    return false;
}

struct MinimizedDialog
{
    std::string name;
    bool restore = false;
};
std::vector<MinimizedDialog>& minimizedDialogs()
{
    static std::vector<MinimizedDialog> dialogs;
    return dialogs;
}
bool lastMinimized = false;

// A "minimize" button left of the close button in the current window's title bar.
bool titleBarMinimizeButton(bool hasCloseButton)
{
    ImGuiContext& g = *ImGui::GetCurrentContext();
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    const ImRect bar = window->TitleBarRect();
    const float size = g.FontSize;
    const float right = bar.Max.x - g.Style.FramePadding.x - (hasCloseButton ? size + g.Style.ItemInnerSpacing.x : 0.0f);
    const ImRect box(ImVec2(right - size, bar.Min.y + g.Style.FramePadding.y), ImVec2(right, bar.Min.y + g.Style.FramePadding.y + size));
    const ImGuiLastItemData saved = g.LastItemData;
    ImGui::PushClipRect(bar.Min, bar.Max, false);
    const ImGuiID id = window->GetID("#MINIMIZE");
    bool pressed = false;
    if (ImGui::ItemAdd(box, id, nullptr, ImGuiItemFlags_NoNav)) {
        bool hovered = false, held = false;
        pressed = ImGui::ButtonBehavior(box, id, &hovered, &held);
        ImDrawList* draw = window->DrawList;
        if (hovered)
            draw->AddRectFilled(box.Min, box.Max, ImGui::GetColorU32(held ? ImGuiCol_ButtonActive : ImGuiCol_ButtonHovered), g.Style.FrameRounding);
        const float inset = size * 0.25f, y = box.Max.y - size * 0.32f;
        draw->AddLine(ImVec2(box.Min.x + inset, y), ImVec2(box.Max.x - inset, y), ImGui::GetColorU32(ImGuiCol_Text), std::max(1.0f, size * 0.09f));
        if (hovered) ImGui::SetTooltip("Minimize: hide the dialog and keep its settings; restore it from the bar above the status bar");
    }
    ImGui::PopClipRect();
    g.LastItemData = saved;
    return pressed;
}
}

bool lastModalMinimized() { return lastMinimized; }

bool beginModal(const char* name, bool* open, ImGuiWindowFlags flags)
{
    auto& dialogs = minimizedDialogs();
    const auto found = std::find_if(dialogs.begin(), dialogs.end(), [&](const MinimizedDialog& d) { return d.name == name; });
    // Restored: reopened here, in the dialog's own ID scope.
    if (found != dialogs.end() && found->restore) {
        ImGui::OpenPopup(name);
        dialogs.erase(found);
    }
    prepare(name);
    flags = (flags & ~ImGuiWindowFlags_NoScrollbar) | ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoCollapse;
    const bool visible = ImGui::BeginPopupModal(name, open, flags);
    lastMinimized = false;
    if (!visible) {
        lastMinimized = std::any_of(dialogs.begin(), dialogs.end(), [&](const MinimizedDialog& d) { return d.name == name; });
        return false;
    }
    // Opened again from its menu while minimized: no longer minimized.
    dialogs.erase(std::remove_if(dialogs.begin(), dialogs.end(), [&](const MinimizedDialog& d) { return d.name == name; }), dialogs.end());
    if (minimizable(name) && titleBarMinimizeButton(open != nullptr)) {
        dialogs.push_back({name, false});
        ImGui::CloseCurrentPopup();
    }
    if (open && *open && escapeCloses()) {
        *open = false;
        ImGui::CloseCurrentPopup();
    }
    return visible;
}

void drawMinimizedDialogs(float statusTop)
{
    auto& dialogs = minimizedDialogs();
    if (dialogs.empty()) return;
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + viewport->WorkSize.x * 0.5f, statusTop - dp(6)), ImGuiCond_Always, ImVec2(0.5f, 1.0f));
    ImGui::SetNextWindowBgAlpha(0.92f);
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
                                   ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoDocking;
    if (ImGui::Begin("##MinimizedDialogs", nullptr, flags)) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Minimized:");
        for (std::size_t i = 0; i < dialogs.size(); ++i) {
            ImGui::SameLine();
            ImGui::PushID(static_cast<int>(i));
            const std::string title = dialogs[i].name.substr(0, dialogs[i].name.find("##"));
            if (button(title.c_str())) dialogs[i].restore = true;
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Restore %s", title.c_str());
            ImGui::PopID();
        }
    }
    ImGui::End();
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
