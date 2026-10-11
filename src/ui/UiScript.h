#pragma once
// Scripted interface actions for the documentation captures
// (ATOMFORGE_TUTORIAL_CAPTURE, ui/TutorialCapture.cpp): the capture requests a
// named action ("bulk.systems" opens the crystal-system list, "csl.build"
// presses Build, ...) and the dialog that owns it performs it the next time it
// is drawn. Without a request every hook is a no-op.
#include "imgui.h"
#include "imgui_internal.h"

#include <set>
#include <string>

namespace uiScript
{
inline std::set<std::string>& pending()
{
    static std::set<std::string> actions;
    return actions;
}
inline void request(const std::string& action) { pending().insert(action); }
inline void clear() { pending().clear(); }
// True once after `action` was requested.
inline bool take(const char* action)
{
    auto& p = pending();
    if (p.empty()) return false;
    const auto it = p.find(action);
    if (it == p.end()) return false;
    p.erase(it);
    return true;
}
// Before BeginCombo(label, ...): opens that combo's list.
inline void openCombo(const char* action, const char* comboLabel)
{
    if (take(action)) ImGui::OpenPopupEx(ImHashStr("##ComboPopup", 0, ImGui::GetID(comboLabel)));
}
// Flags for BeginTabItem: selects the tab.
inline ImGuiTabItemFlags selectTab(const char* action) { return take(action) ? ImGuiTabItemFlags_SetSelected : 0; }
// Before CollapsingHeader/TreeNode: opens it.
inline void openHeader(const char* action) { if (take(action)) ImGui::SetNextItemOpen(true); }
// At the end of a scrolling region's content: scrolls to its bottom.
inline void scrollToBottom(const char* action) { if (take(action)) ImGui::SetScrollY(ImGui::GetScrollMaxY() + 1000.0f); }
}
