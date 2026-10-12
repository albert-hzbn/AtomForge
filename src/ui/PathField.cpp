#include "ui/PathField.h"

#include "imgui.h"
#include "imgui_internal.h"

#include <cstring>
#include <string>

bool pathField(const char* id, char* buffer, std::size_t size)
{
    ImGuiStorage* storage = ImGui::GetStateStorage();
    const ImGuiID key = ImGui::GetID(id);
    const int state = storage->GetInt(key, 0);  // 0 showing, 1 focus requested, 2 waiting for activation, 3 editing
    const float width = ImGui::CalcItemWidth();
    const ImGuiStyle& style = ImGui::GetStyle();
    const float inner = width - 2.0f * style.FramePadding.x;
    if (state == 0 && ImGui::GetActiveID() != key && ImGui::CalcTextSize(buffer).x > inner)
    {
        const ImVec2 min = ImGui::GetCursorScreenPos();
        const ImVec2 max(min.x + width, min.y + ImGui::GetFrameHeight());
        // Its own id: while the mouse is held the label, not the text field, is active.
        ImGui::PushID(id);
        const bool clicked = ImGui::InvisibleButton("label", ImVec2(width, ImGui::GetFrameHeight()));
        ImGui::PopID();
        const bool hovered = ImGui::IsItemHovered();
        ImDrawList* draw = ImGui::GetWindowDrawList();
        draw->AddRectFilled(min, max, ImGui::GetColorU32(hovered ? ImGuiCol_FrameBgHovered : ImGuiCol_FrameBg), style.FrameRounding);
        if (style.FrameBorderSize > 0.0f)
            draw->AddRect(min, max, ImGui::GetColorU32(ImGuiCol_Border), style.FrameRounding, 0, style.FrameBorderSize);
        // The longest tail that fits after an ellipsis, cut at a separator if possible.
        const char* ellipsis = "\xe2\x80\xa6";
        const float ellipsisWidth = ImGui::CalcTextSize(ellipsis).x;
        const std::size_t length = std::strlen(buffer);
        std::size_t start = length;
        while (start > 0 && ImGui::CalcTextSize(buffer + start - 1).x + ellipsisWidth <= inner) --start;
        for (std::size_t k = start; k < length; ++k)
            if (buffer[k] == '/' || buffer[k] == '\\') { start = k; break; }
        const std::string shown = std::string(ellipsis) + (buffer + start);
        draw->AddText(ImVec2(min.x + style.FramePadding.x, min.y + style.FramePadding.y), ImGui::GetColorU32(ImGuiCol_Text), shown.c_str());
        if (hovered) ImGui::SetTooltip("%s", buffer);
        if (clicked) storage->SetInt(key, 1);
        return false;
    }
    if (state == 1)
    {
        ImGui::SetKeyboardFocusHere();
        storage->SetInt(key, 2);
    }
    ImGui::SetNextItemWidth(width);
    const bool entered = ImGui::InputText(id, buffer, size, ImGuiInputTextFlags_EnterReturnsTrue);
    // Focus takes a frame or two to activate the field; editing ends when an
    // active field loses focus.
    if (state == 2 && ImGui::IsItemActive()) storage->SetInt(key, 3);
    if (state == 3 && !ImGui::IsItemActive()) storage->SetInt(key, 0);
    return entered;
}
