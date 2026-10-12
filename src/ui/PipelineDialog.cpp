#include "ui/PipelineDialog.h"
#include "ui/MenuParity.h"
#include "ui/DialogLayout.h"
#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "imgui.h"
#include "imgui_internal.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>

using atomforge::pipeline::Json;
using atomforge::pipeline::ModifierType;

namespace
{

std::string titleOf(const atomforge::pipeline::Modifier& m)
{
    const ModifierType* type = atomforge::pipeline::findModifierType(m.type);
    std::string title = type ? type->title : m.type;
    if (!m.label.empty()) title = m.label + "  (" + title + ")";
    return title;
}

// A row of controls that wraps instead of running past the panel's edge:
// next(width) keeps the next control on the line when it fits.
struct Flow
{
    float right = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
    bool first = true;
    void next(float width)
    {
        if (!first && ImGui::GetItemRectMax().x + ImGui::GetStyle().ItemSpacing.x + width <= right) ImGui::SameLine();
        first = false;
    }
    void button(const char* label) { next(ImGui::CalcTextSize(label, nullptr, true).x + 2 * ImGui::GetStyle().FramePadding.x); }
    void checkbox(const char* label)
    {
        next(ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x + ImGui::CalcTextSize(label, nullptr, true).x);
    }
};
}

void PipelineDialog::requestDialog(std::size_t index)
{
    if (index >= m_editor.size() || menuPathOf(m_editor.modifier(index).type).empty()) return;
    const auto& m = m_editor.modifier(index);
    m_dialogRequest = StepEdit{};
    m_dialogRequest.step = m.type;
    m_dialogRequest.parameters = m.parameters;
    // The dialog sees the structure that enters the step.
    if (m_editor.active()) {
        try { m_editor.evaluate(); } catch (const std::exception&) {}
        m_dialogRequest.input = index == 0 || m_editor.stages().size() < index ? m_editor.input() : m_editor.stageOutput(index - 1).structure;
    }
    const std::string type = m.type;
    m_dialogRequest.commit = [this, index, type](const Json& parameters) { updateStep(index, type, parameters); };
    m_dialogRequested = true;
}

bool PipelineDialog::consumeDialogRequest(StepEdit& request)
{
    if (!m_dialogRequested) return false;
    request = std::move(m_dialogRequest);
    m_dialogRequest = StepEdit{};
    m_dialogRequested = false;
    return true;
}

void PipelineDialog::updateStep(std::size_t index, const std::string& step, const Json& parameters)
{
    // The steps may have been moved while the dialog was open.
    if (index >= m_editor.size() || m_editor.modifier(index).type != step) {
        m_message = "The step edited in the dialog is no longer at position " + std::to_string(index + 1);
        return;
    }
    if (!parameters.isObject()) return;
    for (const auto& [name, value] : parameters.members()) m_editor.setParameter(index, name, value);
    m_selected = static_cast<int>(index);
    m_message.clear();
}

void PipelineDialog::addSteps(const std::vector<atomforge::pipeline::Modifier>& steps, const Structure& structure)
{
    if (!m_editor.active() && !structure.atoms.empty()) m_editor.setInput(structure);
    atomforge::pipeline::Pipeline pipeline = m_editor.pipeline();
    for (const auto& step : steps) pipeline.modifiers.push_back(step);
    m_editor.setPipeline(pipeline);
    m_selected = static_cast<int>(m_editor.size()) - 1;
    m_open = true;
}

void PipelineDialog::drawMenuItem(bool hasStructure)
{
    if (ImGui::MenuItem("Structure Pipeline", nullptr, m_open)) m_open = true;
    if (ImGui::IsItemHovered() && !hasStructure) ImGui::SetTooltip("Load a structure first; the pipeline edits a copy of it.");
}

std::string PipelineDialog::snapshot() const
{
    Json state = m_editor.toJson();
    state["open"] = m_open;
    state["auto_update"] = m_autoUpdate;
    state["colour_selection"] = m_colourSelection;
    state["selected"] = m_selected;
    return state.dump();
}

void PipelineDialog::restore(const std::string& text, const Structure* input)
{
    m_editor.reset();
    m_selected = -1;
    if (text.empty()) return;
    try {
        const Json state = Json::parse(text);
        m_editor.fromJson(state);
        if (input && state.contains("active") && state.at("active").boolean()) m_editor.setInput(*input);
        if (const Json* v = state.find("open")) m_open = v->boolean();
        if (const Json* v = state.find("auto_update")) m_autoUpdate = v->boolean();
        if (const Json* v = state.find("colour_selection")) m_colourSelection = v->boolean();
        if (const Json* v = state.find("selected")) m_selected = std::min(static_cast<int>(v->number()), static_cast<int>(m_editor.size()) - 1);
    } catch (const std::exception& error) {
        m_message = std::string("The saved pipeline could not be restored: ") + error.what();
    }
}

void PipelineDialog::show(Structure& structure, const std::function<void(Structure&)>& update, const std::function<void()>& showAtomProperty)
{
    const auto& output = m_editor.evaluate();
    structure = output.structure;
    if (m_colourSelection && !structure.atoms.empty()) {
        structure.atomProperty.assign(output.selected.begin(), output.selected.end());
        structure.atomPropertyName = "Pipeline selection (1 = selected)";
        if (showAtomProperty) showAtomProperty();
    }
    m_updating = true;
    update(structure);
    m_updating = false;
}

void PipelineDialog::record(const std::string& name)
{
    const ImVec2 a = ImGui::GetItemRectMin(), b = ImGui::GetItemRectMax();
    m_controls.push_back({name, {a.x, a.y, b.x, b.y}});
}

void PipelineDialog::moveStep(std::size_t from, std::size_t to)
{
    if (from == to || from >= m_editor.size() || to >= m_editor.size()) return;
    m_editor.move(from, to);
    // The selection follows the moved step (and the steps it passed).
    if (m_selected == static_cast<int>(from)) m_selected = static_cast<int>(to);
    else if (from < to && m_selected > static_cast<int>(from) && m_selected <= static_cast<int>(to)) --m_selected;
    else if (from > to && m_selected >= static_cast<int>(to) && m_selected < static_cast<int>(from)) ++m_selected;
}

void PipelineDialog::removeStep(std::size_t index)
{
    if (index >= m_editor.size()) return;
    m_editor.remove(index);
    if (m_selected > static_cast<int>(index)) --m_selected;
    m_selected = std::min(m_selected, static_cast<int>(m_editor.size()) - 1);
}

void PipelineDialog::duplicateStep(std::size_t index)
{
    if (index < m_editor.size()) m_selected = static_cast<int>(m_editor.duplicate(index));
}

void PipelineDialog::drawAddMenu(const Structure& structure)
{
    if (responsive::button("Add Modifier")) ImGui::OpenPopup("Add modifier");
    if (ImGui::BeginPopup("Add modifier")) {
        const std::string type = drawModifierMenu("", "add");
        if (!type.empty()) {
            // The first modifier starts the pipeline on the active structure.
            if (!m_editor.active() && !structure.atoms.empty()) m_editor.setInput(structure);
            // New modifiers go after the selected one, else at the end.
            const std::size_t at = m_selected >= 0 ? static_cast<std::size_t>(m_selected) + 1 : m_editor.size();
            m_selected = static_cast<int>(m_editor.add(type, at));
        }
        ImGui::EndPopup();
    }
}

std::string PipelineDialog::drawModifierMenu(const std::string& exclude, const std::string& prefix)
{
    // Every modifier by category (the menus of "Add modifier", "Replace with", "Insert ...").
    std::string chosen, category;
    bool menuOpen = false;
    for (const auto& type : atomforge::pipeline::modifierTypes()) {
        if (category != type.category) {
            if (menuOpen) ImGui::EndMenu();
            category = type.category;
            menuOpen = ImGui::BeginMenu(type.category);
            record(prefix + " " + category);
        }
        if (!menuOpen) continue;
        if (ImGui::MenuItem(type.title, nullptr, exclude == type.id, exclude != type.id)) chosen = type.id;
        record(prefix + " " + type.id);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", modifierDescription(type).c_str());
    }
    if (menuOpen) ImGui::EndMenu();
    return chosen;
}

void PipelineDialog::drawList()
{
    const auto& stages = m_editor.stages();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const float rail = ImGui::GetCursorScreenPos().x + responsive::dp(7);
    const float indent = responsive::dp(20);
    const ImU32 lineColour = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    const ImU32 accent = ImGui::GetColorU32(ImGuiCol_CheckMark);
    std::vector<float> nodes;  // screen y of each node, input first
    const auto node = [&] { nodes.push_back(ImGui::GetCursorScreenPos().y + ImGui::GetFrameHeight() * 0.5f); };
    node();
    ImGui::Indent(indent);
    ImGui::AlignTextToFramePadding();
    if (m_editor.active()) {
        const Structure& input = m_editor.input();
        ImGui::Text("Input: %zu atoms%s", input.atoms.size(), input.hasUnitCell ? ", periodic cell" : "");
    } else ImGui::TextDisabled("Input: not started (the active structure is used when it starts)");
    ImGui::Unindent(indent);

    // The step being dragged (from an earlier frame), shown faded in place.
    int dragging = -1;
    if (const ImGuiPayload* payload = ImGui::GetDragDropPayload(); payload && payload->IsDataType("PIPELINE_STEP"))
        dragging = *static_cast<const int*>(payload->Data);
    // Esc cancels a drag (the step stays where it was).
    if (dragging >= 0 && ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
        ImGui::ClearDragDrop();
        ImGui::ClearActiveID();
        dragging = -1;
        m_dropGap = -1;
    }

    int moveFrom = -1, moveTo = -1, remove = -1, duplicate = -1;
    std::string replaceWith, insertType;
    std::size_t replaceAt = 0, insertAt = 0;
    std::vector<float> rowTop, rowBottom;
    float left = 0, right = 0;
    for (std::size_t i = 0; i < m_editor.size(); ++i) {
        const auto& m = m_editor.modifier(i);
        ImGui::PushID(static_cast<int>(i));
        rowTop.push_back(ImGui::GetCursorScreenPos().y);
        node();
        ImGui::Indent(indent);
        if (dragging == static_cast<int>(i)) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.35f);
        // Drag either the grip or the step's name.
        const auto dragSource = [&] {
            if (ImGui::BeginDragDropSource()) {
                const int index = static_cast<int>(i);
                ImGui::SetDragDropPayload("PIPELINE_STEP", &index, sizeof(index));
                m_selected = index;
                ImGui::Text("Move %s", titleOf(m).c_str());
                ImGui::TextDisabled("Drop between steps; Esc cancels");
                ImGui::EndDragDropSource();
            }
        };
        const ImVec2 grip = ImGui::GetCursorScreenPos();
        const float gripWidth = responsive::dp(12);
        ImGui::InvisibleButton("##grip", ImVec2(gripWidth, ImGui::GetFrameHeight()));
        const bool gripHot = ImGui::IsItemHovered() || ImGui::IsItemActive();
        if (gripHot) ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeAll);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Drag to move this step");
        dragSource();
        record("grip " + std::to_string(i));
        // Six dots: the grip.
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 2; ++c)
                draw->AddCircleFilled(ImVec2(grip.x + responsive::dp(3.5f + 5.0f * c), grip.y + ImGui::GetFrameHeight() * (0.3f + 0.2f * r)),
                                      responsive::dp(1.3f), gripHot ? accent : lineColour);
        ImGui::SameLine();
        bool enabled = m.enabled;
        if (ImGui::Checkbox("##enabled", &enabled)) m_editor.setEnabled(i, enabled);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(enabled ? "Applied; untick to bypass this step" : "Bypassed");
        ImGui::SameLine();
        const std::string title = std::to_string(i + 1) + "  " + titleOf(m);
        // The step's result after the last evaluation: on the same line when it
        // fits beside the name, else on a line of its own below.
        std::string status;
        bool failed = false;
        if (i < stages.size()) {
            const auto& stage = stages[i];
            if (!stage.error.empty()) failed = true;
            else if (stage.skipped) status = m.enabled ? "not run (an earlier step failed)" : "bypassed";
            else {
                char text[96];
                std::snprintf(text, sizeof(text), "%zu atoms, %zu selected, %.1f ms", stage.atoms, stage.selected, stage.milliseconds);
                status = text;
            }
        }
        // A bypassed step is greyed but still selectable, draggable and right-clickable.
        if (!m.enabled) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        const ImGuiStyle& style = ImGui::GetStyle();
        // Up, down and remove: three square buttons.
        const float buttons = 3 * ImGui::GetFrameHeight() + 3 * style.ItemSpacing.x;
        const float nameWidth = std::max(responsive::dp(80), ImGui::GetContentRegionAvail().x - buttons);
        // Beside the name in full, else shortened ("512 atoms, 192 sel."), else below.
        const float titleWidth = ImGui::CalcTextSize(title.c_str()).x + style.FramePadding.x * 2;
        const auto fits = [&](const std::string& text) {
            return ImGui::CalcTextSize(text.c_str()).x + 2 * style.ItemSpacing.x + titleWidth <= nameWidth;
        };
        if (!status.empty() && !fits(status) && i < stages.size() && stages[i].error.empty() && !stages[i].skipped) {
            const std::string shorter = std::to_string(stages[i].atoms) + " atoms, " + std::to_string(stages[i].selected) + " sel.";
            if (fits(shorter)) status = shorter;
        }
        const float statusWidth = status.empty() ? 0.0f : ImGui::CalcTextSize(status.c_str()).x + 2 * style.ItemSpacing.x;
        const bool statusBeside = !status.empty() && fits(status);
        // The checkbox's text baseline would push the name below the row's centre.
        ImGui::GetCurrentWindow()->DC.CurrLineTextBaseOffset = 0.0f;
        const ImVec2 nameMin = ImGui::GetCursorScreenPos();
        ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.0f, 0.5f));
        const bool clicked = ImGui::Selectable(title.c_str(), m_selected == static_cast<int>(i), ImGuiSelectableFlags_AllowOverlap,
                                               ImVec2(nameWidth, ImGui::GetFrameHeight()));
        ImGui::PopStyleVar();
        if (clicked) m_selected = static_cast<int>(i);
        if (statusBeside) {
            // Right-aligned in the name's row, clipped to it.
            const ImVec2 at(nameMin.x + nameWidth - statusWidth + style.ItemSpacing.x,
                            nameMin.y + 0.5f * (ImGui::GetFrameHeight() - ImGui::GetTextLineHeight()));
            draw->AddText(at, ImGui::GetColorU32(ImGuiCol_TextDisabled), status.c_str());
        }
        if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) requestDialog(i);
        if (ImGui::IsItemHovered() && !ImGui::IsMouseDown(ImGuiMouseButton_Left) && dragging < 0)
            ImGui::SetTooltip("Drag to move; right-click to replace, insert or remove; double-click to edit in the %s dialog", menuPathOf(m.type).c_str());
        if (!m.enabled) ImGui::PopStyleColor();
        record("step " + std::to_string(i));
        left = ImGui::GetItemRectMin().x;
        dragSource();
        if (ImGui::BeginPopupContextItem("step menu")) {
            m_selected = static_cast<int>(i);
            ImGui::TextDisabled("%s", titleOf(m).c_str());
            ImGui::Separator();
            if (ImGui::MenuItem("Edit in Dialog")) requestDialog(i);
            if (ImGui::BeginMenu("Replace With")) {
                record("replace with");
                const std::string type = drawModifierMenu(m.type, "replace");
                if (!type.empty()) { replaceWith = type; replaceAt = i; }
                ImGui::EndMenu();
            } else record("replace with");
            if (ImGui::BeginMenu("Insert Above")) {
                record("insert above");
                const std::string type = drawModifierMenu("", "above");
                if (!type.empty()) { insertType = type; insertAt = i; }
                ImGui::EndMenu();
            } else record("insert above");
            if (ImGui::BeginMenu("Insert Below")) {
                record("insert below");
                const std::string type = drawModifierMenu("", "below");
                if (!type.empty()) { insertType = type; insertAt = i + 1; }
                ImGui::EndMenu();
            } else record("insert below");
            ImGui::Separator();
            if (ImGui::MenuItem("Move Up", "Alt+Up", false, i > 0)) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) - 1; }
            if (ImGui::MenuItem("Move Down", "Alt+Down", false, i + 1 < m_editor.size())) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) + 1; }
            if (ImGui::MenuItem("Duplicate")) duplicate = static_cast<int>(i);
            if (ImGui::MenuItem(m.enabled ? "Disable" : "Enable")) m_editor.setEnabled(i, !m.enabled);
            ImGui::Separator();
            if (ImGui::MenuItem("Delete", "Del")) remove = static_cast<int>(i);
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        ImGui::BeginDisabled(i == 0);
        if (ImGui::ArrowButton("##up", ImGuiDir_Up)) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) - 1; }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Move up");
        record("up " + std::to_string(i));
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(i + 1 >= m_editor.size());
        if (ImGui::ArrowButton("##down", ImGuiDir_Down)) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) + 1; }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Move down");
        record("down " + std::to_string(i));
        ImGui::EndDisabled();
        ImGui::SameLine();
        // A square button with a drawn cross (the same in every font).
        if (ImGui::Button("##remove", ImVec2(ImGui::GetFrameHeight(), ImGui::GetFrameHeight()))) remove = static_cast<int>(i);
        {
            const ImVec2 a = ImGui::GetItemRectMin(), b = ImGui::GetItemRectMax();
            const ImVec2 c((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f);
            const float r = (b.x - a.x) * 0.2f;
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 colour = ImGui::GetColorU32(ImGuiCol_Text);
            dl->AddLine(ImVec2(c.x - r, c.y - r), ImVec2(c.x + r, c.y + r), colour, 1.5f);
            dl->AddLine(ImVec2(c.x - r, c.y + r), ImVec2(c.x + r, c.y - r), colour, 1.5f);
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Remove this step");
        record("remove " + std::to_string(i));
        right = ImGui::GetItemRectMax().x;
        // An error, or a status too long for the name's row, below it.
        if (failed || (!status.empty() && !statusBeside)) {
            const float under = gripWidth + ImGui::GetFrameHeight() + 2 * style.ItemSpacing.x;
            ImGui::Indent(under);
            ImGui::PushTextWrapPos(right - ImGui::GetWindowPos().x + ImGui::GetScrollX());
            if (failed) {
                ImGui::PushStyleColor(ImGuiCol_Text, themeStatusError());
                ImGui::TextUnformatted(stages[i].error.c_str());
                ImGui::PopStyleColor();
            } else ImGui::TextDisabled("%s", status.c_str());
            ImGui::PopTextWrapPos();
            ImGui::Unindent(under);
        }
        if (dragging == static_cast<int>(i)) ImGui::PopStyleVar();
        ImGui::Unindent(indent);
        rowBottom.push_back(ImGui::GetCursorScreenPos().y - ImGui::GetStyle().ItemSpacing.y);
        ImGui::PopID();
    }
    // The hint sits beside the pipe, not over it.
    ImGui::Indent(indent);
    ImGui::PushTextWrapPos(0.0f);
    if (m_editor.size() == 0) ImGui::TextDisabled("No modifiers yet: use Add Modifier, or right-click here. The output equals the input.");
    else ImGui::TextDisabled("Drag a step by its grip or name to reorder; right-click it for more.");
    ImGui::PopTextWrapPos();
    ImGui::Unindent(indent);
    node();
    ImGui::Indent(indent);
    ImGui::AlignTextToFramePadding();
    if (!stages.empty() && stages.size() == m_editor.size())
        ImGui::Text("Output: %zu atoms, %zu selected", stages.back().atoms, stages.back().selected);
    else ImGui::TextUnformatted("Output");
    ImGui::Unindent(indent);
    // The pipe: one line through every node, a dot per step.
    if (nodes.size() >= 2) draw->AddLine(ImVec2(rail, nodes.front()), ImVec2(rail, nodes.back()), lineColour, responsive::dp(2));
    for (std::size_t k = 0; k < nodes.size(); ++k) {
        const bool ends = k == 0 || k + 1 == nodes.size();
        const bool off = !ends && !m_editor.modifier(k - 1).enabled;
        const ImU32 fill = ends ? accent : off ? ImGui::GetColorU32(ImGuiCol_FrameBg) : ImGui::GetColorU32(ImGuiCol_ButtonActive);
        draw->AddCircleFilled(ImVec2(rail, nodes[k]), responsive::dp(ends ? 5.5f : 4.5f), fill);
        draw->AddCircle(ImVec2(rail, nodes[k]), responsive::dp(ends ? 5.5f : 4.5f), lineColour, 0, responsive::dp(1.2f));
    }

    // Dropping a dragged step: anywhere over the list, into the gap nearest the mouse,
    // which an accent line marks while dragging.
    if (!rowTop.empty()) {
        const std::size_t n = rowTop.size();
        const float mouse = ImGui::GetIO().MousePos.y;
        std::size_t gap = n;
        for (std::size_t k = 0; k < n; ++k)
            if (mouse < 0.5f * (rowTop[k] + rowBottom[k])) { gap = k; break; }
        // The list, and a margin above and below it so that the first and last gaps
        // can be reached when the list is scrolled.
        const float margin = responsive::dp(28);
        const ImVec2 position = ImGui::GetWindowPos(), size = ImGui::GetWindowSize();
        const ImRect area(ImVec2(position.x, position.y - margin), ImVec2(position.x + size.x, position.y + size.y + margin));
        // Near the top or bottom edge while dragging, the list scrolls.
        if (dragging >= 0 && area.Contains(ImGui::GetIO().MousePos)) {
            float speed = 0;
            if (mouse < position.y + margin) speed = mouse - (position.y + margin);
            else if (mouse > position.y + size.y - margin) speed = mouse - (position.y + size.y - margin);
            if (speed != 0) ImGui::SetScrollY(std::clamp(ImGui::GetScrollY() + 0.5f * speed, 0.0f, ImGui::GetScrollMaxY()));
        }
        // The target's hit test is clipped to the list; widen the clip to its margin.
        ImGui::PushClipRect(area.Min, area.Max, false);
        const bool target = ImGui::BeginDragDropTargetCustom(area, ImGui::GetID("pipeline list"));
        ImGui::PopClipRect();
        if (target) {
            const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("PIPELINE_STEP", ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
            if (payload) {
                const int from = *static_cast<const int*>(payload->Data);
                const std::size_t to = gap > static_cast<std::size_t>(from) ? gap - 1 : gap;
                // The insertion line (not where the step already is).
                if (to != static_cast<std::size_t>(from)) {
                    const float y = gap == 0 ? rowTop[0] - responsive::dp(2) : gap == n ? rowBottom[n - 1] + responsive::dp(2) : 0.5f * (rowBottom[gap - 1] + rowTop[gap]);
                    draw->AddLine(ImVec2(rail, y), ImVec2(right, y), accent, responsive::dp(3));
                    draw->AddTriangleFilled(ImVec2(left - responsive::dp(8), y - responsive::dp(5)), ImVec2(left - responsive::dp(8), y + responsive::dp(5)),
                                            ImVec2(left - responsive::dp(1), y), accent);
                    m_dropGap = static_cast<int>(gap);
                } else m_dropGap = -1;
                if (payload->IsDelivery() && to != static_cast<std::size_t>(from)) { moveFrom = from; moveTo = static_cast<int>(to); }
            }
            ImGui::EndDragDropTarget();
        }
        if (dragging < 0) m_dropGap = -1;
    }
    // Right-click on empty space: add a modifier there.
    if (ImGui::BeginPopupContextWindow("list menu", ImGuiPopupFlags_MouseButtonRight | ImGuiPopupFlags_NoOpenOverItems)) {
        if (ImGui::BeginMenu("Add modifier")) {
            record("add modifier");
            const std::string type = drawModifierMenu("", "add");
            if (!type.empty()) { insertType = type; insertAt = m_editor.size(); }
            ImGui::EndMenu();
        } else record("add modifier");
        ImGui::EndPopup();
    }
    // Keys for the selected step while the list has focus.
    if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && !ImGui::GetIO().WantTextInput &&
        m_selected >= 0 && static_cast<std::size_t>(m_selected) < m_editor.size()) {
        const auto& io = ImGui::GetIO();
        if (ImGui::IsKeyPressed(ImGuiKey_Delete, false)) remove = m_selected;
        else if (io.KeyAlt && ImGui::IsKeyPressed(ImGuiKey_UpArrow) && m_selected > 0) { moveFrom = m_selected; moveTo = m_selected - 1; }
        else if (io.KeyAlt && ImGui::IsKeyPressed(ImGuiKey_DownArrow) && static_cast<std::size_t>(m_selected) + 1 < m_editor.size()) { moveFrom = m_selected; moveTo = m_selected + 1; }
    }

    // Apply list edits after drawing, so indices stay valid during the loop.
    if (moveFrom >= 0 && moveTo >= 0) moveStep(static_cast<std::size_t>(moveFrom), static_cast<std::size_t>(moveTo));
    if (duplicate >= 0) duplicateStep(static_cast<std::size_t>(duplicate));
    if (remove >= 0) removeStep(static_cast<std::size_t>(remove));
    if (!replaceWith.empty()) {
        m_editor.replace(replaceAt, replaceWith);
        m_selected = static_cast<int>(replaceAt);
    }
    if (!insertType.empty()) m_selected = static_cast<int>(m_editor.add(insertType, insertAt));
}

void PipelineDialog::drawEditor()
{
    if (m_selected < 0 || static_cast<std::size_t>(m_selected) >= m_editor.size()) {
        ImGui::TextDisabled("Select a step to edit its parameters.");
        return;
    }
    const std::size_t index = static_cast<std::size_t>(m_selected);
    const auto& m = m_editor.modifier(index);
    const ModifierType* type = atomforge::pipeline::findModifierType(m.type);
    ImGui::SeparatorText(("Step " + std::to_string(index + 1) + ": " + type->title).c_str());
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", modifierDescription(*type).c_str());
    ImGui::PopTextWrapPos();
    // Parameters first: they are what a step is selected for.
    char label[128];
    std::snprintf(label, sizeof(label), "%s", m.label.c_str());
    if (ImGui::BeginTable("##name", 2, ImGuiTableFlags_SizingStretchProp)) {
        // The same columns as the parameter form below.
        ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.42f);
        ImGui::TableSetupColumn("field", ImGuiTableColumnFlags_WidthStretch, 0.58f);
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Display name");
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::InputTextWithHint("##label", "optional", label, sizeof(label))) m_editor.setLabel(index, label);
        ImGui::EndTable();
    }
    drawModifierParameters(*type, m.parameters, [&](const std::string& name, const Json& value) { m_editor.setParameter(index, name, value); });
    ImGui::Spacing();
    Flow flow;
    flow.button("Edit in Dialog");
    // The same dialog as the Build or Edit menu, editing this step.
    if (responsive::button("Edit in Dialog")) requestDialog(index);
    record("edit in dialog");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Open %s on this step", menuPathOf(m.type).c_str());
    flow.checkbox("Enabled");
    bool enabled = m.enabled;
    if (ImGui::Checkbox("Enabled", &enabled)) m_editor.setEnabled(index, enabled);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Untick to bypass this step");
    flow.button("Move Up");
    ImGui::BeginDisabled(index == 0);
    if (responsive::button("Move Up")) { moveStep(index, index - 1); ImGui::EndDisabled(); return; }
    record("move up");
    ImGui::EndDisabled();
    flow.button("Move Down");
    ImGui::BeginDisabled(index + 1 >= m_editor.size());
    if (responsive::button("Move Down")) { moveStep(index, index + 1); ImGui::EndDisabled(); return; }
    record("move down");
    ImGui::EndDisabled();
    flow.button("Duplicate");
    if (responsive::button("Duplicate")) { duplicateStep(index); return; }
    record("duplicate");
    flow.button("Remove");
    if (dialogLayout::dangerButton("Remove")) { removeStep(index); return; }
    record("remove");
}

void PipelineDialog::drawText()
{
    // The pipeline in the shell-pipe syntax, usable with AtomForge --pipe.
    const std::string text = m_editor.pipeline().toText();
    ImGui::SeparatorText("As a command");
    ImGui::PushTextWrapPos(0.0f);
    if (text.empty()) ImGui::TextDisabled("(empty pipeline)");
    else ImGui::TextUnformatted(text.c_str());
    ImGui::PopTextWrapPos();
    Flow flow;
    flow.button("Copy as Command");
    if (responsive::button("Copy as Command"))
        ImGui::SetClipboardText(("AtomForge --pipe \"" + text + "\" --input INPUT --output OUTPUT").c_str());
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Copy an AtomForge --pipe command line that runs these steps");
    flow.button(m_showText ? "Cancel Text Edit" : "Edit as Text");
    if (responsive::button(m_showText ? "Cancel Text Edit" : "Edit as Text")) {
        m_showText = !m_showText;
        std::snprintf(m_textBuffer.data(), m_textBuffer.size(), "%s", text.c_str());
    }
    if (m_showText) {
        ImGui::InputTextMultiline("##pipelinetext", m_textBuffer.data(), m_textBuffer.size(), ImVec2(-1, ImGui::GetTextLineHeight() * 4));
        if (responsive::button("Apply Text")) {
            try {
                m_editor.setPipeline(atomforge::pipeline::Pipeline::parse(m_textBuffer.data()));
                m_selected = m_editor.size() ? 0 : -1;
                m_showText = false;
                m_message.clear();
            } catch (const std::exception& error) { m_message = error.what(); }
        }
    }
}

void PipelineDialog::draw(Structure& structure, const std::function<void(Structure&)>& update, const std::function<void()>& showAtomProperty)
{
    // The output follows every edit, also while the panel is closed (e.g. after loading a project).
    if (m_editor.active() && m_editor.dirty() && m_autoUpdate) {
        try { show(structure, update, showAtomProperty); }
        catch (const std::exception& error) { m_message = error.what(); }
    }
    if (!m_open) return;
    m_controls.clear();
    responsive::windowSize(ImVec2(560, 760), ImGuiCond_FirstUseEver);
    responsive::windowConstraints(ImVec2(420, 480), ImVec2(FLT_MAX, FLT_MAX));
    if (responsive::begin("Structure Pipeline###Structure pipeline", &m_open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Modifiers change a copy of the structure in order; the input is kept, so steps can be edited, bypassed, reordered or removed at any time.");
        ImGui::PopTextWrapPos();
        if (!m_editor.active()) {
            ImGui::BeginDisabled(structure.atoms.empty());
            if (responsive::button("Start from the Active Structure", ImVec2(-1, responsive::dp(34)))) {
                m_editor.setInput(structure);
                m_message.clear();
            }
            ImGui::EndDisabled();
            if (structure.atoms.empty()) ImGui::TextDisabled("Load or build a structure first.");
            else ImGui::TextDisabled("Adding a modifier also starts the pipeline from the active structure.");
        } else {
            // Edits made outside the pipeline would be overwritten by its next update.
            const bool changedOutside = !m_editor.dirty() && !m_editor.isCurrentOutput(structure) && !m_colourSelection;
            if (changedOutside) {
                ImGui::PushStyleColor(ImGuiCol_Text, themeStatusError());
                ImGui::TextWrapped("The structure was changed outside the pipeline; its next update would replace those changes.");
                ImGui::PopStyleColor();
                Flow choices;
                choices.button("Use It as the Pipeline Input");
                if (responsive::button("Use It as the Pipeline Input")) m_editor.setInput(structure);
                choices.button("Restore the Pipeline Output");
                if (responsive::button("Restore the Pipeline Output")) show(structure, update, showAtomProperty);
            }
            Flow flow;
            flow.button("Bake");
            if (dialogLayout::primaryButton("Bake")) {
                structure = m_editor.bake();
                m_selected = -1;
                update(structure);  // a normal edit: it enters undo history
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Make the output the structure and end the pipeline.");
            flow.button("Discard");
            if (dialogLayout::dangerButton("Discard")) {
                structure = m_editor.input();
                m_editor.reset();
                m_selected = -1;
                m_updating = true;
                update(structure);
                m_updating = false;
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("End the pipeline and show the input structure again.");
            flow.button("Use the Active Structure as Input");
            if (responsive::button("Use the Active Structure as Input")) m_editor.setInput(structure);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Restart the pipeline from the structure now in the view.");
        }
        {
            ImGui::Separator();
            // Exactly as tall as the input, the steps and the hint, between
            // three rows and 45% of the panel; it scrolls beyond that.
            const float row = ImGui::GetFrameHeightWithSpacing();
            const float maximum = std::max(4.0f * row, ImGui::GetContentRegionAvail().y * 0.45f);
            ImGui::SetNextWindowSizeConstraints(ImVec2(0, std::min(4.0f * row, maximum)), ImVec2(FLT_MAX, maximum));
            responsive::beginChild("Pipeline steps", ImVec2(0, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_AutoResizeY);
            drawList();
            ImGui::EndChild();
            Flow flow;
            flow.button("Add Modifier");
            drawAddMenu(structure);
            flow.checkbox("Update automatically");
            ImGui::Checkbox("Update automatically", &m_autoUpdate);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Recompute the output after every edit");
            if (!m_autoUpdate) {
                flow.button("Update");
                ImGui::BeginDisabled(!m_editor.dirty());
                if (responsive::button("Update")) show(structure, update, showAtomProperty);
                ImGui::EndDisabled();
            }
            flow.checkbox("Colour atoms by selection");
            if (ImGui::Checkbox("Colour atoms by selection", &m_colourSelection) && m_editor.active()) show(structure, update, showAtomProperty);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Colour the output's selected atoms (1) apart from the others (0)");
            // The selected step's form, the command and the file: the rest of the panel.
            const float messageHeight = m_message.empty() ? 0.0f : 2.0f * ImGui::GetTextLineHeightWithSpacing();
            responsive::beginChild("Step editor", ImVec2(0, -messageHeight), ImGuiChildFlags_Borders);
            drawEditor();
            drawText();
            ImGui::SeparatorText("Pipeline file");
            Flow file;
            file.button("Save Pipeline");
            if (responsive::button("Save Pipeline")) { m_pickerAction = 2; m_picker.open("Save pipeline", true, "pipeline.json"); }
            file.button("Load Pipeline");
            if (responsive::button("Load Pipeline")) { m_pickerAction = 1; m_picker.open("Load pipeline (JSON or text)", false, "pipeline.json"); }
            ImGui::EndChild();
        }
        if (!m_message.empty()) {
            ImGui::PushStyleColor(ImGuiCol_Text, themeStatusError());
            ImGui::TextWrapped("%s", m_message.c_str());
            ImGui::PopStyleColor();
        }
        if (auto path = m_picker.draw()) {
            try {
                const auto file = std::filesystem::u8path(*path);
                if (m_pickerAction == 2) {
                    std::ofstream out(file);
                    out << m_editor.pipeline().toJson().dump(2) << "\n";
                    if (!out) throw std::runtime_error("Cannot write " + *path);
                    m_message.clear();
                } else {
                    std::ifstream in(file);
                    if (!in) throw std::runtime_error("Cannot read " + *path);
                    std::stringstream content;
                    content << in.rdbuf();
                    const std::string body = content.str();
                    const auto start = body.find_first_not_of(" \t\r\n");
                    m_editor.setPipeline(start != std::string::npos && body[start] == '{'
                        ? atomforge::pipeline::Pipeline::fromJson(Json::parse(body)) : atomforge::pipeline::Pipeline::parse(body));
                    m_selected = m_editor.size() ? 0 : -1;
                    m_message.clear();
                }
            } catch (const std::exception& error) { m_message = error.what(); }
        }
    }
    ImGui::End();
}
