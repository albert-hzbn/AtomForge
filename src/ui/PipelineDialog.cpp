#include "ui/PipelineDialog.h"
#include "ui/MenuParity.h"
#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "imgui.h"

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
const ImVec4 kError(0.95f, 0.42f, 0.38f, 1.0f);

std::string titleOf(const atomforge::pipeline::Modifier& m)
{
    const ModifierType* type = atomforge::pipeline::findModifierType(m.type);
    std::string title = type ? type->title : m.type;
    if (!m.label.empty()) title = m.label + "  (" + title + ")";
    return title;
}
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
    if (ImGui::MenuItem("Structure pipeline", nullptr, m_open)) m_open = true;
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
    if (responsive::button("Add modifier")) ImGui::OpenPopup("Add modifier");
    if (ImGui::BeginPopup("Add modifier")) {
        std::string category;
        bool menuOpen = false;
        for (const auto& type : atomforge::pipeline::modifierTypes()) {
            if (category != type.category) {
                if (menuOpen) ImGui::EndMenu();
                category = type.category;
                menuOpen = ImGui::BeginMenu(type.category);
            }
            if (menuOpen && ImGui::MenuItem(type.title)) {
                // The first modifier starts the pipeline on the active structure.
                if (!m_editor.active() && !structure.atoms.empty()) m_editor.setInput(structure);
                // New modifiers go after the selected one, else at the end.
                const std::size_t at = m_selected >= 0 ? static_cast<std::size_t>(m_selected) + 1 : m_editor.size();
                m_selected = static_cast<int>(m_editor.add(type.id, at));
            }
            if (menuOpen && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", type.help);
        }
        if (menuOpen) ImGui::EndMenu();
        ImGui::EndPopup();
    }
}

void PipelineDialog::drawList()
{
    const auto& stages = m_editor.stages();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const float rail = ImGui::GetCursorScreenPos().x + responsive::dp(7);
    const float indent = responsive::dp(20);
    const ImU32 lineColour = ImGui::GetColorU32(ImGuiCol_TextDisabled);
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
    int moveFrom = -1, moveTo = -1, remove = -1, duplicate = -1;
    for (std::size_t i = 0; i < m_editor.size(); ++i) {
        const auto& m = m_editor.modifier(i);
        ImGui::PushID(static_cast<int>(i));
        node();
        ImGui::Indent(indent);
        bool enabled = m.enabled;
        if (ImGui::Checkbox("##enabled", &enabled)) m_editor.setEnabled(i, enabled);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(enabled ? "Applied; untick to bypass this step" : "Bypassed");
        ImGui::SameLine();
        const std::string title = std::to_string(i + 1) + "  " + titleOf(m);
        ImGui::BeginDisabled(!m.enabled);
        const ImGuiStyle& style = ImGui::GetStyle();
        const float buttons = 2 * ImGui::GetFrameHeight() + ImGui::CalcTextSize("X").x + 2 * style.FramePadding.x + 3 * style.ItemSpacing.x;
        if (ImGui::Selectable(title.c_str(), m_selected == static_cast<int>(i), ImGuiSelectableFlags_AllowOverlap,
                              ImVec2(std::max(responsive::dp(80), ImGui::GetContentRegionAvail().x - buttons), ImGui::GetFrameHeight())))
            m_selected = static_cast<int>(i);
        if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) requestDialog(i);
        if (ImGui::IsItemHovered() && !ImGui::IsMouseDown(ImGuiMouseButton_Left)) ImGui::SetTooltip("Double-click to edit in the %s dialog", menuPathOf(m.type).c_str());
        ImGui::EndDisabled();
        record("step " + std::to_string(i));
        // Drag a step onto another to move it there.
        if (ImGui::BeginDragDropSource()) {
            const int index = static_cast<int>(i);
            ImGui::SetDragDropPayload("PIPELINE_STEP", &index, sizeof(index));
            m_selected = index;
            ImGui::Text("Move %s", titleOf(m).c_str());
            ImGui::EndDragDropSource();
        }
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("PIPELINE_STEP")) {
                moveFrom = *static_cast<const int*>(payload->Data);
                moveTo = static_cast<int>(i);
            }
            ImGui::EndDragDropTarget();
        }
        if (ImGui::BeginPopupContextItem("step menu")) {
            if (ImGui::MenuItem("Move up", nullptr, false, i > 0)) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) - 1; }
            if (ImGui::MenuItem("Move down", nullptr, false, i + 1 < m_editor.size())) { moveFrom = static_cast<int>(i); moveTo = static_cast<int>(i) + 1; }
            if (ImGui::MenuItem("Edit in dialog")) requestDialog(i);
            if (ImGui::MenuItem("Duplicate")) duplicate = static_cast<int>(i);
            if (ImGui::MenuItem(m.enabled ? "Disable" : "Enable")) m_editor.setEnabled(i, !m.enabled);
            if (ImGui::MenuItem("Delete")) remove = static_cast<int>(i);
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
        if (ImGui::Button("X")) remove = static_cast<int>(i);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Remove this step (right-click a step to duplicate or disable it; drag it to move it)");
        record("remove " + std::to_string(i));
        // Status of the step after the last evaluation.
        if (i < stages.size()) {
            const auto& stage = stages[i];
            ImGui::Indent(responsive::dp(34));
            if (!stage.error.empty()) {
                ImGui::PushStyleColor(ImGuiCol_Text, kError);
                ImGui::TextWrapped("%s", stage.error.c_str());
                ImGui::PopStyleColor();
            } else if (stage.skipped) ImGui::TextDisabled(m.enabled ? "not run (an earlier step failed)" : "bypassed");
            else ImGui::TextDisabled("%zu atoms, %zu selected, %.1f ms", stage.atoms, stage.selected, stage.milliseconds);
            ImGui::Unindent(responsive::dp(34));
        }
        ImGui::Unindent(indent);
        ImGui::PopID();
    }
    if (m_editor.size() == 0) ImGui::TextDisabled("No modifiers yet: add one below. The output equals the input.");
    else ImGui::TextDisabled("Drag a step onto another to move it there.");
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
        const ImU32 fill = ends ? ImGui::GetColorU32(ImGuiCol_CheckMark) : off ? ImGui::GetColorU32(ImGuiCol_FrameBg) : ImGui::GetColorU32(ImGuiCol_ButtonActive);
        draw->AddCircleFilled(ImVec2(rail, nodes[k]), responsive::dp(ends ? 5.5f : 4.5f), fill);
        draw->AddCircle(ImVec2(rail, nodes[k]), responsive::dp(ends ? 5.5f : 4.5f), lineColour, 0, responsive::dp(1.2f));
    }
    // Apply list edits after drawing, so indices stay valid during the loop.
    if (moveFrom >= 0 && moveTo >= 0) moveStep(static_cast<std::size_t>(moveFrom), static_cast<std::size_t>(moveTo));
    if (duplicate >= 0) duplicateStep(static_cast<std::size_t>(duplicate));
    if (remove >= 0) removeStep(static_cast<std::size_t>(remove));
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
    ImGui::SeparatorText((std::to_string(index + 1) + "  " + type->title).c_str());
    ImGui::BeginDisabled(index == 0);
    if (responsive::button("Move up")) { moveStep(index, index - 1); ImGui::EndDisabled(); return; }
    record("move up");
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(index + 1 >= m_editor.size());
    if (responsive::button("Move down")) { moveStep(index, index + 1); ImGui::EndDisabled(); return; }
    record("move down");
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (responsive::button("Duplicate")) { duplicateStep(index); return; }
    record("duplicate");
    ImGui::SameLine();
    if (responsive::button("Remove")) { removeStep(index); return; }
    record("remove");
    ImGui::SameLine();
    bool enabled = m.enabled;
    if (ImGui::Checkbox("Enabled", &enabled)) m_editor.setEnabled(index, enabled);
    // The same dialog as the Build or Edit menu, editing this step.
    if (responsive::button("Edit in dialog")) requestDialog(index);
    record("edit in dialog");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Open %s on this step", menuPathOf(m.type).c_str());
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", type->help);
    ImGui::PopTextWrapPos();
    char label[128];
    std::snprintf(label, sizeof(label), "%s", m.label.c_str());
    ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x, responsive::dp(300)));
    if (ImGui::InputTextWithHint("Name##label", "optional display name", label, sizeof(label))) m_editor.setLabel(index, label);
    drawModifierParameters(*type, m.parameters, [&](const std::string& name, const Json& value) { m_editor.setParameter(index, name, value); });
}

void PipelineDialog::drawText()
{
    // The pipeline in the shell-pipe syntax, usable with AtomForge --pipe.
    const std::string text = m_editor.pipeline().toText();
    ImGui::SeparatorText("As a command");
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(text.empty() ? "(empty pipeline)" : text.c_str());
    ImGui::PopTextWrapPos();
    if (responsive::button("Copy as command"))
        ImGui::SetClipboardText(("AtomForge --pipe \"" + text + "\" --input INPUT --output OUTPUT").c_str());
    ImGui::SameLine();
    if (responsive::button(m_showText ? "Cancel text edit" : "Edit as text")) {
        m_showText = !m_showText;
        std::snprintf(m_textBuffer.data(), m_textBuffer.size(), "%s", text.c_str());
    }
    if (m_showText) {
        ImGui::InputTextMultiline("##pipelinetext", m_textBuffer.data(), m_textBuffer.size(), ImVec2(-1, ImGui::GetTextLineHeight() * 4));
        if (responsive::button("Apply text")) {
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
    if (responsive::begin("Structure pipeline", &m_open, ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Modifiers transform a copy of the structure in order, like a shell pipe; the input is kept, so steps can be edited, bypassed, reordered or removed at any time.");
        ImGui::PopTextWrapPos();
        if (!m_editor.active()) {
            ImGui::BeginDisabled(structure.atoms.empty());
            if (responsive::button("Start from the active structure", ImVec2(-1, responsive::dp(34)))) {
                m_editor.setInput(structure);
                m_message.clear();
            }
            ImGui::EndDisabled();
            if (structure.atoms.empty()) ImGui::TextDisabled("Load or build a structure first.");
            else ImGui::TextDisabled("Adding a modifier also starts the pipeline from the active structure.");
        } else {
            // Edits made outside the pipeline would be overwritten by its next update.
            if (!m_editor.dirty() && !m_editor.isCurrentOutput(structure) && !m_colourSelection) {
                ImGui::PushStyleColor(ImGuiCol_Text, kError);
                ImGui::TextWrapped("The structure was changed outside the pipeline.");
                ImGui::PopStyleColor();
                if (responsive::button("Use it as the pipeline input")) m_editor.setInput(structure);
                ImGui::SameLine();
                if (responsive::button("Restore the pipeline output")) show(structure, update, showAtomProperty);
            }
            if (responsive::button("Use the active structure as input")) m_editor.setInput(structure);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Restart the pipeline from the structure now in the view.");
            ImGui::SameLine();
            if (responsive::button("Bake")) {
                structure = m_editor.bake();
                m_selected = -1;
                update(structure);  // a normal edit: it enters undo history
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Make the output the structure and end the pipeline.");
            ImGui::SameLine();
            if (responsive::button("Discard")) {
                structure = m_editor.input();
                m_editor.reset();
                m_selected = -1;
                m_updating = true;
                update(structure);
                m_updating = false;
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("End the pipeline and show the input structure again.");
        }
        {
            ImGui::Separator();
            const float listHeight = std::max(responsive::dp(160), ImGui::GetContentRegionAvail().y * 0.42f);
            responsive::beginChild("Pipeline steps", ImVec2(0, listHeight), true);
            drawList();
            ImGui::EndChild();
            drawAddMenu(structure);
            ImGui::SameLine();
            ImGui::Checkbox("Update automatically", &m_autoUpdate);
            if (!m_autoUpdate) {
                ImGui::SameLine();
                ImGui::BeginDisabled(!m_editor.dirty());
                if (responsive::button("Update")) show(structure, update, showAtomProperty);
                ImGui::EndDisabled();
            }
            if (ImGui::Checkbox("Colour atoms by selection", &m_colourSelection) && m_editor.active()) show(structure, update, showAtomProperty);
            responsive::beginChild("Step editor", ImVec2(0, 0), false);
            drawEditor();
            drawText();
            ImGui::SeparatorText("Pipeline file");
            if (responsive::button("Save pipeline")) { m_pickerAction = 2; m_picker.open("Save pipeline", true, "pipeline.json"); }
            ImGui::SameLine();
            if (responsive::button("Load pipeline")) { m_pickerAction = 1; m_picker.open("Load pipeline (JSON or text)", false, "pipeline.json"); }
            ImGui::EndChild();
        }
        if (!m_message.empty()) {
            ImGui::PushStyleColor(ImGuiCol_Text, kError);
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
