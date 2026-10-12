#include "ui/OperationDialog.h"
#include "ui/DialogLayout.h"
#include "ui/MenuParity.h"
#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "imgui.h"

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <stdexcept>

using atomforge::pipeline::Json;
using atomforge::pipeline::ModifierType;

void OperationDialog::drawMenuItems(const char* menu, const char* submenu, bool hasStructure)
{
    for (const auto& location : menuLocations()) {
        if (location.dialog || std::strcmp(location.menu, menu) != 0 || std::strcmp(location.submenu, submenu) != 0) continue;
        const ModifierType* type = atomforge::pipeline::findModifierType(location.step);
        if (!type) continue;  // e.g. Build steps in a build without the command-line builders
        // Builders that make a structure need nothing loaded; the rest work on one.
        const bool needsStructure = std::string(location.step) != "build-bulk";
        if (ImGui::MenuItem(location.label, nullptr, m_open && m_step == location.step, hasStructure || !needsStructure)) open(location.step);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
            const std::string help = type->help;
            ImGui::SetTooltip("%s\nAlso a pipeline step: %s", help.substr(0, help.find('\n')).c_str(), location.step);
        }
    }
}

bool OperationDialog::open(const std::string& step)
{
    if (!atomforge::pipeline::findModifierType(step)) return false;
    if (m_edit.active()) { m_edit.finish(); m_step.clear(); }  // from the menu: no longer editing a step
    if (m_step != step) {
        m_step = step;
        m_parameters = atomforge::pipeline::makeModifier(step).parameters;
        m_message.clear();
    }
    m_open = true;
    return true;
}

bool OperationDialog::editStep(StepEdit edit)
{
    if (!open(edit.step)) return false;
    m_edit = std::move(edit);
    // The step's own parameters over the defaults (older steps may lack some).
    for (const auto& [name, value] : m_edit.parameters.members()) m_parameters[name] = value;
    m_message.clear();
    return true;
}

bool OperationDialog::actsOnAtoms(const std::string& step)
{
    if (step == "delete-selected") return true;
    const ModifierType* type = atomforge::pipeline::findModifierType(step);
    if (!type) return false;
    for (const auto& p : type->parameters)
        if (std::string(p.name) == "target") return true;
    return false;
}

bool OperationDialog::isSelectionStep(const std::string& step)
{
    const ModifierType* type = atomforge::pipeline::findModifierType(step);
    return type && std::string(type->category) == "Selection";
}

void OperationDialog::setTarget(Target target, const std::string& condition)
{
    m_target = target;
    std::snprintf(m_condition.data(), m_condition.size(), "%s", condition.c_str());
}

std::vector<int> OperationDialog::apply(Structure& structure, const std::vector<int>& viewSelection)
{
    const ModifierType* type = atomforge::pipeline::findModifierType(m_step);
    if (!type) throw std::runtime_error("Unknown operation " + m_step);
    atomforge::pipeline::PipelineData data;
    data.structure = structure;
    data.clearSelection();
    // Selection steps start from the view's selection, so add/subtract/intersect combine with it.
    if (isSelectionStep(m_step) || (actsOnAtoms(m_step) && m_target == ViewSelection))
        for (int i : viewSelection)
            if (i >= 0 && static_cast<std::size_t>(i) < data.selected.size()) data.selected[static_cast<std::size_t>(i)] = 1;
    Json parameters = m_parameters;
    if (actsOnAtoms(m_step)) {
        if (m_target == AllAtoms) std::fill(data.selected.begin(), data.selected.end(), 1);
        else if (m_target == Condition) {
            Json condition = atomforge::pipeline::makeModifier("select-expression").parameters;
            condition["expression"] = std::string(m_condition.data());
            atomforge::pipeline::findModifierType("select-expression")->apply(condition, data);
        }
        if (m_target != AllAtoms && data.selectedCount() == 0) throw std::runtime_error("No atoms match; nothing was changed");
        if (parameters.contains("target")) parameters["target"] = std::string("selected");
    }
    data.notes.clear();
    type->apply(parameters, data);
    m_message = data.notes.empty() ? std::string(type->title) + " applied" : data.notes.back();
    if (isSelectionStep(m_step)) {
        std::vector<int> selected;
        for (std::size_t i = 0; i < data.selected.size(); ++i)
            if (data.selected[i]) selected.push_back(static_cast<int>(i));
        return selected;
    }
    structure = std::move(data.structure);
    return {};
}

void OperationDialog::draw(Structure& structure, const std::vector<int>& viewSelection, const Callbacks& callbacks)
{
    if (!m_open) return;
    const ModifierType* type = atomforge::pipeline::findModifierType(m_step);
    if (!type) { m_open = false; return; }
    const std::string path = menuPathOf(m_step);
    const bool editing = m_edit.active();
    const bool selecting = isSelectionStep(m_step);
    // A fixed width and a height that follows the content.
    const float width = responsive::dp(500);
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0), ImVec2(width, FLT_MAX));
    // Titled like the menu item that opens it (the step's own title otherwise).
    const auto separator = path.rfind(" > ");
    const std::string windowName = separator == std::string::npos ? std::string(type->title) : path.substr(separator + 3);
    const std::string title = windowName + (editing ? " (pipeline step)" : "") + "###Structure operation";
    if (responsive::begin(title.c_str(), &m_open, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + width - 2 * ImGui::GetStyle().WindowPadding.x);
        ImGui::TextDisabled("%s", path.c_str());
        if (editing) dialogLayout::stepBanner(m_step.c_str());
        ImGui::TextUnformatted(modifierDescription(*type).c_str());
        ImGui::PopTextWrapPos();
        if (!type->parameters.empty()) {
            dialogLayout::section("Settings");
            drawModifierParameters(*type, m_parameters, [&](const std::string& name, const Json& value) { m_parameters[name] = value; });
        }
        if (!editing && actsOnAtoms(m_step)) {
            dialogLayout::section("Apply to");
            ImGui::RadioButton("All atoms", &m_target, AllAtoms);
            ImGui::SameLine();
            ImGui::RadioButton(("Selected in the view (" + std::to_string(viewSelection.size()) + ")").c_str(), &m_target, ViewSelection);
            ImGui::SameLine();
            ImGui::RadioButton("Atoms where", &m_target, Condition);
            if (m_target == Condition) {
                ImGui::SetNextItemWidth(-1);
                ImGui::InputTextWithHint("##condition", "e.g. fz > 0.5 && element == O", m_condition.data(), m_condition.size());
            }
        } else if (!editing && selecting) {
            ImGui::Spacing();
            const bool combines = std::any_of(type->parameters.begin(), type->parameters.end(),
                                              [](const auto& p) { return std::string(p.name) == "mode"; });
            if (combines)
                dialogLayout::notef(dialogLayout::Note::Info, "Sets the selection in the view (%zu atoms selected now), combined with it by the mode above.", viewSelection.size());
            else
                dialogLayout::notef(dialogLayout::Note::Info, "Changes the selection in the view (%zu atoms selected now).", viewSelection.size());
        }
        if (!m_message.empty()) {
            ImGui::Spacing();
            dialogLayout::note(m_error ? dialogLayout::Note::Error : dialogLayout::Note::Success, m_message.c_str());
        }
        if (editing) {
            // A pipeline step acts on the selection made by the steps above it.
            const int pressed = dialogLayout::footer("Update Step", "Close");
            m_updateButton = {ImGui::GetItemRectMin().x, ImGui::GetItemRectMin().y, ImGui::GetItemRectMax().x, ImGui::GetItemRectMax().y};
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Write these settings into the pipeline step");
            if (pressed == 1) {
                m_edit.commitParameters(m_parameters);
                m_message = "Step updated";
                m_error = false;
            } else if (pressed == 2) m_open = false;
        } else {
            const char* action = selecting ? "Select" : "Apply";
            const float minimum = ImGui::GetFontSize() * 5.5f;
            dialogLayout::beginFooter(dialogLayout::buttonRowWidth({"Add to Pipeline", action}, minimum));
            if (dialogLayout::secondaryButton("Add to Pipeline") && callbacks.addToPipeline) {
                std::vector<atomforge::pipeline::Modifier> steps;
                if (actsOnAtoms(m_step) && m_target == Condition) {
                    auto condition = atomforge::pipeline::makeModifier("select-expression");
                    condition.parameters["expression"] = std::string(m_condition.data());
                    steps.push_back(condition);
                }
                auto step = atomforge::pipeline::makeModifier(m_step);
                step.parameters = m_parameters;
                if (actsOnAtoms(m_step) && step.parameters.contains("target"))
                    step.parameters["target"] = std::string(m_target == AllAtoms ? "all" : "selected");
                if (m_step == "delete-selected" && m_target == AllAtoms) {
                    m_message = "Choose a condition: deleting every atom is not a useful pipeline step";
                    m_error = true;
                } else {
                    steps.push_back(step);
                    callbacks.addToPipeline(steps);
                    m_message = "Added to the structure pipeline" + std::string(m_target == ViewSelection && actsOnAtoms(m_step)
                        ? " (add a selection step there: the view's selection is not part of a pipeline)" : "");
                    m_error = false;
                }
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Append this step to Edit > Structure pipeline instead of changing the structure now");
            ImGui::SameLine();
            if (dialogLayout::primaryButton(action, ImVec2(std::max(minimum, ImGui::CalcTextSize(action).x + 2 * ImGui::GetStyle().FramePadding.x), 0))) {
                try {
                    const auto selected = apply(structure, viewSelection);
                    m_error = false;
                    if (selecting) { if (callbacks.select) callbacks.select(selected); }
                    else {
                        if (m_step == "compute-property" && callbacks.showAtomProperty) callbacks.showAtomProperty();
                        if (callbacks.update) callbacks.update(structure);
                    }
                } catch (const std::exception& error) {
                    m_message = error.what();
                    m_error = true;
                }
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(selecting ? "Change the selection in the view" : "Change the structure (Edit > Undo reverts it)");
        }
    }
    ImGui::End();
    if (!m_open && m_edit.active()) { m_edit.finish(); m_step.clear(); }
}
