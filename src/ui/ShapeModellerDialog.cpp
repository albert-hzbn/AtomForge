#include "ui/ShapeModellerDialog.h"
#include "ui/DialogLayout.h"
#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "util/ElementData.h"
#include "imgui.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>

using atomforge::pipeline::Json;
namespace shape = atomforge::shape;

namespace
{
constexpr const char* kStep = "build-shape";

int atomicNumberOf(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    return 0;
}

ImU32 elementColour(const std::string& symbol, float alpha)
{
    const int z = atomicNumberOf(symbol);
    float r = 0.62f, g = 0.66f, b = 0.72f;
    if (z > 0) getDefaultElementColor(z, r, g, b);
    return ImGui::GetColorU32(ImVec4(r, g, b, alpha));
}

double cellVolume(const Structure& s)
{
    if (!s.hasUnitCell) return 0;
    const auto& c = s.cellVectors;
    return std::abs(c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) +
                    c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]));
}
}  // namespace

void ShapeModellerDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Shape Modeller (CAD)", nullptr, m_open, enabled)) open();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
        ImGui::SetTooltip("Carve the open crystal into any shape built from primitives, booleans and transforms, as in CAD\n"
                          "Also a pipeline step: %s", kStep);
}

void ShapeModellerDialog::setScript(const std::string& script)
{
    std::snprintf(m_script.data(), m_script.size(), "%s", script.c_str());
}

bool ShapeModellerDialog::editStep(StepEdit edit)
{
    m_edit = std::move(edit);
    m_parameters = atomforge::pipeline::makeModifier(kStep).parameters;
    for (const auto& [name, value] : m_edit.parameters.members()) m_parameters[name] = value;
    if (m_parameters.contains("script") && m_parameters.at("script").isString()) setScript(m_parameters.at("script").string());
    m_open = true;
    m_message.clear();
    return true;
}

Json ShapeModellerDialog::stepParameters() const
{
    Json p = m_parameters.isObject() && !m_parameters.members().empty() ? m_parameters : atomforge::pipeline::makeModifier(kStep).parameters;
    p["script"] = std::string(m_script.data());
    return p;
}

void ShapeModellerDialog::insert(const std::string& snippet)
{
    std::string text = m_script.data();
    if (!text.empty() && text.back() != '\n') text += '\n';
    text += snippet;
    // Operations take the shapes that follow (or a { } block).
    if (snippet.find('{') == std::string::npos && snippet.find(';') == std::string::npos) text += " ";
    else text += "\n";
    setScript(text);
}

void ShapeModellerDialog::reparse(const Structure& structure)
{
    std::array<double, 9> cell{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) cell[static_cast<std::size_t>(3 * r + c)] = structure.hasUnitCell ? structure.cellVectors[r][c] : 0.0;
    if (m_parsedScript == m_script.data() && cell == m_parsedCell && m_parsedFor == &structure) return;
    m_parsedScript = m_script.data();
    m_parsedCell = cell;
    m_parsedFor = &structure;
    m_model.reset();
    m_parseError.clear();
    m_errorLine = 0;
    m_volume = 0;
    try {
        glm::dmat3 lattice(1.0);
        if (structure.hasUnitCell)
            for (int r = 0; r < 3; ++r) lattice[r] = glm::dvec3(structure.cellVectors[r][0], structure.cellVectors[r][1], structure.cellVectors[r][2]);
        m_model = shape::Model::parse(m_parsedScript, structure.hasUnitCell ? &lattice : nullptr);
        const shape::Bounds b = m_model->bounds();
        if (b.finite) {
            const glm::dvec3 extent = b.high - b.low;
            const double spacing = std::max(0.25, std::cbrt(extent.x * extent.y * extent.z) / 60);
            m_volume = m_model->volume(spacing);
        }
    } catch (const shape::ShapeError& error) {
        m_parseError = error.what();
        m_errorLine = error.line;
    } catch (const std::exception& error) {
        m_parseError = error.what();
    }
    m_rasterKey.clear();
}

void ShapeModellerDialog::drawSection(const Structure& structure, float size)
{
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 p1(p0.x + size, p0.y + size);
    draw->AddRectFilled(p0, p1, ImGui::GetColorU32(ImGuiCol_FrameBg), ImGui::GetStyle().FrameRounding);
    ImGui::InvisibleButton("##section", ImVec2(size, size));
    if (!m_model) {
        draw->AddText(ImVec2(p0.x + 8, p0.y + 8), ImGui::GetColorU32(ImGuiCol_TextDisabled), "No valid shape");
        return;
    }
    const shape::Bounds b = m_model->bounds();
    glm::dvec3 centre(0.0), extent(60.0);
    if (b.finite) {
        centre = 0.5 * (b.low + b.high);
        extent = b.high - b.low;
    }
    // In-plane axes (u, v) and the normal axis w.
    const int axes[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
    const int u = axes[m_plane][0], v = axes[m_plane][1], w = axes[m_plane][2];
    const double span = 1.1 * std::max({extent[u], extent[v], 1.0});
    const double depth = centre[w] + (m_sectionAt - 0.5) * extent[w];
    const int n = 96;
    char key[256];
    std::snprintf(key, sizeof key, "%d %.4f %d", m_plane, depth, n);
    std::string fullKey = std::string(key) + m_parsedScript;
    if (fullKey != m_rasterKey) {
        m_rasterKey = fullKey;
        m_rasterSize = n;
        m_raster.assign(static_cast<std::size_t>(n * n), 0);
        std::string body = structure.atoms.empty() ? std::string() : structure.atoms.front().symbol;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                glm::dvec3 q(0.0);
                q[u] = centre[u] + ((i + 0.5) / n - 0.5) * span;
                q[v] = centre[v] + (0.5 - (j + 0.5) / n) * span;
                q[w] = depth;
                if (!m_model->contains(q)) continue;
                const std::string* region = m_model->material(q);
                ImU32 colour = region ? (*region == "none" ? ImGui::GetColorU32(ImVec4(0.9f, 0.2f, 0.2f, 0.35f)) : elementColour(*region, 0.95f))
                                      : elementColour(body, 0.85f);
                m_raster[static_cast<std::size_t>(j * n + i)] = colour;
            }
    }
    const float cellSize = size / static_cast<float>(m_rasterSize);
    for (int j = 0; j < m_rasterSize; ++j)
        for (int i = 0; i < m_rasterSize;) {
            const ImU32 colour = m_raster[static_cast<std::size_t>(j * m_rasterSize + i)];
            int run = i + 1;
            while (run < m_rasterSize && m_raster[static_cast<std::size_t>(j * m_rasterSize + run)] == colour) ++run;
            if (colour)
                draw->AddRectFilled(ImVec2(p0.x + i * cellSize, p0.y + j * cellSize), ImVec2(p0.x + run * cellSize, p0.y + (j + 1) * cellSize), colour);
            i = run;
        }
    // The shape's origin and a scale bar.
    const auto toScreen = [&](double a, double c) {
        return ImVec2(p0.x + static_cast<float>((a - centre[u]) / span + 0.5) * size, p0.y + static_cast<float>(0.5 - (c - centre[v]) / span) * size);
    };
    const ImVec2 o = toScreen(0, 0);
    const ImU32 axisColour = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    draw->AddLine(ImVec2(o.x - 6, o.y), ImVec2(o.x + 6, o.y), axisColour);
    draw->AddLine(ImVec2(o.x, o.y - 6), ImVec2(o.x, o.y + 6), axisColour);
    const char* names = "xyz";
    char label[96];
    std::snprintf(label, sizeof label, "%c-%c section at %c = %.1f A", names[u], names[v], names[w], depth);
    draw->AddText(ImVec2(p0.x + 6, p0.y + 4), ImGui::GetColorU32(ImGuiCol_Text), label);
    double bar = std::pow(10.0, std::floor(std::log10(span / 3)));
    if (span / bar > 8) bar *= 2;
    const float barLength = static_cast<float>(bar / span) * size;
    const ImVec2 b0(p1.x - barLength - 8, p1.y - 10);
    draw->AddLine(b0, ImVec2(b0.x + barLength, b0.y), ImGui::GetColorU32(ImGuiCol_Text), 2.0f);
    std::snprintf(label, sizeof label, "%g A", bar);
    draw->AddText(ImVec2(b0.x, b0.y - ImGui::GetTextLineHeight() - 2), ImGui::GetColorU32(ImGuiCol_Text), label);
}

void ShapeModellerDialog::draw(Structure& structure, const Callbacks& callbacks)
{
    if (!m_open) return;
    if (m_parameters.members().empty()) {
        m_parameters = atomforge::pipeline::makeModifier(kStep).parameters;
        if (!m_script[0]) setScript(shape::templates()[1].second);
    }
    const bool editing = m_edit.active();
    reparse(editing ? m_edit.input : structure);
    const Structure& crystal = editing ? m_edit.input : structure;
    ImGui::SetNextWindowSize(ImVec2(responsive::dp(980), responsive::dp(700)), ImGuiCond_FirstUseEver);
    const std::string title = std::string("Shape Modeller") + (editing ? " (pipeline step)" : "") + "###Shape Modeller";
    if (responsive::begin(title.c_str(), &m_open, ImGuiWindowFlags_NoCollapse)) {
        if (editing) dialogLayout::stepBanner(kStep);
        const float footerHeight = ImGui::GetFrameHeightWithSpacing() * 2.2f;
        const float leftWidth = ImGui::GetContentRegionAvail().x * 0.56f;
        ImGui::BeginChild("##shape-script", ImVec2(leftWidth, -footerHeight), false);
        dialogLayout::section("Shape");
        // Templates and commands to insert.
        if (ImGui::Button("Templates")) ImGui::OpenPopup("##shape-templates");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Replace the script with a ready-made shape");
        if (ImGui::BeginPopup("##shape-templates")) {
            for (const auto& [name, script] : shape::templates())
                if (ImGui::MenuItem(name.c_str())) setScript(script);
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Insert")) ImGui::OpenPopup("##shape-insert");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Add a primitive, boolean, transform, modifier, pattern or region to the end of the script");
        if (ImGui::BeginPopup("##shape-insert")) {
            std::vector<std::string> groups;
            for (const auto& c : shape::commands())
                if (std::find(groups.begin(), groups.end(), c.group) == groups.end()) groups.push_back(c.group);
            for (const std::string& group : groups)
                if (ImGui::BeginMenu(group.c_str())) {
                    for (const auto& c : shape::commands()) {
                        if (group != c.group) continue;
                        if (ImGui::MenuItem(c.name)) insert(c.snippet);
                        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s\n%s", c.help, c.snippet);
                    }
                    ImGui::EndMenu();
                }
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        dialogLayout::caption("Operations act on the shape after them or on a { } block.");
        ImGui::InputTextMultiline("##shape-source", m_script.data(), m_script.size(), ImVec2(-1, responsive::dp(230)), ImGuiInputTextFlags_AllowTabInput);
        reparse(crystal);
        if (!m_parseError.empty())
            dialogLayout::note(dialogLayout::Note::Error, m_parseError.c_str());
        else if (m_model) {
            const shape::Bounds b = m_model->bounds();
            char status[512];
            if (b.finite) {
                const glm::dvec3 e = b.high - b.low;
                const double volume = cellVolume(crystal);
                if (volume > 0 && !crystal.atoms.empty())
                    std::snprintf(status, sizeof status, "%d primitive%s; %.1f x %.1f x %.1f A; volume %.0f A^3; about %.0f atoms of this crystal",
                                  m_model->primitives(), m_model->primitives() == 1 ? "" : "s", e.x, e.y, e.z, m_volume,
                                  m_volume / volume * crystal.atoms.size());
                else
                    std::snprintf(status, sizeof status, "%d primitive%s; %.1f x %.1f x %.1f A; volume %.0f A^3",
                                  m_model->primitives(), m_model->primitives() == 1 ? "" : "s", e.x, e.y, e.z, m_volume);
                dialogLayout::note(dialogLayout::Note::Success, status);
            } else
                dialogLayout::note(dialogLayout::Note::Warning,
                                   "The shape is unbounded (a half-space or slab not intersected with a finite shape): it can delete or "
                                   "select atoms, but carving needs a finite shape.");
            const auto regions = m_model->regionElements();
            if (!regions.empty()) {
                std::string list;
                for (const auto& r : regions) list += (list.empty() ? "" : ", ") + r;
                dialogLayout::caption(("Regions: " + list).c_str());
            }
        }
        dialogLayout::section("Settings");
        // The step's settings except the script (edited above).
        const auto* type = atomforge::pipeline::findModifierType(kStep);
        if (type) {
            atomforge::pipeline::ModifierType settings = *type;
            settings.parameters.erase(std::remove_if(settings.parameters.begin(), settings.parameters.end(),
                                                     [](const auto& p) {
                                                         const std::string name = p.name;
                                                         return name == "script" || name == "mesh_file";
                                                     }),
                                      settings.parameters.end());
            drawModifierParameters(settings, m_parameters, [&](const std::string& name, const Json& value) { m_parameters[name] = value; });
        }
        if (ImGui::CollapsingHeader("Command reference")) {
            if (ImGui::BeginTable("##shape-commands", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
                for (const auto& c : shape::commands()) {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    if (ImGui::Selectable(c.name, false)) insert(c.snippet);
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", c.snippet);
                    ImGui::TableNextColumn();
                    ImGui::TextWrapped("%s", c.help);
                }
                ImGui::EndTable();
            }
            ImGui::TextWrapped("Variables: R = 20; then sphere(r = R / 2). Arithmetic + - * / ^, functions sqrt, sin, cos, tan (degrees), "
                               "abs, min, max, pow. Operations act on the shape after them or on a { } block. Comments: // or #.");
        }
        ImGui::EndChild();
        ImGui::SameLine();
        ImGui::BeginChild("##shape-side", ImVec2(0, -footerHeight), false);
        dialogLayout::section("Section view");
        ImGui::RadioButton("x-y", &m_plane, 0);
        ImGui::SameLine();
        ImGui::RadioButton("x-z", &m_plane, 1);
        ImGui::SameLine();
        ImGui::RadioButton("y-z", &m_plane, 2);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(-1);
        ImGui::SliderFloat("##section-at", &m_sectionAt, 0.0f, 1.0f, "depth %.2f");
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Where the section cuts through the shape (0: one side of its bounds, 1: the other)");
        const float side = std::min(ImGui::GetContentRegionAvail().x, responsive::dp(380));
        drawSection(crystal, side);
        dialogLayout::section("Surface mesh");
        ImGui::SetNextItemWidth(-ImGui::CalcTextSize("Browse...").x - 2 * ImGui::GetStyle().FramePadding.x - ImGui::GetStyle().ItemSpacing.x);
        ImGui::InputTextWithHint("##shape-mesh", "shape.obj", m_meshPath, sizeof m_meshPath);
        ImGui::SameLine();
        if (ImGui::Button("Browse...")) m_meshPicker.open("Save the shape surface", true, m_meshPath[0] ? m_meshPath : "shape.obj");
        if (auto chosen = m_meshPicker.draw()) std::snprintf(m_meshPath, sizeof m_meshPath, "%s", chosen->c_str());
        if (ImGui::Button("Export Surface (OBJ)") && m_model) {
            try {
                const double spacing = m_parameters.contains("mesh_spacing") && m_parameters.at("mesh_spacing").isNumber()
                                           ? m_parameters.at("mesh_spacing").number() : 0.5;
                const auto mesh = m_model->mesh(spacing);
                const std::string path = m_meshPath[0] ? m_meshPath : "shape.obj";
                shape::saveObj(mesh, path);
                m_message = "Wrote " + std::to_string(mesh.vertices.size() / 3) + " triangles to " + path;
                m_error = false;
            } catch (const std::exception& error) {
                m_message = error.what();
                m_error = true;
            }
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("The shape's surface as a triangle mesh, in the shape's own coordinates (Angstrom)");
        if (!m_message.empty()) dialogLayout::note(m_error ? dialogLayout::Note::Error : dialogLayout::Note::Success, m_message.c_str());
        ImGui::EndChild();

        const Json parameters = stepParameters();
        if (editing) {
            const int pressed = dialogLayout::footer("Update Step", "Close", m_model.has_value());
            if (pressed == 1) {
                m_edit.commitParameters(parameters);
                m_message = "Step updated";
                m_error = false;
            } else if (pressed == 2)
                m_open = false;
        } else {
            const std::string mode = parameters.contains("mode") && parameters.at("mode").isString() ? parameters.at("mode").string() : "carve";
            const char* action = mode == "carve" ? "Carve" : mode == "select" ? "Select" : "Apply";
            dialogLayout::beginFooter(dialogLayout::buttonRowWidth({"Add to Pipeline", action}, ImGui::GetFontSize() * 5.5f));
            if (dialogLayout::secondaryButton("Add to Pipeline") && callbacks.addToPipeline) {
                auto step = atomforge::pipeline::makeModifier(kStep);
                step.parameters = parameters;
                callbacks.addToPipeline({step});
                m_message = "Added to the structure pipeline";
                m_error = false;
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(!m_model.has_value());
            if (dialogLayout::primaryButton(action, ImVec2(std::max(ImGui::GetFontSize() * 5.5f, ImGui::CalcTextSize(action).x + 2 * ImGui::GetStyle().FramePadding.x), 0))) {
                try {
                    atomforge::pipeline::PipelineData data;
                    data.structure = structure;
                    data.clearSelection();
                    atomforge::pipeline::findModifierType(kStep)->apply(parameters, data);
                    m_message = data.notes.empty() ? "Done" : data.notes.front();
                    m_error = false;
                    if (mode == "carve") {
                        if (callbacks.build) callbacks.build(data.structure);
                    } else if (mode == "select") {
                        std::vector<int> selected;
                        for (std::size_t i = 0; i < data.selected.size(); ++i)
                            if (data.selected[i]) selected.push_back(static_cast<int>(i));
                        if (callbacks.select) callbacks.select(selected);
                    } else {
                        structure = std::move(data.structure);
                        if (callbacks.update) callbacks.update(structure);
                    }
                } catch (const std::exception& error) {
                    m_message = error.what();
                    m_error = true;
                }
            }
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip(mode == "carve" ? "Fill the shape with this crystal (a new tab)" : "Change the structure (Edit > Undo reverts it)");
        }
    }
    ImGui::End();
    if (!m_open && m_edit.active()) m_edit.finish();
}
