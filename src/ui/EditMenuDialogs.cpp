#include "ui/ResponsiveLayout.h"
#include "ui/UiScript.h"
#include "ui/DialogLayout.h"
#include "EditMenuDialogs.h"

#include "ElementData.h"
#include "PeriodicTableDialog.h"
#include "io/StructureLoader.h"
#include "math/StructureMath.h"
#include "science/ScienceData.h"
#include "ui/ThemeUtils.h"

#include "imgui.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>

namespace
{
constexpr int kMinElementZ = 1;
constexpr int kMaxElementZ = 118;

bool isValidElementNumber(int z)
{
    return z >= kMinElementZ && z <= kMaxElementZ;
}

std::vector<float> makeDefaultElementShininess()
{
    // Moderate baseline specular exponent for all elements.
    return std::vector<float>(119, 32.0f);
}
}

EditMenuDialogs::EditMenuDialogs()
    : elementRadii(makeLiteratureCovalentRadii())
    , elementColors(makeDefaultElementColors())
    , elementShininess(makeDefaultElementShininess())
{}

void EditMenuDialogs::drawMenuItems()
{
    if (ImGui::MenuItem("Edit Structure"))
    {
        m_step.finish();
        m_openEditStructure = true;
    }
}

bool EditMenuDialogs::editStep(StepEdit edit)
{
    if (edit.step != "set-cell" && edit.step != "add-atom")
        return false;
    m_step = std::move(edit);
    using atomforge::pipeline::Json;
    using atomforge::pipeline::parameter;
    const atomforge::pipeline::ModifierType* type = atomforge::pipeline::findModifierType(m_step.step);
    const auto vector = [](const Json& v, double (&out)[3]) {
        for (int k = 0; k < 3 && k < (int)v.items().size(); ++k)
            out[k] = v.items()[(size_t)k].number();
    };
    try
    {
        if (m_step.step == "set-cell")
        {
            // Default vectors stand for "not set yet": show the cell entering the step.
            bool defaults = true;
            const char* names[3] = {"a", "b", "c"};
            for (int r = 0; r < 3; ++r)
            {
                const Json value = parameter(m_step.parameters, *type, names[r]);
                defaults = defaults && value.dump() == parameter(Json::object(), *type, names[r]).dump();
                vector(value, m_stepCell[r]);
            }
            if (defaults && m_step.input.hasUnitCell)
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        m_stepCell[r][c] = m_step.input.cellVectors[r][c];
            m_stepScaleAtoms = parameter(m_step.parameters, *type, "scale_atoms").boolean();
        }
        else
        {
            const int z = atomforge::science::atomicNumber(parameter(m_step.parameters, *type, "element").string());
            if (isValidElementNumber(z))
                m_stepElement = z;
            vector(parameter(m_step.parameters, *type, "position"), m_stepPosition);
            m_stepFractional = parameter(m_step.parameters, *type, "fractional").boolean();
            m_stepSelect = parameter(m_step.parameters, *type, "select").boolean();
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "[Edit Structure] Step parameters not read: " << e.what() << std::endl;
    }
    m_openEditStructure = true;
    return true;
}

atomforge::pipeline::Json EditMenuDialogs::stepParameters() const
{
    using atomforge::pipeline::Json;
    const auto vector = [](const double (&v)[3]) { return Json::array({v[0], v[1], v[2]}); };
    Json p = m_step.parameters;
    if (m_step.step == "set-cell")
    {
        p["a"] = vector(m_stepCell[0]);
        p["b"] = vector(m_stepCell[1]);
        p["c"] = vector(m_stepCell[2]);
        p["scale_atoms"] = m_stepScaleAtoms;
    }
    else
    {
        p["element"] = std::string(elementSymbol(m_stepElement));
        p["position"] = vector(m_stepPosition);
        p["fractional"] = m_stepFractional;
        p["select"] = m_stepSelect;
    }
    return p;
}

// The Edit Structure dialog on a pipeline step: only the part the step can express.
void EditMenuDialogs::drawStepEditor()
{
    dialogLayout::stepBanner(m_step.step.c_str());
    ImGui::Separator();

    const Structure& input = m_step.input;
    if (m_step.step == "set-cell")
    {
        ImGui::Text("Lattice Vectors (Cartesian)");
        ImGui::DragScalarN("a", ImGuiDataType_Double, m_stepCell[0], 3, 0.01f, nullptr, nullptr, "%.6f");
        ImGui::DragScalarN("b", ImGuiDataType_Double, m_stepCell[1], 3, 0.01f, nullptr, nullptr, "%.6f");
        ImGui::DragScalarN("c", ImGuiDataType_Double, m_stepCell[2], 3, 0.01f, nullptr, nullptr, "%.6f");
        if (input.hasUnitCell && responsive::button("Use Cell Entering The Step"))
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    m_stepCell[r][c] = input.cellVectors[r][c];
        ImGui::Checkbox("Scale atoms with the cell (keep fractional coordinates)", &m_stepScaleAtoms);
        if (m_stepScaleAtoms && !input.hasUnitCell)
            ImGui::TextColored(themeStatusWarn(), "Scaling atoms needs a structure with a unit cell.");
    }
    else
    {
        ImGui::Text("Add Atom");
        drawPeriodicTableInlineSelector(m_stepElement);
        if (!isValidElementNumber(m_stepElement))
            m_stepElement = 29;
        ImGui::Text("Element: %s", elementSymbol(m_stepElement));

        const char* modeLabels[] = { "Cartesian", "Direct" };
        int mode = m_stepFractional ? 1 : 0;
        if (ImGui::Combo("Position Mode", &mode, modeLabels, 2) && (mode == 1) != m_stepFractional)
        {
            // Keep the point where it is when the cell entering the step allows converting.
            const glm::vec3 from((float)m_stepPosition[0], (float)m_stepPosition[1], (float)m_stepPosition[2]);
            glm::vec3 to(0.0f);
            const bool converted = mode == 1 ? tryCartesianToFractional(input, from, to)
                                             : tryFractionalToCartesian(input, from, to);
            if (converted)
                for (int k = 0; k < 3; ++k) m_stepPosition[k] = to[k];
            m_stepFractional = mode == 1;
        }
        ImGui::DragScalarN(m_stepFractional ? "Position (fractional)" : "Position (\xc3\x85)", ImGuiDataType_Double,
                           m_stepPosition, 3, 0.005f, nullptr, nullptr, "%.6f");
        if (m_stepFractional && !input.hasUnitCell)
            ImGui::TextColored(themeStatusWarn(), "A fractional position needs a structure with a unit cell.");
        ImGui::Checkbox("Select the new atom", &m_stepSelect);
    }

    ImGui::Separator();
    if (responsive::button("Update Step"))
    {
        m_step.commitParameters(stepParameters());
        ImGui::CloseCurrentPopup();
        m_step.finish();
    }
    ImGui::SameLine();
    if (responsive::button("Close"))
    {
        ImGui::CloseCurrentPopup();
        m_step.finish();
    }
}

void EditMenuDialogs::drawSettingsMenuItems()
{
    if (ImGui::MenuItem("Atomic Sizes"))
        m_openAtomicSize = true;

    if (ImGui::MenuItem("Display Settings"))
        m_openElementColor = true;
}

void EditMenuDialogs::drawPopups(Structure& structure,
                                 const std::function<void(Structure&)>& updateBuffers)
{
    if (m_openAtomicSize)   { ImGui::OpenPopup("Atomic Sizes##edit");   m_openAtomicSize   = false; }
    if (m_openElementColor) { ImGui::OpenPopup("Display Settings##edit"); m_openElementColor = false; }
    if (m_openEditStructure) { ImGui::OpenPopup("Edit Structure##edit"); m_openEditStructure = false; }

    // ------------------------------------------------------------------
    // Atomic Sizes modal
    // ------------------------------------------------------------------

    bool atomicSizesOpen = true;
    if (responsive::beginModal("Atomic Sizes##edit", &atomicSizesOpen,
                               ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::TextDisabled("Pick an element and set the radius its atoms are drawn with. Defaults: covalent radii of Cordero et al., Dalton Trans. 2008.");
        ImGui::Spacing();
        drawPeriodicTableInlineSelector(m_selectedRadiusElement);

        if (isValidElementNumber(m_selectedRadiusElement))
        {
            dialogLayout::section((std::string(elementSymbol(m_selectedRadiusElement)) + " (Z = " + std::to_string(m_selectedRadiusElement) + ")").c_str());

            float radius = elementRadii[m_selectedRadiusElement];
            ImGui::SetNextItemWidth(responsive::dp(240.0f));
            if (ImGui::DragFloat("Radius", &radius,
                                 0.01f, 0.20f, 3.50f, "%.2f \xc3\x85"))
            {
                elementRadii[m_selectedRadiusElement] = radius;
                updateBuffers(structure);
            }

            ImGui::SameLine();
            if (dialogLayout::secondaryButton("Reset Element"))
            {
                const std::vector<float> defaults = makeLiteratureCovalentRadii();
                elementRadii[m_selectedRadiusElement] =
                    defaults[m_selectedRadiusElement];
                updateBuffers(structure);
            }
        }
        else
        {
            ImGui::TextDisabled("No element selected.");
        }

        const int sizesPressed = dialogLayout::footer("Close", "Reset All Elements");
        if (sizesPressed == 2)
        {
            elementRadii = makeLiteratureCovalentRadii();
            updateBuffers(structure);
        }
        else if (sizesPressed == 1)
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
    if (!atomicSizesOpen)
        ImGui::CloseCurrentPopup();

    // ------------------------------------------------------------------
    // Display Settings modal
    // ------------------------------------------------------------------

    bool elementColorsOpen = true;
    if (responsive::beginModal("Display Settings##edit", &elementColorsOpen,
                               ImGuiWindowFlags_AlwaysAutoResize))
    {
        constexpr float kSliderW = 340.0f;

        if (ImGui::BeginTabBar("##displayTabs"))
        {
            // ==============================================================
            // TAB: Lighting
            // ==============================================================
            if (ImGui::BeginTabItem("Lighting", nullptr, uiScript::selectTab("display.lighting")))
            {
                ImGui::Spacing();
                ImGui::TextDisabled("Controls the scene illumination model.");
                ImGui::Spacing();

                auto lightRow = [&](const char* label, float* v,
                                    float mn, float mx, const char* fmt,
                                    const char* tip)
                {
                    ImGui::SetNextItemWidth(kSliderW);
                    ImGui::SliderFloat(label, v, mn, mx, fmt);
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tip);
                };

                lightRow("Ambient##L",        &lightAmbient,        0.00f, 0.60f, "%.2f",
                         "Shadow-side base brightness (default 0.18).");
                lightRow("Saturation##L",     &lightSaturation,     0.50f, 3.00f, "%.2f",
                         "Colour vividness — 1.0 = neutral, >1 = more vivid (default 1.55).");
                lightRow("Contrast##L",       &lightContrast,       0.50f, 2.50f, "%.2f",
                         "Dark-to-light stretch — 1.0 = none, >1 = more contrast (default 1.25).");
                lightRow("Shadow Strength##L",&lightShadowStrength, 0.00f, 1.00f, "%.2f",
                         "Shadow darkness — 0 = no shadow, 1 = fully black (default 0.75).");

                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
                if (responsive::button("Reset Lighting Defaults"))
                {
                    lightAmbient        = 0.18f;
                    lightSaturation     = 1.55f;
                    lightContrast       = 1.25f;
                    lightShadowStrength = 0.75f;
                }
                ImGui::Spacing();
                ImGui::EndTabItem();
            }

            // ==============================================================
            // TAB: Material
            // ==============================================================
            if (ImGui::BeginTabItem("Material", nullptr, uiScript::selectTab("display.material")))
            {
                ImGui::Spacing();
                ImGui::TextDisabled("Global material properties applied to all atoms and bonds.");
                ImGui::Spacing();

                auto matRow = [&](const char* label, float* v,
                                  float mn, float mx, const char* fmt,
                                  const char* tip)
                {
                    ImGui::SetNextItemWidth(kSliderW);
                    ImGui::SliderFloat(label, v, mn, mx, fmt);
                    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", tip);
                };

                matRow("Specular Intensity##M", &materialSpecularIntensity, 0.00f, 2.00f, "%.2f",
                       "Highlight brightness on all surfaces (default 0.65).");
                matRow("Shininess Scale##M",    &materialShininessScale,    0.10f, 5.00f, "%.2f",
                       "Multiplier on per-element shininess — higher = tighter highlights (default 1.50).");
                matRow("Shininess Floor##M",    &materialShininessFloor,    4.00f, 256.0f, "%.1f",
                       "Minimum specular exponent — prevents very diffuse highlights (default 32).");

                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
                if (responsive::button("Reset Material Defaults"))
                {
                    materialSpecularIntensity = 0.65f;
                    materialShininessScale    = 1.5f;
                    materialShininessFloor    = 32.0f;
                }
                ImGui::Spacing();
                ImGui::EndTabItem();
            }

            // ==============================================================
            // TAB: Colors
            // ==============================================================
            if (ImGui::BeginTabItem("Colours", nullptr, uiScript::selectTab("display.colors")))
            {
                ImGui::Spacing();
                ImGui::TextDisabled("The colour each element is drawn with.");
                ImGui::Spacing();

                drawPeriodicTableInlineSelector(m_selectedColorElement);

                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();

                if (isValidElementNumber(m_selectedColorElement))
                {
                    ImGui::Text("%s (Z = %d)", elementSymbol(m_selectedColorElement), m_selectedColorElement);

                    float color[3] = {
                        elementColors[m_selectedColorElement].r,
                        elementColors[m_selectedColorElement].g,
                        elementColors[m_selectedColorElement].b
                    };
                    ImGui::SetNextItemWidth(kSliderW);
                    if (ImGui::ColorEdit3("Colour##el", color))
                    {
                        elementColors[m_selectedColorElement] =
                            glm::vec3(color[0], color[1], color[2]);
                        updateBuffers(structure);
                    }

                    ImGui::SameLine();
                    if (dialogLayout::secondaryButton("Reset Element##sel"))
                    {
                        float r, g, b;
                        getDefaultElementColor(m_selectedColorElement, r, g, b);
                        elementColors[m_selectedColorElement] = glm::vec3(r, g, b);
                        updateBuffers(structure);
                    }
                }
                else
                {
                    ImGui::TextDisabled("No element selected.");
                }

                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
                if (dialogLayout::secondaryButton("Reset All Colours"))
                {
                    elementColors = makeDefaultElementColors();
                    updateBuffers(structure);
                }
                ImGui::Spacing();
                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }

        if (dialogLayout::footer("Close", nullptr) == 1)
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (!elementColorsOpen)
        ImGui::CloseCurrentPopup();

    // ------------------------------------------------------------------
    // Edit Structure modal
    // ------------------------------------------------------------------

    responsive::windowSize(ImVec2(980.0f, 640.0f), ImGuiCond_FirstUseEver);
    bool editStructureOpen = true;
    if (m_step.active())
    {
        if (responsive::beginModal("Edit Structure##edit", &editStructureOpen,
                                   ImGuiWindowFlags_NoResize))
        {
            drawStepEditor();
            ImGui::EndPopup();
        }
        else
        {
            // Closed (X): the menu path edits the structure again.
            m_step.finish();
        }
    }
    else if (responsive::beginModal("Edit Structure##edit", &editStructureOpen,
                                    ImGuiWindowFlags_NoResize))
    {
        ImGui::TextDisabled("Edit the cell and each atom; click an element to substitute it. Edit > Undo reverts changes.");

        // Lattice editing
        dialogLayout::section("Cell Vectors (\xc3\x85)");
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("This structure has no unit cell.");
        }
        else
        {
            bool latticeChanged = false;
            const char* names[3] = {"a", "b", "c"};
            const double step = 0.01;
            for (int r = 0; r < 3; ++r)
            {
                ImGui::SetNextItemWidth(responsive::dp(520.0f));
                latticeChanged |= ImGui::DragScalarN(names[r], ImGuiDataType_Double, structure.cellVectors[r].data(), 3,
                                                     static_cast<float>(step), nullptr, nullptr, "%.5f");
            }
            if (latticeChanged)
                updateBuffers(structure);
        }

        const std::string atomsTitle = "Atoms (" + std::to_string(structure.atoms.size()) + ")";
        dialogLayout::section(atomsTitle.c_str());

        const char* modeLabels[] = { "Cartesian (\xc3\x85)", "Direct (fractional)" };
        int mode = m_useDirectCoords ? 1 : 0;
        if (!structure.hasUnitCell) ImGui::BeginDisabled();
        ImGui::SetNextItemWidth(responsive::dp(220.0f));
        if (ImGui::Combo("Coordinates", &mode, modeLabels, 2))
            m_useDirectCoords = (mode == 1);
        if (!structure.hasUnitCell) ImGui::EndDisabled();
        if (!structure.hasUnitCell)
            m_useDirectCoords = false;

        ImGui::SameLine();
        if (dialogLayout::secondaryButton("Add Atom"))
        {
            AtomSite newAtom;
            newAtom.atomicNumber = 1;
            newAtom.symbol = "H";
            newAtom.x = 0.0;
            newAtom.y = 0.0;
            newAtom.z = 0.0;
            if (1 < (int)elementColors.size())
            {
                newAtom.r = elementColors[1].r;
                newAtom.g = elementColors[1].g;
                newAtom.b = elementColors[1].b;
            }
            else
            {
                getDefaultElementColor(1, newAtom.r, newAtom.g, newAtom.b);
            }
            structure.atoms.push_back(newAtom);
            m_selectedEditAtom = (int)structure.atoms.size() - 1;
            m_scrollEditRowsToBottom = true;
            std::cout << "[Operation] Added atom (edit structure): "
                      << newAtom.symbol << "(" << newAtom.atomicNumber << ")"
                      << " at [" << newAtom.x << ", " << newAtom.y << ", " << newAtom.z << "]"
                      << std::endl;
            updateBuffers(structure);
        }
        dialogLayout::hint("Adds a hydrogen atom at the origin; set its element and position in its row");

        const float footerHeight = ImGui::GetFrameHeight() + 3.0f * ImGui::GetStyle().ItemSpacing.y + responsive::dp(2) + 2.0f;
        if (structure.atoms.empty())
        {
            ImGui::TextDisabled("No atoms. Add Atom places one at the origin.");
        }
        else
        {
            int pendingDelete = -1;
            bool anyPositionChanged = false;
            const std::string positionHeader = m_useDirectCoords ? "Position (u, v, w)" : "Position (x, y, z in \xc3\x85)";

            const ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersInnerH;
            if (ImGui::BeginTable("##atom-edit-rows", 4, flags, ImVec2(0.0f, std::max(responsive::dp(160.0f), ImGui::GetContentRegionAvail().y - footerHeight))))
            {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("000000").x);
                ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("Element").x + ImGui::GetStyle().FramePadding.x * 2);
                ImGui::TableSetupColumn(positionHeader.c_str(), ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("##delete", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("Delete").x + ImGui::GetStyle().FramePadding.x * 2);
                ImGui::TableHeadersRow();

                ImGuiListClipper clipper;
                clipper.Begin(static_cast<int>(structure.atoms.size()));
                while (clipper.Step())
                {
                    for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; ++i)
                    {
                        AtomSite& atom = structure.atoms[static_cast<std::size_t>(i)];
                        ImGui::PushID(i);
                        ImGui::TableNextRow();

                        ImGui::TableSetColumnIndex(0);
                        ImGui::AlignTextToFramePadding();
                        ImGui::TextDisabled("%d", i);

                        ImGui::TableSetColumnIndex(1);
                        std::string elemBtn = std::string(atom.symbol) + "##edit-atom-element";
                        if (ImGui::Button(elemBtn.c_str(), ImVec2(-FLT_MIN, 0)))
                        {
                            m_selectedEditAtom = i;
                            m_editStructureElementTargetAtom = i;
                            int z = atom.atomicNumber;
                            if (isValidElementNumber(z))
                                m_selectedEditElement = z;
                            m_showEditStructureElementPicker = true;
                            m_restoreEditStructureAfterElementPicker = true;
                        }
                        dialogLayout::hint("Substitute this atom's element");

                        ImGui::TableSetColumnIndex(2);
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        if (m_useDirectCoords && structure.hasUnitCell)
                        {
                            glm::vec3 frac(0.0f);
                            const glm::vec3 cart((float)atom.x, (float)atom.y, (float)atom.z);
                            tryCartesianToFractional(structure, cart, frac);
                            float editPos[3] = { frac.x, frac.y, frac.z };
                            if (ImGui::DragFloat3("##DirectPos", editPos, 0.001f, -1000.0f, 1000.0f, "%.5f"))
                            {
                                glm::vec3 newCart(0.0f);
                                if (tryFractionalToCartesian(structure, glm::vec3(editPos[0], editPos[1], editPos[2]), newCart))
                                {
                                    atom.x = newCart.x;
                                    atom.y = newCart.y;
                                    atom.z = newCart.z;
                                    anyPositionChanged = true;
                                }
                            }
                        }
                        else
                        {
                            double editPos[3] = { atom.x, atom.y, atom.z };
                            if (ImGui::DragScalarN("##CartesianPos", ImGuiDataType_Double, editPos, 3, 0.005f, nullptr, nullptr, "%.5f"))
                            {
                                atom.x = editPos[0];
                                atom.y = editPos[1];
                                atom.z = editPos[2];
                                anyPositionChanged = true;
                            }
                        }

                        ImGui::TableSetColumnIndex(3);
                        if (dialogLayout::dangerButton("Delete"))
                            pendingDelete = i;

                        ImGui::PopID();
                    }
                }

                if (m_scrollEditRowsToBottom)
                {
                    ImGui::SetScrollY(ImGui::GetScrollMaxY());
                    m_scrollEditRowsToBottom = false;
                }
                ImGui::EndTable();
            }

            if (pendingDelete >= 0 && pendingDelete < (int)structure.atoms.size())
            {
                const AtomSite deletedAtom = structure.atoms[pendingDelete];
                structure.eraseAtom(static_cast<std::size_t>(pendingDelete));
                if (m_selectedEditAtom >= (int)structure.atoms.size())
                    m_selectedEditAtom = std::max(0, (int)structure.atoms.size() - 1);
                std::cout << "[Operation] Deleted atom (edit structure row): "
                          << deletedAtom.symbol << "(" << deletedAtom.atomicNumber << ")"
                          << " row=" << pendingDelete << std::endl;
                updateBuffers(structure);
            }
            else if (anyPositionChanged)
            {
                updateBuffers(structure);
            }

        }
        if (dialogLayout::footer("Close", nullptr) == 1)
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (!editStructureOpen)
        ImGui::CloseCurrentPopup();

    if (m_showEditStructureElementPicker)
    {
        ImGui::OpenPopup("Edit Atom Element##edit-structure");
        m_showEditStructureElementPicker = false;
    }

    responsive::windowSize(ImVec2(940.0f, 560.0f), ImGuiCond_Appearing);
    bool editAtomElementOpen = true;
    auto clearEditStructureElementPicker = [&]() {
        m_editStructureElementTargetAtom = -1;
        if (m_restoreEditStructureAfterElementPicker)
            m_openEditStructure = true;
        m_restoreEditStructureAfterElementPicker = false;
    };

    if (responsive::beginModal("Edit Atom Element##edit-structure", &editAtomElementOpen,
                               ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::TextDisabled("Click an element or type its symbol; double-click or Enter applies it.");
        bool activated = false;
        drawPeriodicTableInlineSelector(m_selectedEditElement, &activated, true);
        const bool enter = ImGui::IsKeyPressed(ImGuiKey_Enter, false) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter, false);
        const int pressed = dialogLayout::footer("Apply", "Cancel");
        if (pressed == 2)
            editAtomElementOpen = false;
        if (pressed == 1 || activated || enter)
        {
            if (m_editStructureElementTargetAtom >= 0 &&
                m_editStructureElementTargetAtom < (int)structure.atoms.size())
            {
                AtomSite& target = structure.atoms[m_editStructureElementTargetAtom];
                const int oldAtomicNumber = target.atomicNumber;
                const std::string oldSymbol = target.symbol;
                target.atomicNumber = m_selectedEditElement;
                target.symbol = elementSymbol(m_selectedEditElement);
                if (m_selectedEditElement >= 0 &&
                    m_selectedEditElement < (int)elementColors.size())
                {
                    target.r = elementColors[m_selectedEditElement].r;
                    target.g = elementColors[m_selectedEditElement].g;
                    target.b = elementColors[m_selectedEditElement].b;
                }
                else
                {
                    getDefaultElementColor(m_selectedEditElement, target.r, target.g, target.b);
                }
                std::cout << "[Operation] Substituted atom (edit structure): "
                          << "row=" << m_editStructureElementTargetAtom << " "
                          << oldSymbol << "(" << oldAtomicNumber << ") -> "
                          << target.symbol << "(" << target.atomicNumber << ")"
                          << std::endl;
                updateBuffers(structure);
            }

            clearEditStructureElementPicker();
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }
    if (!editAtomElementOpen)
    {
        clearEditStructureElementPicker();
        ImGui::CloseCurrentPopup();
    }
}
