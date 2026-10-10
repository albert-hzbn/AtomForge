#include "ui/ResponsiveLayout.h"
#include "EditMenuDialogs.h"

#include "ElementData.h"
#include "PeriodicTableDialog.h"
#include "io/StructureLoader.h"
#include "math/StructureMath.h"
#include "science/ScienceData.h"
#include "ui/ThemeUtils.h"

#include "imgui.h"

#include <algorithm>
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
    ImGui::TextColored(ImVec4(0.45f, 0.75f, 1.0f, 1.0f), "Editing pipeline step %s: Update step writes these settings into the pipeline.", m_step.step.c_str());
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
        ImGui::DragScalarN(m_stepFractional ? "Position (fractional)" : "Position (A)", ImGuiDataType_Double,
                           m_stepPosition, 3, 0.005f, nullptr, nullptr, "%.6f");
        if (m_stepFractional && !input.hasUnitCell)
            ImGui::TextColored(themeStatusWarn(), "A fractional position needs a structure with a unit cell.");
        ImGui::Checkbox("Select the new atom", &m_stepSelect);
    }

    ImGui::Separator();
    if (responsive::button("Update step"))
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
        ImGui::Text("Select an element in the periodic table and edit its radius.");
        ImGui::Text("Defaults: Cordero et al., Dalton Trans. 2008.");

        if (responsive::button("Reset to Literature Defaults"))
        {
            elementRadii = makeLiteratureCovalentRadii();
            updateBuffers(structure);
        }

        ImGui::Separator();
        drawPeriodicTableInlineSelector(m_selectedRadiusElement);

        if (isValidElementNumber(m_selectedRadiusElement))
        {
            ImGui::Text("Selected: Z=%d (%s)",
                        m_selectedRadiusElement,
                        elementSymbol(m_selectedRadiusElement));

            float radius = elementRadii[m_selectedRadiusElement];
            if (ImGui::DragFloat("Atomic Radius (A)", &radius,
                                 0.01f, 0.20f, 3.50f, "%.2f A"))
            {
                elementRadii[m_selectedRadiusElement] = radius;
                updateBuffers(structure);
            }

            if (responsive::button("Apply Literature Radius To Selected Element"))
            {
                const std::vector<float> defaults = makeLiteratureCovalentRadii();
                elementRadii[m_selectedRadiusElement] =
                    defaults[m_selectedRadiusElement];
                updateBuffers(structure);
            }
        }
        else
        {
            ImGui::Text("Selected: None");
        }

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
            if (ImGui::BeginTabItem("Lighting"))
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
                         "Color vividness — 1.0 = neutral, >1 = more vivid (default 1.55).");
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
            if (ImGui::BeginTabItem("Material"))
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
            if (ImGui::BeginTabItem("Colors"))
            {
                ImGui::Spacing();
                ImGui::TextDisabled("Per-element display color.");
                ImGui::Spacing();

                drawPeriodicTableInlineSelector(m_selectedColorElement);

                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();

                if (isValidElementNumber(m_selectedColorElement))
                {
                    ImGui::Text("Selected: Z=%d  %s",
                                m_selectedColorElement,
                                elementSymbol(m_selectedColorElement));
                    ImGui::Spacing();

                    float color[3] = {
                        elementColors[m_selectedColorElement].r,
                        elementColors[m_selectedColorElement].g,
                        elementColors[m_selectedColorElement].b
                    };
                    ImGui::SetNextItemWidth(kSliderW);
                    if (ImGui::ColorEdit3("Color##el", color))
                    {
                        elementColors[m_selectedColorElement] =
                            glm::vec3(color[0], color[1], color[2]);
                        updateBuffers(structure);
                    }

                    ImGui::Spacing();
                    if (responsive::button("Reset Color##sel"))
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
                if (responsive::button("Reset All Colors"))
                {
                    elementColors = makeDefaultElementColors();
                    updateBuffers(structure);
                }
                ImGui::Spacing();
                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }

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
        ImGui::Text("Modify lattice vectors and atom list (add/edit/delete).\n"
                "Click the element name next to position to substitute atom.");
        ImGui::Separator();

        // Lattice editing
        ImGui::Text("Lattice Vectors (Cartesian)");
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("No unit cell present in this structure.");
        }
        else
        {
            float aVec[3] = {
                (float)structure.cellVectors[0][0],
                (float)structure.cellVectors[0][1],
                (float)structure.cellVectors[0][2]
            };
            float bVec[3] = {
                (float)structure.cellVectors[1][0],
                (float)structure.cellVectors[1][1],
                (float)structure.cellVectors[1][2]
            };
            float cVec[3] = {
                (float)structure.cellVectors[2][0],
                (float)structure.cellVectors[2][1],
                (float)structure.cellVectors[2][2]
            };

            bool latticeChanged = false;
            latticeChanged |= ImGui::DragFloat3("a", aVec, 0.01f, -1000.0f, 1000.0f, "%.6f");
            latticeChanged |= ImGui::DragFloat3("b", bVec, 0.01f, -1000.0f, 1000.0f, "%.6f");
            latticeChanged |= ImGui::DragFloat3("c", cVec, 0.01f, -1000.0f, 1000.0f, "%.6f");

            if (latticeChanged)
            {
                structure.cellVectors[0][0] = aVec[0];
                structure.cellVectors[0][1] = aVec[1];
                structure.cellVectors[0][2] = aVec[2];
                structure.cellVectors[1][0] = bVec[0];
                structure.cellVectors[1][1] = bVec[1];
                structure.cellVectors[1][2] = bVec[2];
                structure.cellVectors[2][0] = cVec[0];
                structure.cellVectors[2][1] = cVec[1];
                structure.cellVectors[2][2] = cVec[2];
                updateBuffers(structure);
            }
        }

        ImGui::Separator();
        ImGui::Text("Atom Editing");

        const char* modeLabels[] = { "Cartesian", "Direct" };
        int mode = m_useDirectCoords ? 1 : 0;
        if (!structure.hasUnitCell) ImGui::BeginDisabled();
        if (ImGui::Combo("Position Mode", &mode, modeLabels, 2))
            m_useDirectCoords = (mode == 1);
        if (!structure.hasUnitCell) ImGui::EndDisabled();
        if (!structure.hasUnitCell)
            m_useDirectCoords = false;

        if (responsive::button("Add Atom"))
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

        ImGui::SameLine();
        ImGui::TextDisabled("Rows update automatically as atoms are added/removed.");

        if (structure.atoms.empty())
        {
            ImGui::TextDisabled("No atoms available.");
        }
        else
        {
            ImGui::Separator();
            ImGui::Text("Per-atom rows: Delete | Position | Element");

            int pendingDelete = -1;
            bool anyPositionChanged = false;

            if (responsive::beginChild("##atom-edit-rows", responsive::size(930.0f,260.0f), true))
            {
                for (int i = 0; i < (int)structure.atoms.size(); ++i)
                {
                    AtomSite& atom = structure.atoms[i];

                    ImGui::PushID(i);
                    if (ImGui::SmallButton("Delete"))
                        pendingDelete = i;

                    ImGui::SameLine();
                    ImGui::Text("#%d", i);

                    ImGui::SameLine();

                    glm::vec3 cart((float)atom.x, (float)atom.y, (float)atom.z);
                    float editPos[3] = { cart.x, cart.y, cart.z };

                    if (m_useDirectCoords && structure.hasUnitCell)
                    {
                        glm::vec3 frac(0.0f);
                        if (tryCartesianToFractional(structure, cart, frac))
                        {
                            editPos[0] = frac.x;
                            editPos[1] = frac.y;
                            editPos[2] = frac.z;
                        }
                    }

                    ImGui::PushItemWidth(responsive::dp(420.0f));
                    bool rowPosChanged = ImGui::DragFloat3(m_useDirectCoords ? "##DirectPos" : "##CartesianPos",
                                                           editPos, 0.005f, -1000.0f, 1000.0f, "%.6f");
                    ImGui::PopItemWidth();

                    if (rowPosChanged)
                    {
                        if (m_useDirectCoords && structure.hasUnitCell)
                        {
                            glm::vec3 frac(editPos[0], editPos[1], editPos[2]);
                            glm::vec3 newCart(0.0f);
                            if (tryFractionalToCartesian(structure, frac, newCart))
                            {
                                atom.x = newCart.x;
                                atom.y = newCart.y;
                                atom.z = newCart.z;
                            }
                        }
                        else
                        {
                            atom.x = editPos[0];
                            atom.y = editPos[1];
                            atom.z = editPos[2];
                        }
                        anyPositionChanged = true;
                    }

                    ImGui::SameLine();
                    std::string elemBtn = std::string(atom.symbol) + "##edit-atom-element";
                    if (ImGui::SmallButton(elemBtn.c_str()))
                    {
                        m_selectedEditAtom = i;
                        m_editStructureElementTargetAtom = i;
                        int z = atom.atomicNumber;
                        if (isValidElementNumber(z))
                            m_selectedEditElement = z;
                        m_showEditStructureElementPicker = true;
                        m_restoreEditStructureAfterElementPicker = true;
                    }
                    if (ImGui::IsItemHovered())
                        ImGui::SetTooltip("Click to substitute element");

                    ImGui::PopID();
                }

                if (m_scrollEditRowsToBottom)
                {
                    ImGui::SetScrollY(ImGui::GetScrollMaxY());
                    m_scrollEditRowsToBottom = false;
                }
                ImGui::EndChild();
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
                               ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select replacement element from the periodic table.");
        ImGui::Separator();

        drawPeriodicTableInlineSelector(m_selectedEditElement);

        ImGui::Separator();
        if (responsive::button("Apply Element##edit-structure-window"))
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
