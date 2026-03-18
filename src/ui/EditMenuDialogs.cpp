#include "EditMenuDialogs.h"

#include "ElementData.h"
#include "PeriodicTableDialog.h"
#include "io/StructureLoader.h"
#include "math/StructureMath.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>

namespace
{
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
    if (ImGui::MenuItem("Atomic Sizes..."))
        m_openAtomicSize = true;

    if (ImGui::MenuItem("Element Colors..."))
        m_openElementColor = true;

    if (ImGui::MenuItem("Edit Structure..."))
        m_openEditStructure = true;
}

void EditMenuDialogs::drawPopups(Structure& structure,
                                 const std::function<void(Structure&)>& updateBuffers)
{
    if (m_openAtomicSize)   { ImGui::OpenPopup("Atomic Sizes##edit");   m_openAtomicSize   = false; }
    if (m_openElementColor) { ImGui::OpenPopup("Element Colors##edit"); m_openElementColor = false; }
    if (m_openEditStructure) { ImGui::OpenPopup("Edit Structure##edit"); m_openEditStructure = false; }

    // ------------------------------------------------------------------
    // Atomic Sizes modal
    // ------------------------------------------------------------------

    bool atomicSizesOpen = true;
    if (ImGui::BeginPopupModal("Atomic Sizes##edit", &atomicSizesOpen,
                               ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Select an element in the periodic table and edit its radius.");
        ImGui::Text("Defaults: Cordero et al., Dalton Trans. 2008.");

        if (ImGui::Button("Reset to Literature Defaults"))
        {
            elementRadii = makeLiteratureCovalentRadii();
            updateBuffers(structure);
        }

        ImGui::Separator();
        drawPeriodicTableInlineSelector(m_selectedRadiusElement);

        if (m_selectedRadiusElement >= 1 && m_selectedRadiusElement <= 118)
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

            if (ImGui::Button("Apply Literature Radius To Selected Element"))
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
    // Element Colors modal
    // ------------------------------------------------------------------

    bool elementColorsOpen = true;
    if (ImGui::BeginPopupModal("Element Colors##edit", &elementColorsOpen,
                               ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Select an element in the periodic table and edit its color.");
        ImGui::Separator();

        drawPeriodicTableInlineSelector(m_selectedColorElement);

        if (m_selectedColorElement >= 1 && m_selectedColorElement <= 118)
        {
            ImGui::Text("Selected: Z=%d (%s)",
                        m_selectedColorElement,
                        elementSymbol(m_selectedColorElement));

            float color[3] = {
                elementColors[m_selectedColorElement].r,
                elementColors[m_selectedColorElement].g,
                elementColors[m_selectedColorElement].b
            };

            if (ImGui::ColorEdit3("Element Color", color))
            {
                elementColors[m_selectedColorElement] =
                    glm::vec3(color[0], color[1], color[2]);
                updateBuffers(structure);
            }

            float shininess = elementShininess[m_selectedColorElement];
            if (ImGui::SliderFloat("Material Shininess", &shininess, 4.0f, 128.0f, "%.1f"))
            {
                elementShininess[m_selectedColorElement] = shininess;
                updateBuffers(structure);
            }

            if (ImGui::Button("Reset Selected Element Color"))
            {
                float r, g, b;
                getDefaultElementColor(m_selectedColorElement, r, g, b);
                elementColors[m_selectedColorElement] = glm::vec3(r, g, b);
                updateBuffers(structure);
            }

            ImGui::SameLine();
            if (ImGui::Button("Reset Selected Shininess"))
            {
                elementShininess[m_selectedColorElement] = 32.0f;
                updateBuffers(structure);
            }
        }
        else
        {
            ImGui::Text("Selected: None");
        }

        ImGui::SameLine();
        if (ImGui::Button("Reset All Colors"))
        {
            elementColors = makeDefaultElementColors();
            updateBuffers(structure);
        }

        ImGui::SameLine();
        if (ImGui::Button("Reset All Shininess"))
        {
            elementShininess = makeDefaultElementShininess();
            updateBuffers(structure);
        }

        ImGui::EndPopup();
    }
    if (!elementColorsOpen)
        ImGui::CloseCurrentPopup();

    // ------------------------------------------------------------------
    // Edit Structure modal
    // ------------------------------------------------------------------

    ImGui::SetNextWindowSize(ImVec2(980.0f, 640.0f), ImGuiCond_FirstUseEver);
    bool editStructureOpen = true;
    if (ImGui::BeginPopupModal("Edit Structure##edit", &editStructureOpen,
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

        if (ImGui::Button("Add Atom"))
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

            if (ImGui::BeginChild("##atom-edit-rows", ImVec2(930.0f, 260.0f), true))
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

                    ImGui::PushItemWidth(420.0f);
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
                        if (z >= 1 && z <= 118)
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
                structure.atoms.erase(structure.atoms.begin() + pendingDelete);
                if (m_selectedEditAtom >= (int)structure.atoms.size())
                    m_selectedEditAtom = std::max(0, (int)structure.atoms.size() - 1);
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

    ImGui::SetNextWindowSize(ImVec2(940.0f, 560.0f), ImGuiCond_Appearing);
    bool editAtomElementOpen = true;
    if (ImGui::BeginPopupModal("Edit Atom Element##edit-structure", &editAtomElementOpen,
                               ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select replacement element from the periodic table.");
        ImGui::Separator();

        drawPeriodicTableInlineSelector(m_selectedEditElement);

        ImGui::Separator();
        if (ImGui::Button("Apply Element##edit-structure-window"))
        {
            if (m_editStructureElementTargetAtom >= 0 &&
                m_editStructureElementTargetAtom < (int)structure.atoms.size())
            {
                AtomSite& target = structure.atoms[m_editStructureElementTargetAtom];
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
                updateBuffers(structure);
            }

            m_editStructureElementTargetAtom = -1;
            if (m_restoreEditStructureAfterElementPicker)
                m_openEditStructure = true;
            m_restoreEditStructureAfterElementPicker = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }
    if (!editAtomElementOpen)
    {
        m_editStructureElementTargetAtom = -1;
        if (m_restoreEditStructureAfterElementPicker)
            m_openEditStructure = true;
        m_restoreEditStructureAfterElementPicker = false;
        ImGui::CloseCurrentPopup();
    }
}
