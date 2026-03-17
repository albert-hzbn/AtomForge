#include "AtomContextMenu.h"

#include "PeriodicTableDialog.h"
#include "io/StructureLoader.h"

#include "imgui.h"

void AtomContextMenu::open()
{
    m_openRequested = true;
}

void AtomContextMenu::draw(Structure& structure,
                           SceneBuffers& sceneBuffers,
                           const std::vector<glm::vec3>& elementColors,
                           std::vector<int>& selectedInstanceIndices,
                           AtomRequests& requests,
                           const std::function<void(Structure&)>& updateBuffers)
{
    bool doOpenPeriodicTable = false;

    if (m_openRequested)
    {
        ImGui::OpenPopup("##atomCtx");
        m_openRequested = false;
    }

    if (ImGui::BeginPopup("##atomCtx"))
    {
        if (ImGui::MenuItem("Substitute Atom..."))
        {
            m_pendingAction  = PeriodicAction::Substitute;
            doOpenPeriodicTable = true;
        }

        bool canInsert = selectedInstanceIndices.size() >= 2;
        if (!canInsert) ImGui::BeginDisabled();
        if (ImGui::MenuItem("Insert Atom at Midpoint..."))
        {
            m_pendingAction  = PeriodicAction::InsertMidpoint;
            doOpenPeriodicTable = true;
        }
        if (!canInsert) ImGui::EndDisabled();

        bool canMeasureDistance = selectedInstanceIndices.size() == 2;
        if (!canMeasureDistance) ImGui::BeginDisabled();
        if (ImGui::MenuItem("Measure Distance"))
            requests.measureDistance = true;
        if (!canMeasureDistance) ImGui::EndDisabled();

        bool canMeasureAngle = selectedInstanceIndices.size() == 3;
        if (!canMeasureAngle) ImGui::BeginDisabled();
        if (ImGui::MenuItem("Measure Angle"))
            requests.measureAngle = true;
        if (!canMeasureAngle) ImGui::EndDisabled();

        bool canShowInfo = selectedInstanceIndices.size() == 1;
        if (!canShowInfo) ImGui::BeginDisabled();
        if (ImGui::MenuItem("Atom Info"))
            requests.atomInfo = true;
        if (!canShowInfo) ImGui::EndDisabled();

        if (ImGui::MenuItem("Delete"))
            requests.doDelete = true;

        if (ImGui::MenuItem("Deselect"))
        {
            for (int idx : selectedInstanceIndices)
                sceneBuffers.restoreAtomColor(idx);
            selectedInstanceIndices.clear();
        }

        ImGui::EndPopup();
    }

    if (doOpenPeriodicTable)
    {
        openPeriodicTable();
        m_ownsPeriodicPopup = true;
    }

    // ------------------------------------------------------------------
    // Periodic table picker – process confirmed selection
    // ------------------------------------------------------------------

    std::vector<ElementSelection> selections;
    if (m_ownsPeriodicPopup && drawPeriodicTable(selections))
    {
        if (!selections.empty())
        {
            const auto& sel = selections[0];

            if (m_pendingAction == PeriodicAction::Substitute &&
                !selectedInstanceIndices.empty())
            {
                for (int idx : selectedInstanceIndices)
                {
                    if (idx < 0 || idx >= (int)sceneBuffers.atomIndices.size())
                        continue;
                    int baseIdx = sceneBuffers.atomIndices[idx];
                    if (baseIdx < 0 || baseIdx >= (int)structure.atoms.size())
                        continue;

                    structure.atoms[baseIdx].symbol      = sel.symbol;
                    structure.atoms[baseIdx].atomicNumber = sel.atomicNumber;
                    if (sel.atomicNumber >= 0 &&
                        sel.atomicNumber < (int)elementColors.size())
                    {
                        structure.atoms[baseIdx].r = elementColors[sel.atomicNumber].r;
                        structure.atoms[baseIdx].g = elementColors[sel.atomicNumber].g;
                        structure.atoms[baseIdx].b = elementColors[sel.atomicNumber].b;
                    }
                }
                updateBuffers(structure);
            }

            if (m_pendingAction == PeriodicAction::InsertMidpoint &&
                selectedInstanceIndices.size() >= 2)
            {
                double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
                int validCount = 0;

                for (int idx : selectedInstanceIndices)
                {
                        if (idx < 0 || idx >= (int)sceneBuffers.atomPositions.size())
                        continue;
                        // Use the actual rendered world position so that picking a
                        // supercell copy gives the correct Cartesian coordinates.
                        sumX += sceneBuffers.atomPositions[idx].x;
                        sumY += sceneBuffers.atomPositions[idx].y;
                        sumZ += sceneBuffers.atomPositions[idx].z;
                    ++validCount;
                }

                if (validCount >= 2)
                {
                    AtomSite newAtom;
                    newAtom.symbol      = sel.symbol;
                    newAtom.atomicNumber = sel.atomicNumber;
                    newAtom.x = sumX / (double)validCount;
                    newAtom.y = sumY / (double)validCount;
                    newAtom.z = sumZ / (double)validCount;

                    if (sel.atomicNumber >= 0 &&
                        sel.atomicNumber < (int)elementColors.size())
                    {
                        newAtom.r = elementColors[sel.atomicNumber].r;
                        newAtom.g = elementColors[sel.atomicNumber].g;
                        newAtom.b = elementColors[sel.atomicNumber].b;
                    }
                    else
                    {
                        getDefaultElementColor(sel.atomicNumber,
                                               newAtom.r, newAtom.g, newAtom.b);
                    }

                    structure.atoms.push_back(newAtom);
                    updateBuffers(structure);
                }
            }
        }

        m_pendingAction = PeriodicAction::None;
        m_ownsPeriodicPopup = false;
    }

    // If the picker was dismissed without a selection, clear the pending action.
    if (m_ownsPeriodicPopup && !ImGui::IsPopupOpen("Periodic Table##picker"))
    {
        m_pendingAction = PeriodicAction::None;
        m_ownsPeriodicPopup = false;
    }
}
