#include "AtomContextMenu.h"

#include "ElementData.h"
#include "PeriodicTableDialog.h"
#include "io/StructureLoader.h"

#include "imgui.h"

#include <iostream>
#include <map>
#include <sstream>

namespace
{
bool isValidIndex(int idx, std::size_t size)
{
    return idx >= 0 && static_cast<std::size_t>(idx) < size;
}

void setAtomColorFromElement(AtomSite& atom,
                             int atomicNumber,
                             const std::vector<glm::vec3>& elementColors)
{
    if (isValidIndex(atomicNumber, elementColors.size()))
    {
        atom.r = elementColors[atomicNumber].r;
        atom.g = elementColors[atomicNumber].g;
        atom.b = elementColors[atomicNumber].b;
        return;
    }

    getDefaultElementColor(atomicNumber, atom.r, atom.g, atom.b);
}

void deselectAtoms(SceneBuffers& sceneBuffers, std::vector<int>& selectedInstanceIndices)
{
    for (int idx : selectedInstanceIndices)
        sceneBuffers.restoreAtomColor(idx);
    selectedInstanceIndices.clear();
}

std::string summarizeFromElements(const std::map<int, int>& fromElementCounts)
{
    std::ostringstream fromSummary;
    bool first = true;
    for (std::map<int, int>::const_iterator it = fromElementCounts.begin(); it != fromElementCounts.end(); ++it)
    {
        if (!first)
            fromSummary << ", ";
        first = false;
        fromSummary << elementSymbol(it->first) << "(" << it->first << "):" << it->second;
    }
    return fromSummary.str();
}

void substituteSelectedAtoms(Structure& structure,
                            SceneBuffers& sceneBuffers,
                            const std::vector<glm::vec3>& elementColors,
                            const std::vector<int>& selectedInstanceIndices,
                            const ElementSelection& selection,
                            const std::function<void(Structure&)>& updateBuffers)
{
    std::map<int, int> fromElementCounts;
    int replacedCount = 0;

    for (int idx : selectedInstanceIndices)
    {
        if (!isValidIndex(idx, sceneBuffers.atomIndices.size()))
            continue;

        const int baseIdx = sceneBuffers.atomIndices[idx];
        if (!isValidIndex(baseIdx, structure.atoms.size()))
            continue;

        AtomSite& atom = structure.atoms[baseIdx];
        fromElementCounts[atom.atomicNumber]++;
        replacedCount++;

        atom.symbol = selection.symbol;
        atom.atomicNumber = selection.atomicNumber;
        setAtomColorFromElement(atom, selection.atomicNumber, elementColors);
    }

    std::cout << "[Operation] Substituted atoms (context menu): "
              << "count=" << replacedCount
              << ", to=" << selection.symbol << "(" << selection.atomicNumber << ")"
              << ", from={" << summarizeFromElements(fromElementCounts) << "}"
              << std::endl;
    updateBuffers(structure);
}

void insertMidpointAtom(Structure& structure,
                        SceneBuffers& sceneBuffers,
                        const std::vector<glm::vec3>& elementColors,
                        const std::vector<int>& selectedInstanceIndices,
                        const ElementSelection& selection,
                        const std::function<void(Structure&)>& updateBuffers)
{
    double sumX = 0.0;
    double sumY = 0.0;
    double sumZ = 0.0;
    int validCount = 0;

    for (int idx : selectedInstanceIndices)
    {
        if (!isValidIndex(idx, sceneBuffers.atomPositions.size()))
            continue;

        // Use rendered world positions so supercell picks keep Cartesian coordinates.
        sumX += sceneBuffers.atomPositions[idx].x;
        sumY += sceneBuffers.atomPositions[idx].y;
        sumZ += sceneBuffers.atomPositions[idx].z;
        validCount++;
    }

    if (validCount < 2)
        return;

    AtomSite newAtom;
    newAtom.symbol = selection.symbol;
    newAtom.atomicNumber = selection.atomicNumber;
    newAtom.x = sumX / static_cast<double>(validCount);
    newAtom.y = sumY / static_cast<double>(validCount);
    newAtom.z = sumZ / static_cast<double>(validCount);
    setAtomColorFromElement(newAtom, selection.atomicNumber, elementColors);

    structure.atoms.push_back(newAtom);
    std::cout << "[Operation] Inserted atom at midpoint: "
              << newAtom.symbol << "(" << newAtom.atomicNumber << ")"
              << " at [" << newAtom.x << ", " << newAtom.y << ", " << newAtom.z << "]"
              << " from_selected=" << validCount
              << std::endl;
    updateBuffers(structure);
}

void applyPeriodicSelection(PeriodicAction pendingAction,
                            Structure& structure,
                            SceneBuffers& sceneBuffers,
                            const std::vector<glm::vec3>& elementColors,
                            const std::vector<int>& selectedInstanceIndices,
                            const ElementSelection& selection,
                            const std::function<void(Structure&)>& updateBuffers)
{
    if (pendingAction == PeriodicAction::Substitute && !selectedInstanceIndices.empty())
    {
        substituteSelectedAtoms(
            structure,
            sceneBuffers,
            elementColors,
            selectedInstanceIndices,
            selection,
            updateBuffers);
        return;
    }

    if (pendingAction == PeriodicAction::InsertMidpoint && selectedInstanceIndices.size() >= 2)
    {
        insertMidpointAtom(
            structure,
            sceneBuffers,
            elementColors,
            selectedInstanceIndices,
            selection,
            updateBuffers);
    }
}
} // namespace

void AtomContextMenu::open()
{
    m_openRequested = true;
}

void AtomContextMenu::openSubstitute()
{
    m_pendingAction = PeriodicAction::Substitute;
    m_ownsPeriodicPopup = true;
    openPeriodicTable();
}

void AtomContextMenu::draw(Structure& structure,
                           SceneBuffers& sceneBuffers,
                           const std::vector<glm::vec3>& elementColors,
                           std::vector<int>& selectedInstanceIndices,
                           AtomRequests& requests,
                           const std::function<void(Structure&)>& updateBuffers)
{
    bool doOpenPeriodicTable = false;

    // Context menu disabled for large structures (CPU caches not available)
    if (sceneBuffers.cpuCachesDisabled)
    {
        m_openRequested = false;
        return;
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("##atomCtx");
        m_openRequested = false;
    }

    if (ImGui::BeginPopup("##atomCtx"))
    {
        if (ImGui::MenuItem("Substitute Atom"))
        {
            m_pendingAction  = PeriodicAction::Substitute;
            doOpenPeriodicTable = true;
        }

        bool canInsert = selectedInstanceIndices.size() >= 2;
        if (!canInsert) ImGui::BeginDisabled();
        if (ImGui::MenuItem("Insert Atom at Midpoint"))
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
            deselectAtoms(sceneBuffers, selectedInstanceIndices);

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
            applyPeriodicSelection(
                m_pendingAction,
                structure,
                sceneBuffers,
                elementColors,
                selectedInstanceIndices,
                selections[0],
                updateBuffers);
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
