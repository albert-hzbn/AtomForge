#include "ui/BulkCrystalBuilderDialog.h"

#include "algorithms/BulkCrystalBuilder.h"
#include "ElementData.h"
#include "imgui.h"
#include "ui/PeriodicTableDialog.h"
#include "ui/DialogLayout.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
// UI-only helper that remains local to the dialog.
void drawLatticeParameterInputs(CrystalSystem system, LatticeParameters& latticeParams)
{
    const float fieldWidth = 100.0f;
    const float spacing = 8.0f;
    const float hintGap = 20.0f;

    auto inputField = [&](const char* label, double* val) {
        ImGui::PushItemWidth(fieldWidth);
        ImGui::InputDouble(label, val, 0.0, 0.0, "%.4f");
        ImGui::PopItemWidth();
    };

    switch (system)
    {
        case CrystalSystem::Triclinic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            inputField("alpha", &latticeParams.alpha);
            ImGui::SameLine(0.0f, spacing); inputField("beta", &latticeParams.beta);
            ImGui::SameLine(0.0f, spacing); inputField("gamma", &latticeParams.gamma);
            break;
        case CrystalSystem::Monoclinic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            inputField("beta", &latticeParams.beta);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("alpha=gamma=90");
            latticeParams.alpha = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Orthorhombic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("alpha=beta=gamma=90");
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Tetragonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=a, angles=90");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Trigonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("hex: b=a, gamma=120");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Hexagonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=a, gamma=120");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Cubic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=c=a, angles=90");
            latticeParams.b = latticeParams.a;
            latticeParams.c = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
    }
}
} // namespace

BulkCrystalBuilderDialog::BulkCrystalBuilderDialog() = default;

void BulkCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Bulk Crystal", NULL, false, enabled))
        m_openRequested = true;
}

void BulkCrystalBuilderDialog::drawDialog(Structure& structure,
                                          const std::vector<glm::vec3>& elementColors,
                                          const std::function<void(Structure&)>& updateBuffers)
{
    static int crystalSystemIndex = (int)CrystalSystem::Cubic;
    static int selectedSpaceGroup = 225;
    static LatticeParameters latticeParams;
    static std::vector<AtomSite> asymmetricAtoms;
    static BulkBuildResult lastResult;
    static int selectedEditElement = 6;
    static int targetAtomIndex = -1;
    static bool showElementPicker = false;
    static bool reopenBulkDialogAfterPicker = false;
    static bool scrollRowsToBottom = false;
    static int lastCrystalSystemIndex = crystalSystemIndex;

    if (asymmetricAtoms.empty())
        addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);

    if (m_openRequested)
    {
        ImGui::OpenPopup("Build Bulk Crystal");
        m_openRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(720.0f, 540.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(640,520),ImVec2(FLT_MAX,FLT_MAX));
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("Build Bulk Crystal", &dialogOpen, ImGuiWindowFlags_None))
    {
#ifndef ATOMS_ENABLE_SPGLIB
        ImGui::TextWrapped("spglib is not available in this build, so symmetry expansion cannot be generated.");
#else
        const char* systemLabels[] = {
            "Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic"
        };

        const float halfWidth = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;
        dialogLayout::section("Crystal symmetry");
        ImGui::TextUnformatted("Crystal system");
        ImGui::SameLine(ImGui::GetStyle().WindowPadding.x+halfWidth+ImGui::GetStyle().ItemSpacing.x);
        ImGui::TextUnformatted("Space group");
        ImGui::PushItemWidth(halfWidth);
        if (ImGui::Combo("##system", &crystalSystemIndex, systemLabels, 7))
        {
            const SpaceGroupRange& range = currentRange((CrystalSystem)crystalSystemIndex);
            selectedSpaceGroup = range.first;
        }
        ImGui::PopItemWidth();
        if (crystalSystemIndex != lastCrystalSystemIndex)
        {
            applySystemConstraints((CrystalSystem)crystalSystemIndex, latticeParams);
            lastCrystalSystemIndex = crystalSystemIndex;
        }

        const SpaceGroupRange& range = currentRange((CrystalSystem)crystalSystemIndex);
        if (selectedSpaceGroup < range.first || selectedSpaceGroup > range.last)
            selectedSpaceGroup = range.first;
        int sgIndex = selectedSpaceGroup - range.first;
        auto sgGetter = [](void* data, int idx) -> const char* {
            static char label[64];
            const SpaceGroupRange* sgRange = static_cast<const SpaceGroupRange*>(data);
            int sgNumber = sgRange->first + idx;
            std::snprintf(label, sizeof(label), "%s", spaceGroupLabel(sgNumber).c_str());
            return label;
        };
        ImGui::SameLine();
        ImGui::PushItemWidth(halfWidth);
        if (ImGui::Combo("##spacegroup", &sgIndex, sgGetter, (void*)&range, range.last - range.first + 1))
            selectedSpaceGroup = range.first + sgIndex;
        ImGui::PopItemWidth();

        dialogLayout::section("Lattice parameters (A / degrees)");
        drawLatticeParameterInputs((CrystalSystem)crystalSystemIndex, latticeParams);

        ImGui::Separator();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Asymmetric unit");
        ImGui::SameLine();
        if (ImGui::Button("+ Add"))
        {
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            scrollRowsToBottom = true;
            const AtomSite& atom = asymmetricAtoms.back();
            std::cout << "[Operation] Added asymmetric atom (bulk builder): "
                      << atom.symbol << "(" << atom.atomicNumber << ")"
                      << " at [" << atom.x << ", " << atom.y << ", " << atom.z << "]"
                      << std::endl;
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear"))
        {
            asymmetricAtoms.clear();
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            std::cout << "[Operation] Cleared asymmetric atoms (bulk builder); reset to default H(1)" << std::endl;
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(fractional coords)");

        int pendingDelete = -1;
        float listHeight = ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 2 - ImGui::GetStyle().ItemSpacing.y * 2;
        if (listHeight < 60.0f) listHeight = 60.0f;
        if (ImGui::BeginChild("##bulk-atom-rows", ImVec2(-1.0f, listHeight), true))
        {
            ImGui::TextDisabled("Sites / fractional coordinates (x, y, z) / element");
            ImGui::Separator();
            for (int i = 0; i < (int)asymmetricAtoms.size(); ++i)
            {
                AtomSite& atom = asymmetricAtoms[i];
                ImGui::PushID(i);

                if (ImGui::Button("Delete"))
                    pendingDelete = i;
                ImGui::SameLine();
                ImGui::Text("#%d", i);
                ImGui::SameLine();

                float editPos[3] = { (float)atom.x, (float)atom.y, (float)atom.z };

                ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 50.0f);
                if (ImGui::DragFloat3("##bulk-direct", editPos, 0.005f, -10.0f, 10.0f, "%.4f"))
                {
                    atom.x = editPos[0];
                    atom.y = editPos[1];
                    atom.z = editPos[2];
                }
                ImGui::PopItemWidth();

                ImGui::SameLine();
                std::string elemButton = std::string(atom.symbol) + "##bulk-element";
                if (ImGui::Button(elemButton.c_str()))
                {
                    targetAtomIndex = i;
                    selectedEditElement = atom.atomicNumber;
                    showElementPicker = true;
                    reopenBulkDialogAfterPicker = true;
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Click to substitute element");

                ImGui::PopID();
            }

            if (scrollRowsToBottom)
            {
                ImGui::SetScrollY(ImGui::GetScrollMaxY());
                scrollRowsToBottom = false;
            }
        }
        ImGui::EndChild();

        if (pendingDelete >= 0 && pendingDelete < (int)asymmetricAtoms.size())
        {
            const AtomSite deletedAtom = asymmetricAtoms[pendingDelete];
            asymmetricAtoms.erase(asymmetricAtoms.begin() + pendingDelete);
            if (asymmetricAtoms.empty())
                addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            std::cout << "[Operation] Deleted asymmetric atom (bulk builder): "
                      << deletedAtom.symbol << "(" << deletedAtom.atomicNumber << ")"
                      << " row=" << pendingDelete << std::endl;
        }

        ImGui::Separator();
        if (dialogLayout::primaryButton("Build",dialogLayout::actionSize()))
        {
            applySystemConstraints((CrystalSystem)crystalSystemIndex, latticeParams);
            lastResult = buildBulkCrystal(structure,
                                          (CrystalSystem)crystalSystemIndex,
                                          selectedSpaceGroup,
                                          latticeParams,
                                          asymmetricAtoms,
                                          elementColors);
            if (lastResult.success)
            {
                updateBuffers(structure);
                std::cout << "[Operation] Built bulk crystal: system="
                          << crystalSystemLabel((CrystalSystem)crystalSystemIndex)
                          << ", SG=" << lastResult.spaceGroupNumber
                          << " (" << lastResult.spaceGroupSymbol << ")"
                          << ", asym_atoms=" << asymmetricAtoms.size()
                          << ", generated_atoms=" << lastResult.generatedAtoms
                          << std::endl;
            }
            else
            {
                std::cout << "[Operation] Bulk crystal build failed: " << lastResult.message << std::endl;
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Close",dialogLayout::actionSize()))
        {
            dialogOpen = false;
            ImGui::CloseCurrentPopup();
        }

        if (!lastResult.message.empty())
        {
            if (lastResult.success)
                ImGui::Text("OK: %d atoms, SG %d (%s)", lastResult.generatedAtoms, lastResult.spaceGroupNumber, lastResult.spaceGroupSymbol.c_str());
            else
                ImGui::TextWrapped("Error: %s", lastResult.message.c_str());
        }
#endif

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();

    if (showElementPicker)
    {
        ImGui::OpenPopup("Bulk Crystal Element Picker");
        showElementPicker = false;
    }

    ImGui::SetNextWindowSize(ImVec2(940.0f, 560.0f), ImGuiCond_Appearing);
    bool pickerOpen = true;
    if (ImGui::BeginPopupModal("Bulk Crystal Element Picker", &pickerOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select replacement element from the periodic table.");
        ImGui::Separator();
        drawPeriodicTableInlineSelector(selectedEditElement);
        ImGui::Separator();
        if (ImGui::Button("Apply Element##bulk-builder"))
        {
            if (targetAtomIndex >= 0 && targetAtomIndex < (int)asymmetricAtoms.size())
            {
                AtomSite& atom = asymmetricAtoms[targetAtomIndex];
                const int oldAtomicNumber = atom.atomicNumber;
                const std::string oldSymbol = atom.symbol;
                applyElementToAtom(atom, selectedEditElement, elementColors);
                std::cout << "[Operation] Substituted asymmetric atom (bulk builder): "
                          << "row=" << targetAtomIndex << " "
                          << oldSymbol << "(" << oldAtomicNumber << ") -> "
                          << atom.symbol << "(" << atom.atomicNumber << ")"
                          << std::endl;
            }
            targetAtomIndex = -1;
            if (reopenBulkDialogAfterPicker)
                m_openRequested = true;
            reopenBulkDialogAfterPicker = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }
    if (!pickerOpen)
    {
        targetAtomIndex = -1;
        if (reopenBulkDialogAfterPicker)
            m_openRequested = true;
        reopenBulkDialogAfterPicker = false;
        ImGui::CloseCurrentPopup();
    }
}
