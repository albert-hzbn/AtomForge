#include "ui/ResponsiveLayout.h"
#include "ui/UiScript.h"
#include "ui/BulkCrystalBuilderDialog.h"

#include "algorithms/BulkCrystalBuilder.h"
#include "ElementData.h"
#include "imgui.h"
#include "ui/PeriodicTableDialog.h"
#include "ui/DialogLayout.h"
#include "util/PathUtils.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cctype>
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
    const float fieldWidth = responsive::dp(100.0f);
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
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            inputField("\xce\xb1##alpha", &latticeParams.alpha);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("\xce\xb2##beta", &latticeParams.beta);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("\xce\xb3##gamma", &latticeParams.gamma);
            break;
        case CrystalSystem::Monoclinic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            inputField("\xce\xb2##beta", &latticeParams.beta);
            ImGui::SameLine(responsive::dp(0.0f), hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("\xce\xb1 = \xce\xb3 = 90\xc2\xb0");
            latticeParams.alpha = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Orthorhombic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("\xce\xb1 = \xce\xb2 = \xce\xb3 = 90\xc2\xb0");
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Tetragonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(responsive::dp(0.0f), hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b = a, all angles 90\xc2\xb0");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Trigonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(responsive::dp(0.0f), hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("hexagonal axes: b = a, \xce\xb3 = 120\xc2\xb0");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Hexagonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(responsive::dp(0.0f), hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b = a, \xce\xb3 = 120\xc2\xb0");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Cubic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(responsive::dp(0.0f), hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b = c = a, all angles 90\xc2\xb0");
            latticeParams.b = latticeParams.a;
            latticeParams.c = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
    }
}

// Names of the crystal systems as --system takes them, in CrystalSystem order.
const char* const kSystemNames[] = {
    "triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"
};
} // namespace

BulkCrystalBuilderDialog::BulkCrystalBuilderDialog() = default;

void BulkCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Bulk Crystal", NULL, false, enabled))
    {
        m_step.finish();
        m_openRequested = true;
    }
}

bool BulkCrystalBuilderDialog::editStep(StepEdit edit)
{
    m_step = std::move(edit);
    const auto reader = m_step.reader();

    // Missing flags take the defaults of AtomForge --build bulk.
    m_crystalSystemIndex = (int)CrystalSystem::Cubic;
    std::string system = reader.text("--system", "cubic");
    for (char& c : system) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (int i = 0; i < 7; ++i)
        if (system == kSystemNames[i]) m_crystalSystemIndex = i;
    m_lastCrystalSystemIndex = m_crystalSystemIndex;
    m_selectedSpaceGroup = reader.integer("--spacegroup", 225);
    m_originChoice = reader.integer("--origin", 1) == 2 ? 2 : 1;

    m_latticeParams.a = reader.number("--a", 4.0);
    m_latticeParams.b = reader.number("--b", m_latticeParams.a);
    m_latticeParams.c = reader.number("--c", m_latticeParams.a);
    m_latticeParams.alpha = reader.number("--alpha", 90.0);
    m_latticeParams.beta = reader.number("--beta", 90.0);
    m_latticeParams.gamma = reader.number("--gamma", 90.0);
    applySystemConstraints((CrystalSystem)m_crystalSystemIndex, m_latticeParams);

    // Each --atom is "SYMBOL fx fy fz"; colours are applied in drawDialog.
    m_asymmetricAtoms.clear();
    for (const auto& value : reader.values("--atom"))
    {
        std::istringstream in(value);
        std::string symbol;
        AtomSite site;
        if (!(in >> symbol >> site.x >> site.y >> site.z)) continue;
        const int z = atomicNumberFromSymbol(symbol);
        if (z <= 0) continue;
        site.atomicNumber = z;
        site.symbol = elementSymbol(z);
        m_asymmetricAtoms.push_back(site);
    }
    m_recolorAtoms = true;
    m_openRequested = true;
    return true;
}

std::string BulkCrystalBuilderDialog::stepOptions() const
{
    LatticeParameters lp = m_latticeParams;
    applySystemConstraints((CrystalSystem)m_crystalSystemIndex, lp);
    atomforge::pipeline::options::Writer out;
    out.add("--system", kSystemNames[std::clamp(m_crystalSystemIndex, 0, 6)])
       .add("--spacegroup", m_selectedSpaceGroup)
       .add("--a", lp.a).add("--b", lp.b).add("--c", lp.c)
       .add("--alpha", lp.alpha).add("--beta", lp.beta).add("--gamma", lp.gamma);
    if (hasTwoOrigins(m_selectedSpaceGroup) && m_originChoice == 2) out.add("--origin", 2);
    for (const AtomSite& atom : m_asymmetricAtoms)
        out.add("--atom", atom.symbol + " " + out.number(atom.x) + " " + out.number(atom.y) + " " + out.number(atom.z));
    return out.str();
}

void BulkCrystalBuilderDialog::drawDialog(Structure& structure,
                                          const std::vector<glm::vec3>& elementColors,
                                          const std::function<void(Structure&)>& updateBuffers)
{
    int& crystalSystemIndex = m_crystalSystemIndex;
    int& selectedSpaceGroup = m_selectedSpaceGroup;
    LatticeParameters& latticeParams = m_latticeParams;
    std::vector<AtomSite>& asymmetricAtoms = m_asymmetricAtoms;
    static BulkBuildResult lastResult;
    static int selectedEditElement = 6;
    static int targetAtomIndex = -1;
    static bool showElementPicker = false;
    static bool reopenBulkDialogAfterPicker = false;
    static bool scrollRowsToBottom = false;
    int& lastCrystalSystemIndex = m_lastCrystalSystemIndex;

    if (m_recolorAtoms)
    {
        for (AtomSite& atom : asymmetricAtoms)
            applyElementToAtom(atom, atom.atomicNumber, elementColors);
        m_recolorAtoms = false;
    }
    if (asymmetricAtoms.empty())
        addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);

    if (m_openRequested)
    {
        ImGui::OpenPopup("Build Bulk Crystal");
        m_openRequested = false;
    }

    responsive::windowSize(ImVec2(720.0f, 450.0f), ImGuiCond_FirstUseEver);
    responsive::windowConstraints(ImVec2(640,520),ImVec2(FLT_MAX,FLT_MAX));
    bool dialogOpen = true;
    m_isOpen = responsive::beginModal("Build Bulk Crystal", &dialogOpen, ImGuiWindowFlags_None);
    if (m_isOpen)
    {
        if (m_step.active())
            dialogLayout::stepBanner(m_step.step.c_str());
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
        uiScript::openCombo("bulk.systems", "##system");
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
        if (hasTwoOrigins(selectedSpaceGroup))
        {
            ImGui::TextUnformatted("Origin choice");
            ImGui::SameLine();
            ImGui::RadioButton("1##origin", &m_originChoice, 1);
            ImGui::SameLine();
            ImGui::RadioButton("2##origin", &m_originChoice, 2);
            dialogLayout::hint("The International Tables give this group two origins; site coordinates differ between them\n"
                               "(Si in Fd-3m: 0 0 0 in choice 1, 1/8 1/8 1/8 in choice 2). Match the source of your coordinates.");
        }

        dialogLayout::section("Lattice parameters (\xc3\x85, \xc2\xb0)");
        drawLatticeParameterInputs((CrystalSystem)crystalSystemIndex, latticeParams);

        dialogLayout::section("Asymmetric unit");
        ImGui::TextDisabled("Sites in fractional coordinates; the space group generates the rest.");
        ImGui::SameLine();
        const float unitButtons = dialogLayout::buttonRowWidth({"Add Site", "Reset"});
        ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - unitButtons));
        if (dialogLayout::secondaryButton("Add Site"))
        {
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            scrollRowsToBottom = true;
            const AtomSite& atom = asymmetricAtoms.back();
            std::cout << "[Operation] Added asymmetric atom (bulk builder): "
                      << atom.symbol << "(" << atom.atomicNumber << ")"
                      << " at [" << atom.x << ", " << atom.y << ", " << atom.z << "]"
                      << std::endl;
        }
        dialogLayout::hint("Adds a hydrogen site at the origin; set its element and position in its row");
        ImGui::SameLine();
        if (dialogLayout::secondaryButton("Reset"))
        {
            asymmetricAtoms.clear();
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            std::cout << "[Operation] Cleared asymmetric atoms (bulk builder); reset to default H(1)" << std::endl;
        }
        dialogLayout::hint("Replaces every site with a single hydrogen site at the origin");

        int pendingDelete = -1;
        const float resultHeight = lastResult.message.empty() ? 0.0f : ImGui::GetFrameHeightWithSpacing() * 1.6f;
        const float footerHeight = ImGui::GetFrameHeight() + 3.0f * ImGui::GetStyle().ItemSpacing.y + responsive::dp(2) + 2.0f;
        const float listHeight = std::max(responsive::dp(110.0f), ImGui::GetContentRegionAvail().y - footerHeight - resultHeight);
        const ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersInnerH;
        if (ImGui::BeginTable("##bulk-atom-rows", 4, flags, ImVec2(0.0f, listHeight)))
        {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("0000").x);
            ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("Element").x + ImGui::GetStyle().FramePadding.x * 2);
            ImGui::TableSetupColumn("Position (u, v, w)", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("##delete", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("Delete").x + ImGui::GetStyle().FramePadding.x * 2);
            ImGui::TableHeadersRow();
            for (int i = 0; i < (int)asymmetricAtoms.size(); ++i)
            {
                AtomSite& atom = asymmetricAtoms[i];
                ImGui::PushID(i);
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                ImGui::AlignTextToFramePadding();
                ImGui::TextDisabled("%d", i);

                ImGui::TableSetColumnIndex(1);
                std::string elemButton = std::string(atom.symbol) + "##bulk-element";
                if (ImGui::Button(elemButton.c_str(), ImVec2(-FLT_MIN, 0)))
                {
                    targetAtomIndex = i;
                    selectedEditElement = atom.atomicNumber;
                    showElementPicker = true;
                    reopenBulkDialogAfterPicker = true;
                }
                dialogLayout::hint("Substitute this site's element");

                ImGui::TableSetColumnIndex(2);
                ImGui::SetNextItemWidth(-FLT_MIN);
                double editPos[3] = { atom.x, atom.y, atom.z };
                if (ImGui::DragScalarN("##bulk-direct", ImGuiDataType_Double, editPos, 3, 0.005f, nullptr, nullptr, "%.4f"))
                {
                    atom.x = editPos[0];
                    atom.y = editPos[1];
                    atom.z = editPos[2];
                }

                ImGui::TableSetColumnIndex(3);
                ImGui::BeginDisabled(asymmetricAtoms.size() == 1);
                if (dialogLayout::dangerButton("Delete"))
                    pendingDelete = i;
                ImGui::EndDisabled();

                ImGui::PopID();
            }

            if (scrollRowsToBottom)
            {
                ImGui::SetScrollY(ImGui::GetScrollMaxY());
                scrollRowsToBottom = false;
            }
            ImGui::EndTable();
        }

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

        if (!lastResult.message.empty())
        {
            ImGui::Spacing();
            if (lastResult.success)
                dialogLayout::notef(dialogLayout::Note::Success, "Built %d atoms in space group %d (%s).",
                                    lastResult.generatedAtoms, lastResult.spaceGroupNumber, lastResult.spaceGroupSymbol.c_str());
            else
                dialogLayout::note(dialogLayout::Note::Error, lastResult.message.c_str());
        }

        const int pressed = dialogLayout::footer(m_step.active() ? "Update Step" : "Build", "Close");
        if (pressed == 1 && m_step.active())
        {
            m_step.commitOptions(stepOptions());
            m_step.finish();
            dialogOpen = false;
            ImGui::CloseCurrentPopup();
        }
        else if (pressed == 1)
        {
            applySystemConstraints((CrystalSystem)crystalSystemIndex, latticeParams);
            lastResult = buildBulkCrystal(structure,
                                          (CrystalSystem)crystalSystemIndex,
                                          selectedSpaceGroup,
                                          latticeParams,
                                          asymmetricAtoms,
                                          elementColors,
                                          hasTwoOrigins(selectedSpaceGroup) ? m_originChoice : 1);
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
        else if (pressed == 2)
        {
            m_step.finish();
            dialogOpen = false;
            ImGui::CloseCurrentPopup();
        }
#endif

        ImGui::EndPopup();
    }
    else if (m_step.active() && !m_openRequested && !reopenBulkDialogAfterPicker)
    {
        // Closed with the X (or by another popup): the next menu build is normal.
        m_step.finish();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();

    if (showElementPicker)
    {
        ImGui::OpenPopup("Bulk Crystal Element Picker");
        showElementPicker = false;
    }

    responsive::windowSize(ImVec2(940.0f, 560.0f), ImGuiCond_Appearing);
    bool pickerOpen = true;
    if (responsive::beginModal("Bulk Crystal Element Picker", &pickerOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::TextDisabled("Click an element or type its symbol; double-click or Enter applies it.");
        bool activated = false;
        drawPeriodicTableInlineSelector(selectedEditElement, &activated, true);
        const bool enter = ImGui::IsKeyPressed(ImGuiKey_Enter, false) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter, false);
        const int pressed = dialogLayout::footer("Apply", "Cancel");
        if (pressed == 2)
            pickerOpen = false;
        if (pressed == 1 || activated || enter)
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
