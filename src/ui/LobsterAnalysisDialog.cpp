#include "ui/LobsterAnalysisDialog.h"
#include "ui/ResponsiveLayout.h"
#include "ui/DialogLayout.h"
#include "ui/LinePlot.h"
#include "imgui.h"

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <filesystem>
#include <stdexcept>

using namespace atomforge::electronic;

void LobsterAnalysisDialog::drawMenuItem()
{
    if (ImGui::MenuItem("LOBSTER Bonding Analysis")) m_openRequested = true;
}

void LobsterAnalysisDialog::loadCohpcar(const std::string& path)
{
    m_cohpcarError.clear();
    try { m_cohpcar = readCohpcar(path); m_cohpcarPath = path; m_selectedBond = -1; m_spinChannel = 0; }
    catch (const std::exception& e) { m_cohpcarError = e.what(); m_cohpcar = {}; m_cohpcarPath.clear(); }
}

void LobsterAnalysisDialog::loadIcohplist(const std::string& path)
{
    m_icohplistError.clear();
    try { m_icohplist = readIcohplist(path); m_icohplistPath = path; }
    catch (const std::exception& e) { m_icohplistError = e.what(); m_icohplist = {}; m_icohplistPath.clear(); }
}

void LobsterAnalysisDialog::drawDialog()
{
    if (m_openRequested) { m_open = true; m_openRequested = false; }
    if (!m_open) return;

    responsive::windowSize(ImVec2(1180,820),ImGuiCond_FirstUseEver);
    if (!responsive::begin("LOBSTER Bonding Analysis",&m_open,ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    const bool stackPanels = responsive::stacked(1000);
    const float sidebar = stackPanels ? ImGui::GetContentRegionAvail().x : std::clamp(ImGui::GetContentRegionAvail().x*.38f,responsive::dp(340),responsive::dp(440));
    responsive::beginChild("Lobster controls",ImVec2(sidebar,stackPanels?responsive::dp(420):0),true);

    dialogLayout::section("COHP / COOP / COBI curves");
    if (responsive::button("Open COHPCAR...")) m_cohpcarPicker.open("Open LOBSTER COHPCAR",false,m_cohpcarPath.empty()?"COHPCAR.lobster":m_cohpcarPath);
    if (!m_cohpcarPath.empty())
        ImGui::TextWrapped("Loaded: %s (%zu bonds, %s)",std::filesystem::u8path(m_cohpcarPath).filename().u8string().c_str(),
            m_cohpcar.bonds.size(),m_cohpcar.spinPolarized?"spin-polarized":"non-spin-polarized");
    if (!m_cohpcarError.empty()) ImGui::TextWrapped("%s",m_cohpcarError.c_str());
    if (!m_cohpcar.energies.empty())
    {
        ImGui::Text("Fermi energy: %.4f eV",m_cohpcar.fermiEnergy);
        if (m_cohpcar.spinPolarized)
        {
            ImGui::SetNextItemWidth(responsive::dp(-1));
            ImGui::Combo("Spin channel",&m_spinChannel,"Up\0Down\0");
        }
        m_spinChannel = std::clamp(m_spinChannel,0,m_cohpcar.spinPolarized?1:0);
        responsive::beginChild("Bond list",responsive::size(0,160),true);
        if (ImGui::Selectable("Average (all bonds)",m_selectedBond==-1)) m_selectedBond=-1;
        for (int i=0;i<static_cast<int>(m_cohpcar.bonds.size());++i)
        {
            const auto& bond = m_cohpcar.bonds[static_cast<std::size_t>(i)];
            char label[128];
            std::snprintf(label,sizeof(label),"Bond %d: atom %d - atom %d (%.3f A)",i,bond.atom1,bond.atom2,bond.length);
            ImGui::PushID(i);
            if (ImGui::Selectable(label,m_selectedBond==i)) m_selectedBond=i;
            ImGui::PopID();
        }
        ImGui::EndChild();
    }

    dialogLayout::section("ICOHPLIST / ICOOPLIST / ICOBILIST summary");
    if (responsive::button("Open ICOHPLIST...")) m_icohplistPicker.open("Open LOBSTER ICOHPLIST",false,m_icohplistPath.empty()?"ICOHPLIST.lobster":m_icohplistPath);
    if (!m_icohplistPath.empty())
        ImGui::TextWrapped("Loaded: %s (%zu entries, %s)",std::filesystem::u8path(m_icohplistPath).filename().u8string().c_str(),
            m_icohplist.entries.size(),m_icohplist.spinPolarized?"spin-polarized":"non-spin-polarized");
    if (!m_icohplistError.empty()) ImGui::TextWrapped("%s",m_icohplistError.c_str());
    ImGui::EndChild();

    responsive::nextPanel(stackPanels);
    responsive::beginChild("Lobster results",responsive::size(0,stackPanels?520:0),false);
    if (!m_cohpcar.energies.empty())
    {
        const std::size_t spin = static_cast<std::size_t>(m_spinChannel);
        const auto& cohp = m_selectedBond==-1 ? m_cohpcar.averageCohp.at(spin) : m_cohpcar.cohp.at(spin).at(static_cast<std::size_t>(m_selectedBond));
        const auto& icohp = m_selectedBond==-1 ? m_cohpcar.averageIcohp.at(spin) : m_cohpcar.icohp.at(spin).at(static_cast<std::size_t>(m_selectedBond));
        dialogLayout::section(m_selectedBond==-1 ? "Average pCOHP / IpCOHP" : "Selected-bond pCOHP / IpCOHP");
        float low=FLT_MAX, high=-FLT_MAX;
        for (double v : cohp) { low=std::min(low,static_cast<float>(v)); high=std::max(high,static_cast<float>(v)); }
        for (double v : icohp) { low=std::min(low,static_cast<float>(v)); high=std::max(high,static_cast<float>(v)); }
        ImGui::TextWrapped("Series 1: pCOHP. Series 2: IpCOHP (running integral). Fermi level at 0 eV, range %.3f to %.3f eV.",
            m_cohpcar.energies.front(),m_cohpcar.energies.back());
        uiPlot::drawColoredLines("cohpplot",2,static_cast<int>(m_cohpcar.energies.size()),
            [&](int series,int sample){ return series==0 ? cohp[static_cast<std::size_t>(sample)] : icohp[static_cast<std::size_t>(sample)]; },
            low,high,"Energy - Fermi (eV)");
    }
    if (!m_icohplist.entries.empty())
    {
        dialogLayout::section("ICOHPLIST entries");
        responsive::beginChild("icohptable",responsive::size(0,220),true);
        ImGui::TextUnformatted(m_icohplist.spinPolarized ? "atom1  atom2  length_A  numBonds  icohp_up  icohp_down" : "atom1  atom2  length_A  numBonds  icohp");
        for (const auto& entry : m_icohplist.entries)
        {
            if (m_icohplist.spinPolarized)
                ImGui::Text("%5d  %5d  %8.4f  %8d  %8.5f  %8.5f",entry.atom1,entry.atom2,entry.length,entry.numBonds,entry.icohp[0],entry.icohp[1]);
            else
                ImGui::Text("%5d  %5d  %8.4f  %8d  %8.5f",entry.atom1,entry.atom2,entry.length,entry.numBonds,entry.icohp[0]);
        }
        ImGui::EndChild();
    }
    ImGui::EndChild();

    if (auto path=m_cohpcarPicker.draw()) loadCohpcar(*path);
    if (auto path=m_icohplistPicker.draw()) loadIcohplist(*path);
    ImGui::End();
}
