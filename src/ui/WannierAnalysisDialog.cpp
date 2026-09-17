#include "ui/WannierAnalysisDialog.h"
#include "ui/ResponsiveLayout.h"
#include "ui/DialogLayout.h"
#include "ui/LinePlot.h"
#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <glm/glm.hpp>
#include <sstream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
// One k-point (fractional reciprocal) per line; a path needs at least two,
// interpolated linearly with `samplesPerSegment` points between each pair
// (plus the final anchor), matching a conventional band-structure k-path.
std::vector<glm::dvec3> interpolatePath(const std::string& text, int samplesPerSegment)
{
    std::vector<glm::dvec3> anchors;
    std::istringstream input(text);
    for (std::string line; std::getline(input,line);)
    {
        std::istringstream row(line);
        glm::dvec3 point;
        if (row >> point.x >> point.y >> point.z) anchors.push_back(point);
    }
    if (anchors.size() < 2) throw std::invalid_argument("Supply at least two k-points (one per line, fractional reciprocal)");
    std::vector<glm::dvec3> path;
    for (std::size_t segment = 0; segment + 1 < anchors.size(); ++segment)
        for (int s = 0; s < samplesPerSegment; ++s)
        {
            const double t = static_cast<double>(s) / samplesPerSegment;
            path.push_back(anchors[segment] * (1 - t) + anchors[segment + 1] * t);
        }
    path.push_back(anchors.back());
    return path;
}

}

void WannierAnalysisDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Wannier Analysis")) m_openRequested = true;
}

void WannierAnalysisDialog::load(const std::string& path)
{
    m_loadError.clear();
    try { m_model = readWannierHr(path); m_loadedPath = path; }
    catch (const std::exception& e) { m_loadError = e.what(); m_model = {}; m_loadedPath.clear(); }
}

void WannierAnalysisDialog::drawDialog()
{
    if (m_openRequested) { m_open = true; m_openRequested = false; }

    if (m_bandsTask.poll()) { if (m_bandsTask.result()) m_bands = *m_bandsTask.result(); }
    if (m_curvatureTask.poll()) { if (m_curvatureTask.result()) m_curvature = *m_curvatureTask.result(); }
    if (m_chernTask.poll()) { if (m_chernTask.result()) { m_chernNumber = *m_chernTask.result(); m_chernComputed = true; } }
    if (!m_open) return;

    responsive::windowSize(ImVec2(1180,820),ImGuiCond_FirstUseEver);
    if (!responsive::begin("Wannier Analysis",&m_open,ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    const bool stackPanels = responsive::stacked(1000);
    const float sidebar = stackPanels ? ImGui::GetContentRegionAvail().x : std::clamp(ImGui::GetContentRegionAvail().x*.38f,responsive::dp(340),responsive::dp(440));
    responsive::beginChild("Wannier controls",ImVec2(sidebar,stackPanels?responsive::dp(360):0),true);

    if (responsive::button("Open seedname_hr.dat...")) m_picker.open("Open Wannier90 Hamiltonian",false,m_loadedPath.empty()?"wannier90_hr.dat":m_loadedPath);
    if (!m_loadedPath.empty()) ImGui::TextWrapped("Loaded: %s (%d Wannier functions)",std::filesystem::u8path(m_loadedPath).filename().u8string().c_str(),m_model.numWann);
    if (!m_loadError.empty()) ImGui::TextWrapped("%s",m_loadError.c_str());
    ImGui::TextWrapped("H(k) = sum_R exp(2*pi*i*k.R) H(R) / ndegen(R), k fractional reciprocal.");

    ImGui::BeginDisabled(m_model.numWann==0);
    dialogLayout::section("Bands");
    ImGui::TextWrapped("One fractional k-point (kx ky kz) per line; at least two, forming a piecewise path.");
    ImGui::InputTextMultiline("##kpath",m_kpath,sizeof(m_kpath),ImVec2(-1,90));
    ImGui::SetNextItemWidth(responsive::dp(-1));
    ImGui::InputInt("Samples per segment",&m_samplesPerSegment);
    ImGui::BeginDisabled(m_bandsTask.running());
    if (dialogLayout::primaryButton("Compute bands",ImVec2(-FLT_MIN,0)))
    {
        const std::string text(m_kpath); const int samples=std::clamp(m_samplesPerSegment,1,2000); const auto model=m_model;
        m_bandsTask.start([=]
        {
            const auto path = interpolatePath(text,samples);
            const auto energies = wannierBands(model,path);
            return BandsResult{energies,model.numWann,static_cast<int>(path.size())};
        });
    }
    ImGui::EndDisabled();
    if (m_bandsTask.running()) ImGui::TextUnformatted("Computing...");
    if (!m_bandsTask.error().empty()) ImGui::TextWrapped("%s",m_bandsTask.error().c_str());

    dialogLayout::section("Berry curvature (Hamiltonian-gauge/Kubo term)");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::InputFloat3("k-point",m_curvatureK);
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::Combo("Plane",&m_curvaturePlane,"kx-ky\0kx-kz\0ky-kz\0");
    ImGui::BeginDisabled(m_curvatureTask.running());
    if (dialogLayout::primaryButton("Compute curvature",ImVec2(-FLT_MIN,0)))
    {
        const glm::dvec3 k(m_curvatureK[0],m_curvatureK[1],m_curvatureK[2]);
        const int planeA=std::array<int,3>{0,0,1}[m_curvaturePlane], planeB=std::array<int,3>{1,2,2}[m_curvaturePlane];
        const auto model=m_model;
        m_curvatureTask.start([=] { return berryCurvature(model,k,planeA,planeB); });
    }
    ImGui::EndDisabled();
    if (!m_curvatureTask.error().empty()) ImGui::TextWrapped("%s",m_curvatureTask.error().c_str());
    ImGui::TextWrapped("Diverges at exact band degeneracies; use Chern number below for those k-points instead.");

    dialogLayout::section("Chern number (Fukui-Hatsugai-Suzuki)");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::InputInt("Band (0-based)",&m_chernBand);
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::Combo("Plane##chern",&m_chernPlane,"kx-ky\0kx-kz\0ky-kz\0");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::SliderInt("Mesh resolution",&m_chernGrid,4,120);
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::InputFloat("Fixed fractional coordinate",&m_chernFixed);
    ImGui::BeginDisabled(m_chernTask.running());
    if (dialogLayout::primaryButton("Compute Chern number",ImVec2(-FLT_MIN,0)))
    {
        const int planeA=std::array<int,3>{0,0,1}[m_chernPlane], planeB=std::array<int,3>{1,2,2}[m_chernPlane];
        const int band=m_chernBand, grid=std::clamp(m_chernGrid,2,2000); const double fixed=m_chernFixed; const auto model=m_model;
        m_chernComputed=false;
        m_chernTask.start([=] { return chernNumber(model,band,planeA,planeB,grid,fixed); });
    }
    ImGui::EndDisabled();
    if (!m_chernTask.error().empty()) ImGui::TextWrapped("%s",m_chernTask.error().c_str());
    ImGui::EndDisabled();
    ImGui::EndChild();

    responsive::nextPanel(stackPanels);
    responsive::beginChild("Wannier results",responsive::size(0,stackPanels?520:0),false);
    if (m_bands.numWann>0)
    {
        ImGui::TextWrapped("%d bands, %d k-points along the path (index axis; no reciprocal-lattice distance is stored in seedname_hr.dat).",m_bands.numWann,m_bands.numK);
        float low=FLT_MAX, high=-FLT_MAX;
        for (double v : m_bands.energies) { low=std::min(low,static_cast<float>(v)); high=std::max(high,static_cast<float>(v)); }
        uiPlot::drawColoredLines("bandplot",m_bands.numWann,m_bands.numK,
            [&](int series,int sample){ return m_bands.energies[static_cast<std::size_t>(sample)*m_bands.numWann+series]; },
            low,high,"k-point index");
    }
    if (!m_curvature.empty())
    {
        dialogLayout::section("Berry curvature per band");
        for (std::size_t i=0;i<m_curvature.size();++i) ImGui::Text("Band %zu: %.6g",i,m_curvature[i]);
    }
    if (m_chernComputed)
    {
        dialogLayout::section("Chern number result");
        ImGui::Text("C = %.6f",m_chernNumber);
        if (std::abs(m_chernNumber-std::round(m_chernNumber))>0.05)
            ImGui::TextWrapped("Far from an integer: refine the mesh resolution or check the band/plane for a degeneracy.");
    }
    ImGui::EndChild();

    if (auto path=m_picker.draw()) load(*path);
    ImGui::End();
}
