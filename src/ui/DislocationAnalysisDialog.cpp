#include "ui/DislocationAnalysisDialog.h"
#include "ui/ResponsiveLayout.h"
#include "ui/DialogLayout.h"
#include "io/StructureLoader.h"
#include "imgui.h"

#include <algorithm>
#include <filesystem>
#include <fstream>

using namespace atomforge;

void DislocationAnalysisDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Dislocation Post-Processing")) m_openRequested = true;
}

void DislocationAnalysisDialog::loadReference(const std::string& path)
{
    m_referenceError.clear();
    m_referenceLoaded = loadStructureFromFile(path, m_reference, m_referenceError);
    m_referencePath = m_referenceLoaded ? path : std::string();
}

void DislocationAnalysisDialog::loadDeformed(const std::string& path)
{
    m_deformedError.clear();
    m_deformedLoaded = loadStructureFromFile(path, m_deformed, m_deformedError);
    m_deformedPath = m_deformedLoaded ? path : std::string();
}

void DislocationAnalysisDialog::drawDialog()
{
    if (m_openRequested) { m_open = true; m_openRequested = false; }

    if (m_nyeTask.poll() && m_nyeTask.result()) m_nyeResult = *m_nyeTask.result();
    if (m_vitekTask.poll() && m_vitekTask.result()) m_vitekResult = *m_vitekTask.result();
    if (m_patternBuildTask.poll() && m_patternBuildTask.result()) m_pattern = *m_patternBuildTask.result();
    if (m_patternDetectTask.poll() && m_patternDetectTask.result()) m_patternResult = *m_patternDetectTask.result();
    if (m_dragTask.poll() && m_dragTask.result()) m_dragResult = *m_dragTask.result();
    if (m_fitTask.poll() && m_fitTask.result()) m_fitResult = *m_fitTask.result();

    if (!m_open) return;

    responsive::windowSize(ImVec2(1180, 860), ImGuiCond_FirstUseEver);
    if (!responsive::begin("Dislocation Post-Processing", &m_open, ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }

    const bool stackPanels = responsive::stacked(1000);
    const float sidebar = stackPanels ? ImGui::GetContentRegionAvail().x
        : std::clamp(ImGui::GetContentRegionAvail().x * .40f, responsive::dp(360), responsive::dp(460));
    responsive::beginChild("Dislocation analysis controls", ImVec2(sidebar, stackPanels ? responsive::dp(420) : 0), true);

    dialogLayout::section("Structures");
    ImGui::TextWrapped("Both structures need the same atom count and ordering (e.g. a perfect crystal and the same crystal after inserting a dislocation).");
    if (responsive::button("Open reference (undeformed)...", ImVec2(-FLT_MIN, 0))) m_referencePicker.open("Open reference structure", false, m_referencePath);
    if (m_referenceLoaded)
        ImGui::TextWrapped("Reference: %s (%zu atoms)", std::filesystem::u8path(m_referencePath).filename().u8string().c_str(), m_reference.atoms.size());
    if (!m_referenceError.empty()) ImGui::TextWrapped("%s", m_referenceError.c_str());
    if (responsive::button("Open deformed...", ImVec2(-FLT_MIN, 0))) m_deformedPicker.open("Open deformed structure", false, m_deformedPath);
    if (m_deformedLoaded)
        ImGui::TextWrapped("Deformed: %s (%zu atoms)", std::filesystem::u8path(m_deformedPath).filename().u8string().c_str(), m_deformed.atoms.size());
    if (!m_deformedError.empty()) ImGui::TextWrapped("%s", m_deformedError.c_str());

    const bool haveBoth = m_referenceLoaded && m_deformedLoaded;
    ImGui::BeginDisabled(!haveBoth);

    dialogLayout::section("Nye tensor");
    ImGui::TextWrapped("Dislocation density tensor per atom (Hartley & Mishin lattice-correspondence method, matching BABEL's nyeTensor.f90).");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Cutoff (A)##nye", &m_nyeCutoff, 0.05f, 0.1f, 20.0f, "%.2f");
    ImGui::Checkbox("Ignore periodicity##nye", &m_nyeNoPbc);
    ImGui::BeginDisabled(m_nyeTask.running());
    if (dialogLayout::primaryButton("Compute Nye tensor", ImVec2(-FLT_MIN, 0)))
    {
        const Structure reference = m_reference, deformed = m_deformed;
        const float cutoff = m_nyeCutoff; const bool usePbc = !m_nyeNoPbc;
        m_nyeTask.start([=] { return computeNyeTensor(reference, deformed, cutoff, usePbc); });
    }
    ImGui::EndDisabled();
    if (m_nyeTask.running()) ImGui::TextUnformatted("Computing...");
    if (!m_nyeTask.error().empty()) ImGui::TextWrapped("%s", m_nyeTask.error().c_str());

    dialogLayout::section("Vitek differential-displacement map");
    ImGui::TextWrapped("Displacement difference between neighbor pairs, split into screw (along line) and edge components (matching BABEL's vitek.f90).");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::InputFloat3("Line direction##vitek", m_vitekLine);
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Burgers magnitude (A)##vitek", &m_vitekBurgers, 0.01f, 0.0f, 20.0f, "%.3f");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Cutoff (A, 0=auto)##vitek", &m_vitekCutoff, 0.05f, 0.0f, 20.0f, "%.2f");
    ImGui::Checkbox("Ignore periodicity##vitek", &m_vitekNoPbc);
    ImGui::BeginDisabled(m_vitekTask.running());
    if (dialogLayout::primaryButton("Compute Vitek map", ImVec2(-FLT_MIN, 0)))
    {
        const Structure reference = m_reference, deformed = m_deformed;
        const glm::dvec3 line(m_vitekLine[0], m_vitekLine[1], m_vitekLine[2]);
        const double burgers = m_vitekBurgers; const float cutoff = m_vitekCutoff; const bool usePbc = !m_vitekNoPbc;
        m_vitekTask.start([=] { return computeVitekMap(reference, deformed, line, burgers, cutoff, usePbc); });
    }
    ImGui::EndDisabled();
    if (m_vitekTask.running()) ImGui::TextUnformatted("Computing...");
    if (!m_vitekTask.error().empty()) ImGui::TextWrapped("%s", m_vitekTask.error().c_str());

    dialogLayout::section("Pattern matching");
    ImGui::TextWrapped("Flags atoms whose local neighbor environment no longer matches the reference lattice (matching BABEL's patternInit/patternDetect).");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Cutoff (A)##pattern", &m_patternCutoff, 0.05f, 0.1f, 20.0f, "%.2f");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Angle threshold (deg)##pattern", &m_patternAngleThreshold, 0.2f, 0.5f, 90.0f, "%.1f");
    ImGui::Checkbox("Ignore periodicity##pattern", &m_patternNoPbc);
    ImGui::BeginDisabled(m_patternBuildTask.running() || m_patternDetectTask.running());
    if (dialogLayout::primaryButton("Build pattern + detect", ImVec2(-FLT_MIN, 0)))
    {
        const Structure reference = m_reference, deformed = m_deformed;
        const float cutoff = m_patternCutoff; const bool usePbc = !m_patternNoPbc;
        const double angleThreshold = m_patternAngleThreshold;
        m_patternResult = PatternDetectResult{};
        m_patternBuildTask.start([=] { return buildPattern(reference, cutoff, usePbc); });
        // Detection is queued once the pattern build result is polled, below.
        (void)deformed; (void)angleThreshold;
    }
    ImGui::EndDisabled();
    if (m_patternBuildTask.running()) ImGui::TextUnformatted("Building pattern...");
    if (!m_patternBuildTask.error().empty()) ImGui::TextWrapped("%s", m_patternBuildTask.error().c_str());
    if (m_pattern.success && !m_patternDetectTask.running() && m_patternResult.rows.empty())
    {
        // Pattern just became available (or is available from a prior run):
        // chain the detection pass automatically.
        const Structure deformed = m_deformed;
        const float cutoff = m_patternCutoff; const bool usePbc = !m_patternNoPbc;
        const double angleThreshold = m_patternAngleThreshold; const CrystalPattern pattern = m_pattern;
        m_patternDetectTask.start([=] { return detectPattern(deformed, pattern, cutoff, angleThreshold, usePbc); });
    }
    if (m_patternDetectTask.running()) ImGui::TextUnformatted("Detecting...");
    if (!m_patternDetectTask.error().empty()) ImGui::TextWrapped("%s", m_patternDetectTask.error().c_str());

    dialogLayout::section("Prepare drag (migration barrier)");
    ImGui::TextWrapped("Interpolates reference -> deformed at reaction coordinate zeta for a constrained-minimization calculation (matching BABEL's prepareDrag). AtomForge does not itself run the constrained minimization.");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::SliderFloat("Zeta##drag", &m_dragZeta, 0.0f, 1.0f, "%.3f");
    ImGui::Checkbox("Clip displacement to minimum image##drag", &m_dragClipDisplacement);
    ImGui::Checkbox("Ignore periodicity##drag", &m_dragNoPbc);
    ImGui::BeginDisabled(m_dragTask.running());
    if (dialogLayout::primaryButton("Interpolate", ImVec2(-FLT_MIN, 0)))
    {
        const Structure reference = m_reference, deformed = m_deformed;
        const double zeta = m_dragZeta; const bool clip = m_dragClipDisplacement; const bool usePbc = !m_dragNoPbc;
        m_dragTask.start([=] { return prepareDrag(reference, deformed, zeta, clip, usePbc); });
    }
    ImGui::EndDisabled();
    if (m_dragTask.running()) ImGui::TextUnformatted("Interpolating...");
    if (!m_dragTask.error().empty()) ImGui::TextWrapped("%s", m_dragTask.error().c_str());
    if (m_dragResult.success)
    {
        if (responsive::button("Save interpolated structure...", ImVec2(-FLT_MIN, 0))) m_dragSavePicker.open("Save interpolated structure", true, "interpolated.vasp");
    }

    dialogLayout::section("Fit dislocation (position/Burgers vector)");
    ImGui::TextWrapped("Recovers a dislocation's position and Burgers vector from the Nye tensor field (weighted centroid/integral), the same purpose as BABEL's displacementFit/vitekFit.");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::InputFloat3("Line direction##fit", m_fitLine);
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Cutoff (A)##fit", &m_fitCutoff, 0.05f, 0.1f, 20.0f, "%.2f");
    ImGui::SetNextItemWidth(responsive::dp(-1)); ImGui::DragFloat("Area per atom (A^2)##fit", &m_fitArea, 0.05f, 0.001f, 1000.0f, "%.3f");
    ImGui::Checkbox("Ignore periodicity##fit", &m_fitNoPbc);
    ImGui::BeginDisabled(m_fitTask.running());
    if (dialogLayout::primaryButton("Fit dislocation", ImVec2(-FLT_MIN, 0)))
    {
        const Structure reference = m_reference, deformed = m_deformed;
        const glm::dvec3 line(m_fitLine[0], m_fitLine[1], m_fitLine[2]);
        const float cutoff = m_fitCutoff; const double area = m_fitArea; const bool usePbc = !m_fitNoPbc;
        m_fitTask.start([=]() -> DislocationFitResult
        {
            const NyeTensorResult nye = computeNyeTensor(reference, deformed, cutoff, usePbc);
            if (!nye.success)
            {
                DislocationFitResult failed;
                failed.message = "Nye tensor computation failed: " + nye.message;
                return failed;
            }
            std::vector<glm::dvec3> positions(reference.atoms.size());
            for (std::size_t i = 0; i < positions.size(); ++i)
                positions[i] = glm::dvec3(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);
            return fitDislocationFromNye(positions, nye, line, area);
        });
    }
    ImGui::EndDisabled();
    if (m_fitTask.running()) ImGui::TextUnformatted("Fitting (computes the Nye tensor internally)...");
    if (!m_fitTask.error().empty()) ImGui::TextWrapped("%s", m_fitTask.error().c_str());

    ImGui::EndDisabled(); // haveBoth
    if (!haveBoth) ImGui::TextWrapped("Load both a reference and a deformed structure to enable these analyses.");
    ImGui::EndChild();

    responsive::nextPanel(stackPanels);
    responsive::beginChild("Dislocation analysis results", responsive::size(0, stackPanels ? 520 : 0), false);

    if (m_nyeResult.success)
    {
        dialogLayout::section("Nye tensor result");
        ImGui::TextWrapped("%s", m_nyeResult.message.c_str());
        std::vector<std::size_t> order(m_nyeResult.alphaNorm.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::partial_sort(order.begin(), order.begin() + std::min<std::size_t>(10, order.size()), order.end(),
            [&](std::size_t a, std::size_t b) { return m_nyeResult.alphaNorm[a] > m_nyeResult.alphaNorm[b]; });
        ImGui::TextUnformatted("Top atoms by |alpha| (dislocation content):");
        for (std::size_t k = 0; k < std::min<std::size_t>(10, order.size()); ++k)
            ImGui::Text("  atom %zu: |alpha| = %.4g", order[k], m_nyeResult.alphaNorm[order[k]]);
        if (responsive::button("Export full Nye tensor CSV...", ImVec2(-FLT_MIN, 0))) m_nyeExportPicker.open("Export Nye tensor", true, "nye_tensor.csv");
    }

    if (m_vitekResult.success)
    {
        dialogLayout::section("Vitek map result");
        ImGui::TextWrapped("%s", m_vitekResult.message.c_str());
        if (responsive::button("Export full Vitek map CSV...", ImVec2(-FLT_MIN, 0))) m_vitekExportPicker.open("Export Vitek map", true, "vitek_map.csv");
    }

    if (m_patternResult.success)
    {
        dialogLayout::section("Pattern matching result");
        ImGui::TextWrapped("Pattern: %s", m_pattern.message.c_str());
        ImGui::TextWrapped("Detection: %s", m_patternResult.message.c_str());
        if (responsive::button("Export full pattern-match CSV...", ImVec2(-FLT_MIN, 0))) m_patternExportPicker.open("Export pattern match", true, "pattern_match.csv");
    }

    if (m_dragResult.success)
    {
        dialogLayout::section("Drag preparation result");
        ImGui::TextWrapped("%s", m_dragResult.message.c_str());
        ImGui::TextWrapped("Interpolated structure ready (%zu atoms); use \"Save interpolated structure...\" above.", m_dragResult.interpolated.atoms.size());
    }

    if (m_fitResult.success)
    {
        dialogLayout::section("Dislocation fit result");
        ImGui::Text("Line position: (%.4f, %.4f, %.4f) A", m_fitResult.linePosition.x, m_fitResult.linePosition.y, m_fitResult.linePosition.z);
        ImGui::Text("Burgers vector (uncalibrated by area per atom): (%.4f, %.4f, %.4f)", m_fitResult.burgersVector.x, m_fitResult.burgersVector.y, m_fitResult.burgersVector.z);
        const double mag = glm::length(m_fitResult.burgersVector);
        ImGui::Text("Magnitude: %.4f", mag);
    }
    else if (!m_fitResult.message.empty())
    {
        ImGui::TextWrapped("%s", m_fitResult.message.c_str());
    }

    ImGui::EndChild();

    if (auto path = m_referencePicker.draw()) loadReference(*path);
    if (auto path = m_deformedPicker.draw()) loadDeformed(*path);

    if (auto path = m_nyeExportPicker.draw())
    {
        std::ofstream file(std::filesystem::u8path(*path));
        file << "index,symbol,alpha_xx,alpha_xy,alpha_xz,alpha_yx,alpha_yy,alpha_yz,alpha_zx,alpha_zy,alpha_zz,norm\n";
        for (std::size_t i = 0; i < m_nyeResult.alpha.size(); ++i)
        {
            file << i << ',' << (i < m_deformed.atoms.size() ? m_deformed.atoms[i].symbol : std::string());
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) file << ',' << m_nyeResult.alpha[i][a][b];
            file << ',' << m_nyeResult.alphaNorm[i] << '\n';
        }
    }
    if (auto path = m_vitekExportPicker.draw())
    {
        std::ofstream file(std::filesystem::u8path(*path));
        file << "index_i,index_j,xi_A,yi_A,zi_A,xj_A,yj_A,zj_A,screw_A,edge_1_A,edge_2_A\n";
        for (const auto& p : m_vitekResult.pairs)
            file << p.indexI << ',' << p.indexJ << ','
                 << p.positionI.x << ',' << p.positionI.y << ',' << p.positionI.z << ','
                 << p.positionJ.x << ',' << p.positionJ.y << ',' << p.positionJ.z << ','
                 << p.screwComponent << ',' << p.edgeComponent.x << ',' << p.edgeComponent.y << '\n';
    }
    if (auto path = m_patternExportPicker.draw())
    {
        std::ofstream file(std::filesystem::u8path(*path));
        file << "index,matched,neighbor_count,max_angle_deviation_deg\n";
        for (const auto& row : m_patternResult.rows)
            file << row.index << ',' << (row.matched ? 1 : 0) << ',' << row.neighborCount << ',' << row.maxAngleDeviationDeg << '\n';
    }
    if (auto path = m_dragSavePicker.draw())
    {
        const std::string ext = std::filesystem::u8path(*path).extension().u8string();
        const std::string format = ext.size() > 1 ? ext.substr(1) : "vasp";
        saveStructure(m_dragResult.interpolated, *path, format);
    }

    ImGui::End();
}
