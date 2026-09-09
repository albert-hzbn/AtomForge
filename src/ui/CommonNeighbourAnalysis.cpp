#include "ui/CommonNeighbourAnalysis.h"

#include "ui/ThemeUtils.h"
#include "imgui.h"


#include <cstdio>
#include <string>

using atomforge::analysis::CnaResult;
using atomforge::analysis::Signature;

namespace
{
std::string signatureToString(const Signature& s)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "1-%d-%d-%d", s.common, s.bonds, s.chain);
    return std::string(buffer);
}

void drawSignatureTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-signatures", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Sortable))
    {
        ImGui::TableSetupColumn("Signature", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Pair Count", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Fraction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableHeadersRow();

        for (const auto& [signature, count] : result.signatureCounts)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", signatureToString(signature).c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", count);
            ImGui::TableSetColumnIndex(2);
            float f = (result.pairCount > 0) ? (float)count / (float)result.pairCount : 0.0f;
            ImGui::Text("%.4f", f);
        }

        ImGui::EndTable();
    }
}

void drawEnvironmentTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-environments", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
    {
        ImGui::TableSetupColumn("Environment", ImGuiTableColumnFlags_WidthFixed, 160.0f);
        ImGui::TableSetupColumn("Atom Count", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Fraction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableHeadersRow();

        for (const auto& [environment, count] : result.environmentCounts)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", environment.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", count);
            ImGui::TableSetColumnIndex(2);
            float f = (result.atomCount > 0) ? (float)count / (float)result.atomCount : 0.0f;
            ImGui::Text("%.4f", f);
        }

        ImGui::EndTable();
    }
}

void drawAtomTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-per-atom", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY))
    {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 52.0f);
        ImGui::TableSetupColumn("El", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Z", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("CN", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Dominant Signature", ImGuiTableColumnFlags_WidthFixed, 160.0f);
        ImGui::TableSetupColumn("Sig Count", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Environment", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (const auto& row : result.atomRows)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", row.index);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s", row.symbol.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", row.atomicNumber);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", row.coordination);
            ImGui::TableSetColumnIndex(4);
            if (row.dominantSignatureCount > 0)
                ImGui::Text("%s", signatureToString(row.dominantSignature).c_str());
            else
                ImGui::Text("N/A");
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%d", row.dominantSignatureCount);
            ImGui::TableSetColumnIndex(6);
            ImGui::Text("%s", row.environment.c_str());
        }

        ImGui::EndTable();
    }
}

void drawCnaSummary(const CnaResult& result)
{
    ImGui::Separator();
    ImGui::Text("Status: %s", result.message.c_str());
    ImGui::Text("Atoms: %d", result.atomCount);
    ImGui::Text("Bonded pairs analyzed: %d", result.pairCount);
    ImGui::Text("PBC used in analysis: %s", result.pbcUsed ? "Yes" : "No");
}

void drawCnaDetails(const CnaResult& result)
{
    if (!result.valid)
        return;

    ImGui::Separator();
    ImGui::Text("Pair Signature Distribution (Honeycutt-Andersen form: 1-j-k-l)");
    ImGui::BeginChild("##cna-sig-child", ImVec2(0.0f, 180.0f), true);
    drawSignatureTable(result);
    ImGui::EndChild();

    ImGui::Separator();
    ImGui::Text("Per-Atom Environment Summary");
    ImGui::BeginChild("##cna-env-child", ImVec2(0.0f, 130.0f), true);
    drawEnvironmentTable(result);
    ImGui::EndChild();

    ImGui::Separator();
    ImGui::Text("Per-Atom CNA Details");
    ImGui::BeginChild("##cna-atom-child", ImVec2(0.0f, 240.0f), true);
    drawAtomTable(result);
    ImGui::EndChild();
}

} // namespace

void CommonNeighbourAnalysisDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Common Neighbour Analysis", nullptr, false, enabled))
        m_openRequested = true;
}

void CommonNeighbourAnalysisDialog::drawDialog(const Structure& structure)
{
    if (m_task.poll())
        m_computeCompleted = m_task.result().has_value();

    bool& usePbc = m_params.usePbc;
    float& cutoffScale = m_params.cutoffScale;

    if (m_openRequested)
    {
        ImGui::OpenPopup("Common Neighbour Analysis");
        m_openRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(1100.0f, 760.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("Common Neighbour Analysis", &dialogOpen, ImGuiWindowFlags_NoResize))
    {
        bool changed = false;
        changed |= ImGui::Checkbox("Use PBC when unit cell is available", &usePbc);
        changed |= ImGui::SliderFloat("Bond cutoff scale", &cutoffScale, 1.00f, 1.60f, "%.2f");
        ImGui::SameLine();
        
        bool computeRequested = false;
        if (ImGui::Button("Run CNA") && !m_task.running())
            computeRequested = true;
        
        ImGui::SameLine();
        if (m_task.running())
            ImGui::TextColored(themeStatusComputing(), "Computing...");
        else if (m_computeCompleted)
            ImGui::TextColored(themeStatusGood(), "Done");
        if (!m_task.error().empty())
            ImGui::TextWrapped("%s", m_task.error().c_str());

        // Launch computation in background thread
        if (computeRequested && !m_task.running())
        {
            m_computeCompleted = false;
            m_task.start([snapshot = structure, params = m_params]() {
                return atomforge::analysis::computeCna(snapshot, params);
            });
        }

        if (changed && !m_task.running())
        {
            // Clear results if parameters change while not computing
            m_task.clearResult();
            m_computeCompleted = false;
        }

        // Display results if available
        if (m_task.result())
        {
            const auto* resultPtr = &*m_task.result();
            drawCnaSummary(*resultPtr);
            drawCnaDetails(*resultPtr);
        }

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();
}
