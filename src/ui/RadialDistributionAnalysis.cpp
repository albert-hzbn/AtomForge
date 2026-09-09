#include "ui/RadialDistributionAnalysis.h"

#include "ElementData.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"


#include <algorithm>
#include <cstdio>
#include <set>
#include <string>
#include <vector>
#include <utility>

using atomforge::analysis::RdfResult;

namespace
{
void drawPlot(const RdfResult& result,
              bool showRawCounts,
              bool showCumulative)
{
    const ImVec2 canvasSize(-1.0f, 260.0f);
    ImGui::BeginChild("##rdf-plot-child", canvasSize, true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (avail.x < 50.0f) avail.x = 50.0f;
    if (avail.y < 50.0f) avail.y = 50.0f;
    ImGui::InvisibleButton("##rdf-canvas", avail);

    ImVec2 pMin = origin;
    ImVec2 pMax(origin.x + avail.x, origin.y + avail.y);
    const bool lightRdf = isLightTheme();
    drawList->AddRectFilled(pMin, pMax,
                            lightRdf ? IM_COL32(245, 246, 248, 255) : IM_COL32(20, 24, 31, 255));
    drawList->AddRect(pMin, pMax,
                      lightRdf ? IM_COL32(140, 145, 160, 255) : IM_COL32(110, 125, 150, 255));

    const ImU32 noDataCol = lightRdf ? IM_COL32(80, 80, 100, 255)  : IM_COL32(220, 220, 220, 255);
    const ImU32 gridCol   = lightRdf ? IM_COL32(170, 175, 190, 180) : IM_COL32(50, 60, 75, 255);
    const ImU32 labelCol  = lightRdf ? IM_COL32(50, 50, 65, 255)   : IM_COL32(200, 200, 200, 255);
    const ImU32 rdfLineCol   = lightRdf ? IM_COL32(20,  120, 210, 255) : IM_COL32(80,  220, 255, 255);
    const ImU32 rawLineCol   = lightRdf ? IM_COL32(200, 110,  20, 220) : IM_COL32(255, 180,  80, 220);
    const ImU32 cumLineCol   = lightRdf ? IM_COL32( 30, 150,  40, 220) : IM_COL32(140, 255, 140, 220);
    const ImU32 rdfLabelCol  = lightRdf ? IM_COL32(20,  120, 210, 255) : IM_COL32(80,  220, 255, 255);
    const ImU32 rawLabelCol  = lightRdf ? IM_COL32(200, 110,  20, 255) : IM_COL32(255, 180,  80, 255);
    const ImU32 cumLabelCol  = lightRdf ? IM_COL32( 30, 150,  40, 255) : IM_COL32(140, 255, 140, 255);

    if (!result.valid || result.bins.empty())
    {
        drawList->AddText(ImVec2(pMin.x + 10.0f, pMin.y + 10.0f), noDataCol, "No RDF data.");
        ImGui::EndChild();
        return;
    }

    const float leftPad = 52.0f;
    const float rightPad = 16.0f;
    const float topPad = 16.0f;
    const float bottomPad = 28.0f;
    ImVec2 gMin(pMin.x + leftPad, pMin.y + topPad);
    ImVec2 gMax(pMax.x - rightPad, pMax.y - bottomPad);

    drawList->AddRect(gMin, gMax, lightRdf ? IM_COL32(150, 155, 170, 255) : IM_COL32(80, 90, 110, 255));

    float maxY = 0.0f;
    for (int i = 0; i < (int)result.bins.size(); ++i)
    {
        maxY = std::max(maxY, result.bins[i].g);
        if (showRawCounts)
            maxY = std::max(maxY, result.bins[i].rawCount);
        if (showCumulative)
            maxY = std::max(maxY, result.bins[i].cumulative);
    }
    maxY = std::max(maxY, 1.0f);

    for (int grid = 0; grid <= 4; ++grid)
    {
        float t = (float)grid / 4.0f;
        float y = gMax.y + (gMin.y - gMax.y) * t;
        drawList->AddLine(ImVec2(gMin.x, y), ImVec2(gMax.x, y), gridCol);
        char label[32];
        std::snprintf(label, sizeof(label), "%.2f", maxY * t);
        drawList->AddText(ImVec2(pMin.x + 4.0f, y - 7.0f), labelCol, label);
    }

    for (int grid = 0; grid <= 4; ++grid)
    {
        float t = (float)grid / 4.0f;
        float x = gMin.x + (gMax.x - gMin.x) * t;
        drawList->AddLine(ImVec2(x, gMin.y), ImVec2(x, gMax.y), gridCol);
        float r = result.rMin + (result.rMax - result.rMin) * t;
        char label[32];
        std::snprintf(label, sizeof(label), "%.2f", r);
        drawList->AddText(ImVec2(x - 10.0f, gMax.y + 6.0f), labelCol, label);
    }

    std::vector<ImVec2> rdfPoints;
    std::vector<ImVec2> rawPoints;
    std::vector<ImVec2> cumulativePoints;
    rdfPoints.reserve(result.bins.size());
    rawPoints.reserve(result.bins.size());
    cumulativePoints.reserve(result.bins.size());

    for (int i = 0; i < (int)result.bins.size(); ++i)
    {
        float tx = (result.bins[i].rCenter - result.rMin) / std::max(result.rMax - result.rMin, 1e-6f);
        float x = gMin.x + (gMax.x - gMin.x) * tx;

        float yRdf = gMax.y - (gMax.y - gMin.y) * (result.bins[i].g / maxY);
        rdfPoints.push_back(ImVec2(x, yRdf));

        float yRaw = gMax.y - (gMax.y - gMin.y) * (result.bins[i].rawCount / maxY);
        rawPoints.push_back(ImVec2(x, yRaw));

        float yCum = gMax.y - (gMax.y - gMin.y) * (result.bins[i].cumulative / maxY);
        cumulativePoints.push_back(ImVec2(x, yCum));
    }

    if (rdfPoints.size() >= 2)
        drawList->AddPolyline(&rdfPoints[0], (int)rdfPoints.size(), rdfLineCol, 0, 2.0f);
    if (showRawCounts && rawPoints.size() >= 2)
        drawList->AddPolyline(&rawPoints[0], (int)rawPoints.size(), rawLineCol, 0, 2.0f);
    if (showCumulative && cumulativePoints.size() >= 2)
        drawList->AddPolyline(&cumulativePoints[0], (int)cumulativePoints.size(), cumLineCol, 0, 2.0f);

    drawList->AddText(ImVec2(gMin.x + 8.0f, gMin.y + 6.0f), rdfLabelCol, result.normalized ? "g(r)" : "Histogram");
    if (showRawCounts)
        drawList->AddText(ImVec2(gMin.x + 64.0f, gMin.y + 6.0f), rawLabelCol, "Counts");
    if (showCumulative)
        drawList->AddText(ImVec2(gMin.x + 130.0f, gMin.y + 6.0f), cumLabelCol, "Cumulative CN");

    ImGui::EndChild();
}

const char* speciesLabelGetter(void* userData, int idx)
{
    const std::vector<std::pair<int, std::string>>* entries = static_cast<const std::vector<std::pair<int, std::string>>*>(userData);
    return (*entries)[idx].second.c_str();
}

std::vector<std::pair<int, std::string>> buildSpeciesOptions(const Structure& structure)
{
    std::vector<std::pair<int, std::string>> options;
    options.push_back(std::make_pair(0, std::string("All elements")));

    std::set<int> seen;
    for (const auto& atom : structure.atoms)
    {
        const int z = atom.atomicNumber;
        if (!seen.insert(z).second)
            continue;

        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%s (%d)", elementSymbol(z), z);
        options.push_back(std::make_pair(z, std::string(buffer)));
    }

    return options;
}

void clampSpeciesSelectionIndices(int& refSpeciesIndex,
                                  int& targetSpeciesIndex,
                                  int optionCount)
{
    if (refSpeciesIndex >= optionCount)
        refSpeciesIndex = 0;
    if (targetSpeciesIndex >= optionCount)
        targetSpeciesIndex = 0;
}

void drawRdfSummary(const RdfResult& result)
{
    ImGui::SeparatorText("RDF Summary");
    ImGui::Text("Status: %s", result.message.c_str());
    ImGui::Text("Atoms: %d    Ref: %d    Target: %d", result.atomCount, result.refCount, result.targetCount);
    ImGui::Text("PBC used: %s", result.pbcUsed ? "Yes" : "No");
    ImGui::Text("Volume: %.6f A^3    Target density: %.6f A^-3", result.volume, result.density);
    ImGui::Text("Bin width: %.5f A", result.binWidth);

    if (result.hasFirstPeak)
        ImGui::Text("First peak: r = %.4f A, value = %.4f", result.firstPeakR, result.firstPeakValue);
    else
        ImGui::Text("First peak: not detected");

    if (result.hasFirstMinimum)
        ImGui::Text("First minimum after peak: r = %.4f A, value = %.4f", result.firstMinimumR, result.firstMinimumValue);
    else
        ImGui::Text("First minimum after peak: not detected");

    if (result.hasFirstPeakFwhm)
        ImGui::Text("First-peak FWHM: %.4f A", result.firstPeakFwhm);
}

void drawDistortionSummary(const RdfResult& result)
{
    ImGui::SeparatorText("Alloy/Compound Distortion");
    if (!result.valid)
    {
        ImGui::TextDisabled("Run RDF to compute distortion metrics.");
        return;
    }

    if (!result.hasDistortionMetrics)
    {
        ImGui::TextWrapped("No valid shell window was found for distortion analysis. Use manual shell bounds or adjust radius/bin settings.");
        return;
    }

    ImGui::Text("Shell window: [%.4f, %.4f] A (%s)",
                result.distortionWindowMin,
                result.distortionWindowMax,
                result.distortionWindowAuto ? "auto" : "manual");
    ImGui::Text("Shell mean bond length: %.5f A", result.shellMeanDistance);
    ImGui::Text("Shell sigma (bond spread): %.5f A", result.shellStdDev);
    ImGui::Text("Relative distortion sigma/r_mean: %.3f %%", result.shellRelativeDistortion);
    ImGui::Text("Shell coordination (pair-resolved): %.4f", result.shellCoordination);
}

} // namespace

void RadialDistributionAnalysisDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Radial Distribution Function", nullptr, false, enabled))
        m_openRequested = true;
}

void RadialDistributionAnalysisDialog::drawDialog(const Structure& structure)
{
    if (m_task.poll())
        m_computeCompleted = m_task.result().has_value();

    bool& usePbc = m_params.usePbc;
    bool& normalize = m_params.normalize;
    bool& showRawCounts = m_showRawCounts;
    bool& showCumulative = m_showCumulative;
    bool& enableDistortionAnalysis = m_params.enableDistortionAnalysis;
    bool& autoDistortionWindow = m_params.autoDistortionWindow;
    float& manualDistortionMin = m_params.manualDistortionMin;
    float& manualDistortionMax = m_params.manualDistortionMax;
    float& rMin = m_params.rMin;
    float& rMax = m_params.rMax;
    int& binCount = m_params.binCount;
    int& smoothingPasses = m_params.smoothingPasses;
    int& refSpeciesIndex = m_refSpeciesIndex;
    int& targetSpeciesIndex = m_targetSpeciesIndex;

    if (m_openRequested)
    {
        ImGui::OpenPopup("Radial Distribution Function");
        m_openRequested = false;
    }

    std::vector<std::pair<int, std::string>> speciesOptions = buildSpeciesOptions(structure);
    clampSpeciesSelectionIndices(refSpeciesIndex, targetSpeciesIndex, (int)speciesOptions.size());

    ImGui::SetNextWindowSize(ImVec2(1180.0f, 820.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("Radial Distribution Function", &dialogOpen, ImGuiWindowFlags_None))
    {
        bool changed = false;
        if (ImGui::BeginTable("##rdf-layout", 2, ImGuiTableFlags_SizingStretchProp))
        {
            ImGui::TableSetupColumn("##rdf-controls", ImGuiTableColumnFlags_WidthStretch, 0.43f);
            ImGui::TableSetupColumn("##rdf-results", ImGuiTableColumnFlags_WidthStretch, 0.57f);
            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::BeginChild("##rdf-controls-child", ImVec2(0.0f, 0.0f), true);
            ImGui::SeparatorText("Scope and Pair Selection");
            changed |= ImGui::Combo("Reference species", &refSpeciesIndex, speciesLabelGetter, &speciesOptions, (int)speciesOptions.size());
            changed |= ImGui::Combo("Target species", &targetSpeciesIndex, speciesLabelGetter, &speciesOptions, (int)speciesOptions.size());
            changed |= ImGui::Checkbox("Use PBC when unit cell is available", &usePbc);
            changed |= ImGui::Checkbox("Normalize to g(r)", &normalize);

            ImGui::SeparatorText("Sampling");
            changed |= ImGui::DragFloatRange2("Radius range", &rMin, &rMax, 0.01f, 0.0f, 50.0f, "r_min = %.2f", "r_max = %.2f");
            changed |= ImGui::SliderInt("Bins", &binCount, 32, 2000);
            changed |= ImGui::SliderInt("Smoothing passes", &smoothingPasses, 0, 8);

            ImGui::SeparatorText("Plot Overlays");
            changed |= ImGui::Checkbox("Overlay raw counts", &showRawCounts);
            changed |= ImGui::Checkbox("Overlay cumulative coordination", &showCumulative);

            ImGui::SeparatorText("Distortion Metrics (Alloy/Compound)");
            changed |= ImGui::Checkbox("Enable distortion analysis", &enableDistortionAnalysis);
            if (enableDistortionAnalysis)
            {
                changed |= ImGui::Checkbox("Auto shell window from first-shell minima", &autoDistortionWindow);
                if (!autoDistortionWindow)
                    changed |= ImGui::DragFloatRange2("Manual shell window", &manualDistortionMin, &manualDistortionMax, 0.005f, 0.0f, 50.0f, "r_min = %.3f", "r_max = %.3f");
                ImGui::TextWrapped("Use pair-resolved RDF (for example A-B) to track bond-length disorder and local distortion in multicomponent systems.");
            }

            bool computeRequested = false;
            if (ImGui::Button("Run RDF", ImVec2(140.0f, 0.0f)) && !m_task.running())
                computeRequested = true;
            ImGui::SameLine();
            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
                dialogOpen = false;
            ImGui::SameLine();
            if (m_task.running())
                ImGui::TextColored(themeStatusComputing(), "Computing...");
            else if (m_computeCompleted)
                ImGui::TextColored(themeStatusGood(), "Done");
            if (!m_task.error().empty())
                ImGui::TextWrapped("%s", m_task.error().c_str());
            
            if (changed)
                m_computeCompleted = false;
            
            ImGui::EndChild();

            // Launch computation in background thread
            if (computeRequested && !m_task.running())
            {
                auto params = m_params;
                params.refAtomicNumberFilter = speciesOptions[refSpeciesIndex].first;
                params.targetAtomicNumberFilter = speciesOptions[targetSpeciesIndex].first;
                m_computeCompleted = false;
                m_task.start([snapshot = structure, params]() {
                    return atomforge::analysis::computeRdf(snapshot, params);
                });
            }

            ImGui::TableSetColumnIndex(1);
            ImGui::BeginChild("##rdf-results-child", ImVec2(0.0f, 0.0f), true);
            
            // Display results if available
            if (m_task.result())
            {
                const auto* resultPtr = &*m_task.result();
                drawRdfSummary(*resultPtr);
                drawPlot(*resultPtr, showRawCounts, showCumulative);
                drawDistortionSummary(*resultPtr);
            }
            
            ImGui::EndChild();

            ImGui::EndTable();
        }

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();
}
