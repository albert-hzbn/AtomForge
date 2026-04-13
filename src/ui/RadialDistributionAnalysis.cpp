#include "ui/RadialDistributionAnalysis.h"

#include "ElementData.h"
#include "math/StructureMath.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <memory>
#include <atomic>

namespace
{
constexpr float kPi = 3.14159265358979323846f;
constexpr float kMinDistance = 1e-6f;

struct RdfBin
{
    float rCenter = 0.0f;
    float g = 0.0f;
    float rawCount = 0.0f;
    float cumulative = 0.0f;
};

struct RdfResult
{
    bool valid = false;
    std::string message;

    bool pbcUsed = false;
    bool normalized = false;
    int atomCount = 0;
    int refCount = 0;
    int targetCount = 0;
    float volume = 0.0f;
    float density = 0.0f;
    float rMin = 0.0f;
    float rMax = 0.0f;
    int binCount = 0;
    float binWidth = 0.0f;

    float firstPeakR = 0.0f;
    float firstPeakValue = 0.0f;
    float firstMinimumR = 0.0f;
    float firstMinimumValue = 0.0f;
    bool hasFirstPeak = false;
    bool hasFirstMinimum = false;

    bool hasDistortionMetrics = false;
    bool distortionWindowAuto = true;
    float distortionWindowMin = 0.0f;
    float distortionWindowMax = 0.0f;
    float shellMeanDistance = 0.0f;
    float shellStdDev = 0.0f;
    float shellRelativeDistortion = 0.0f;
    float shellCoordination = 0.0f;
    float firstPeakFwhm = 0.0f;
    bool hasFirstPeakFwhm = false;

    std::vector<RdfBin> bins;
};

float computeBoundingVolume(const Structure& structure)
{
    if (structure.atoms.empty())
        return 0.0f;

    glm::vec3 minP((float)structure.atoms[0].x, (float)structure.atoms[0].y, (float)structure.atoms[0].z);
    glm::vec3 maxP = minP;
    for (int i = 1; i < (int)structure.atoms.size(); ++i)
    {
        glm::vec3 p((float)structure.atoms[i].x, (float)structure.atoms[i].y, (float)structure.atoms[i].z);
        minP = glm::min(minP, p);
        maxP = glm::max(maxP, p);
    }

    glm::vec3 extents = glm::max(maxP - minP, glm::vec3(1.0f));
    return extents.x * extents.y * extents.z;
}

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
    drawList->AddRectFilled(pMin, pMax, IM_COL32(20, 24, 31, 255));
    drawList->AddRect(pMin, pMax, IM_COL32(110, 125, 150, 255));

    if (!result.valid || result.bins.empty())
    {
        drawList->AddText(ImVec2(pMin.x + 10.0f, pMin.y + 10.0f), IM_COL32(220, 220, 220, 255), "No RDF data.");
        ImGui::EndChild();
        return;
    }

    const float leftPad = 52.0f;
    const float rightPad = 16.0f;
    const float topPad = 16.0f;
    const float bottomPad = 28.0f;
    ImVec2 gMin(pMin.x + leftPad, pMin.y + topPad);
    ImVec2 gMax(pMax.x - rightPad, pMax.y - bottomPad);

    drawList->AddRect(gMin, gMax, IM_COL32(80, 90, 110, 255));

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
        drawList->AddLine(ImVec2(gMin.x, y), ImVec2(gMax.x, y), IM_COL32(50, 60, 75, 255));
        char label[32];
        std::snprintf(label, sizeof(label), "%.2f", maxY * t);
        drawList->AddText(ImVec2(pMin.x + 4.0f, y - 7.0f), IM_COL32(200, 200, 200, 255), label);
    }

    for (int grid = 0; grid <= 4; ++grid)
    {
        float t = (float)grid / 4.0f;
        float x = gMin.x + (gMax.x - gMin.x) * t;
        drawList->AddLine(ImVec2(x, gMin.y), ImVec2(x, gMax.y), IM_COL32(50, 60, 75, 255));
        float r = result.rMin + (result.rMax - result.rMin) * t;
        char label[32];
        std::snprintf(label, sizeof(label), "%.2f", r);
        drawList->AddText(ImVec2(x - 10.0f, gMax.y + 6.0f), IM_COL32(200, 200, 200, 255), label);
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
        drawList->AddPolyline(&rdfPoints[0], (int)rdfPoints.size(), IM_COL32(80, 220, 255, 255), 0, 2.0f);
    if (showRawCounts && rawPoints.size() >= 2)
        drawList->AddPolyline(&rawPoints[0], (int)rawPoints.size(), IM_COL32(255, 180, 80, 220), 0, 2.0f);
    if (showCumulative && cumulativePoints.size() >= 2)
        drawList->AddPolyline(&cumulativePoints[0], (int)cumulativePoints.size(), IM_COL32(140, 255, 140, 220), 0, 2.0f);

    drawList->AddText(ImVec2(gMin.x + 8.0f, gMin.y + 6.0f), IM_COL32(80, 220, 255, 255), result.normalized ? "g(r)" : "Histogram");
    if (showRawCounts)
        drawList->AddText(ImVec2(gMin.x + 64.0f, gMin.y + 6.0f), IM_COL32(255, 180, 80, 255), "Counts");
    if (showCumulative)
        drawList->AddText(ImVec2(gMin.x + 130.0f, gMin.y + 6.0f), IM_COL32(140, 255, 140, 255), "Cumulative CN");

    ImGui::EndChild();
}

RdfResult runRdf(const Structure& structure,
                 int refAtomicNumberFilter,
                 int targetAtomicNumberFilter,
                 bool usePbcRequest,
                 bool normalize,
                 float rMin,
                 float rMax,
                 int binCount,
                 int smoothingPasses,
                 bool enableDistortionAnalysis,
                 bool autoDistortionWindow,
                 float manualDistortionMin,
                 float manualDistortionMax)
{
    RdfResult result;
    result.atomCount = (int)structure.atoms.size();
    result.normalized = normalize;
    result.rMin = rMin;
    result.rMax = rMax;
    result.binCount = binCount;

    if (structure.atoms.empty())
    {
        result.message = "No atoms available.";
        return result;
    }
    if (rMax <= rMin)
    {
        result.message = "r_max must be greater than r_min.";
        return result;
    }
    if (binCount < 8)
    {
        result.message = "Bin count too small.";
        return result;
    }

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    bool usePbc = false;
    float volume = 0.0f;
    if (usePbcRequest && structure.hasUnitCell)
    {
        if (tryMakeCellMatrices(structure, cell, invCell))
        {
            volume = std::abs(glm::determinant(cell));
            usePbc = true;
        }
    }
    if (!usePbc)
        volume = computeBoundingVolume(structure);

    result.pbcUsed = usePbc;
    result.volume = volume;
    result.binWidth = (rMax - rMin) / (float)binCount;

    std::vector<int> refIndices;
    std::vector<int> targetIndices;
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        int z = structure.atoms[i].atomicNumber;
        if (refAtomicNumberFilter <= 0 || z == refAtomicNumberFilter)
            refIndices.push_back(i);
        if (targetAtomicNumberFilter <= 0 || z == targetAtomicNumberFilter)
            targetIndices.push_back(i);
    }

    result.refCount = (int)refIndices.size();
    result.targetCount = (int)targetIndices.size();
    if (refIndices.empty() || targetIndices.empty())
    {
        result.message = "Reference or target species selection is empty.";
        return result;
    }
    if (volume <= 0.0f)
    {
        result.message = "Unable to determine analysis volume.";
        return result;
    }

    std::vector<float> histogram(binCount, 0.0f);
    std::vector<glm::vec3> positions(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
        positions[i] = glm::vec3((float)structure.atoms[i].x, (float)structure.atoms[i].y, (float)structure.atoms[i].z);

    for (int a = 0; a < (int)refIndices.size(); ++a)
    {
        int i = refIndices[a];
        for (int b = 0; b < (int)targetIndices.size(); ++b)
        {
            int j = targetIndices[b];
            if (i == j)
                continue;

            glm::vec3 delta = minimumImageDelta(positions[j] - positions[i], usePbc, cell, invCell);
            float r = glm::length(delta);
            if (r < rMin || r >= rMax || r <= kMinDistance)
                continue;

            int bin = (int)((r - rMin) / result.binWidth);
            if (bin >= 0 && bin < binCount)
                histogram[bin] += 1.0f;
        }
    }

    if (smoothingPasses > 0)
    {
        for (int pass = 0; pass < smoothingPasses; ++pass)
        {
            std::vector<float> smoothed = histogram;
            for (int i = 1; i + 1 < binCount; ++i)
                smoothed[i] = 0.25f * histogram[i - 1] + 0.5f * histogram[i] + 0.25f * histogram[i + 1];
            histogram.swap(smoothed);
        }
    }

    const float rhoTarget = (float)targetIndices.size() / volume;
    result.density = rhoTarget;
    result.bins.resize(binCount);

    float cumulative = 0.0f;
    for (int i = 0; i < binCount; ++i)
    {
        float r0 = rMin + i * result.binWidth;
        float r1 = r0 + result.binWidth;
        float rc = 0.5f * (r0 + r1);
        float shellVolume = (4.0f / 3.0f) * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
        float expected = (float)refIndices.size() * rhoTarget * shellVolume;
        float g = normalize && expected > 1e-12f ? histogram[i] / expected : histogram[i];

        cumulative += histogram[i] / (float)refIndices.size();

        result.bins[i].rCenter = rc;
        result.bins[i].rawCount = histogram[i];
        result.bins[i].g = g;
        result.bins[i].cumulative = cumulative;
    }

    float bestPeak = -std::numeric_limits<float>::max();
    int peakIndex = -1;
    for (int i = 1; i + 1 < binCount; ++i)
    {
        float y = result.bins[i].g;
        if (y >= result.bins[i - 1].g && y >= result.bins[i + 1].g && y > bestPeak)
        {
            bestPeak = y;
            peakIndex = i;
        }
    }
    if (peakIndex >= 0)
    {
        result.hasFirstPeak = true;
        result.firstPeakR = result.bins[peakIndex].rCenter;
        result.firstPeakValue = result.bins[peakIndex].g;

        int minIndex = -1;
        for (int i = peakIndex + 1; i + 1 < binCount; ++i)
        {
            float y = result.bins[i].g;
            if (y <= result.bins[i - 1].g && y <= result.bins[i + 1].g)
            {
                minIndex = i;
                break;
            }
        }
        if (minIndex >= 0)
        {
            result.hasFirstMinimum = true;
            result.firstMinimumR = result.bins[minIndex].rCenter;
            result.firstMinimumValue = result.bins[minIndex].g;
        }

        float halfHeight = 0.5f * result.firstPeakValue;
        float leftCrossR = result.bins[peakIndex].rCenter;
        float rightCrossR = result.bins[peakIndex].rCenter;

        for (int i = peakIndex; i > 0; --i)
        {
            float y0 = result.bins[i - 1].g;
            float y1 = result.bins[i].g;
            if (y0 <= halfHeight && y1 >= halfHeight)
            {
                float t = std::abs(y1 - y0) > 1e-12f ? (halfHeight - y0) / (y1 - y0) : 0.0f;
                float r0 = result.bins[i - 1].rCenter;
                float r1 = result.bins[i].rCenter;
                leftCrossR = r0 + t * (r1 - r0);
                break;
            }
        }
        for (int i = peakIndex; i + 1 < binCount; ++i)
        {
            float y0 = result.bins[i].g;
            float y1 = result.bins[i + 1].g;
            if (y0 >= halfHeight && y1 <= halfHeight)
            {
                float t = std::abs(y1 - y0) > 1e-12f ? (halfHeight - y0) / (y1 - y0) : 0.0f;
                float r0 = result.bins[i].rCenter;
                float r1 = result.bins[i + 1].rCenter;
                rightCrossR = r0 + t * (r1 - r0);
                break;
            }
        }

        if (rightCrossR > leftCrossR)
        {
            result.hasFirstPeakFwhm = true;
            result.firstPeakFwhm = rightCrossR - leftCrossR;
        }
    }

    if (enableDistortionAnalysis)
    {
        float shellMin = manualDistortionMin;
        float shellMax = manualDistortionMax;
        bool hasWindow = false;

        if (autoDistortionWindow)
        {
            if (result.hasFirstPeak)
            {
                int leftMinIndex = -1;
                int rightMinIndex = -1;

                for (int i = peakIndex - 1; i >= 1; --i)
                {
                    const float y = result.bins[i].g;
                    if (y <= result.bins[i - 1].g && y <= result.bins[i + 1].g)
                    {
                        leftMinIndex = i;
                        break;
                    }
                }
                for (int i = peakIndex + 1; i + 1 < binCount; ++i)
                {
                    const float y = result.bins[i].g;
                    if (y <= result.bins[i - 1].g && y <= result.bins[i + 1].g)
                    {
                        rightMinIndex = i;
                        break;
                    }
                }

                if (rightMinIndex >= 0)
                {
                    shellMin = (leftMinIndex >= 0) ? result.bins[leftMinIndex].rCenter : rMin;
                    shellMax = result.bins[rightMinIndex].rCenter;
                    hasWindow = shellMax > shellMin;
                }
            }
        }
        else
        {
            hasWindow = shellMax > shellMin;
        }

        if (hasWindow)
        {
            float sumW = 0.0f;
            float sumWR = 0.0f;
            float shellRawSum = 0.0f;

            for (int i = 0; i < (int)result.bins.size(); ++i)
            {
                const float r = result.bins[i].rCenter;
                if (r < shellMin || r > shellMax)
                    continue;

                const float w = result.bins[i].rawCount;
                sumW += w;
                sumWR += w * r;
                shellRawSum += w;
            }

            if (sumW > 1e-8f && result.refCount > 0)
            {
                const float mean = sumWR / sumW;
                float var = 0.0f;
                for (int i = 0; i < (int)result.bins.size(); ++i)
                {
                    const float r = result.bins[i].rCenter;
                    if (r < shellMin || r > shellMax)
                        continue;

                    const float w = result.bins[i].rawCount;
                    const float d = r - mean;
                    var += w * d * d;
                }
                var /= sumW;

                result.hasDistortionMetrics = true;
                result.distortionWindowAuto = autoDistortionWindow;
                result.distortionWindowMin = shellMin;
                result.distortionWindowMax = shellMax;
                result.shellMeanDistance = mean;
                result.shellStdDev = std::sqrt(std::max(var, 0.0f));
                result.shellRelativeDistortion = mean > 1e-8f ? (result.shellStdDev / mean) * 100.0f : 0.0f;
                result.shellCoordination = shellRawSum / (float)result.refCount;
            }
        }
    }

    result.valid = true;
    result.message = normalize ? "RDF computed." : "Radial histogram computed.";
    return result;
}

const char* speciesLabelGetter(void* userData, int idx)
{
    const std::vector<std::pair<int, std::string> >* entries = static_cast<const std::vector<std::pair<int, std::string> >*>(userData);
    return (*entries)[idx].second.c_str();
}

std::vector<std::pair<int, std::string> > buildSpeciesOptions(const Structure& structure)
{
    std::vector<std::pair<int, std::string> > options;
    options.push_back(std::make_pair(0, std::string("All elements")));

    std::set<int> seen;
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const int z = structure.atoms[i].atomicNumber;
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

// Deleter function for the opaque result pointer
void deleteRdfResult(void* ptr)
{
    if (ptr)
    {
        delete static_cast<RdfResult*>(ptr);
    }
}

RadialDistributionAnalysisDialog::~RadialDistributionAnalysisDialog()
{
    if (m_computeThread && m_computeThread->joinable())
        m_computeThread->join();
    deleteRdfResult(m_result);
}

void RadialDistributionAnalysisDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Radial Distribution Function", NULL, false, enabled))
        m_openRequested = true;
}

void RadialDistributionAnalysisDialog::drawDialog(const Structure& structure)
{
    static bool usePbc = true;
    static bool normalize = true;
    static bool showRawCounts = false;
    static bool showCumulative = false;
    static bool enableDistortionAnalysis = true;
    static bool autoDistortionWindow = true;
    static float manualDistortionMin = 1.8f;
    static float manualDistortionMax = 3.2f;
    static float rMin = 0.0f;
    static float rMax = 8.0f;
    static int binCount = 200;
    static int smoothingPasses = 0;
    static int refSpeciesIndex = 0;
    static int targetSpeciesIndex = 0;

    if (m_openRequested)
    {
        ImGui::OpenPopup("Radial Distribution Function");
        m_openRequested = false;
    }

    std::vector<std::pair<int, std::string> > speciesOptions = buildSpeciesOptions(structure);
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
            if (ImGui::Button("Run RDF", ImVec2(140.0f, 0.0f)) && !m_isComputing)
                computeRequested = true;
            ImGui::SameLine();
            if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
                dialogOpen = false;
            ImGui::SameLine();
            if (m_isComputing)
                ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Computing...");
            else if (m_computeCompleted)
                ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "Done");
            
            if (changed)
                m_computeCompleted = false;
            
            ImGui::EndChild();

            // Launch computation in background thread
            if (computeRequested && !m_isComputing)
            {
                // Store current parameters
                m_lastStructure = structure;
                m_lastRefZ = speciesOptions[refSpeciesIndex].first;
                m_lastTargetZ = speciesOptions[targetSpeciesIndex].first;
                m_lastUsePbc = usePbc;
                m_lastNormalize = normalize;
                m_lastRmin = rMin;
                m_lastRmax = rMax;
                m_lastBinCount = binCount;
                m_lastSmoothingPasses = smoothingPasses;
                m_lastEnableDistortionAnalysis = enableDistortionAnalysis;
                m_lastAutoDistortionWindow = autoDistortionWindow;
                m_lastManualDistortionMin = manualDistortionMin;
                m_lastManualDistortionMax = manualDistortionMax;
                
                // Join previous thread if it exists
                if (m_computeThread && m_computeThread->joinable())
                    m_computeThread->join();
                
                // Clean up old result
                deleteRdfResult(m_result);
                m_result = nullptr;
                
                // Start new computation thread
                m_isComputing = true;
                m_computeCompleted = false;
                m_computeThread = std::make_unique<std::thread>([this]() {
                    RdfResult* result = new RdfResult();
                    *result = runRdf(m_lastStructure,
                                     m_lastRefZ,
                                     m_lastTargetZ,
                                     m_lastUsePbc,
                                     m_lastNormalize,
                                     m_lastRmin,
                                     m_lastRmax,
                                     m_lastBinCount,
                                     m_lastSmoothingPasses,
                                     m_lastEnableDistortionAnalysis,
                                     m_lastAutoDistortionWindow,
                                     m_lastManualDistortionMin,
                                     m_lastManualDistortionMax);
                    m_result = result;
                    m_isComputing = false;
                    m_computeCompleted = true;
                });
            }

            ImGui::TableSetColumnIndex(1);
            ImGui::BeginChild("##rdf-results-child", ImVec2(0.0f, 0.0f), true);
            
            // Display results if available
            if (m_result)
            {
                RdfResult* resultPtr = static_cast<RdfResult*>(m_result);
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
