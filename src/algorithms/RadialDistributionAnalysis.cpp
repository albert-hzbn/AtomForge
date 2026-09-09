#include "algorithms/RadialDistributionAnalysis.h"
#include "math/StructureMath.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace atomforge::analysis
{
namespace
{
constexpr float kPi = 3.14159265358979323846f;
constexpr float kMinDistance = 1e-6f;

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

} // namespace

RdfResult computeRdf(const Structure& structure, const RdfParams& p)
{
    if (!std::isfinite(p.rMin) || !std::isfinite(p.rMax) || p.rMin < 0.0f
        || !std::isfinite(p.manualDistortionMin) || !std::isfinite(p.manualDistortionMax))
    {
        RdfResult result;
        result.message = "Radius parameters must be finite and r_min non-negative.";
        return result;
    }
    return runRdf(structure, p.refAtomicNumberFilter, p.targetAtomicNumberFilter,
                  p.usePbc, p.normalize, p.rMin, p.rMax, p.binCount, p.smoothingPasses,
                  p.enableDistortionAnalysis, p.autoDistortionWindow,
                  p.manualDistortionMin, p.manualDistortionMax);
}
} // namespace atomforge::analysis
