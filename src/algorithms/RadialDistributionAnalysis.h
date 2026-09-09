#pragma once

#include "model/Structure.h"
#include <string>
#include <vector>

namespace atomforge::analysis
{
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

struct RdfParams
{
    int refAtomicNumberFilter = 0;
    int targetAtomicNumberFilter = 0;
    bool usePbc = true;
    bool normalize = true;
    float rMin = 0.0f;
    float rMax = 8.0f;
    int binCount = 200;
    int smoothingPasses = 0;
    bool enableDistortionAnalysis = true;
    bool autoDistortionWindow = true;
    float manualDistortionMin = 1.8f;
    float manualDistortionMax = 3.2f;
};

[[nodiscard]] RdfResult computeRdf(const Structure& structure, const RdfParams& params = {});
} // namespace atomforge::analysis
