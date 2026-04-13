#pragma once

#include "io/StructureLoader.h"

#include <thread>
#include <memory>
#include <atomic>

struct RadialDistributionAnalysisDialog
{
    RadialDistributionAnalysisDialog() = default;
    ~RadialDistributionAnalysisDialog();

    void drawMenuItem(bool enabled);
    void drawDialog(const Structure& structure);

private:
    bool m_openRequested = false;
    std::unique_ptr<std::thread> m_computeThread;
    void* m_result = nullptr;  // Opaque pointer to internal RdfResult
    std::atomic<bool> m_isComputing{false};
    std::atomic<bool> m_computeCompleted{false};
    Structure m_lastStructure;
    int m_lastRefZ = 0;
    int m_lastTargetZ = 0;
    bool m_lastUsePbc = true;
    bool m_lastNormalize = true;
    float m_lastRmin = 0.0f;
    float m_lastRmax = 8.0f;
    int m_lastBinCount = 200;
    int m_lastSmoothingPasses = 0;
    bool m_lastEnableDistortionAnalysis = true;
    bool m_lastAutoDistortionWindow = true;
    float m_lastManualDistortionMin = 1.8f;
    float m_lastManualDistortionMax = 3.2f;
    
    friend void deleteRdfResult(void* ptr);  // Forward declare deleter
};
