#pragma once

#include "io/StructureLoader.h"

#include <thread>
#include <memory>
#include <atomic>

struct CommonNeighbourAnalysisDialog
{
    CommonNeighbourAnalysisDialog() = default;
    ~CommonNeighbourAnalysisDialog();

    void drawMenuItem(bool enabled);
    void drawDialog(const Structure& structure);

private:
    bool m_openRequested = false;
    std::unique_ptr<std::thread> m_computeThread;
    void* m_result = nullptr;  // Opaque pointer to internal CnaResult
    std::atomic<bool> m_isComputing{false};
    std::atomic<bool> m_computeCompleted{false};
    Structure m_lastStructure;
    float m_lastCutoffScale = 1.18f;
    bool m_lastUsePbc = true;
    
    friend void deleteAnalysisResult(void* ptr);  // Forward declare deleter
};
