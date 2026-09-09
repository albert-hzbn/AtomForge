#pragma once

#include "algorithms/ShortRangeOrderAnalysis.h"
#include "model/Structure.h"
#include <memory>
#include <thread>
#include <atomic>
#include <vector>
#include <string>

struct ShortRangeOrderDialogState {
    ShortRangeOrderDialogState()  = default;
    ~ShortRangeOrderDialogState();
    // Non-copyable due to thread + atomics
    ShortRangeOrderDialogState(const ShortRangeOrderDialogState&)            = delete;
    ShortRangeOrderDialogState& operator=(const ShortRangeOrderDialogState&) = delete;

    bool openRequested = false;
    bool isOpen        = false;

    // Parameters
    int    nShells        = 2;
    double shellTolerance = 0.05;

    // Results (main-thread only after swap)
    bool              hasResults = false;
    SroAnalysisResult results;
    int               selectedMethod = 0;   // 0=Warren-Cowley  1=Rao-Curtin
    int               selectedShell  = 0;

    // Background thread state
    std::unique_ptr<std::thread> workerThread;
    std::atomic<bool>            isComputing{false};
    std::atomic<bool>            computeCompleted{false};
    Structure                    workerStructure;   // copy handed to worker
    SroAnalysisResult            workerResult;      // written by worker
};

/**
 * Draw the menu item for SRO analysis dialog.
 * Should be called within an open ImGui menu context.
 */
void drawShortRangeOrderMenuItem(bool enabled, ShortRangeOrderDialogState& state);

/**
 * Draw the SRO analysis dialog.
 * Should be called once per frame outside of any menu context.
 */
void drawShortRangeOrderDialog(ShortRangeOrderDialogState& state,
                               const Structure& structure);
