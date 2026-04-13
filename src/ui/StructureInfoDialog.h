#pragma once

#include "io/StructureLoader.h"

#include <future>
#include <map>
#include <string>

struct StructureInfoComputedData
{
    bool ready = false;
    std::size_t atomCount = 0;
    std::map<int, int> elementCounts;
    std::string formula;
    bool symmetryAttempted = false;
    bool symmetrySuccess = false;
    int symmetrySpaceGroupNumber = 0;
    std::string symmetryInternationalSymbol;
    std::string symmetryHallSymbol;
    std::string symmetryPointGroup;
    std::string symmetryError;
};

struct StructureInfoDialogState
{
    bool openRequested = false;
    bool workerRunning = false;
    std::size_t pendingFingerprint = 0;
    std::size_t cachedFingerprint = 0;
    std::future<StructureInfoComputedData> workerFuture;
    StructureInfoComputedData cachedData;
};

void drawStructureInfoDialog(StructureInfoDialogState& state,
                             bool requestOpen,
                             const Structure& structure);