#pragma once

#include "model/Structure.h"
#include <string>
#include <vector>
#include <map>

// ---------------------------------------------------------------------------
// Angular Distribution Function (ADF) — i-j-k triplet bond angles
//
// For every ordered triplet (j, i, k) where both j and i, and k and i are
// within rCutoff, the angle at vertex i is accumulated into a histogram.
// Three centre-selection modes are supported:
//   All        — every atom can be a centre
//   ByElement  — only atoms of centreSymbol are centres
//   ByPair     — centre of type centreSymbol, neighbours limited to
//                neighSymbol1 and neighSymbol2
// ---------------------------------------------------------------------------

enum class AdfCentreMode {
    All,        // centre = any atom
    ByElement,  // centre = centreSymbol only
    ByPair      // explicit i, j, k types
};

struct AdfParams {
    float        rCutoff       = 3.5f;   // neighbour cutoff (Å)
    int          binCount      = 180;    // histogram bins in [0, 180] deg
    bool         usePbc        = true;
    AdfCentreMode centreMode   = AdfCentreMode::All;
    std::string  centreSymbol;           // used for ByElement / ByPair
    std::string  neighSymbol1;           // used for ByPair (j type)
    std::string  neighSymbol2;           // used for ByPair (k type)
    bool         normalize     = true;   // normalise to peak = 1
    int          smoothPasses  = 2;      // Gaussian smoothing passes
};

struct AdfBin {
    float angleDeg = 0.0f;   // bin centre (°)
    float count    = 0.0f;   // raw count
    float value    = 0.0f;   // after normalisation / smoothing
};

struct AdfPeakInfo {
    float angleDeg  = 0.0f;
    float value     = 0.0f;
    std::string label;        // e.g. "tetrahedral 109.5°"
};

struct AdfResult {
    bool        valid         = false;
    std::string message;

    // Parameters echoed back
    float       rCutoff       = 0.0f;
    int         binCount      = 0;
    float       binWidth      = 0.0f;
    bool        pbcUsed       = false;
    bool        normalized    = false;
    int         nAtoms        = 0;
    int         nCentreAtoms  = 0;
    long long   nTriplets     = 0;

    // Element counts for display
    std::map<std::string, int> elementCounts;

    // Histogram
    std::vector<AdfBin> bins;

    // Auto-detected peaks
    std::vector<AdfPeakInfo> peaks;

    // Coordination statistics per centre element
    struct CoordStat { float mean = 0.0f; float stddev = 0.0f; int count = 0; };
    std::map<std::string, CoordStat> coordStats;
};

// Main entry point – pure compute, no OpenGL / ImGui
AdfResult computeADF(const Structure& structure, const AdfParams& params);
