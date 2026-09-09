#include "AngularDistributionAnalysis.h"
#include "math/StructureMath.h"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <cmath>
#include <set>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static constexpr float kPiAdf = 3.14159265358979323846f;

static glm::vec3 adfMinImage(glm::vec3 delta, const glm::mat3& cell,
                              const glm::mat3& invCell, bool usePbc)
{
    if (!usePbc) return delta;
    glm::vec3 frac = invCell * delta;
    frac.x -= std::round(frac.x);
    frac.y -= std::round(frac.y);
    frac.z -= std::round(frac.z);
    return cell * frac;
}

static float angleDeg(const glm::vec3& rij, const glm::vec3& rik)
{
    float d1 = glm::length(rij);
    float d2 = glm::length(rik);
    if (d1 < 1e-6f || d2 < 1e-6f) return 0.0f;
    float cosA = glm::dot(rij, rik) / (d1 * d2);
    cosA = std::max(-1.0f, std::min(1.0f, cosA));
    return std::acos(cosA) * (180.0f / kPiAdf);
}

// One-pass Gaussian-like box smooth
static void smoothHistogram(std::vector<float>& v, int passes)
{
    if (v.size() < 3) return;
    for (int p = 0; p < passes; ++p) {
        std::vector<float> tmp(v.size(), 0.0f);
        for (int i = 0; i < (int)v.size(); ++i) {
            float s = v[i] * 0.5f;
            if (i > 0)              s += v[i-1] * 0.25f;
            if (i < (int)v.size()-1) s += v[i+1] * 0.25f;
            tmp[i] = s;
        }
        v = tmp;
    }
}

// Simple peak finder: local max with value > threshold * globalMax
static std::vector<int> findPeaks(const std::vector<float>& v, float threshold = 0.08f)
{
    float gmax = *std::max_element(v.begin(), v.end());
    if (gmax < 1e-10f) return {};
    std::vector<int> peaks;
    for (int i = 1; i + 1 < (int)v.size(); ++i) {
        if (v[i] > v[i-1] && v[i] >= v[i+1] && v[i] > threshold * gmax)
            peaks.push_back(i);
    }
    return peaks;
}

// Label common coordination geometries from angle
static std::string geometryLabel(float deg)
{
    struct Known { float angle; const char* label; };
    static const Known kKnown[] = {
        { 60.0f,   "triangular 60\xc2\xb0" },
        { 70.5f,   "FCC/HCP 70.5\xc2\xb0" },
        { 90.0f,   "octahedral 90\xc2\xb0" },
        { 109.47f, "tetrahedral 109.5\xc2\xb0" },
        { 120.0f,  "trigonal planar 120\xc2\xb0" },
        { 135.0f,  "square antiprism 135\xc2\xb0" },
        { 150.0f,  "BCC 150\xc2\xb0" },
        { 180.0f,  "linear 180\xc2\xb0" },
    };
    for (const auto& k : kKnown)
        if (std::abs(deg - k.angle) < 3.5f)
            return k.label;
    return "";
}

// ---------------------------------------------------------------------------
// Main compute
// ---------------------------------------------------------------------------

AdfResult computeADF(const Structure& structure, const AdfParams& params)
{
    AdfResult result;
    result.rCutoff   = params.rCutoff;
    result.binCount  = params.binCount;
    result.binWidth  = 180.0f / static_cast<float>(params.binCount);
    result.normalized = params.normalize;
    result.nAtoms    = static_cast<int>(structure.atoms.size());

    if (structure.atoms.empty()) {
        result.message = "No atoms in structure.";
        return result;
    }
    if (params.rCutoff < 0.1f) {
        result.message = "Cutoff radius too small.";
        return result;
    }

    // Count elements
    for (const auto& a : structure.atoms)
        result.elementCounts[a.symbol]++;

    // Cell matrices
    glm::mat3 cell(1.0f), invCell(1.0f);
    bool usePbc = false;
    if (params.usePbc && structure.hasUnitCell) {
        if (tryMakeCellMatrices(structure, cell, invCell))
            usePbc = true;
    }
    result.pbcUsed = usePbc;

    // Determine which atoms can be centres / neighbours
    auto isCentre = [&](const AtomSite& a) -> bool {
        if (params.centreMode == AdfCentreMode::All) return true;
        return a.symbol == params.centreSymbol;
    };
    auto isNeigh1 = [&](const AtomSite& a) -> bool {
        if (params.centreMode != AdfCentreMode::ByPair) return true;
        return params.neighSymbol1.empty() || a.symbol == params.neighSymbol1;
    };
    auto isNeigh2 = [&](const AtomSite& a) -> bool {
        if (params.centreMode != AdfCentreMode::ByPair) return true;
        return params.neighSymbol2.empty() || a.symbol == params.neighSymbol2;
    };

    const float rc2 = params.rCutoff * params.rCutoff;
    int N = result.nAtoms;

    // Histogram bins [0, 180°]
    std::vector<float> hist(params.binCount, 0.0f);

    // Coordination counts for stats
    std::map<std::string, std::vector<int>> coordCounts;

    long long nTriplets = 0;

    for (int i = 0; i < N; ++i)
    {
        const AtomSite& ai = structure.atoms[i];
        if (!isCentre(ai)) continue;

        glm::vec3 pi(ai.x, ai.y, ai.z);

        // Collect neighbours of i
        std::vector<std::pair<int, glm::vec3>> neighs;
        for (int j = 0; j < N; ++j)
        {
            if (j == i) continue;
            const AtomSite& aj = structure.atoms[j];
            glm::vec3 pj(aj.x, aj.y, aj.z);
            glm::vec3 rij = adfMinImage(pj - pi, cell, invCell, usePbc);
            float d2 = glm::dot(rij, rij);
            if (d2 > 1e-12f && d2 <= rc2)
                neighs.push_back({j, rij});
        }

        coordCounts[ai.symbol].push_back(static_cast<int>(neighs.size()));

        // All ordered (j, k) pairs → angle at i
        for (int a = 0; a < (int)neighs.size(); ++a)
        {
            if (!isNeigh1(structure.atoms[neighs[a].first])) continue;
            for (int b = a + 1; b < (int)neighs.size(); ++b)
            {
                if (!isNeigh2(structure.atoms[neighs[b].first])) continue;
                float ang = angleDeg(neighs[a].second, neighs[b].second);
                int bin = static_cast<int>(ang / result.binWidth);
                if (bin >= params.binCount) bin = params.binCount - 1;
                hist[bin] += 1.0f;
                ++nTriplets;
            }
        }
    }

    result.nTriplets    = nTriplets;
    result.nCentreAtoms = 0;
    for (int i = 0; i < N; ++i)
        if (isCentre(structure.atoms[i])) ++result.nCentreAtoms;

    if (nTriplets == 0) {
        result.message = "No triplets found — try a larger cutoff radius.";
        return result;
    }

    // Smooth
    smoothHistogram(hist, params.smoothPasses);

    // Build bin structs
    result.bins.resize(params.binCount);
    float gmax = *std::max_element(hist.begin(), hist.end());
    for (int b = 0; b < params.binCount; ++b) {
        result.bins[b].angleDeg = (b + 0.5f) * result.binWidth;
        result.bins[b].count    = hist[b];
        result.bins[b].value    = (params.normalize && gmax > 1e-10f)
                                  ? hist[b] / gmax : hist[b];
    }

    // Detect peaks
    std::vector<float> plotVals(params.binCount);
    for (int b = 0; b < params.binCount; ++b)
        plotVals[b] = result.bins[b].value;
    for (int idx : findPeaks(plotVals)) {
        AdfPeakInfo pk;
        pk.angleDeg = result.bins[idx].angleDeg;
        pk.value    = result.bins[idx].value;
        pk.label    = geometryLabel(pk.angleDeg);
        result.peaks.push_back(pk);
    }

    // Coordination stats
    for (auto& [sym, cnts] : coordCounts) {
        float sum = 0.0f, sum2 = 0.0f;
        for (int c : cnts) { sum += c; sum2 += c * c; }
        float mean = sum / cnts.size();
        float var  = sum2 / cnts.size() - mean * mean;
        result.coordStats[sym] = { mean, std::sqrt(std::max(0.0f, var)),
                                   static_cast<int>(cnts.size()) };
    }

    result.valid   = true;
    result.message = "ADF computed successfully.";
    return result;
}
