#pragma once

#include "io/StructureLoader.h"
#include "algorithms/AmorphousBuilder.h"

#include <functional>
#include <string>
#include <vector>

#include <glm/glm.hpp>

// ---------------------------------------------------------------------------
// AmorphousBuilderDialog
//
// Modal ImGui dialog for generating amorphous structures via Random Sequential
// Addition (RSA).  No 3-D preview — the result is sent to a new tab.
//
// Features exposed in the UI:
//   • Composition table  (element picker, per-element count, derived at%)
//   • Box mode toggle    (manual a/b/c vs. auto-from-target-density)
//   • Cell scale factor  (>1 = expand box to ease packing, like MD pre-melt)
//   • Advanced section   (element-pair min-distance overrides, RNG seed,
//                          max placement attempts, periodic toggle)
//   • Live summary       (total atoms, box volume, estimated density)
// ---------------------------------------------------------------------------
struct AmorphousBuilderDialog
{
    AmorphousBuilderDialog();

    void drawMenuItem(bool enabled);

    // elementColors  – size-119 CPK colors (from EditMenuDialogs::elementColors).
    // covalentRadii  – size-119 covalent radii (from makeLiteratureCovalentRadii).
    // updateBuffers  – callback that opens the result in a new tab.
    void drawDialog(const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>&     covalentRadii,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }

private:
    // -----------------------------------------------------------------------
    // Dialog open/close state
    // -----------------------------------------------------------------------
    bool m_openRequested = false;
    bool m_isOpen        = false;

    // -----------------------------------------------------------------------
    // Composition
    // -----------------------------------------------------------------------
    struct ElementRow
    {
        int atomicNumber = 14;  // default: Si
        int count        = 100;
    };
    std::vector<ElementRow> m_elements;

    // Index of the row currently requesting the element picker (-1 = none).
    int m_pickerTarget = -1;

    // -----------------------------------------------------------------------
    // Box settings
    // -----------------------------------------------------------------------
    int   m_boxMode        = 1;    // 0 = Manual, 1 = AutoFromDensity
    float m_boxA           = 20.0f;
    float m_boxB           = 20.0f;
    float m_boxC           = 20.0f;
    float m_targetDensity  = 2.0f;
    float m_cellScaleFactor= 1.0f;
    bool  m_periodic       = true;

    // -----------------------------------------------------------------------
    // Advanced: pair-distance overrides
    // -----------------------------------------------------------------------
    // Stored as parallel arrays for compact UI iteration.
    struct PairRow
    {
        int   z1      = 14;
        int   z2      = 8;
        float minDist = 1.60f;
        bool  enabled = true;  // when false, row is ignored (uses default)
    };
    std::vector<PairRow> m_pairRows;
    bool m_showAdvanced = false;

    // -----------------------------------------------------------------------
    // Advanced: placement
    // -----------------------------------------------------------------------
    int          m_seed        = 42;
    int          m_maxAttempts = 1000;
    float        m_covTolerance= 0.75f;

    // -----------------------------------------------------------------------
    // Status
    // -----------------------------------------------------------------------
    std::string m_statusMsg;
    bool        m_statusIsError = false;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------
    int  totalAtomCount() const;
    void rebuildPairRows(const std::vector<float>& covalentRadii);
    float defaultMinDist(int z1, int z2, const std::vector<float>& covalentRadii) const;
};
