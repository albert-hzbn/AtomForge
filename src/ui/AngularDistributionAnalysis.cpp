#include "AngularDistributionAnalysis.h"
#include "ThemeUtils.h"

#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

// ============================================================
// Destructor
// ============================================================

AngularDistributionAnalysisDialog::~AngularDistributionAnalysisDialog()
{
    if (m_workerThread && m_workerThread->joinable())
        m_workerThread->join();
}

// ============================================================
// Menu item
// ============================================================

void AngularDistributionAnalysisDialog::drawMenuItem(bool enabled)
{
    if (!enabled) ImGui::BeginDisabled();
    if (ImGui::MenuItem("Angular Distribution Function"))
        m_openRequested = true;
    if (!enabled)
    {
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Load a structure first.");
    }
}

// ============================================================
// Background thread
// ============================================================

void AngularDistributionAnalysisDialog::startCompute(const Structure& structure)
{
    if (m_isComputing) return;

    m_workerStructure  = structure;
    m_isComputing      = true;
    m_computeCompleted = false;
    m_paramsDirty      = false;

    AdfParams p;
    p.rCutoff      = m_rCutoff;
    p.binCount     = m_binCount;
    p.usePbc       = m_usePbc;
    p.normalize    = m_normalize;
    p.smoothPasses = m_smoothPasses;
    p.centreSymbol = m_centreSymbol;
    p.neighSymbol1 = m_neighSym1;
    p.neighSymbol2 = m_neighSym2;

    switch (m_centreMode)
    {
        case 0:  p.centreMode = AdfCentreMode::All;       break;
        case 1:  p.centreMode = AdfCentreMode::ByElement; break;
        default: p.centreMode = AdfCentreMode::ByPair;    break;
    }

    if (m_workerThread && m_workerThread->joinable())
        m_workerThread->join();

    m_workerThread = std::make_unique<std::thread>([this, p]() {
        m_workerResult     = computeADF(m_workerStructure, p);
        m_computeCompleted = true;
        m_isComputing      = false;
    });
}

void AngularDistributionAnalysisDialog::pollWorker()
{
    if (m_computeCompleted.load())
    {
        m_result           = m_workerResult;
        m_hasResult        = m_result.valid;
        m_computeCompleted = false;
        m_xMin             = 0.0f;
        m_xMax             = 180.0f;
    }
}

void AngularDistributionAnalysisDialog::updateElementList(const Structure& structure)
{
    m_elements.clear();
    for (const auto& a : structure.atoms)
    {
        bool found = false;
        for (const auto& s : m_elements)
            if (s == a.symbol) { found = true; break; }
        if (!found) m_elements.push_back(a.symbol);
    }
}

// ============================================================
// Plot
// ============================================================

void AngularDistributionAnalysisDialog::drawPlot()
{
    const bool lt = isLightTheme();

    const ImU32 bgCol      = lt ? IM_COL32(245, 246, 248, 255) : IM_COL32(20, 24, 31, 255);
    const ImU32 borderCol  = lt ? IM_COL32(140, 145, 160, 255) : IM_COL32(110, 125, 150, 255);
    const ImU32 gridCol    = lt ? IM_COL32(170, 175, 190, 180) : IM_COL32(50, 60, 75, 255);
    const ImU32 labelCol   = lt ? IM_COL32(50,  50,  65, 255)  : IM_COL32(200, 200, 200, 255);
    const ImU32 lineCol    = lt ? IM_COL32(20, 120, 210, 255)  : IM_COL32(80, 220, 255, 255);
    const ImU32 fillCol    = lt ? IM_COL32(20, 120, 210,  40)  : IM_COL32(80, 220, 255,  28);
    const ImU32 peakCol    = lt ? IM_COL32(210,  35,  35, 220) : IM_COL32(255, 100, 100, 230);
    const ImU32 peakLineC  = lt ? IM_COL32(210,  35,  35,  70) : IM_COL32(255, 100, 100,  55);
    const ImU32 noDataCol  = lt ? IM_COL32( 80,  80, 100, 255) : IM_COL32(200, 200, 200, 255);
    const ImU32 refLineCol = lt ? IM_COL32(130, 142, 160, 130) : IM_COL32( 65,  75,  95, 200);
    const ImU32 refLblCol  = lt ? IM_COL32(110, 122, 145, 210) : IM_COL32(120, 138, 165, 220);
    const ImU32 crosshairC = lt ? IM_COL32( 80,  88, 108, 160) : IM_COL32(185, 195, 215, 150);

    ImGui::BeginChild("##adf-plot", ImVec2(-1.0f, 270.0f), true,
                      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 origin  = ImGui::GetCursorScreenPos();
    ImVec2 avail   = ImGui::GetContentRegionAvail();
    if (avail.x < 60.0f) avail.x = 60.0f;
    if (avail.y < 60.0f) avail.y = 60.0f;
    ImGui::InvisibleButton("##adf-canvas", avail);
    const bool plotHovered = ImGui::IsItemHovered();

    ImVec2 pMin = origin;
    ImVec2 pMax(origin.x + avail.x, origin.y + avail.y);
    dl->AddRectFilled(pMin, pMax, bgCol);
    dl->AddRect(pMin, pMax, borderCol);

    if (!m_hasResult || m_result.bins.empty())
    {
        const char* msg = m_isComputing.load()
            ? "Computing \xe2\x80\x94 please wait..."
            : "No data.  Configure parameters and click  Compute ADF.";
        ImVec2 ts = ImGui::CalcTextSize(msg);
        dl->AddText(ImVec2(pMin.x + (avail.x - ts.x) * 0.5f,
                           pMin.y + (avail.y - ts.y) * 0.5f),
                    noDataCol, msg);
        ImGui::EndChild();
        return;
    }

    const float leftPad   = 54.0f;
    const float rightPad  = 16.0f;
    const float topPad    = 22.0f;
    const float bottomPad = 34.0f;
    ImVec2 gMin(pMin.x + leftPad, pMin.y + topPad);
    ImVec2 gMax(pMax.x - rightPad, pMax.y - bottomPad);
    dl->AddRect(gMin, gMax, lt ? IM_COL32(150,155,170,255) : IM_COL32(80,90,110,255));

    const float xRange = std::max(m_xMax - m_xMin, 1.0f);

    // Max Y in visible range
    float maxY = 0.0f;
    for (const auto& b : m_result.bins)
    {
        if (b.angleDeg < m_xMin || b.angleDeg > m_xMax) continue;
        maxY = std::max(maxY, m_showRawCounts ? b.count : b.value);
    }
    maxY = std::max(maxY, 1e-6f);

    // Y grid
    for (int g = 0; g <= 4; ++g)
    {
        float t = (float)g / 4.0f;
        float y = gMax.y + (gMin.y - gMax.y) * t;
        dl->AddLine(ImVec2(gMin.x, y), ImVec2(gMax.x, y), gridCol);
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", maxY * t);
        dl->AddText(ImVec2(pMin.x + 3.0f, y - 7.0f), labelCol, buf);
    }

    // Reference geometry lines with abbreviated labels
    struct RefAngle { float angle; const char* abbr; };
    static const RefAngle kRefs[] = {
        {  60.00f, "tri"  }, {  70.50f, "fcc"  }, {  90.00f, "oct"  },
        { 109.47f, "tet"  }, { 120.00f, "trpl" }, { 150.00f, "bcc"  },
        { 180.00f, "lin"  },
    };
    // Track last placed abbr label x to avoid overlaps
    float lastAbbrRight = -9999.0f;
    for (const auto& ref : kRefs)
    {
        if (ref.angle < m_xMin || ref.angle > m_xMax) continue;
        float tx = (ref.angle - m_xMin) / xRange;
        float x  = gMin.x + (gMax.x - gMin.x) * tx;
        // Dashed vertical line
        const float gH = gMax.y - gMin.y;
        for (int s = 0; s < 5; s += 2)
        {
            float y0 = gMin.y + gH * (float)s     / 5.0f;
            float y1 = gMin.y + gH * (float)(s+1) / 5.0f;
            dl->AddLine(ImVec2(x, y0), ImVec2(x, y1), refLineCol, 1.0f);
        }
        // Abbr label at top — only draw if it won't overlap previous label
        float alw = ImGui::CalcTextSize(ref.abbr).x;
        float lblX = x - alw * 0.5f;
        if (lblX > lastAbbrRight + 2.0f)
        {
            dl->AddText(ImVec2(lblX, gMin.y + 2.0f), refLblCol, ref.abbr);
            lastAbbrRight = lblX + alw;
        }
    }

    // Build polyline
    std::vector<ImVec2> pts;
    pts.reserve(m_result.bins.size());
    for (const auto& b : m_result.bins)
    {
        if (b.angleDeg < m_xMin || b.angleDeg > m_xMax) continue;
        float tx = (b.angleDeg - m_xMin) / xRange;
        float x  = gMin.x + (gMax.x - gMin.x) * tx;
        float v  = m_showRawCounts ? b.count : b.value;
        float fy = std::max(0.0f, std::min(1.0f, v / maxY));
        float y  = gMax.y - (gMax.y - gMin.y) * fy;
        pts.push_back(ImVec2(x, y));
    }

    // Filled area under the curve (one trapezoid per bin segment)
    if (pts.size() >= 2)
    {
        for (int i = 0; i + 1 < (int)pts.size(); ++i)
        {
            dl->AddQuadFilled(
                ImVec2(pts[i].x,   gMax.y),
                ImVec2(pts[i+1].x, gMax.y),
                pts[i+1], pts[i], fillCol);
        }
    }

    // Curve line on top
    if ((int)pts.size() >= 2)
        dl->AddPolyline(pts.data(), (int)pts.size(), lineCol, 0, 2.0f);

    // Peak markers
    if (m_showPeaks)
    {
        for (const auto& pk : m_result.peaks)
        {
            if (pk.angleDeg < m_xMin || pk.angleDeg > m_xMax) continue;
            float tx  = (pk.angleDeg - m_xMin) / xRange;
            float x   = gMin.x + (gMax.x - gMin.x) * tx;
            float fy  = std::max(0.0f, std::min(1.0f, pk.value / maxY));
            float y   = gMax.y - (gMax.y - gMin.y) * fy;
            dl->AddLine(ImVec2(x, gMin.y + 2.0f), ImVec2(x, y - 6.0f), peakLineC, 1.5f);
            dl->AddCircleFilled(ImVec2(x, y), 5.0f, peakCol);
            dl->AddCircle(ImVec2(x, y), 5.0f,
                          lt ? IM_COL32(170,20,20,200) : IM_COL32(255,170,170,200));
            char buf[20];
            std::snprintf(buf, sizeof(buf), "%.1f\xc2\xb0", pk.angleDeg);
            float bw2 = ImGui::CalcTextSize(buf).x;
            float fh  = ImGui::GetFontSize();
            float lblY = std::max(gMin.y + 14.0f, y - fh - 6.0f);
            dl->AddText(ImVec2(x - bw2 * 0.5f, lblY), peakCol, buf);
        }
    }

    // Y-axis label
    const char* yLbl = m_showRawCounts ? "Counts"
                     : (m_result.normalized ? "g(\xce\xb8) [norm]" : "g(\xce\xb8)");
    dl->AddText(ImVec2(gMin.x + 5.0f, pMin.y + 4.0f), lineCol, yLbl);

    // X-axis label
    const char* xLbl = "Angle (\xc2\xb0)";
    float xlw = ImGui::CalcTextSize(xLbl).x;
    dl->AddText(ImVec2(gMin.x + (gMax.x - gMin.x) * 0.5f - xlw * 0.5f, gMax.y + 20.0f),
                labelCol, xLbl);

    // Hover crosshair + tooltip
    if (plotHovered)
    {
        ImVec2 mp = ImGui::GetMousePos();
        if (mp.x >= gMin.x && mp.x <= gMax.x)
        {
            float tx    = (mp.x - gMin.x) / (gMax.x - gMin.x);
            float angle = m_xMin + tx * xRange;
            float bwBin = std::max(m_result.binWidth, 1e-4f);
            int binIdx  = static_cast<int>(angle / bwBin);
            binIdx = std::max(0, std::min(binIdx, (int)m_result.bins.size() - 1));
            const auto& bin = m_result.bins[binIdx];
            float val = m_showRawCounts ? bin.count : bin.value;

            dl->AddLine(ImVec2(mp.x, gMin.y), ImVec2(mp.x, gMax.y), crosshairC, 1.0f);

            float fy = std::max(0.0f, std::min(1.0f, val / maxY));
            float cy = gMax.y - (gMax.y - gMin.y) * fy;
            dl->AddCircleFilled(ImVec2(mp.x, cy), 3.5f, crosshairC);

            ImGui::BeginTooltip();
            ImGui::Text("\xce\xb8 = %.2f\xc2\xb0", bin.angleDeg);
            ImGui::Text("g(\xce\xb8) = %.5f", val);
            if (m_result.normalized && bin.count > 0.0f)
                ImGui::TextDisabled("Raw count: %.0f", bin.count);
            ImGui::EndTooltip();
        }
    }

    ImGui::EndChild();
}

// ============================================================
// Peak table
// ============================================================

void AngularDistributionAnalysisDialog::drawPeakTable()
{
    if (!m_hasResult || m_result.peaks.empty())
    {
        ImGui::Spacing();
        ImGui::TextDisabled("No peaks detected. Try reducing smooth passes or recompute with finer bins.");
        return;
    }

    ImGui::Spacing();
    ImGui::TextDisabled("%d peak(s) detected:", (int)m_result.peaks.size());
    ImGui::Spacing();

    if (ImGui::BeginTable("##adf-peaks", 3,
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
        ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY))
    {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Angle",     ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableSetupColumn("Intensity", ImGuiTableColumnFlags_WidthStretch, 0.38f);
        ImGui::TableSetupColumn("Geometry",  ImGuiTableColumnFlags_WidthStretch, 0.40f);
        ImGui::TableHeadersRow();

        for (const auto& pk : m_result.peaks)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%.2f\xc2\xb0", pk.angleDeg);

            ImGui::TableSetColumnIndex(1);
            {
                char valBuf[20];
                std::snprintf(valBuf, sizeof(valBuf), "%.4f", pk.value);
                ImGui::TextUnformatted(valBuf);
                ImGui::SameLine();
                float cw = ImGui::GetContentRegionAvail().x - 4.0f;
                if (cw > 8.0f)
                {
                    float barW = cw * std::min(pk.value, 1.0f);
                    ImVec2 pos = ImGui::GetCursorScreenPos();
                    const bool lt2 = isLightTheme();
                    ImGui::GetWindowDrawList()->AddRectFilled(
                        pos, ImVec2(pos.x + barW, pos.y + ImGui::GetFrameHeight() - 4.0f),
                        lt2 ? IM_COL32(20,120,210,180) : IM_COL32(80,220,255,160));
                    ImGui::Dummy(ImVec2(cw, ImGui::GetFrameHeight() - 4.0f));
                }
            }

            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(pk.label.empty() ? "\xe2\x80\x94" : pk.label.c_str());
        }
        ImGui::EndTable();
    }
}

// ============================================================
// Coordination stats
// ============================================================

void AngularDistributionAnalysisDialog::drawCoordStats()
{
    if (!m_hasResult || m_result.coordStats.empty())
    {
        ImGui::Spacing();
        ImGui::TextDisabled("No coordination data. Run a computation first.");
        return;
    }

    ImGui::Spacing();
    ImGui::TextDisabled("Mean coordination within cutoff %.2f \xc3\x85:", m_result.rCutoff);
    ImGui::Spacing();

    float maxCN = 1.0f;
    for (const auto& [s, cs] : m_result.coordStats)
        maxCN = std::max(maxCN, cs.mean);

    if (ImGui::BeginTable("##adf-cn", 5,
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
        ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Elem",    ImGuiTableColumnFlags_WidthFixed, 44.0f);
        ImGui::TableSetupColumn("Mean CN", ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableSetupColumn("",        ImGuiTableColumnFlags_WidthStretch, 0.30f);
        ImGui::TableSetupColumn("\xc2\xb1 Std",    ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableSetupColumn("Centres", ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableHeadersRow();

        for (const auto& [sym, cs] : m_result.coordStats)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(sym.c_str());

            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.2f", cs.mean);

            ImGui::TableSetColumnIndex(2);
            {
                float cw = ImGui::GetContentRegionAvail().x;
                float barW = cw * (cs.mean / maxCN);
                ImVec2 pos = ImGui::GetCursorScreenPos();
                const bool lt2 = isLightTheme();
                ImGui::GetWindowDrawList()->AddRectFilled(
                    ImVec2(pos.x, pos.y + 1.0f),
                    ImVec2(pos.x + barW, pos.y + ImGui::GetFrameHeight() - 3.0f),
                    lt2 ? IM_COL32(30,140,60,160) : IM_COL32(80,230,120,140));
                ImGui::Dummy(ImVec2(cw, ImGui::GetFrameHeight() - 4.0f));
            }

            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.2f", cs.stddev);

            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%d", cs.count);
        }
        ImGui::EndTable();
    }
}

// ============================================================
// Settings panel (scrollable part)
// ============================================================

void AngularDistributionAnalysisDialog::drawSettings(const Structure& structure)
{
    const bool computing = m_isComputing.load();

    // Label column wide enough for "Normalise" — the longest label used below.
    // Using a BeginTable instead of SameLine(offset) avoids all clipping issues
    // caused by SameLine's window-origin-relative coordinate system.
    const float labelColW = std::ceil(ImGui::CalcTextSize("Normalise").x)
                          + ImGui::GetStyle().ItemSpacing.x;

    // Element-combo helper — returns true when selection changed
    auto elemCombo = [&](const char* id, char* buf, int bufSz) -> bool
    {
        bool changed = false;
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo(id, buf[0] ? buf : "(any)", ImGuiComboFlags_HeightSmall))
        {
            if (ImGui::Selectable("(any)", buf[0] == '\0'))
                { buf[0] = '\0'; changed = true; }
            for (const auto& s : m_elements)
                if (ImGui::Selectable(s.c_str(), s == buf))
                    { std::snprintf(buf, bufSz, "%s", s.c_str()); changed = true; }
            ImGui::EndCombo();
        }
        return changed;
    };

    constexpr ImGuiTableFlags kRowFlags =
        ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_NoBordersInBody;

    if (computing) ImGui::BeginDisabled();

    // ── Computation ───────────────────────────────────────────
    ImGui::SeparatorText("Computation");

    if (ImGui::BeginTable("##adf-comp", 2, kRowFlags))
    {
        ImGui::TableSetupColumn("##lc", ImGuiTableColumnFlags_WidthFixed, labelColW);
        ImGui::TableSetupColumn("##rc", ImGuiTableColumnFlags_WidthStretch);

        // Cutoff
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Cutoff");
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::DragFloat("##cutoff", &m_rCutoff, 0.05f, 0.5f, 15.0f, "%.2f \xc3\x85"))
            m_paramsDirty = true;
        if (m_rCutoff < 0.5f)  m_rCutoff = 0.5f;
        if (m_rCutoff > 15.0f) m_rCutoff = 15.0f;
        ImGui::SetItemTooltip(
            "Maximum neighbour distance (\xc3\x85) for\n"
            "finding bonded triplets.\n"
            "Drag or Ctrl+click to type a value.");

        // Bins
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Bins");
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::DragInt("##bins", &m_binCount, 1.0f, 30, 360))
            m_paramsDirty = true;
        if (m_binCount < 30)  m_binCount = 30;
        if (m_binCount > 360) m_binCount = 360;
        {
            float bwDeg = 180.0f / (float)std::max(1, m_binCount);
            ImGui::SetItemTooltip("%d bins \xe2\x80\x94 %.2f\xc2\xb0 per bin.", m_binCount, bwDeg);
        }

        // Use PBC (nested disable when no unit cell)
        {
            bool pbcAvail = structure.hasUnitCell;
            if (!pbcAvail) ImGui::BeginDisabled();
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted("Use PBC");
            ImGui::TableSetColumnIndex(1);
            if (ImGui::Checkbox("##pbc", &m_usePbc))
                m_paramsDirty = true;
            if (!pbcAvail)
            {
                ImGui::EndDisabled();
                ImGui::SetItemTooltip("No unit cell \xe2\x80\x94 PBC unavailable.");
            }
            else
            {
                ImGui::SetItemTooltip("Apply minimum-image PBC when searching for neighbours.");
            }
        }

        ImGui::EndTable();
    }

    // ── Atom Selection ────────────────────────────────────────
    ImGui::SeparatorText("Atom Selection");

    if (ImGui::RadioButton("All",     &m_centreMode, 0)) m_paramsDirty = true;
    ImGui::SetItemTooltip("Every atom is a valid triplet centre.");
    if (ImGui::RadioButton("Element", &m_centreMode, 1)) m_paramsDirty = true;
    ImGui::SetItemTooltip("Only atoms of the chosen element type are centres.");
    if (ImGui::RadioButton("Triplet", &m_centreMode, 2)) m_paramsDirty = true;
    ImGui::SetItemTooltip("Specify centre and both neighbour element\ntypes for j\xe2\x80\x93i\xe2\x80\x93k triplets.");

    if (m_centreMode >= 1)
    {
        ImGui::Spacing();
        const float filtLabelW = std::ceil(ImGui::CalcTextSize("Neigh j:").x)
                               + ImGui::GetStyle().ItemSpacing.x;
        if (ImGui::BeginTable("##adf-filt", 2, kRowFlags))
        {
            ImGui::TableSetupColumn("##fl", ImGuiTableColumnFlags_WidthFixed, filtLabelW);
            ImGui::TableSetupColumn("##fc", ImGuiTableColumnFlags_WidthStretch);

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted("Centre:");
            ImGui::TableSetColumnIndex(1);
            if (elemCombo("##c", m_centreSymbol, sizeof(m_centreSymbol)))
                m_paramsDirty = true;

            if (m_centreMode == 2)
            {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::AlignTextToFramePadding();
                ImGui::TextUnformatted("Neigh j:");
                ImGui::TableSetColumnIndex(1);
                if (elemCombo("##n1", m_neighSym1, sizeof(m_neighSym1)))
                    m_paramsDirty = true;

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::AlignTextToFramePadding();
                ImGui::TextUnformatted("Neigh k:");
                ImGui::TableSetColumnIndex(1);
                if (elemCombo("##n2", m_neighSym2, sizeof(m_neighSym2)))
                    m_paramsDirty = true;
            }
            ImGui::EndTable();
        }
    }

    // ── Post-processing ───────────────────────────────────────
    ImGui::SeparatorText("Post-processing");

    if (ImGui::BeginTable("##adf-post", 2, kRowFlags))
    {
        ImGui::TableSetupColumn("##lp", ImGuiTableColumnFlags_WidthFixed, labelColW);
        ImGui::TableSetupColumn("##rp", ImGuiTableColumnFlags_WidthStretch);

        // Smooth
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Smooth");
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::SliderInt("##smooth", &m_smoothPasses, 0, 8, "%d pass"))
            m_paramsDirty = true;
        ImGui::SetItemTooltip("3-point box smoothing passes (0 = off).\n2" "\xe2\x80\x93" "4 recommended.");

        // Normalise
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Normalise");
        ImGui::TableSetColumnIndex(1);
        if (ImGui::Checkbox("##norm", &m_normalize))
            m_paramsDirty = true;
        ImGui::SetItemTooltip("Scale histogram so tallest bin = 1.\nDisable to compare raw counts across runs.");

        ImGui::EndTable();
    }

    if (computing) ImGui::EndDisabled();

    // ── Display ───────────────────────────────────────────────
    ImGui::SeparatorText("Display");

    ImGui::Checkbox("Peak markers", &m_showPeaks);
    ImGui::Checkbox("Raw counts",   &m_showRawCounts);

    ImGui::Spacing();
    ImGui::TextDisabled("Angle range:");
    {
        const float btnColW = ImGui::CalcTextSize("All").x
                            + ImGui::GetStyle().FramePadding.x * 2.0f;
        if (ImGui::BeginTable("##adf-angrange", 2,
                ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_NoBordersInBody))
        {
            ImGui::TableSetupColumn("##ar-range", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("##ar-btn",   ImGuiTableColumnFlags_WidthFixed, btnColW);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            float rangeVals[2] = { m_xMin, m_xMax };
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::DragFloatRange2("##angrange", &rangeVals[0], &rangeVals[1],
                                        0.5f, 0.0f, 180.0f, "%.0f\xc2\xb0", "%.0f\xc2\xb0",
                                        ImGuiSliderFlags_AlwaysClamp))
            {
                m_xMin = rangeVals[0];
                m_xMax = rangeVals[1];
                if (m_xMax - m_xMin < 5.0f)
                    m_xMax = std::min(180.0f, m_xMin + 5.0f);
            }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip))
                ImGui::SetTooltip("Zoom the plot to an angle sub-range.\nDoes not affect the computation.");
            ImGui::TableSetColumnIndex(1);
            if (ImGui::Button("All"))
                { m_xMin = 0.0f; m_xMax = 180.0f; }
            ImGui::SetItemTooltip("Reset to full 0" "\xe2\x80\x93" "180\xc2\xb0 view.");
            ImGui::EndTable();
        }
    }
}

// ============================================================
// Result summary (below plot)
// ============================================================

void AngularDistributionAnalysisDialog::drawResultSummary()
{
    if (!m_hasResult) return;
    ImGui::Separator();
    ImGui::TextDisabled(
        "Atoms: %d  |  Centres: %d  |  Triplets: %lld  |  Bins: %d  (%.2f\xc2\xb0 / bin)",
        m_result.nAtoms, m_result.nCentreAtoms, m_result.nTriplets,
        m_result.binCount, m_result.binWidth);
    if (m_result.pbcUsed)
        ImGui::TextDisabled("PBC applied  |  Cutoff: %.2f \xc3\x85", m_result.rCutoff);
    else
        ImGui::TextDisabled("No PBC  |  Cutoff: %.2f \xc3\x85", m_result.rCutoff);
}

// ============================================================
// Main dialog
// ============================================================

void AngularDistributionAnalysisDialog::drawDialog(const Structure& structure)
{
    pollWorker();

    if (m_openRequested)
    {
        ImGui::OpenPopup("Angular Distribution Function");
        updateElementList(structure);
        m_openRequested = false;
    }

    ImVec2 centre = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(centre, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(920.0f, 660.0f), ImGuiCond_Appearing);
    ImGui::SetNextWindowSizeConstraints(ImVec2(640.0f, 420.0f), ImVec2(FLT_MAX, FLT_MAX));

    if (!ImGui::BeginPopupModal("Angular Distribution Function", nullptr,
                                ImGuiWindowFlags_NoScrollbar))
        return;

    const bool computing = m_isComputing.load();

    // ── Left panel ────────────────────────────────────────────
    const float leftW        = 300.0f;
    const float computeAreaH = ImGui::GetFrameHeightWithSpacing() * 2.4f;

    ImGui::BeginChild("##adf-left", ImVec2(leftW, -36.0f), true);
    {
        ImGui::SeparatorText("ADF Settings");

        // Scrollable settings (leaves fixed space at bottom for compute button)
        float settingsH = ImGui::GetContentRegionAvail().y
                        - computeAreaH
                        - ImGui::GetStyle().ItemSpacing.y * 2.0f;
        if (settingsH < 60.0f) settingsH = 60.0f;

        ImGui::BeginChild("##adf-scroll", ImVec2(-1.0f, settingsH), false);
        {
            drawSettings(structure);

            // Composition bars (in scrollable area, below settings)
            if (m_hasResult && !m_result.elementCounts.empty())
            {
                ImGui::Spacing();
                ImGui::SeparatorText("Composition");
                const int total = std::max(1, m_result.nAtoms);
                for (const auto& [sym, cnt] : m_result.elementCounts)
                {
                    float frac = (float)cnt / (float)total;
                    char lbl[32];
                    std::snprintf(lbl, sizeof(lbl), "%-2s  %.1f%%", sym.c_str(), frac * 100.0f);
                    ImGui::TextUnformatted(lbl);
                    float bw = ImGui::GetContentRegionAvail().x;
                    ImVec2 pos = ImGui::GetCursorScreenPos();
                    const bool lt2 = isLightTheme();
                    ImGui::GetWindowDrawList()->AddRectFilled(
                        pos, ImVec2(pos.x + bw * frac, pos.y + 6.0f),
                        lt2 ? IM_COL32(30,130,220,200) : IM_COL32(80,200,255,200));
                    ImGui::GetWindowDrawList()->AddRect(
                        pos, ImVec2(pos.x + bw, pos.y + 6.0f),
                        lt2 ? IM_COL32(120,130,150,180) : IM_COL32(90,100,120,180));
                    ImGui::Dummy(ImVec2(bw, 8.0f));
                }
            }
        }
        ImGui::EndChild();

        // ── Fixed compute area ─────────────────────────────────
        ImGui::Separator();

        const bool wantHighlight = m_paramsDirty && !computing;
        if (wantHighlight)
        {
            ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0.15f, 0.47f, 0.78f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.20f, 0.56f, 0.92f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(0.11f, 0.38f, 0.64f, 1.0f));
        }

        if (computing)
        {
            ImGui::BeginDisabled();
            float t = (float)ImGui::GetTime();
            char spinBtn[40];
            std::snprintf(spinBtn, sizeof(spinBtn), "Computing %c###computebtn",
                          "|/-\\"[(int)(t * 4.0f) & 3]);
            ImGui::Button(spinBtn, ImVec2(-FLT_MIN, 0.0f));
            ImGui::EndDisabled();
        }
        else
        {
            if (ImGui::Button("Compute ADF###computebtn", ImVec2(-FLT_MIN, 0.0f)))
                startCompute(structure);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip(m_paramsDirty
                    ? "Parameters changed \xe2\x80\x94 click to recompute."
                    : "Run the Angular Distribution Function analysis.");
        }

        if (wantHighlight)
            ImGui::PopStyleColor(3);

        // Status line
        if (computing)
            ImGui::TextDisabled("Running on background thread...");
        else if (m_hasResult)
        {
            if (m_result.valid)
            {
                ImGui::TextColored(themeStatusGood(), "Complete");
                ImGui::SameLine();
                ImGui::TextDisabled("%lld triplets", m_result.nTriplets);
            }
            else
            {
                ImGui::TextColored(ImVec4(0.90f, 0.30f, 0.30f, 1.0f), "Failed:");
                ImGui::SameLine();
                ImGui::TextDisabled("%s", m_result.message.c_str());
            }
        }
        else
        {
            ImGui::TextDisabled("Click Compute ADF to start.");
        }
    }
    ImGui::EndChild();

    ImGui::SameLine();

    // ── Right panel ───────────────────────────────────────────
    ImGui::BeginChild("##adf-right", ImVec2(-1.0f, -36.0f));
    {
        if (ImGui::BeginTabBar("##adf-tabs"))
        {
            if (ImGui::BeginTabItem("Distribution"))
            {
                drawPlot();
                drawResultSummary();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Peaks"))
            {
                drawPeakTable();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Coordination"))
            {
                drawCoordStats();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Guide"))
            {
                ImGui::Spacing();
                ImGui::SeparatorText("Method");
                ImGui::TextWrapped(
                    "The ADF g(\xce\xb8) accumulates the frequency of bond angles \xce\xb8 "
                    "formed by atom triplets (j\xe2\x80\x93i\xe2\x80\x93k), where i is the "
                    "centre atom and both j and k lie within the cutoff radius.");
                ImGui::Spacing();
                ImGui::SeparatorText("Reference Geometries");
                if (ImGui::BeginTable("##adf-guide", 3,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_SizingStretchProp))
                {
                    ImGui::TableSetupColumn("Angle",    ImGuiTableColumnFlags_WidthFixed, 68.0f);
                    ImGui::TableSetupColumn("Abbrev.",  ImGuiTableColumnFlags_WidthFixed, 52.0f);
                    ImGui::TableSetupColumn("Geometry / Example", ImGuiTableColumnFlags_WidthStretch);
                    ImGui::TableHeadersRow();
                    struct GRow { const char* ang; const char* abbr; const char* desc; };
                    static const GRow kRows[] = {
                        { "60\xc2\xb0",     "tri",  "Close-packed triangular layer" },
                        { "70.5\xc2\xb0",   "fcc",  "FCC / HCP second peak" },
                        { "90\xc2\xb0",     "oct",  "Octahedral (NaCl, MgO)" },
                        { "109.5\xc2\xb0",  "tet",  "Tetrahedral (diamond, ZnS)" },
                        { "120\xc2\xb0",    "trpl", "Trigonal planar (graphene)" },
                        { "150\xc2\xb0",    "bcc",  "BCC second-shell" },
                        { "180\xc2\xb0",    "lin",  "Linear (CO\xe2\x82\x82, triple bonds)" },
                    };
                    for (const auto& r : kRows)
                    {
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(r.ang);
                        ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(r.abbr);
                        ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(r.desc);
                    }
                    ImGui::EndTable();
                }
                ImGui::Spacing();
                ImGui::SeparatorText("Tips");
                ImGui::TextWrapped(
                    "Cutoff: set to the first-shell peak of the RDF for best results.\n\n"
                    "Bins: 90" "\xe2\x80\x93" "180 is a good default. More bins give finer "
                    "resolution but noisier histograms.\n\n"
                    "Smooth passes: 2" "\xe2\x80\x93" "4 passes reduce noise without hiding real "
                    "peaks. Set to 0 to see the raw histogram.\n\n"
                    "Atom Selection: use Element or Triplet mode to isolate the ADF "
                    "around a specific species in multi-element structures.\n\n"
                    "Angle range: drag the range control to zoom into a region of "
                    "interest without recomputing.");
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
    }
    ImGui::EndChild();

    // ── Footer ────────────────────────────────────────────────
    ImGui::Separator();
    if (ImGui::Button("Close", ImVec2(110.0f, 0.0f)))
        ImGui::CloseCurrentPopup();

    ImGui::EndPopup();
}