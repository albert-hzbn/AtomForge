#include "ShortRangeOrderDialog.h"
#include "ThemeUtils.h"
#include <imgui.h>
#include <imgui_internal.h>
#include <algorithm>
#include <cmath>
#include <cstdio>

// ---------------------------------------------------------------------------
// State destructor — must join the worker before the object is destroyed
// ---------------------------------------------------------------------------
ShortRangeOrderDialogState::~ShortRangeOrderDialogState()
{
    if (workerThread && workerThread->joinable())
        workerThread->join();
}

// ===========================================================================
// Internal helpers
// ===========================================================================
namespace {

// Diverging colormap: blue (-1) -> center (0) -> red (+1), theme-aware center
static ImVec4 sroAlphaColor(double alpha)
{
    float t = static_cast<float>(std::max(-1.0, std::min(1.0, alpha)));
    if (isLightTheme()) {
        // On light background, use stronger saturated colors and mid-gray center
        if (t < 0.0f) {
            float s = -t;
            return ImVec4(0.60f - s * 0.40f, 0.60f - s * 0.20f, 0.60f + s * 0.40f, 1.0f);
        }
        return ImVec4(0.60f + t * 0.40f, 0.60f - t * 0.40f, 0.60f - t * 0.50f, 1.0f);
    }
    if (t < 0.0f) {
        float s = -t;
        return ImVec4(1.0f - s * 0.55f, 1.0f - s * 0.45f, 1.0f, 1.0f);
    }
    return ImVec4(1.0f, 1.0f - t * 0.72f, 1.0f - t * 0.88f, 1.0f);
}

// Needle bar showing where alpha sits in [-1, +1]
static void drawAlphaBar(double alpha, float bw = 88.0f, float bh = 9.0f)
{
    float t        = static_cast<float>(std::max(-1.0, std::min(1.0, alpha)));
    float fraction = (t + 1.0f) * 0.5f;

    ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* dl = ImGui::GetWindowDrawList();

    // Gradient background
    dl->AddRectFilledMultiColor(
        pos, { pos.x + bw * 0.5f, pos.y + bh },
        IM_COL32(55, 95, 215, 170), IM_COL32(195, 195, 195, 110),
        IM_COL32(195, 195, 195, 110), IM_COL32(55, 95, 215, 170));
    dl->AddRectFilledMultiColor(
        { pos.x + bw * 0.5f, pos.y }, { pos.x + bw, pos.y + bh },
        IM_COL32(195, 195, 195, 110), IM_COL32(205, 55, 55, 170),
        IM_COL32(205, 55, 55, 170), IM_COL32(195, 195, 195, 110));
    dl->AddRect(pos, { pos.x + bw, pos.y + bh }, IM_COL32(110, 110, 110, 150), 1.5f);

    // Centre tick
    dl->AddLine({ pos.x + bw * 0.5f, pos.y }, { pos.x + bw * 0.5f, pos.y + bh },
                IM_COL32(130, 130, 130, 160), 1.0f);

    // Needle
    const ImU32 needleCol  = isLightTheme() ? IM_COL32(30,  30,  30,  230) : IM_COL32(255, 255, 255, 240);
    const ImU32 triangleCol= isLightTheme() ? IM_COL32(30,  30,  30,  210) : IM_COL32(255, 255, 255, 230);
    float mx = pos.x + fraction * bw;
    dl->AddLine({ mx, pos.y - 1 }, { mx, pos.y + bh + 1 },
                needleCol, 2.0f);
    // Small triangle pointer
    dl->AddTriangleFilled(
        { mx - 3, pos.y - 4 }, { mx + 3, pos.y - 4 }, { mx, pos.y - 1 },
        triangleCol);

    ImGui::Dummy({ bw, bh });
}

// Full-width colormap legend strip with ticks
static void drawColorLegend(float width)
{
    ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    float h = 13.0f;

    dl->AddRectFilledMultiColor(
        pos, { pos.x + width * 0.5f, pos.y + h },
        IM_COL32(55, 95, 215, 235), IM_COL32(225, 225, 225, 235),
        IM_COL32(225, 225, 225, 235), IM_COL32(55, 95, 215, 235));
    dl->AddRectFilledMultiColor(
        { pos.x + width * 0.5f, pos.y }, { pos.x + width, pos.y + h },
        IM_COL32(225, 225, 225, 235), IM_COL32(205, 55, 55, 235),
        IM_COL32(205, 55, 55, 235), IM_COL32(225, 225, 225, 235));
    dl->AddRect(pos, { pos.x + width, pos.y + h }, IM_COL32(110, 110, 110, 200), 2.0f);

    bool light = isLightTheme();
    ImU32 tickCol = light ? IM_COL32(30, 30, 30, 210) : IM_COL32(215, 215, 215, 210);
    const float tfs[] = { 0.0f, 0.25f, 0.5f, 0.75f, 1.0f };
    const char* tls[] = { "-1", "-0.5", "0", "+0.5", "+1" };
    for (int i = 0; i < 5; ++i) {
        float x = pos.x + tfs[i] * width;
        dl->AddLine({ x, pos.y + h }, { x, pos.y + h + 4 }, tickCol);
        dl->AddText({ x - 6.0f, pos.y + h + 5 }, tickCol, tls[i]);
    }
    ImGui::Dummy({ width, h + 20.0f });

    ImGui::TextColored({ 0.25f, 0.45f, 1.0f, 1.0f }, "Ordering");
    ImGui::SameLine(width * 0.38f);
    ImGui::TextDisabled("Random");
    ImGui::SameLine(width * 0.75f);
    ImGui::TextColored({ 1.0f, 0.28f, 0.28f, 1.0f }, "Clustering");
}

// Animated spinner string (call each frame while running)
static const char* sroSpinner()
{
    static const char* frames[] = { "|", "/", "-", "\\" };
    return frames[static_cast<int>(ImGui::GetTime() * 9.0) % 4];
}

// ---------------------------------------------------------------------------
// Warren-Cowley table (grouped by shell, unique table IDs per shell)
// ---------------------------------------------------------------------------
static void drawWarrenCowleyResults(const SroReport& report)
{
    ImGui::TextDisabled(
        u8"\u03b1\u1d62\u2c7c(s) = 1 \u2212 P(j|i,s) / c\u2c7c"
        "   [Cowley, Phys. Rev. 138 (1965) A1384]");
    ImGui::Spacing();

    int  prevShell = -1;
    bool tableOpen = false;

    for (const auto& e : report.entries)
    {
        if (e.shell != prevShell)
        {
            if (tableOpen) { ImGui::EndTable(); tableOpen = false; ImGui::Spacing(); }

            prevShell = e.shell;
            char hdr[80];
            std::snprintf(hdr, sizeof(hdr),
                          " Shell %d  \u2014  %.4f \u00c5##WCsh%d",
                          e.shell + 1, e.shellDistance, e.shell);

            if (ImGui::CollapsingHeader(hdr, ImGuiTreeNodeFlags_DefaultOpen))
            {
                char tblId[32];
                std::snprintf(tblId, sizeof(tblId), "SROWCtbl%d", e.shell);

                if (ImGui::BeginTable(tblId, 6,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_PadOuterX))
                {
                    ImGui::TableSetupColumn("Center",   ImGuiTableColumnFlags_WidthFixed, 62.0f);
                    ImGui::TableSetupColumn("Neighbor", ImGuiTableColumnFlags_WidthFixed, 62.0f);
                    ImGui::TableSetupColumn(u8"N\u1d62\u2c7c",  ImGuiTableColumnFlags_WidthFixed, 58.0f);
                    ImGui::TableSetupColumn(u8"Z\u1d62",        ImGuiTableColumnFlags_WidthFixed, 52.0f);
                    ImGui::TableSetupColumn("P(j|i)",           ImGuiTableColumnFlags_WidthFixed, 70.0f);
                    ImGui::TableSetupColumn(u8"\u03b1\u1d62\u2c7c", ImGuiTableColumnFlags_WidthStretch);
                    ImGui::TableHeadersRow();
                    tableOpen = true;
                }
            }
        }

        if (!tableOpen) continue;

        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(e.centerElement.c_str());
        ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(e.neighborElement.c_str());
        ImGui::TableSetColumnIndex(2); ImGui::Text("%d", e.nij);
        ImGui::TableSetColumnIndex(3); ImGui::Text("%d", e.zi);
        ImGui::TableSetColumnIndex(4); ImGui::Text("%.4f", e.pij);
        ImGui::TableSetColumnIndex(5);
        ImGui::TextColored(sroAlphaColor(e.alpha), "% .4f", e.alpha);
        ImGui::SameLine(0, 6);
        drawAlphaBar(e.alpha, 84.0f, 8.0f);
    }
    if (tableOpen) ImGui::EndTable();
}

// ---------------------------------------------------------------------------
// Rao-Curtin table (grouped by shell, symmetric pairs, clear columns)
// ---------------------------------------------------------------------------
static void drawRaoCurtinResults(const SroReportRC& report)
{
    ImGui::TextDisabled(
        u8"\u03b1\u1d62\u2c7c(s) = 1 \u2212 P\u1d62\u2c7c(s) / (2\u00b7c\u1d62\u00b7c\u2c7c)"
        "   [Rao & Curtin, Acta Mater. 226 (2022) 117621]");
    ImGui::Spacing();

    int  prevShell = -1;
    bool tableOpen = false;

    for (const auto& e : report.entries)
    {
        if (e.shell != prevShell)
        {
            if (tableOpen) { ImGui::EndTable(); tableOpen = false; ImGui::Spacing(); }

            prevShell = e.shell;
            double zs = (e.shell < static_cast<int>(report.shellCoordination.size()))
                        ? report.shellCoordination[e.shell] : 0.0;
            char hdr[96];
            std::snprintf(hdr, sizeof(hdr),
                          " Shell %d  \u2014  %.4f \u00c5"
                          "   (avg. coord. Z\u209b = %.2f)##RCsh%d",
                          e.shell + 1, e.shellDistance, zs, e.shell);

            if (ImGui::CollapsingHeader(hdr, ImGuiTreeNodeFlags_DefaultOpen))
            {
                char tblId[32];
                std::snprintf(tblId, sizeof(tblId), "SRORCtbl%d", e.shell);

                if (ImGui::BeginTable(tblId, 8,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_PadOuterX))
                {
                    ImGui::TableSetupColumn("Pair i\u2013j",           ImGuiTableColumnFlags_WidthFixed, 62.0f);
                    ImGui::TableSetupColumn(u8"c\u1d62",               ImGuiTableColumnFlags_WidthFixed, 58.0f);
                    ImGui::TableSetupColumn(u8"c\u2c7c",               ImGuiTableColumnFlags_WidthFixed, 58.0f);
                    ImGui::TableSetupColumn(u8"N\u1d62\u2c7c",         ImGuiTableColumnFlags_WidthFixed, 52.0f);
                    ImGui::TableSetupColumn(u8"N\u2c7c\u1d62",         ImGuiTableColumnFlags_WidthFixed, 52.0f);
                    ImGui::TableSetupColumn(u8"P\u1d62\u2c7c (obs.)",  ImGuiTableColumnFlags_WidthFixed, 78.0f);
                    ImGui::TableSetupColumn(u8"2c\u1d62c\u2c7c (rand)",ImGuiTableColumnFlags_WidthFixed, 88.0f);
                    ImGui::TableSetupColumn(u8"\u03b1\u1d62\u2c7c",    ImGuiTableColumnFlags_WidthStretch);
                    ImGui::TableHeadersRow();
                    tableOpen = true;
                }
            }
        }

        if (!tableOpen) continue;

        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%s\u2013%s", e.elem_i.c_str(), e.elem_j.c_str());
        ImGui::TableSetColumnIndex(1); ImGui::Text("%.4f", e.ci);
        ImGui::TableSetColumnIndex(2); ImGui::Text("%.4f", e.cj);
        ImGui::TableSetColumnIndex(3); ImGui::Text("%d", e.n_ij);
        ImGui::TableSetColumnIndex(4); ImGui::Text("%d", e.n_ji);
        ImGui::TableSetColumnIndex(5); ImGui::Text("%.5f", e.p_ij);
        ImGui::TableSetColumnIndex(6); ImGui::Text("%.5f", 2.0 * e.ci * e.cj);
        ImGui::TableSetColumnIndex(7);
        ImGui::TextColored(sroAlphaColor(e.alpha), "% .4f", e.alpha);
        ImGui::SameLine(0, 6);
        drawAlphaBar(e.alpha, 84.0f, 8.0f);
    }
    if (tableOpen) ImGui::EndTable();
}

} // anonymous namespace

// ===========================================================================
// Public API
// ===========================================================================

void drawShortRangeOrderMenuItem(bool enabled, ShortRangeOrderDialogState& state)
{
    if (!enabled) ImGui::BeginDisabled();
    if (ImGui::MenuItem("Short-Range Order (SRO)"))
        state.openRequested = true;
    if (!enabled) ImGui::EndDisabled();
}

void drawShortRangeOrderDialog(ShortRangeOrderDialogState& state,
                               const Structure& structure)
{
    // ------------------------------------------------------------------
    // 1. Poll background thread result (every frame, outside the popup)
    // ------------------------------------------------------------------
    if (state.computeCompleted.load(std::memory_order_acquire))
    {
        state.computeCompleted.store(false, std::memory_order_relaxed);
        state.results    = std::move(state.workerResult);
        state.hasResults = state.results.success;
        state.selectedShell = 0;
        if (state.workerThread && state.workerThread->joinable())
            state.workerThread->join();
    }

    // ------------------------------------------------------------------
    // 2. Open popup if requested
    // ------------------------------------------------------------------
    const char* kTitle = "Short-Range Order Analysis";
    if (state.openRequested)
    {
        ImGui::OpenPopup(kTitle);
        state.openRequested = false;
        state.isOpen        = true;
    }

    ImGui::SetNextWindowSizeConstraints({ 720, 460 }, { 1300, 980 });
    ImGui::SetNextWindowSize({ 920, 680 }, ImGuiCond_FirstUseEver);
    bool open = true;
    if (!ImGui::BeginPopupModal(kTitle, &open)) return;
    if (!open)
    {
        ImGui::CloseCurrentPopup();
        state.isOpen = false;
        ImGui::EndPopup();
        return;
    }

    // ------------------------------------------------------------------
    // 3. Theme colours
    // ------------------------------------------------------------------
    const bool   light    = isLightTheme();
    const ImVec4 kAccent  = light ? ImVec4(0.0f,0.44f,0.78f,1.0f) : ImVec4(0.35f,0.75f,1.0f,1.0f);
    const ImVec4 kGood    = light ? ImVec4(0.0f,0.52f,0.12f,1.0f) : ImVec4(0.18f,0.88f,0.28f,1.0f);
    const ImVec4 kWarn    = ImVec4(1.0f, 0.60f, 0.10f, 1.0f);
    const ImVec4 kMuted   = ImVec4(0.52f, 0.52f, 0.52f, 1.0f);

    const float  footerH  = ImGui::GetFrameHeightWithSpacing() + 6.0f;
    const float  totalH   = ImGui::GetContentRegionAvail().y - footerH;
    const float  leftW    = 218.0f;

    bool isComputing  = state.isComputing.load(std::memory_order_relaxed);
    bool hasStructure = !structure.atoms.empty() && structure.hasUnitCell;

    // ==================================================================
    // LEFT PANEL — settings + status + composition
    // ==================================================================
    ImGui::BeginChild("##SROLeft", { leftW, totalH }, true);

    ImGui::TextColored(kAccent, "Settings");
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::PushItemWidth(-1.0f);

    ImGui::TextDisabled("Neighbor shells");
    ImGui::BeginDisabled(isComputing);
    ImGui::SliderInt("##SROShells", &state.nShells, 1, 8);
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
        ImGui::SetTooltip("Number of nearest-neighbor shells to include.\n"
                          "Higher values are slower on large structures.");
    ImGui::EndDisabled();

    ImGui::Spacing();
    ImGui::TextDisabled("Shell tolerance (\xc3\x85)");
    ImGui::BeginDisabled(isComputing);
    float tol_f = static_cast<float>(state.shellTolerance);
    if (ImGui::SliderFloat("##SROTol", &tol_f, 0.001f, 0.20f, "%.3f"))
        state.shellTolerance = static_cast<double>(tol_f);
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
        ImGui::SetTooltip("Distance window for merging atoms into the same shell.\n"
                          "Increase if shells appear split.");
    ImGui::EndDisabled();

    ImGui::PopItemWidth();
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // ---- Run / Cancel button -------------------------------------------
    bool canRun = hasStructure && !isComputing;
    if (!canRun) ImGui::BeginDisabled();
    ImGui::PushStyleColor(ImGuiCol_Button,        { 0.13f, 0.48f, 0.82f, 1.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, { 0.22f, 0.62f, 1.0f,  1.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonActive,  { 0.08f, 0.33f, 0.62f, 1.0f });
    if (ImGui::Button("  Run Analysis  ##SRO", { -1.0f, 32.0f }))
    {
        // Join any previous (already-finished) thread
        if (state.workerThread && state.workerThread->joinable())
            state.workerThread->join();

        // Snapshot structure for the worker
        state.workerStructure = structure;
        state.isComputing.store(true, std::memory_order_relaxed);
        state.computeCompleted.store(false, std::memory_order_relaxed);
        state.hasResults = false;

        int    nsh = state.nShells;
        double tol = state.shellTolerance;

        state.workerThread = std::make_unique<std::thread>([&state, nsh, tol]()
        {
            SroAnalysisResult res = analyzeSRO(state.workerStructure, nsh, tol);
            state.workerResult = std::move(res);
            state.isComputing.store(false, std::memory_order_release);
            state.computeCompleted.store(true, std::memory_order_release);
        });
    }
    ImGui::PopStyleColor(3);
    if (!canRun) ImGui::EndDisabled();

    ImGui::Spacing();

    // ---- Status area ---------------------------------------------------
    if (isComputing)
    {
        ImGui::TextColored(kWarn, "%s  Computing...", sroSpinner());
        ImGui::Spacing();
        ImGui::TextDisabled("Analysis running in background.\nUI remains responsive.");
    }
    else if (!hasStructure)
    {
        ImGui::TextColored(kWarn, "No unit-cell structure\nloaded.");
    }
    else if (state.hasResults)
    {
        ImGui::TextColored(kGood, "Complete");
        ImGui::Spacing();

        ImGui::TextDisabled("Atoms");    ImGui::SameLine(70.0f);
        ImGui::Text("%d", state.results.warrenCowley.nAtoms);
        ImGui::TextDisabled("Shells");   ImGui::SameLine(70.0f);
        ImGui::Text("%d", state.results.warrenCowley.nShells);

        // Per-element fraction bars
        ImGui::Spacing();
        ImGui::TextColored(kAccent, "Composition");
        ImGui::Separator();
        ImGui::Spacing();

        // Label column: "Fe 62.50%" — measure max label width first
        const float labelColW = 80.0f;  // fixed right column for "Xx 100.00%"
        const float barW      = ImGui::GetContentRegionAvail().x - labelColW - 10.0f;
        ImDrawList* dl = ImGui::GetWindowDrawList();

        for (const auto& [elem, c] : state.results.warrenCowley.composition)
        {
            float frac = static_cast<float>(c);

            // Element label + percentage (left-aligned, fixed width column)
            char label[24];
            std::snprintf(label, sizeof(label), "%-2s  %.1f%%", elem.c_str(), c * 100.0);
            ImGui::TextUnformatted(label);

            // Bar on the same line, remaining width
            ImGui::SameLine(0, 6);
            ImVec2 bpos = ImGui::GetCursorScreenPos();
            float barH = ImGui::GetTextLineHeight() * 0.75f;
            float by   = bpos.y + (ImGui::GetTextLineHeight() - barH) * 0.5f;

            dl->AddRectFilled(
                { bpos.x, by }, { bpos.x + barW * frac, by + barH },
                light ? IM_COL32(50,130,215,160) : IM_COL32(55,175,255,140), 2.0f);
            dl->AddRect(
                { bpos.x, by }, { bpos.x + barW, by + barH },
                light ? IM_COL32(140,148,170,140) : IM_COL32(75,85,105,180), 1.5f);

            ImGui::Dummy({ barW, ImGui::GetTextLineHeight() });
            ImGui::Spacing();
        }
    }
    else if (!state.results.message.empty())
    {
        ImGui::TextColored(kWarn, "%s", state.results.message.c_str());
    }

    ImGui::EndChild(); // left panel
    ImGui::SameLine();

    // ==================================================================
    // RIGHT PANEL — results tabs
    // ==================================================================
    ImGui::BeginChild("##SRORight", { 0.0f, totalH }, false);

    if (!state.hasResults && !isComputing)
    {
        float cy = totalH * 0.38f;
        float cx = (ImGui::GetContentRegionAvail().x - 280.0f) * 0.5f;
        ImGui::SetCursorPos({ std::max(cx, 0.0f), cy });
        ImGui::TextColored(kMuted, "Configure parameters and click  Run Analysis.");
        ImGui::SetCursorPos({ std::max(cx + 30.0f, 0.0f),
                              cy + ImGui::GetTextLineHeightWithSpacing() * 1.2f });
        ImGui::TextColored(kMuted, "A periodic unit cell is required.");
    }
    else if (isComputing && !state.hasResults)
    {
        float cy = totalH * 0.40f;
        float cx = (ImGui::GetContentRegionAvail().x - 200.0f) * 0.5f;
        ImGui::SetCursorPos({ std::max(cx, 0.0f), cy });
        ImGui::TextColored(kWarn, "%s  SRO analysis running...", sroSpinner());
    }
    else if (state.hasResults)
    {
        if (ImGui::BeginTabBar("SROTabs"))
        {
            // ----------------------------------------------------------
            // Tab 1: Warren-Cowley
            // ----------------------------------------------------------
            if (ImGui::BeginTabItem(" Warren-Cowley "))
            {
                state.selectedMethod = 0;
                ImGui::Spacing();
                ImGui::BeginChild("##SROWCScroll", { 0,0 }, false,
                                  ImGuiWindowFlags_HorizontalScrollbar);
                drawWarrenCowleyResults(state.results.warrenCowley);
                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
                ImGui::TextColored(kMuted, u8"\u03b1\u1d62\u2c7c colormap:");
                ImGui::SameLine(0, 6);
                drawColorLegend(260.0f);
                ImGui::EndChild();
                ImGui::EndTabItem();
            }

            // ----------------------------------------------------------
            // Tab 2: Rao-Curtin
            // ----------------------------------------------------------
            if (ImGui::BeginTabItem(" Rao-Curtin (2022) "))
            {
                state.selectedMethod = 1;
                ImGui::Spacing();
                ImGui::BeginChild("##SRORCScroll", { 0,0 }, false,
                                  ImGuiWindowFlags_HorizontalScrollbar);
                drawRaoCurtinResults(state.results.raoCurtin);
                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
                ImGui::TextColored(kMuted, u8"\u03b1\u1d62\u2c7c colormap:");
                ImGui::SameLine(0, 6);
                drawColorLegend(260.0f);
                ImGui::EndChild();
                ImGui::EndTabItem();
            }

            // ----------------------------------------------------------
            // Tab 3: Shells overview
            // ----------------------------------------------------------
            if (ImGui::BeginTabItem(" Shells "))
            {
                ImGui::Spacing();
                const auto& wc = state.results.warrenCowley;
                const auto& rc = state.results.raoCurtin;

                if (ImGui::BeginTable("SROShellTbl", 4,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_SizingFixedFit))
                {
                    ImGui::TableSetupColumn("Shell",              ImGuiTableColumnFlags_WidthFixed, 55.0f);
                    ImGui::TableSetupColumn("Distance (\xc3\x85)",ImGuiTableColumnFlags_WidthFixed, 115.0f);
                    ImGui::TableSetupColumn("Avg. coord. Z",      ImGuiTableColumnFlags_WidthFixed, 115.0f);
                    ImGui::TableSetupColumn("WC pairs analyzed",  ImGuiTableColumnFlags_WidthStretch);
                    ImGui::TableHeadersRow();

                    for (int s = 0; s < static_cast<int>(wc.shellDistances.size()); ++s)
                    {
                        double zs = (s < static_cast<int>(rc.shellCoordination.size()))
                                    ? rc.shellCoordination[s] : -1.0;

                        // Count WC pairs in this shell
                        int wcPairs = 0;
                        for (const auto& e : wc.entries)
                            if (e.shell == s) ++wcPairs;

                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0); ImGui::Text("%d", s + 1);
                        ImGui::TableSetColumnIndex(1); ImGui::Text("%.4f", wc.shellDistances[s]);
                        ImGui::TableSetColumnIndex(2);
                        if (zs >= 0.0) ImGui::Text("%.3f", zs);
                        else           ImGui::TextDisabled("\xe2\x80\x94");
                        ImGui::TableSetColumnIndex(3); ImGui::Text("%d", wcPairs);
                    }
                    ImGui::EndTable();
                }
                ImGui::EndTabItem();
            }

            // ----------------------------------------------------------
            // Tab 4: Interpretation guide
            // ----------------------------------------------------------
            if (ImGui::BeginTabItem(" Guide "))
            {
                ImGui::Spacing();
                ImGui::TextColored(kAccent, "Interpreting SRO Parameters");
                ImGui::Separator();
                ImGui::Spacing();

                // Color legend rows
                struct LegendEntry { ImVec4 col; const char* label; const char* desc; };
                const LegendEntry rows[] = {
                    { {0.25f,0.45f,1.0f,1.0f}, u8"\u03b1 < 0  (blue)",
                      "Ordering tendency. Atoms i and j prefer unlike neighbors "
                      "(attractive effective pair interaction, heterogeneous bonding)." },
                    { {0.65f,0.65f,0.65f,1.0f}, u8"\u03b1 \u2248 0  (grey)",
                      "Near-random distribution. Ideal solid solution — no short-range preference." },
                    { {1.0f,0.28f,0.28f,1.0f}, u8"\u03b1 > 0  (red)",
                      "Clustering tendency. Like atoms prefer like neighbors "
                      "(repulsive effective pair interaction, homogeneous bonding)." },
                };
                for (const auto& row : rows)
                {
                    ImVec2 p = ImGui::GetCursorScreenPos();
                    ImGui::GetWindowDrawList()->AddRectFilled(
                        p, { p.x + 14, p.y + 14 },
                        IM_COL32((int)(row.col.x*255),(int)(row.col.y*255),
                                 (int)(row.col.z*255), 220), 3.0f);
                    ImGui::Dummy({ 14, 14 });
                    ImGui::SameLine(0, 6);
                    ImGui::TextColored(row.col, "%s", row.label);
                    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 20.0f);
                    ImGui::TextWrapped("%s", row.desc);
                    ImGui::Spacing();
                }

                ImGui::Separator();
                ImGui::Spacing();
                ImGui::TextColored(kAccent, "Method Comparison");
                ImGui::Separator();
                ImGui::Spacing();

                ImGui::TextWrapped(
                    "Warren-Cowley  —  Directional conditional probability P(j|i,s): "
                    "counts j-type neighbors around a central i atom. "
                    "alpha_ij and alpha_ji can differ for each pair type.");
                ImGui::Spacing();
                ImGui::TextWrapped(
                    "Rao-Curtin (2022)  —  Symmetric joint probability P_ij(s): "
                    "combines both directions so alpha_ij = alpha_ji by construction. "
                    "Better aligned with free-energy perturbation theory and more "
                    "consistent for multicomponent alloys.");
                ImGui::Spacing();
                ImGui::TextDisabled(
                    "Both methods are accurate for weak-to-moderate interactions "
                    "satisfying  0.5 < exp(-V_ij/kT) < 1.5.");
                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }
    }

    ImGui::EndChild(); // right panel

    // ==================================================================
    // FOOTER
    // ==================================================================
    ImGui::Separator();
    ImGui::Spacing();

    if (isComputing)
    {
        ImGui::TextColored(kWarn, "%s  Running SRO analysis in background...", sroSpinner());
        ImGui::SameLine();
    }

    float closeX = ImGui::GetContentRegionAvail().x - 104.0f;
    if (closeX > 0.0f) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + closeX);
    if (ImGui::Button("Close##SROClose", { 100.0f, 0.0f }))
    {
        ImGui::CloseCurrentPopup();
        state.isOpen = false;
    }

    ImGui::EndPopup();
}
