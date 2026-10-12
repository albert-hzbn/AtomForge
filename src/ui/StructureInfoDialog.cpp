#include "ui/ResponsiveLayout.h"
#include "ui/DialogLayout.h"
#include "StructureInfoDialog.h"

#include "ElementData.h"
#include "math/StructureMath.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdarg>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

namespace
{
constexpr float kRadToDeg = 57.2957795130823208768f;

struct SymmetryInfo
{
    bool attempted = false;
    bool success = false;
    int spaceGroupNumber = 0;
    std::string internationalSymbol;
    std::string hallSymbol;
    std::string pointGroup;
    std::string error;
};

std::size_t hashCombine(std::size_t seed, std::size_t value)
{
    return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

std::size_t structureFingerprint(const Structure& structure)
{
    std::size_t fp = 0;
    fp = hashCombine(fp, static_cast<std::size_t>(structure.atoms.size()));
    fp = hashCombine(fp, structure.hasUnitCell ? 1u : 0u);

    if (structure.hasUnitCell)
    {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                fp = hashCombine(fp, std::hash<double>{}(structure.cellVectors[r][c]));
    }

    const std::size_t sampleCount = std::min<std::size_t>(structure.atoms.size(), 64u);
    for (std::size_t i = 0; i < sampleCount; ++i)
    {
        const AtomSite& atom = structure.atoms[i];
        fp = hashCombine(fp, std::hash<int>{}(atom.atomicNumber));
        fp = hashCombine(fp, std::hash<double>{}(atom.x));
        fp = hashCombine(fp, std::hash<double>{}(atom.y));
        fp = hashCombine(fp, std::hash<double>{}(atom.z));
    }

    return fp;
}

std::map<int, int> buildElementCounts(const Structure& structure)
{
    std::map<int, int> counts;
    for (const AtomSite& atom : structure.atoms)
        counts[atom.atomicNumber]++;
    return counts;
}

std::string buildFormula(const std::map<int, int>& counts)
{
    std::ostringstream out;
    bool first = true;
    for (const auto& [z, count] : counts)
    {
        if (z <= 0 || count <= 0)
            continue;

        if (!first)
            out << " ";
        first = false;
        out << elementSymbol(z);
        if (count > 1)
            out << count;
    }

    return out.str();
}

bool buildLatticeMatrix(const Structure& structure,
                        glm::mat3& lattice,
                        glm::vec3& origin)
{
    if (!structure.hasUnitCell)
        return false;

    lattice = makeCellMatrix(structure);
    origin = glm::vec3((float)structure.cellOffset[0],
                       (float)structure.cellOffset[1],
                       (float)structure.cellOffset[2]);

    return std::abs(glm::determinant(lattice)) > 1e-10f;
}

SymmetryInfo analyzeSymmetryWithSpglib(const Structure& structure)
{
    SymmetryInfo info;
    info.attempted = true;

#ifndef ATOMS_ENABLE_SPGLIB
    info.error = "spglib not available at build time (install libspglib-dev).";
    return info;
#else
    if (!structure.hasUnitCell)
    {
        info.error = "No unit cell available.";
        return info;
    }
    if (structure.atoms.empty())
    {
        info.error = "No atoms available.";
        return info;
    }

    glm::mat3 lattice(1.0f);
    glm::vec3 origin(0.0f);
    if (!buildLatticeMatrix(structure, lattice, origin))
    {
        info.error = "Invalid unit cell (zero volume).";
        return info;
    }

    glm::mat3 invLattice = glm::inverse(lattice);

    double latticeData[3][3] = {
        {structure.cellVectors[0][0], structure.cellVectors[0][1], structure.cellVectors[0][2]},
        {structure.cellVectors[1][0], structure.cellVectors[1][1], structure.cellVectors[1][2]},
        {structure.cellVectors[2][0], structure.cellVectors[2][1], structure.cellVectors[2][2]}
    };

    std::vector<std::array<double, 3>> positions(structure.atoms.size());
    std::vector<int> types(structure.atoms.size(), 0);

    for (size_t index = 0; index < structure.atoms.size(); ++index)
    {
        const AtomSite& atom = structure.atoms[index];
        glm::vec3 cart((float)atom.x, (float)atom.y, (float)atom.z);
        glm::vec3 frac = invLattice * (cart - origin);

        frac.x -= std::floor(frac.x);
        frac.y -= std::floor(frac.y);
        frac.z -= std::floor(frac.z);

        positions[index][0] = frac.x;
        positions[index][1] = frac.y;
        positions[index][2] = frac.z;
        types[index] = (atom.atomicNumber > 0) ? atom.atomicNumber : (int)index + 1000;
    }

    const SpglibDataset* dataset = spg_get_dataset(
        latticeData,
        reinterpret_cast<double (*)[3]>(positions.data()),
        types.data(),
        (int)positions.size(),
        1e-5);

    if (!dataset)
    {
        info.error = "spglib could not identify symmetry (try larger symprec).";
        return info;
    }

    info.success = true;
    info.spaceGroupNumber = dataset->spacegroup_number;
    info.internationalSymbol = dataset->international_symbol;
    info.hallSymbol = dataset->hall_symbol;
    info.pointGroup = dataset->pointgroup_symbol;

    spg_free_dataset(const_cast<SpglibDataset*>(dataset));
    return info;
#endif
}

using Rows = std::vector<std::pair<std::string, std::string>>;

std::string formatText(const char* format, ...)
{
    char text[512];
    va_list args;
    va_start(args, format);
    std::vsnprintf(text, sizeof(text), format, args);
    va_end(args);
    return text;
}

// The formula reduced by the common divisor of the counts, e.g. Al2Cu3.
std::string reducedFormula(const std::map<int, int>& counts)
{
    int divisor = 0;
    for (const auto& [z, count] : counts)
        if (z > 0 && count > 0) divisor = std::gcd(divisor, count);
    std::string formula;
    for (const auto& [z, count] : counts)
    {
        if (z <= 0 || count <= 0) continue;
        formula += elementSymbol(z);
        if (count / divisor > 1) formula += std::to_string(count / divisor);
    }
    return formula;
}

float safeAngleDeg(const glm::vec3& u, const glm::vec3& v)
{
    const float lu = glm::length(u);
    const float lv = glm::length(v);
    if (lu <= 1e-10f || lv <= 1e-10f)
        return 0.0f;
    const float cosang = std::max(-1.0f, std::min(1.0f, glm::dot(u, v) / (lu * lv)));
    return kRadToDeg * std::acos(cosang);
}

Rows summaryRows(const Structure& structure, const StructureInfoComputedData* data)
{
    Rows rows;
    rows.push_back({"Atoms", std::to_string(structure.atoms.size())});
    if (data)
    {
        const std::string reduced = reducedFormula(data->elementCounts);
        std::string full = data->formula;
        full.erase(std::remove(full.begin(), full.end(), ' '), full.end());
        rows.push_back({"Formula", reduced == full ? reduced : reduced + "   (cell: " + full + ")"});
    }
    else
        rows.push_back({"Formula", "computing..."});

    if (!structure.hasUnitCell)
    {
        rows.push_back({"Cell", "none (an isolated molecule or cluster)"});
        return rows;
    }
    const glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
    const glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
    const glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);
    const float volume = std::abs(glm::dot(a, glm::cross(b, c)));
    rows.push_back({"Cell lengths", formatText("a = %.4f \xc3\x85    b = %.4f \xc3\x85    c = %.4f \xc3\x85", glm::length(a), glm::length(b), glm::length(c))});
    rows.push_back({"Cell angles", formatText("\xce\xb1 = %.2f\xc2\xb0    \xce\xb2 = %.2f\xc2\xb0    \xce\xb3 = %.2f\xc2\xb0", safeAngleDeg(b, c), safeAngleDeg(a, c), safeAngleDeg(a, b))});
    rows.push_back({"Volume", formatText("%.3f \xc3\x85\xc2\xb3", volume)});
    if (data && volume > 1e-6f)
    {
        double mass = 0.0;  // atomic mass units
        for (const auto& [z, count] : data->elementCounts)
            if (z > 0) mass += count * elementAtomicMass(z);
        constexpr double kAmuPerCubicAngstromToGramsPerCm3 = 1.66053906660;
        rows.push_back({"Density", formatText("%.4f g/cm\xc2\xb3    (%.5f atoms/\xc3\x85\xc2\xb3)",
                                              mass * kAmuPerCubicAngstromToGramsPerCm3 / volume,
                                              static_cast<double>(structure.atoms.size()) / volume)});
    }
    return rows;
}

Rows symmetryRows(const StructureInfoComputedData* data)
{
    if (!data)
        return {{"Space group", "computing..."}};
    if (!data->symmetrySuccess)
        return {{"Space group", "unavailable (" + data->symmetryError + ")"}};
    return {{"Space group", formatText("%s (No. %d)", data->symmetryInternationalSymbol.c_str(), data->symmetrySpaceGroupNumber)},
            {"Point group", data->symmetryPointGroup},
            {"Hall symbol", data->symmetryHallSymbol}};
}

void drawRows(const char* id, const Rows& rows)
{
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingFixedFit))
        return;
    ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 7.5f);
    ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch);
    for (const auto& [label, value] : rows)
    {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("%s", label.c_str());
        ImGui::TableSetColumnIndex(1);
        ImGui::TextUnformatted(value.c_str());
    }
    ImGui::EndTable();
}

void drawElementCountsTable(const std::map<int, int>& counts, std::size_t total)
{
    if (!ImGui::BeginTable("##element-counts", 4, ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchSame))
        return;

    ImGui::TableSetupColumn("Element");
    ImGui::TableSetupColumn("Z");
    ImGui::TableSetupColumn("Count");
    ImGui::TableSetupColumn("Fraction");
    ImGui::TableHeadersRow();

    for (const auto& [z, count] : counts)
    {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%s", elementSymbol(z));
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%d", z);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%d", count);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.2f %%", total ? 100.0 * count / static_cast<double>(total) : 0.0);
    }

    ImGui::EndTable();
}

StructureInfoComputedData computeStructureInfo(const Structure& structure)
{
    StructureInfoComputedData data;
    data.ready = true;
    data.atomCount = structure.atoms.size();
    data.elementCounts = buildElementCounts(structure);
    data.formula = buildFormula(data.elementCounts);

    const SymmetryInfo symmetry = analyzeSymmetryWithSpglib(structure);
    data.symmetryAttempted = symmetry.attempted;
    data.symmetrySuccess = symmetry.success;
    data.symmetrySpaceGroupNumber = symmetry.spaceGroupNumber;
    data.symmetryInternationalSymbol = symmetry.internationalSymbol;
    data.symmetryHallSymbol = symmetry.hallSymbol;
    data.symmetryPointGroup = symmetry.pointGroup;
    data.symmetryError = symmetry.error;
    return data;
}

void drawAtomicPositionsTable(const Structure& structure, bool hasValidLattice)
{
    if (!ImGui::BeginTable("##position-table", 8,
                           ImGuiTableFlags_BordersInnerH |
                           ImGuiTableFlags_BordersOuter |
                           ImGuiTableFlags_RowBg |
                           ImGuiTableFlags_ScrollY |
                           ImGuiTableFlags_SizingStretchProp))
    {
        return;
    }

    ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthStretch, 0.6f);
    ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthStretch, 0.6f);
    ImGui::TableSetupColumn("x (\xc3\x85)", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupColumn("y (\xc3\x85)", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupColumn("z (\xc3\x85)", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupColumn("u", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupColumn("v", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupColumn("w", ImGuiTableColumnFlags_WidthStretch, 1.0f);
    ImGui::TableSetupScrollFreeze(0, 1);
    ImGui::TableHeadersRow();

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const AtomSite& atom = structure.atoms[i];
        const glm::vec3 cart((float)atom.x, (float)atom.y, (float)atom.z);
        glm::vec3 frac(0.0f);
        if (hasValidLattice)
            tryCartesianToFractional(structure, cart, frac);

        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%d", i);
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%s", atom.symbol.c_str());
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.5f", atom.x);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.5f", atom.y);
        ImGui::TableSetColumnIndex(4);
        ImGui::Text("%.5f", atom.z);
        ImGui::TableSetColumnIndex(5);
        ImGui::Text("%.5f", frac.x);
        ImGui::TableSetColumnIndex(6);
        ImGui::Text("%.5f", frac.y);
        ImGui::TableSetColumnIndex(7);
        ImGui::Text("%.5f", frac.z);
    }

    ImGui::EndTable();
}
} // namespace

void drawStructureInfoDialog(StructureInfoDialogState& state,
                             bool requestOpen,
                             const Structure& structure)
{
    const std::size_t currentFingerprint = structureFingerprint(structure);
    const bool dialogRequestedOrOpen = requestOpen || state.openRequested || ImGui::IsPopupOpen("Structure Info");

    if (requestOpen)
        state.openRequested = true;

    if (state.workerRunning)
    {
        const auto status = state.workerFuture.wait_for(std::chrono::milliseconds(0));
        if (status == std::future_status::ready)
        {
            state.cachedData = state.workerFuture.get();
            state.cachedFingerprint = state.pendingFingerprint;
            state.workerRunning = false;
        }
    }

    if (dialogRequestedOrOpen && !state.workerRunning && state.cachedFingerprint != currentFingerprint)
    {
        const Structure structureSnapshot = structure;
        state.pendingFingerprint = currentFingerprint;
        state.workerFuture = std::async(std::launch::async, [structureSnapshot]() {
            return computeStructureInfo(structureSnapshot);
        });
        state.workerRunning = true;
    }

    if (state.openRequested)
    {
        ImGui::OpenPopup("Structure Info");
        state.openRequested = false;
    }

    // A fixed width and a height that follows the content.
    const float width = responsive::dp(760);
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0), ImVec2(width, FLT_MAX));
    bool popupOpen = true;
    if (responsive::beginModal("Structure Info", &popupOpen, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize))
    {
        const bool hasComputedData = state.cachedData.ready && state.cachedFingerprint == currentFingerprint;
        const StructureInfoComputedData* data = hasComputedData ? &state.cachedData : nullptr;
        const Rows summary = summaryRows(structure, data);
        const Rows symmetry = symmetryRows(data);

        dialogLayout::section("Summary");
        drawRows("##summary", summary);

        dialogLayout::section("Composition");
        if (data)
            drawElementCountsTable(data->elementCounts, structure.atoms.size());
        else
            ImGui::TextDisabled("Counting elements...");

        dialogLayout::section("Symmetry");
        drawRows("##symmetry", symmetry);

        dialogLayout::section("Atomic Positions");
        constexpr std::size_t kMaxAtomsForPositionsTable = 500;
        if (structure.atoms.empty())
            ImGui::TextDisabled("No atoms.");
        else if (structure.atoms.size() > kMaxAtomsForPositionsTable)
        {
            dialogLayout::notef(dialogLayout::Note::Info,
                "The position table lists up to %d atoms; this structure has %zu. Save it (File > Save As) as an XYZ or CIF file to read every position.",
                static_cast<int>(kMaxAtomsForPositionsTable), structure.atoms.size());
        }
        else
        {
            glm::mat3 lattice(1.0f);
            glm::vec3 origin(0.0f);
            bool hasValidLattice = buildLatticeMatrix(structure, lattice, origin);
            const float rows = static_cast<float>(std::min<std::size_t>(structure.atoms.size(), 10) + 1);
            const float rowHeight = ImGui::GetTextLineHeight() + 2 * ImGui::GetStyle().CellPadding.y;
            responsive::beginChild("##positions", ImVec2(0, rows * rowHeight + 2.0f), false);
            drawAtomicPositionsTable(structure, hasValidLattice);
            ImGui::EndChild();
        }

        const int pressed = dialogLayout::footer("Close", "Copy Summary");
        if (pressed == 2)
        {
            std::string text;
            for (const Rows* group : {&summary, &symmetry})
                for (const auto& [label, value] : *group)
                    text += label + ": " + value + "\n";
            if (data)
                for (const auto& [z, count] : data->elementCounts)
                    text += std::string(elementSymbol(z)) + ": " + std::to_string(count) + " atoms\n";
            ImGui::SetClipboardText(text.c_str());
        }
        if (pressed == 1)
            popupOpen = false;

        ImGui::EndPopup();
    }

    if (!popupOpen)
        ImGui::CloseCurrentPopup();
}