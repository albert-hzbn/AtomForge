#include "StructureInfoDialog.h"

#include "ElementData.h"
#include "math/StructureMath.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <map>
#include <sstream>
#include <string>
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

void drawLatticeInfo(const Structure& structure)
{
    glm::vec3 a((float)structure.cellVectors[0][0],
                (float)structure.cellVectors[0][1],
                (float)structure.cellVectors[0][2]);
    glm::vec3 b((float)structure.cellVectors[1][0],
                (float)structure.cellVectors[1][1],
                (float)structure.cellVectors[1][2]);
    glm::vec3 c((float)structure.cellVectors[2][0],
                (float)structure.cellVectors[2][1],
                (float)structure.cellVectors[2][2]);

    float lenA = glm::length(a);
    float lenB = glm::length(b);
    float lenC = glm::length(c);

    auto safeAngleDeg = [](const glm::vec3& u, const glm::vec3& v) -> float {
        float lu = glm::length(u);
        float lv = glm::length(v);
        if (lu <= 1e-10f || lv <= 1e-10f)
            return 0.0f;
        float cosang = glm::dot(u, v) / (lu * lv);
        cosang = std::max(-1.0f, std::min(1.0f, cosang));
        return kRadToDeg * std::acos(cosang);
    };

    float alpha = safeAngleDeg(b, c);
    float beta = safeAngleDeg(a, c);
    float gamma = safeAngleDeg(a, b);

    float volume = std::abs(glm::dot(a, glm::cross(b, c)));

    ImGui::Text("Lattice a,b,c (A):  %.5f  %.5f  %.5f", lenA, lenB, lenC);
    ImGui::Text("Angles alpha,beta,gamma (deg):  %.3f  %.3f  %.3f", alpha, beta, gamma);
    ImGui::Text("Cell volume (A^3):  %.5f", volume);
}

void drawElementCountsTable(const std::map<int, int>& counts)
{
    if (!ImGui::BeginTable("##element-counts", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
        return;

    ImGui::TableSetupColumn("Z", ImGuiTableColumnFlags_WidthFixed, 70.0f);
    ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, 120.0f);
    ImGui::TableSetupColumn("Count", ImGuiTableColumnFlags_WidthFixed, 100.0f);
    ImGui::TableHeadersRow();

    for (const auto& [z, count] : counts)
    {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%d", z);
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%s", elementSymbol(z));
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%d", count);
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

void drawSymmetrySection(const StructureInfoComputedData& data)
{
    ImGui::Separator();
    ImGui::Text("Symmetry (spglib)");

    if (!data.ready)
    {
        ImGui::Text("Computing symmetry...");
        return;
    }

    if (!data.symmetrySuccess)
    {
        ImGui::TextWrapped("Space group: unavailable (%s)", data.symmetryError.c_str());
        return;
    }

    ImGui::Text("Space group: %d (%s)", data.symmetrySpaceGroupNumber, data.symmetryInternationalSymbol.c_str());
    ImGui::Text("Hall symbol: %s", data.symmetryHallSymbol.c_str());
    ImGui::Text("Point group: %s", data.symmetryPointGroup.c_str());
}

void drawAtomicPositionsTable(const Structure& structure, bool hasValidLattice)
{
    if (!ImGui::BeginTable("##position-table", 8,
                           ImGuiTableFlags_Borders |
                           ImGuiTableFlags_RowBg |
                           ImGuiTableFlags_ScrollY))
    {
        return;
    }

    ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 55.0f);
    ImGui::TableSetupColumn("El", ImGuiTableColumnFlags_WidthFixed, 55.0f);
    ImGui::TableSetupColumn("X", ImGuiTableColumnFlags_WidthFixed, 110.0f);
    ImGui::TableSetupColumn("Y", ImGuiTableColumnFlags_WidthFixed, 110.0f);
    ImGui::TableSetupColumn("Z", ImGuiTableColumnFlags_WidthFixed, 110.0f);
    ImGui::TableSetupColumn("u", ImGuiTableColumnFlags_WidthFixed, 110.0f);
    ImGui::TableSetupColumn("v", ImGuiTableColumnFlags_WidthFixed, 110.0f);
    ImGui::TableSetupColumn("w", ImGuiTableColumnFlags_WidthFixed, 110.0f);
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
        ImGui::Text("%.6f", atom.x);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.6f", atom.y);
        ImGui::TableSetColumnIndex(4);
        ImGui::Text("%.6f", atom.z);
        ImGui::TableSetColumnIndex(5);
        ImGui::Text("%.6f", frac.x);
        ImGui::TableSetColumnIndex(6);
        ImGui::Text("%.6f", frac.y);
        ImGui::TableSetColumnIndex(7);
        ImGui::Text("%.6f", frac.z);
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

    ImGui::SetNextWindowSize(ImVec2(980.0f, 700.0f), ImGuiCond_FirstUseEver);
    bool popupOpen = true;
    if (ImGui::BeginPopupModal("Structure Info", &popupOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Structure Summary");
        ImGui::Separator();

        const bool hasComputedData = state.cachedData.ready && state.cachedFingerprint == currentFingerprint;

        const int totalAtoms = static_cast<int>(structure.atoms.size());
        const int uniqueElements = hasComputedData
                                       ? static_cast<int>(state.cachedData.elementCounts.size())
                                       : 0;

        ImGui::Text("Total atoms: %d", totalAtoms);
        if (hasComputedData)
        {
            ImGui::Text("Unique elements: %d", uniqueElements);
            ImGui::Text("Formula: %s", state.cachedData.formula.c_str());
        }
        else
        {
            ImGui::Text("Unique elements: computing...");
            ImGui::Text("Formula: computing...");
        }
        ImGui::Text("Has unit cell: %s", structure.hasUnitCell ? "Yes" : "No");

        if (structure.hasUnitCell)
            drawLatticeInfo(structure);

        ImGui::Separator();
        ImGui::Text("Element Counts");
        if (hasComputedData)
            drawElementCountsTable(state.cachedData.elementCounts);
        else
            ImGui::Text("Computing element counts...");

        drawSymmetrySection(hasComputedData ? state.cachedData : StructureInfoComputedData{});

        ImGui::Separator();
        ImGui::Text("Atomic Positions");

        constexpr std::size_t kMaxAtomsForPositionsTable = 500;
        if (structure.atoms.size() > kMaxAtomsForPositionsTable)
        {
            ImGui::TextWrapped(
                "Atomic position table is disabled for structures larger than %d atoms to keep the UI stable.",
                static_cast<int>(kMaxAtomsForPositionsTable));
        }
        else
        {
            glm::mat3 lattice(1.0f);
            glm::vec3 origin(0.0f);
            bool hasValidLattice = buildLatticeMatrix(structure, lattice, origin);

            ImGui::BeginChild("##positions", ImVec2(940.0f, 250.0f), true);
            drawAtomicPositionsTable(structure, hasValidLattice);
            ImGui::EndChild();
        }

        ImGui::EndPopup();
    }

    if (!popupOpen)
        ImGui::CloseCurrentPopup();
}