#include "ui/CommonNeighbourAnalysis.h"

#include "ElementData.h"
#include "math/StructureMath.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <map>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace
{
constexpr float kDefaultCutoffScale = 1.18f;
constexpr float kSpatialHashCellSize = 4.0f;

struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u) ^ ((size_t)value.y * 19349663u) ^ ((size_t)value.z * 83492791u);
    }
};

glm::ivec3 getGridCell(const glm::vec3& position)
{
    return glm::ivec3((int)std::floor(position.x / kSpatialHashCellSize),
                      (int)std::floor(position.y / kSpatialHashCellSize),
                      (int)std::floor(position.z / kSpatialHashCellSize));
}

struct Signature
{
    int common = 0;
    int bonds = 0;
    int chain = 0;

    bool operator<(const Signature& other) const
    {
        if (common != other.common) return common < other.common;
        if (bonds != other.bonds) return bonds < other.bonds;
        return chain < other.chain;
    }

    bool operator==(const Signature& other) const
    {
        return common == other.common && bonds == other.bonds && chain == other.chain;
    }
};

struct AtomRow
{
    int index = -1;
    int atomicNumber = 0;
    std::string symbol;
    int coordination = 0;
    Signature dominantSignature;
    int dominantSignatureCount = 0;
    std::string environment;
};

struct CnaResult
{
    bool valid = false;
    std::string message;

    int atomCount = 0;
    int pairCount = 0;
    bool pbcUsed = false;

    std::map<Signature, int> signatureCounts;
    std::map<std::string, int> environmentCounts;
    std::vector<AtomRow> atomRows;
};

std::string signatureToString(const Signature& s)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "1-%d-%d-%d", s.common, s.bonds, s.chain);
    return std::string(buffer);
}

std::string classifyEnvironment(const Signature& s)
{
    Signature fcc; fcc.common = 4; fcc.bonds = 2; fcc.chain = 1;
    Signature hcp; hcp.common = 4; hcp.bonds = 2; hcp.chain = 2;
    Signature bccA; bccA.common = 4; bccA.bonds = 4; bccA.chain = 1;
    Signature bccB; bccB.common = 6; bccB.bonds = 6; bccB.chain = 1;
    Signature ico; ico.common = 5; ico.bonds = 5; ico.chain = 1;

    if (s == fcc) return "FCC-like";
    if (s == hcp) return "HCP-like";
    if (s == bccA || s == bccB) return "BCC-like";
    if (s == ico) return "ICO-like";
    return "Unknown";
}

uint64_t edgeKey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

int longestChainLength(const std::vector<int>& commonNodes,
                       const std::unordered_set<uint64_t>& edgeSet)
{
    if (commonNodes.empty())
        return 0;
    if (commonNodes.size() == 1)
        return 1;

    std::map<int, int> localIndex;
    for (int i = 0; i < (int)commonNodes.size(); ++i)
        localIndex[commonNodes[i]] = i;

    std::vector<std::vector<int>> adjacency(commonNodes.size());
    for (int i = 0; i < (int)commonNodes.size(); ++i)
    {
        for (int j = i + 1; j < (int)commonNodes.size(); ++j)
        {
            if (edgeSet.find(edgeKey(commonNodes[i], commonNodes[j])) != edgeSet.end())
            {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
        }
    }

    int best = 1;
    for (int src = 0; src < (int)commonNodes.size(); ++src)
    {
        std::vector<int> dist(commonNodes.size(), -1);
        std::queue<int> q;
        dist[src] = 0;
        q.push(src);

        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            for (int v : adjacency[u])
            {
                if (dist[v] >= 0)
                    continue;
                dist[v] = dist[u] + 1;
                best = std::max(best, dist[v] + 1);
                q.push(v);
            }
        }
    }

    return best;
}

CnaResult runCna(const Structure& structure, float cutoffScale, bool usePbcRequest)
{
    CnaResult result;
    result.atomCount = (int)structure.atoms.size();

    if (structure.atoms.empty())
    {
        result.valid = false;
        result.message = "No atoms available.";
        return result;
    }

    std::vector<float> radii = makeLiteratureCovalentRadii();

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    bool usePbc = false;
    if (usePbcRequest && structure.hasUnitCell)
    {
        usePbc = tryMakeCellMatrices(structure, cell, invCell);
    }
    result.pbcUsed = usePbc;

    std::vector<glm::vec3> positions(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        positions[i] = glm::vec3((float)structure.atoms[i].x,
                                 (float)structure.atoms[i].y,
                                 (float)structure.atoms[i].z);
    }

    std::vector<std::vector<int>> neighbors(structure.atoms.size());
    std::unordered_set<uint64_t> edgeSet;

    auto tryAddBond = [&](int i, int j)
    {
        if (j <= i)
            return;

        int zi = structure.atoms[i].atomicNumber;
        float ri = (zi >= 0 && zi < (int)radii.size()) ? radii[zi] : 1.0f;
        int zj = structure.atoms[j].atomicNumber;
        float rj = (zj >= 0 && zj < (int)radii.size()) ? radii[zj] : 1.0f;

        glm::vec3 delta = minimumImageDelta(positions[j] - positions[i], usePbc, cell, invCell);
        float d = glm::length(delta);
        if (d <= kMinBondDistance)
            return;

        float cutoff = (ri + rj) * cutoffScale;
        if (d > cutoff)
            return;

        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        edgeSet.insert(edgeKey(i, j));
    };

    if (!usePbc)
    {
        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
        grid.reserve(positions.size());
        for (int i = 0; i < (int)positions.size(); ++i)
            grid[getGridCell(positions[i])].push_back(i);

        for (int i = 0; i < (int)positions.size(); ++i)
        {
            const glm::ivec3 cellCoord = getGridCell(positions[i]);
            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        const glm::ivec3 neighborCell(cellCoord.x + dx,
                                                      cellCoord.y + dy,
                                                      cellCoord.z + dz);
                        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash>::const_iterator it = grid.find(neighborCell);
                        if (it == grid.end())
                            continue;

                        const std::vector<int>& candidates = it->second;
                        for (int index = 0; index < (int)candidates.size(); ++index)
                            tryAddBond(i, candidates[index]);
                    }
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            for (int j = i + 1; j < (int)structure.atoms.size(); ++j)
                tryAddBond(i, j);
        }
    }

    for (int i = 0; i < (int)neighbors.size(); ++i)
        std::sort(neighbors[i].begin(), neighbors[i].end());

    std::vector<std::map<Signature, int>> atomSignatures(structure.atoms.size());

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        for (int t = 0; t < (int)neighbors[i].size(); ++t)
        {
            int j = neighbors[i][t];
            if (j <= i)
                continue;

            std::vector<int> common;
            common.reserve(std::min(neighbors[i].size(), neighbors[j].size()));
            std::set_intersection(neighbors[i].begin(), neighbors[i].end(),
                                  neighbors[j].begin(), neighbors[j].end(),
                                  std::back_inserter(common));

            int commonCount = (int)common.size();
            int bondCount = 0;
            for (int a = 0; a < (int)common.size(); ++a)
            {
                for (int b = a + 1; b < (int)common.size(); ++b)
                {
                    if (edgeSet.find(edgeKey(common[a], common[b])) != edgeSet.end())
                        ++bondCount;
                }
            }

            int chainLen = longestChainLength(common, edgeSet);

            Signature sig;
            sig.common = commonCount;
            sig.bonds = bondCount;
            sig.chain = chainLen;

            ++result.signatureCounts[sig];
            ++atomSignatures[i][sig];
            ++atomSignatures[j][sig];
            ++result.pairCount;
        }
    }

    result.atomRows.reserve(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        AtomRow row;
        row.index = i;
        row.atomicNumber = structure.atoms[i].atomicNumber;
        row.symbol = structure.atoms[i].symbol;
        row.coordination = (int)neighbors[i].size();

        Signature dominant;
        int dominantCount = 0;
        for (std::map<Signature, int>::const_iterator it = atomSignatures[i].begin(); it != atomSignatures[i].end(); ++it)
        {
            if (it->second > dominantCount)
            {
                dominant = it->first;
                dominantCount = it->second;
            }
        }

        row.dominantSignature = dominant;
        row.dominantSignatureCount = dominantCount;
        row.environment = (dominantCount > 0) ? classifyEnvironment(dominant) : "Unknown";

        ++result.environmentCounts[row.environment];
        result.atomRows.push_back(row);
    }

    result.valid = true;
    result.message = "CNA completed.";
    return result;
}

void drawSignatureTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-signatures", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Sortable))
    {
        ImGui::TableSetupColumn("Signature", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Pair Count", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Fraction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableHeadersRow();

        for (std::map<Signature, int>::const_iterator it = result.signatureCounts.begin(); it != result.signatureCounts.end(); ++it)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", signatureToString(it->first).c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", it->second);
            ImGui::TableSetColumnIndex(2);
            float f = (result.pairCount > 0) ? (float)it->second / (float)result.pairCount : 0.0f;
            ImGui::Text("%.4f", f);
        }

        ImGui::EndTable();
    }
}

void drawEnvironmentTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-environments", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
    {
        ImGui::TableSetupColumn("Environment", ImGuiTableColumnFlags_WidthFixed, 160.0f);
        ImGui::TableSetupColumn("Atom Count", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Fraction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableHeadersRow();

        for (std::map<std::string, int>::const_iterator it = result.environmentCounts.begin(); it != result.environmentCounts.end(); ++it)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", it->first.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", it->second);
            ImGui::TableSetColumnIndex(2);
            float f = (result.atomCount > 0) ? (float)it->second / (float)result.atomCount : 0.0f;
            ImGui::Text("%.4f", f);
        }

        ImGui::EndTable();
    }
}

void drawAtomTable(const CnaResult& result)
{
    if (ImGui::BeginTable("##cna-per-atom", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY))
    {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 52.0f);
        ImGui::TableSetupColumn("El", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Z", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("CN", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Dominant Signature", ImGuiTableColumnFlags_WidthFixed, 160.0f);
        ImGui::TableSetupColumn("Sig Count", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Environment", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (int i = 0; i < (int)result.atomRows.size(); ++i)
        {
            const AtomRow& row = result.atomRows[i];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", row.index);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s", row.symbol.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", row.atomicNumber);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", row.coordination);
            ImGui::TableSetColumnIndex(4);
            if (row.dominantSignatureCount > 0)
                ImGui::Text("%s", signatureToString(row.dominantSignature).c_str());
            else
                ImGui::Text("N/A");
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%d", row.dominantSignatureCount);
            ImGui::TableSetColumnIndex(6);
            ImGui::Text("%s", row.environment.c_str());
        }

        ImGui::EndTable();
    }
}

} // namespace

void CommonNeighbourAnalysisDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Common Neighbour Analysis...", NULL, false, enabled))
        m_openRequested = true;
}

void CommonNeighbourAnalysisDialog::drawDialog(const Structure& structure)
{
    static bool requestRecompute = true;
    static bool usePbc = true;
    static float cutoffScale = kDefaultCutoffScale;
    static CnaResult result;

    if (m_openRequested)
    {
        ImGui::OpenPopup("Common Neighbour Analysis");
        m_openRequested = false;
        requestRecompute = true;
    }

    ImGui::SetNextWindowSize(ImVec2(1100.0f, 760.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("Common Neighbour Analysis", &dialogOpen, ImGuiWindowFlags_NoResize))
    {
        bool changed = false;
        changed |= ImGui::Checkbox("Use PBC when unit cell is available", &usePbc);
        changed |= ImGui::SliderFloat("Bond cutoff scale", &cutoffScale, 1.00f, 1.60f, "%.2f");
        ImGui::SameLine();
        if (ImGui::Button("Run CNA"))
            requestRecompute = true;

        if (changed)
            requestRecompute = true;

        if (requestRecompute)
        {
            result = runCna(structure, cutoffScale, usePbc);
            requestRecompute = false;
        }

        ImGui::Separator();
        ImGui::Text("Status: %s", result.message.c_str());
        ImGui::Text("Atoms: %d", result.atomCount);
        ImGui::Text("Bonded pairs analyzed: %d", result.pairCount);
        ImGui::Text("PBC used in analysis: %s", result.pbcUsed ? "Yes" : "No");

        if (result.valid)
        {
            ImGui::Separator();
            ImGui::Text("Pair Signature Distribution (Honeycutt-Andersen form: 1-j-k-l)");
            ImGui::BeginChild("##cna-sig-child", ImVec2(0.0f, 180.0f), true);
            drawSignatureTable(result);
            ImGui::EndChild();

            ImGui::Separator();
            ImGui::Text("Per-Atom Environment Summary");
            ImGui::BeginChild("##cna-env-child", ImVec2(0.0f, 130.0f), true);
            drawEnvironmentTable(result);
            ImGui::EndChild();

            ImGui::Separator();
            ImGui::Text("Per-Atom CNA Details");
            ImGui::BeginChild("##cna-atom-child", ImVec2(0.0f, 240.0f), true);
            drawAtomTable(result);
            ImGui::EndChild();
        }

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();
}
