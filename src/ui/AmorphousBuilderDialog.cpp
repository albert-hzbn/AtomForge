#include "ui/AmorphousBuilderDialog.h"

#include "algorithms/AmorphousBuilder.h"
#include "io/StructureLoader.h"
#include "ui/PeriodicTableDialog.h"
#include "util/ElementData.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <sstream>

// ===========================================================================
// Construction
// ===========================================================================

AmorphousBuilderDialog::AmorphousBuilderDialog()
{
    // Default composition: 100 Si atoms.
    m_elements.push_back({ 14, 100 });
}

// ===========================================================================
// Menu item
// ===========================================================================

void AmorphousBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Amorphous Structure", nullptr, false, enabled))
        m_openRequested = true;
}

// ===========================================================================
// Helpers
// ===========================================================================

int AmorphousBuilderDialog::totalAtomCount() const
{
    int n = 0;
    for (const auto& e : m_elements)
        n += std::max(0, e.count);
    return n;
}

float AmorphousBuilderDialog::defaultMinDist(int z1, int z2,
    const std::vector<float>& covalentRadii) const
{
    const float r1 = (z1 >= 1 && z1 < (int)covalentRadii.size()) ? covalentRadii[z1] : 1.5f;
    const float r2 = (z2 >= 1 && z2 < (int)covalentRadii.size()) ? covalentRadii[z2] : 1.5f;
    return (r1 + r2) * m_covTolerance;
}

void AmorphousBuilderDialog::rebuildPairRows(const std::vector<float>& covalentRadii)
{
    // Collect unique atomic numbers from current composition.
    std::vector<int> zList;
    for (const auto& e : m_elements)
    {
        if (e.count > 0)
        {
            const int z = e.atomicNumber;
            if (std::find(zList.begin(), zList.end(), z) == zList.end())
                zList.push_back(z);
        }
    }
    std::sort(zList.begin(), zList.end());

    // Build all unique pairs (z1 <= z2).
    std::vector<PairRow> newRows;
    for (int i = 0; i < (int)zList.size(); ++i)
    for (int j = i; j < (int)zList.size(); ++j)
    {
        const int z1 = zList[i];
        const int z2 = zList[j];

        // Preserve existing override if present.
        bool found = false;
        for (const auto& pr : m_pairRows)
        {
            const int pa = std::min(pr.z1, pr.z2);
            const int pb = std::max(pr.z1, pr.z2);
            if (pa == z1 && pb == z2)
            {
                newRows.push_back(pr);
                found = true;
                break;
            }
        }
        if (!found)
        {
            PairRow pr;
            pr.z1      = z1;
            pr.z2      = z2;
            pr.minDist = defaultMinDist(z1, z2, covalentRadii);
            pr.enabled = false; // start disabled (use default)
            newRows.push_back(pr);
        }
    }
    m_pairRows = std::move(newRows);
}

// ===========================================================================
// Main dialog draw
// ===========================================================================

void AmorphousBuilderDialog::drawDialog(
    const std::vector<glm::vec3>& elementColors,
    const std::vector<float>&     covalentRadii,
    const std::function<void(Structure&)>& updateBuffers)
{
    // -----------------------------------------------------------------------
    // Open the modal when requested.
    // -----------------------------------------------------------------------
    if (m_openRequested)
    {
        ImGui::OpenPopup("Amorphous Structure Builder");
        m_openRequested = false;
        m_isOpen        = true;
        rebuildPairRows(covalentRadii);
        m_statusMsg.clear();
        m_statusIsError = false;
    }

    // -----------------------------------------------------------------------
    // Modal layout.
    // -----------------------------------------------------------------------
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowSize(ImVec2(680.0f, 0.0f), ImGuiCond_Appearing);
    ImGui::SetNextWindowPos(
        ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
        ImGuiCond_Appearing,
        ImVec2(0.5f, 0.5f));

    bool open = true;
    if (!ImGui::BeginPopupModal("Amorphous Structure Builder", &open,
                                 ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (!open)
            m_isOpen = false;
        return;
    }

    // =========================================================================
    // SECTION 1 – Composition
    // =========================================================================
    ImGui::SeparatorText("Composition");

    const float col0 = 100.0f;  // element button
    const float col1 = 260.0f;  // element name
    const float col2 = 340.0f;  // count
    const float col3 = 500.0f;  // at%

    ImGui::TextDisabled("Element");
    ImGui::SameLine(col1);
    ImGui::TextDisabled("Name");
    ImGui::SameLine(col2);
    ImGui::TextDisabled("Count");
    ImGui::SameLine(col3);
    ImGui::TextDisabled("at%%");
    ImGui::Separator();

    const int total = totalAtomCount();
    int removeIdx = -1;
    bool openPicker = false;

    for (int i = 0; i < (int)m_elements.size(); ++i)
    {
        auto& e = m_elements[(size_t)i];
        ImGui::PushID(i);

        // Symbol button → open picker.
        char btnLabel[16];
        std::snprintf(btnLabel, sizeof(btnLabel), "%s##btn",
                      elementSymbol(e.atomicNumber));
        if (ImGui::Button(btnLabel, ImVec2(col0 - 8.0f, 0.0f)))
        {
            m_pickerTarget = i;
            openPicker = true;
        }

        // Element name also acts as a button to open the picker.
        ImGui::SameLine(col1);
        char nameLabel[64];
        std::snprintf(nameLabel, sizeof(nameLabel), "%s##name",
                      elementName(e.atomicNumber));
        ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0, 0, 0, 0));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImGui::GetStyle().Colors[ImGuiCol_ButtonHovered]);
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImGui::GetStyle().Colors[ImGuiCol_ButtonActive]);
        ImGui::PushStyleColor(ImGuiCol_Border,        ImVec4(0, 0, 0, 0));
        if (ImGui::Button(nameLabel, ImVec2(col2 - col1 - 8.0f, 0.0f)))
        {
            m_pickerTarget = i;
            openPicker = true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Click to change element");
        ImGui::PopStyleColor(4);

        ImGui::SameLine(col2);
        ImGui::SetNextItemWidth(130.0f);
        ImGui::DragInt("##count", &e.count, 1.0f, 0, 100000, "%d atoms");

        ImGui::SameLine(col3);
        const float pct = (total > 0) ? 100.0f * (float)e.count / (float)total : 0.0f;
        ImGui::Text("%.1f", pct);

        ImGui::SameLine();
        if (ImGui::SmallButton("X"))
            removeIdx = i;

        ImGui::PopID();
    }

    if (removeIdx >= 0 && m_elements.size() > 1)
    {
        m_elements.erase(m_elements.begin() + removeIdx);
        rebuildPairRows(covalentRadii);
    }

    // Handle periodic-table picker.
    // openPeriodicTable() is intentionally called here (outside PushID loop)
    // so its ID stack depth matches drawPeriodicTable()'s BeginPopupModal call.
    if (openPicker)
        openPeriodicTable();
    {
        std::vector<ElementSelection> sel;
        if (drawPeriodicTable(sel))
        {
            if (m_pickerTarget >= 0 &&
                m_pickerTarget < (int)m_elements.size() &&
                !sel.empty())
            {
                m_elements[(size_t)m_pickerTarget].atomicNumber = sel[0].atomicNumber;
                rebuildPairRows(covalentRadii);
            }
            m_pickerTarget = -1;
        }
    }

    if (ImGui::Button("+ Add element"))
    {
        m_elements.push_back({ 14, 10 });
        rebuildPairRows(covalentRadii);
    }
    ImGui::SameLine();
    ImGui::TextDisabled("Total: %d atoms", total);

    // =========================================================================
    // SECTION 2 – Box dimensions
    // =========================================================================
    ImGui::Spacing();
    ImGui::SeparatorText("Simulation Box");

    ImGui::RadioButton("Auto from target density", &m_boxMode, 1);
    ImGui::SameLine();
    ImGui::RadioButton("Manual dimensions", &m_boxMode, 0);

    if (m_boxMode == 1) // AutoFromDensity
    {
        ImGui::SetNextItemWidth(180.0f);
        ImGui::DragFloat("Target density (g/cm\xc2\xb3)", &m_targetDensity, 0.01f, 0.01f, 30.0f, "%.3f");

        // Show the resulting box side length as a hint.
        if (total > 0)
        {
            double totalMass = 0.0;
            for (const auto& e : m_elements)
                if (e.atomicNumber >= 1 && e.atomicNumber <= 118)
                    totalMass += (double)e.count * elementAtomicMass(e.atomicNumber);
            const double dens = std::max(1e-3, (double)m_targetDensity);
            const double vol  = totalMass * 1.66054 / dens;
            const float  side = (float)std::cbrt(vol) * m_cellScaleFactor;
            ImGui::TextDisabled("  => cubic box %.2f \xc3\x85", side);
        }
    }
    else // Manual
    {
        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragFloat("a (\xc3\x85)##boxa", &m_boxA, 0.5f, 1.0f, 5000.0f, "%.2f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragFloat("b (\xc3\x85)##boxb", &m_boxB, 0.5f, 1.0f, 5000.0f, "%.2f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragFloat("c (\xc3\x85)##boxc", &m_boxC, 0.5f, 1.0f, 5000.0f, "%.2f");
    }

    ImGui::Spacing();
    ImGui::SetNextItemWidth(200.0f);
    ImGui::DragFloat("Cell scale factor", &m_cellScaleFactor, 0.01f, 0.5f, 3.0f, "%.3f");
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip(
            "Values > 1 expand the box before packing (lower effective density).\n"
            "This mirrors the 1.1x cell-expansion trick used in MD melt-quench\n"
            "workflows to help the crystal disorder before quenching.");

    ImGui::Checkbox("Periodic boundary conditions (attach unit cell)", &m_periodic);

    // =========================================================================
    // SECTION 3 – Advanced settings (collapsible)
    // =========================================================================
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Advanced Settings"))
    {
        // -- Placement parameters --
        ImGui::SeparatorText("Placement");
        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragFloat("Covalent radius tolerance", &m_covTolerance, 0.01f, 0.3f, 2.0f, "%.2f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(
                "Default min separation = (r_cov1 + r_cov2) * tolerance.\n"
                "0.75 prevents overlaps without enforcing full bonding.\n"
                "Decrease for denser random packings.");

        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputInt("RNG seed (0 = time)", &m_seed);

        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragInt("Max attempts per atom", &m_maxAttempts, 10.0f, 1, 100000);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(
                "Attempts to find a valid position for each atom.\n"
                "Increase if atoms are skipped; decrease to limit build time.");

        // -- Per-pair minimum distance overrides --
        ImGui::SeparatorText("Element-pair minimum distances");
        ImGui::TextDisabled("Check a row to override the covalent-radii default for that pair.");

        if (m_pairRows.empty())
        {
            ImGui::TextDisabled("(no pairs – add elements above)");
        }
        else
        {
            ImGui::TextDisabled("  Pair");
            ImGui::SameLine(120.0f);
            ImGui::TextDisabled("Default (\xc3\x85)");
            ImGui::SameLine(240.0f);
            ImGui::TextDisabled("Override (\xc3\x85)");
            ImGui::Separator();

            for (int pi = 0; pi < (int)m_pairRows.size(); ++pi)
            {
                auto& pr = m_pairRows[(size_t)pi];
                ImGui::PushID(pi + 1000);

                ImGui::Checkbox("##chk", &pr.enabled);
                ImGui::SameLine();
                ImGui::Text("%s–%s",
                            elementSymbol(pr.z1),
                            elementSymbol(pr.z2));
                ImGui::SameLine(120.0f);
                ImGui::TextDisabled("%.3f", defaultMinDist(pr.z1, pr.z2, covalentRadii));
                ImGui::SameLine(240.0f);
                if (pr.enabled)
                {
                    ImGui::SetNextItemWidth(100.0f);
                    ImGui::DragFloat("##dist", &pr.minDist, 0.01f, 0.1f, 10.0f, "%.3f \xc3\x85");
                }
                else
                {
                    ImGui::TextDisabled("(using default)");
                }

                ImGui::PopID();
            }
        }
    }

    // =========================================================================
    // SECTION 4 – Summary
    // =========================================================================
    ImGui::Spacing();
    ImGui::Separator();

    // Compute summary for display (same logic as builder, non-allocating).
    {
        float bA = m_boxA, bB = m_boxB, bC = m_boxC;
        if (m_boxMode == 1 && total > 0)
        {
            double tm = 0.0;
            for (const auto& e : m_elements)
                if (e.atomicNumber >= 1 && e.atomicNumber <= 118)
                    tm += (double)e.count * elementAtomicMass(e.atomicNumber);
            const double vol  = tm * 1.66054 / std::max(1e-3, (double)m_targetDensity);
            const float  side = (float)std::cbrt(vol) * m_cellScaleFactor;
            bA = bB = bC = side;
        }
        else
        {
            bA *= m_cellScaleFactor;
            bB *= m_cellScaleFactor;
            bC *= m_cellScaleFactor;
        }
        const float volume = bA * bB * bC;
        ImGui::TextDisabled("Box: %.2f x %.2f x %.2f \xc3\x85   |   Volume: %.1f \xc3\x85\xc2\xb3   |   Atoms: %d",
                             bA, bB, bC, volume, total);
    }

    // =========================================================================
    // SECTION 5 – Build / Close
    // =========================================================================
    ImGui::Spacing();
    const bool canBuild = (total > 0);
    if (!canBuild)
        ImGui::BeginDisabled();

    if (ImGui::Button("Build", ImVec2(120.0f, 0.0f)))
    {
        // Assemble AmorphousParams from dialog state.
        AmorphousParams params;

        for (const auto& e : m_elements)
            if (e.count > 0)
                params.elements.push_back({ e.atomicNumber, e.count });

        params.boxMode        = (m_boxMode == 0) ? AmorphousBoxMode::Manual
                                                  : AmorphousBoxMode::AutoFromDensity;
        params.boxA           = m_boxA;
        params.boxB           = m_boxB;
        params.boxC           = m_boxC;
        params.targetDensity  = m_targetDensity;
        params.cellScaleFactor= m_cellScaleFactor;
        params.covalentTolerance = m_covTolerance;
        params.seed           = (unsigned int)std::max(0, m_seed);
        params.maxAttempts    = std::max(1, m_maxAttempts);
        params.periodic       = m_periodic;

        for (const auto& pr : m_pairRows)
        {
            if (pr.enabled)
                params.pairDistances.push_back({ pr.z1, pr.z2, pr.minDist });
        }

        AmorphousResult res = buildAmorphousStructure(params, covalentRadii, elementColors);

        m_statusMsg     = res.message;
        m_statusIsError = !res.success;

        if (res.placedAtoms > 0)
        {
            // Append density to the message.
            std::ostringstream extra;
            extra << "  |  density: " << std::fixed;
            extra.precision(3);
            extra << res.actualDensity << " g/cm\xc2\xb3"
                  << "  |  box: " << res.boxA << " x " << res.boxB
                  << " x " << res.boxC << " \xc3\x85";
            m_statusMsg += extra.str();

            updateBuffers(res.output);
        }
    }

    if (!canBuild)
        ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
    {
        ImGui::CloseCurrentPopup();
        m_isOpen = false;
    }

    // Status line.
    if (!m_statusMsg.empty())
    {
        ImGui::Spacing();
        if (m_statusIsError)
        {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));
            ImGui::TextWrapped("Error: %s", m_statusMsg.c_str());
            ImGui::PopStyleColor();
        }
        else
        {
            ImGui::TextWrapped("%s", m_statusMsg.c_str());
        }
    }

    ImGui::EndPopup();
}
