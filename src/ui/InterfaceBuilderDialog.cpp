#include "ui/InterfaceBuilderDialog.h"

#include "ui/DialogLayout.h"
#include "ui/ResponsiveLayout.h"
#include "ui/ThemeUtils.h"
#include "ui/UiScript.h"

#include "algorithms/InterfaceBuilder.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "util/ElementData.h"
#include "pipeline/Options.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// Stradi, Martinez, Blom, Brandbyge and Stokbro, "Method for determining
// optimal supercell representation of interfaces", J. Phys.: Condens. Matter
// 29, 185901 (2017), as QuantumATK's InterfaceBuilder.

namespace iface = atomforge::interfaces;

namespace
{
const char* strainOnName(iface::StrainOn on) { return on == iface::StrainOn::A ? "A" : on == iface::StrainOn::Both ? "both" : "B"; }

std::string fileNameOf(const std::string& path)
{
    const auto slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

// A label column and a field column, aligned across the rows of a panel.
bool beginForm(const char* id)
{
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingStretchProp)) return false;
    ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.46f);
    ImGui::TableSetupColumn("field", ImGuiTableColumnFlags_WidthStretch, 0.54f);
    return true;
}

void formRow(const char* label)
{
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
}
}

InterfaceBuilderDialog::~InterfaceBuilderDialog()
{
    delete m_sphere;
    delete m_cylinder;
    if (m_fbo) glDeleteFramebuffers(1, &m_fbo);
    if (m_texture) glDeleteTextures(1, &m_texture);
    if (m_depth) glDeleteRenderbuffers(1, &m_depth);
}

void InterfaceBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_sphere = new SphereMesh(24, 24);
    m_cylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_sphere->vbo, m_sphere->ebo, m_sphere->indexCount, m_cylinder->vbo, m_cylinder->vertexCount);
    m_shadow = createShadowMap(1, 1);
    m_glReady = true;
}

void InterfaceBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Interface Builder", nullptr, false, enabled))
    {
        m_step.finish();  // from the menu: build normally
        m_stepSearch = false;
        m_openRequested = true;
    }
}

void InterfaceBuilderDialog::feedDroppedFile(const std::string& path, int slotHint)
{
    if (slotHint == 1 || (slotHint == 0 && m_a.crystal.atoms.empty())) m_pendingDropA = path;
    else m_pendingDropB = path;
}

bool InterfaceBuilderDialog::loadFile(const std::string& path, Input& input)
{
    Structure loaded;
    std::string error;
    if (!loadStructureFromFile(path, loaded, error))
    {
        input.status = "Could not load " + fileNameOf(path) + ": " + error;
        return false;
    }
    if (!loaded.hasUnitCell)
    {
        input.status = fileNameOf(path) + " has no unit cell; the builder needs a crystal.";
        return false;
    }
    input.crystal = std::move(loaded);
    input.path = path;
    input.status = fileNameOf(path) + "  (" + std::to_string(input.crystal.atoms.size()) + " atoms)";
    m_stale = true;
    return true;
}

void InterfaceBuilderDialog::useStructure(const Structure& structure, Input& input, const char* what)
{
    input.crystal = structure;
    input.path.clear();
    input.status = std::string(what) + "  (" + std::to_string(structure.atoms.size()) + " atoms)";
    m_stale = true;
}

void InterfaceBuilderDialog::runSearch()
{
    m_matches.clear();
    m_selected = m_suggested = -1;
    m_searched = true;
    m_stale = false;
    m_previewDirty = true;
    m_preview = Structure();
    const auto side = [](const Input& input) {
        iface::Side s;
        s.crystal = input.crystal;
        s.cleave = input.cleave;
        s.h = input.hkl[0];
        s.k = input.hkl[1];
        s.l = input.hkl[2];
        s.layers = input.layers;
        return s;
    };
    std::string error;
    if (!iface::prepareSlab(side(m_a), m_slabA, error))
    {
        m_searchMessage = "Crystal A: " + error;
        return;
    }
    if (!iface::prepareSlab(side(m_b), m_slabB, error))
    {
        m_searchMessage = "Crystal B: " + error;
        return;
    }
    m_matches = iface::findMatches(m_slabA, m_slabB, m_search);
    if (m_matches.empty())
    {
        char text[256];
        std::snprintf(text, sizeof(text), "No matches within %.3g%% strain and %.0f \xc3\x85\xc2\xb2. Raise the strain or area limit, or the repeats.",
                      m_search.maxStrain * 100.0, m_search.maxArea);
        m_searchMessage = text;
        return;
    }
    m_suggested = iface::suggestedMatch(m_matches);
    m_searchMessage = std::to_string(m_matches.size()) + (m_matches.size() == 1 ? " match" : " matches");
    select(m_suggested);
}

void InterfaceBuilderDialog::select(int index)
{
    if (index < 0 || index >= static_cast<int>(m_matches.size())) return;
    if (index != m_selected) m_scrollToSelected = true;
    m_selected = index;
    m_previewDirty = true;
}

void InterfaceBuilderDialog::rebuildPreview(const std::vector<float>& radii, const std::vector<float>& shininess)
{
    m_previewDirty = false;
    if (m_selected < 0 || m_selected >= static_cast<int>(m_matches.size())) return;
    m_preview = iface::buildInterface(m_slabA, m_slabB, m_matches[static_cast<std::size_t>(m_selected)], m_assembly);
    if (m_repeatX > 1 || m_repeatY > 1) m_preview = repeatInterfaceXY(m_preview, m_repeatX, m_repeatY);
    if (!m_glReady || m_preview.atoms.empty()) return;
    static const int identity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const StructureInstanceData data = buildStructureInstanceData(m_preview, false, identity, radii, shininess);
    std::array<bool, 119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    // Far enough to see the whole cell.
    float radius = 1.0f;
    for (std::size_t i = 0; i < m_previewBuffers.atomPositions.size(); ++i)
        radius = std::max(radius, glm::length(m_previewBuffers.atomPositions[i] - m_previewBuffers.orbitCenter) +
                                      (i < m_previewBuffers.atomRadii.size() ? m_previewBuffers.atomRadii[i] : 0.0f));
    m_distance = std::clamp(radius / std::sin(glm::radians(22.5f)) * 0.9f, Camera::kMinDistance, Camera::kMaxDistance);
}

bool InterfaceBuilderDialog::editStep(StepEdit edit)
{
    m_step = std::move(edit);
    const auto options = m_step.reader();
    // Missing flags take the defaults of AtomForge --build interface.
    const auto plane = [&](const char* flag, Input& input) {
        const auto values = options.numbers(flag);
        input.cleave = values.size() == 3;
        if (input.cleave)
            for (int i = 0; i < 3; ++i) input.hkl[i] = static_cast<int>(std::lround(values[static_cast<std::size_t>(i)]));
    };
    m_a = Input();
    m_b = Input();
    plane("--hklA", m_a);
    plane("--hklB", m_b);
    m_a.layers = std::max(1, options.integer("--layersA", 4));
    m_b.layers = std::max(1, options.integer("--layersB", 4));
    m_search = iface::SearchParams();
    m_search.maxStrain = std::max(1e-5, options.number("--max-strain", 0.03));
    m_search.maxArea = std::max(1.0, options.number("--max-area", 150.0));
    m_search.maxRepeats = std::clamp(options.integer("--repeats", options.integer("--nmax", 6)), 1, 30);
    m_search.minCellAngle = std::clamp(options.number("--min-angle", 20.0), 0.0, 89.0);
    m_search.angleMin = options.number("--angle-min", 0.0);
    m_search.angleMax = options.number("--angle-max", 180.0);
    m_search.angleStep = std::max(0.01, options.number("--angle-step", 1.0));
    m_assembly = iface::Assembly();
    const std::string on = options.text("--strain-on");
    m_assembly.strainOn = on == "A" || on == "a" ? iface::StrainOn::A : on == "both" ? iface::StrainOn::Both : iface::StrainOn::B;
    m_assembly.gap = std::max(0.0, options.number("--gap", 2.0));
    m_assembly.vacuum = std::max(0.0, options.number("--vacuum", 0.0));
    if (const auto shift = options.numbers("--shift"); shift.size() == 2)
    {
        m_assembly.shift[0] = shift[0];
        m_assembly.shift[1] = shift[1];
    }
    m_repeatX = std::max(1, options.integer("--repx", 1));
    m_repeatY = std::max(1, options.integer("--repy", 1));

    // Crystal A: a file named in the options, else the structure entering the step.
    m_a.path = options.text("--layerA");
    if (!m_a.path.empty()) m_pendingDropA = m_a.path;
    else if (!m_step.input.atoms.empty()) useStructure(m_step.input, m_a, "Step input");
    m_b.path = options.text("--layerB");
    if (!m_b.path.empty()) m_pendingDropB = m_b.path;

    m_matches.clear();
    m_searched = false;
    m_selected = -1;
    m_stepPick = options.has("--pick") ? options.integer("--pick", -1) : -1;
    m_stepSearch = true;
    m_openRequested = true;
    return true;
}

std::string InterfaceBuilderDialog::stepOptions() const
{
    atomforge::pipeline::options::Writer out;
    if (!m_a.path.empty()) out.add("--layerA", m_a.path);
    if (!m_b.path.empty()) out.add("--layerB", m_b.path);
    const auto hkl = [](const Input& input) {
        return std::vector<double>{static_cast<double>(input.hkl[0]), static_cast<double>(input.hkl[1]), static_cast<double>(input.hkl[2])};
    };
    if (m_a.cleave) out.add("--hklA", hkl(m_a));
    if (m_b.cleave) out.add("--hklB", hkl(m_b));
    out.add("--layersA", m_a.layers).add("--layersB", m_b.layers);
    out.add("--max-strain", m_search.maxStrain).add("--max-area", m_search.maxArea).add("--repeats", m_search.maxRepeats);
    out.add("--min-angle", m_search.minCellAngle);
    out.add("--angle-min", m_search.angleMin).add("--angle-max", m_search.angleMax).add("--angle-step", m_search.angleStep);
    out.add("--strain-on", strainOnName(m_assembly.strainOn)).add("--gap", m_assembly.gap).add("--vacuum", m_assembly.vacuum);
    if (m_assembly.shift[0] != 0.0 || m_assembly.shift[1] != 0.0) out.add("--shift", std::vector<double>{m_assembly.shift[0], m_assembly.shift[1]});
    out.add("--repx", m_repeatX).add("--repy", m_repeatY);
    if (m_selected >= 0) out.add("--pick", m_selected);
    else if (m_stepPick >= 0) out.add("--pick", m_stepPick);
    return out.str();
}

void InterfaceBuilderDialog::openExample()
{
    m_step.finish();
    const auto fcc = [](const char* symbol, int z, double a) {
        Structure s;
        s.hasUnitCell = true;
        s.cellVectors = {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
        for (const auto& f : {std::array<double, 3>{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0}})
        {
            AtomSite atom;
            atom.symbol = symbol;
            atom.atomicNumber = z;
            atom.x = f[0] * a;
            atom.y = f[1] * a;
            atom.z = f[2] * a;
            getDefaultElementColor(z, atom.r, atom.g, atom.b);
            s.atoms.push_back(atom);
        }
        return s;
    };
    useStructure(fcc("Cu", 29, 3.615), m_a, "Example: fcc Cu");
    useStructure(fcc("Ni", 28, 3.524), m_b, "Example: fcc Ni");
    m_a.cleave = m_b.cleave = true;
    m_a.hkl[0] = m_a.hkl[1] = m_a.hkl[2] = 1;
    m_b.hkl[0] = m_b.hkl[1] = m_b.hkl[2] = 1;
    m_search.maxStrain = 0.03;
    m_matches.clear();
    m_searched = false;
    m_selected = -1;
    m_preview = Structure();
    m_openRequested = true;
}

void InterfaceBuilderDialog::drawInput(const char* title, Input& input, int slot, const Structure& active)
{
    ImGui::PushID(slot);
    dialogLayout::section(title);
    ImGui::PushTextWrapPos(0.0f);
    if (input.crystal.atoms.empty() && input.status.empty())
        ImGui::TextDisabled("Load a crystal, drop a file here, or use the structure in the view.");
    else
        ImGui::TextUnformatted(input.status.c_str());
    ImGui::PopTextWrapPos();
    if (responsive::button("Load File..."))
    {
        m_pickerSlot = slot;
        m_picker.open(slot == 1 ? "Load crystal A" : "Load crystal B", false, input.path);
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(active.atoms.empty() || !active.hasUnitCell);
    if (responsive::button("Use Active Structure")) useStructure(active, input, "Active structure");
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
        ImGui::SetTooltip(active.hasUnitCell ? "Use the structure in the view as this crystal" : "The structure in the view has no unit cell");

    if (beginForm("##side"))
    {
        formRow("Interface plane");
        if (ImGui::Checkbox("Cut along (h k l)", &input.cleave)) m_stale = true;
        dialogLayout::hint("Cut the crystal along a plane with the surface builder. Untick for a structure that is already a slab: its a and b vectors span the interface.");
        if (input.cleave)
        {
            formRow("(h k l)");
            if (ImGui::InputInt3("##hkl", input.hkl)) m_stale = true;
        }
        formRow("Layers");
        if (ImGui::InputInt("##layers", &input.layers, 1, 2))
        {
            input.layers = std::clamp(input.layers, 1, 200);
            m_stale = true;
        }
        dialogLayout::hint("Repeats of the crystal's translation along the interface normal (atomic layers for most planes).");
        ImGui::EndTable();
    }
    ImGui::PopID();
}

void InterfaceBuilderDialog::drawSearchSettings()
{
    dialogLayout::section("Matching");
    if (!beginForm("##search")) return;
    float strainPercent = static_cast<float>(m_search.maxStrain * 100.0);
    formRow("Max strain (%)");
    if (ImGui::InputFloat("##strain", &strainPercent, 0.5f, 1.0f, "%.2f"))
    {
        m_search.maxStrain = std::clamp(strainPercent, 0.01f, 50.0f) / 100.0;
        m_stale = true;
    }
    dialogLayout::hint("Largest mean absolute strain (|e11| + |e22| + |e12|) / 3 of the strained crystal.");
    float area = static_cast<float>(m_search.maxArea);
    formRow("Max area (\xc3\x85\xc2\xb2)");
    if (ImGui::InputFloat("##area", &area, 10.0f, 50.0f, "%.0f"))
    {
        m_search.maxArea = std::clamp(area, 1.0f, 100000.0f);
        m_stale = true;
    }
    dialogLayout::hint("Largest interface cell area: the main limit on the number of atoms.");
    formRow("Repeats");
    if (ImGui::InputInt("##repeats", &m_search.maxRepeats, 1, 2))
    {
        m_search.maxRepeats = std::clamp(m_search.maxRepeats, 1, 30);
        m_stale = true;
    }
    dialogLayout::hint("Largest multiple i, j of B's surface vectors (i b1 + j b2) tried for the interface cell.");
    float minAngle = static_cast<float>(m_search.minCellAngle);
    formRow("Min cell angle (\xc2\xb0)");
    if (ImGui::InputFloat("##minangle", &minAngle, 5.0f, 10.0f, "%.0f"))
    {
        m_search.minCellAngle = std::clamp(minAngle, 0.0f, 89.0f);
        m_stale = true;
    }
    dialogLayout::hint("Smallest angle between the interface cell vectors (excludes very slanted cells).");
    float rotation[3] = {static_cast<float>(m_search.angleMin), static_cast<float>(m_search.angleMax), static_cast<float>(m_search.angleStep)};
    formRow("Rotation from, to, step (\xc2\xb0)");
    if (ImGui::InputFloat3("##rotation", rotation, "%.1f"))
    {
        m_search.angleMin = std::clamp(rotation[0], -180.0f, 360.0f);
        m_search.angleMax = std::clamp(rotation[1], rotation[0], 360.0f);
        m_search.angleStep = std::clamp(rotation[2], 0.05f, 90.0f);
        m_stale = true;
    }
    dialogLayout::hint("Rotations of B about the interface normal that are searched.");
    ImGui::EndTable();
}

void InterfaceBuilderDialog::drawResultsTable(float height)
{
    const ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersInnerV |
                                  ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable | ImGuiTableFlags_Resizable;
    if (!ImGui::BeginTable("##matches", 10, flags, ImVec2(0, height))) return;
    ImGui::TableSetupScrollFreeze(0, 1);
    ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_DefaultSort);
    ImGui::TableSetupColumn("Atoms");
    ImGui::TableSetupColumn("Strain %");
    ImGui::TableSetupColumn("e11 %", ImGuiTableColumnFlags_NoSort);
    ImGui::TableSetupColumn("e22 %", ImGuiTableColumnFlags_NoSort);
    ImGui::TableSetupColumn("e12 %", ImGuiTableColumnFlags_NoSort);
    ImGui::TableSetupColumn("Area \xc3\x85\xc2\xb2");
    ImGui::TableSetupColumn("Cell angle", ImGuiTableColumnFlags_NoSort);
    ImGui::TableSetupColumn("Rotation", ImGuiTableColumnFlags_NoSort);
    ImGui::TableSetupColumn("Cells A / B", ImGuiTableColumnFlags_WidthStretch | ImGuiTableColumnFlags_NoSort);
    ImGui::TableHeadersRow();

    // Rows in the order of the clicked column (the list order: atoms, then strain).
    std::vector<int> order(m_matches.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    if (const ImGuiTableSortSpecs* sort = ImGui::TableGetSortSpecs(); sort && sort->SpecsCount > 0)
    {
        const int column = sort->Specs[0].ColumnIndex;
        const bool descending = sort->Specs[0].SortDirection == ImGuiSortDirection_Descending;
        std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
            const auto& a = m_matches[static_cast<std::size_t>(x)];
            const auto& b = m_matches[static_cast<std::size_t>(y)];
            const double ka = column == 1 ? a.atoms : column == 2 ? a.meanStrain : column == 6 ? a.area : x;
            const double kb = column == 1 ? b.atoms : column == 2 ? b.meanStrain : column == 6 ? b.area : y;
            return descending ? ka > kb : ka < kb;
        });
    }
    for (int index : order)
    {
        const auto& m = m_matches[static_cast<std::size_t>(index)];
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::PushID(index);
        char label[32];
        std::snprintf(label, sizeof(label), "%d%s", index, index == m_suggested ? "  \xe2\x80\xa2" : "");
        if (ImGui::Selectable(label, index == m_selected, ImGuiSelectableFlags_SpanAllColumns)) select(index);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s%s", m.pareto ? "Pareto-optimal: no other match has as few atoms and less strain." : "Another match has as few atoms and less strain.",
                              index == m_suggested ? "\nSuggested (\xe2\x80\xa2): lowest strain (%) + area / 25 \xc3\x85\xc2\xb2 of the Pareto-optimal matches." : "");
        if (index == m_selected && m_scrollToSelected)
        {
            ImGui::SetScrollHereY(0.5f);
            m_scrollToSelected = false;
        }
        ImGui::PopID();
        const auto cell = [](int column, const char* format, double value) {
            ImGui::TableSetColumnIndex(column);
            ImGui::Text(format, value);
        };
        ImGui::TableSetColumnIndex(1);
        if (m.pareto) ImGui::Text("%d", m.atoms);
        else ImGui::TextDisabled("%d", m.atoms);
        cell(2, "%.3f", m.meanStrain * 100.0);
        cell(3, "%.3f", m.strain[0] * 100.0);
        cell(4, "%.3f", m.strain[1] * 100.0);
        cell(5, "%.3f", m.strain[2] * 100.0);
        cell(6, "%.1f", m.area);
        cell(7, "%.1f\xc2\xb0", m.cellAngle);
        cell(8, "%.1f\xc2\xb0", m.thetaDeg);
        ImGui::TableSetColumnIndex(9);
        ImGui::Text("%d / %d", m.cellsA, m.cellsB);
    }
    ImGui::EndTable();
}

void InterfaceBuilderDialog::drawPlot(float width, float height)
{
    responsive::beginChild("##plot", ImVec2(width, height), ImGuiChildFlags_Borders, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImGui::TextDisabled("Strain against atoms (click a point)");
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImGuiStyle& style = ImGui::GetStyle();
    const float labelWidth = ImGui::CalcTextSize("00000").x + style.ItemSpacing.x;
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const ImVec2 available = ImGui::GetContentRegionAvail();
    const float bottom = ImGui::GetTextLineHeightWithSpacing() * 2.0f;
    const ImVec2 lo(origin.x + labelWidth, origin.y + 4.0f);
    const ImVec2 hi(origin.x + available.x - 8.0f, origin.y + available.y - bottom);
    if (hi.x - lo.x < 40.0f || hi.y - lo.y < 40.0f || m_matches.empty()) { ImGui::EndChild(); return; }

    double xMax = 0.0, yMin = 1e300, yMax = 0.0;
    for (const auto& m : m_matches)
    {
        xMax = std::max(xMax, m.meanStrain * 100.0);
        yMin = std::min(yMin, static_cast<double>(m.atoms));
        yMax = std::max(yMax, static_cast<double>(m.atoms));
    }
    xMax = std::max(xMax * 1.08, 1e-3);
    yMin = std::max(0.0, yMin - 0.05 * (yMax - yMin) - 1.0);
    yMax = yMax + 0.05 * (yMax - yMin) + 1.0;
    const auto X = [&](double x) { return lo.x + static_cast<float>(x / xMax) * (hi.x - lo.x); };
    const auto Y = [&](double y) { return hi.y - static_cast<float>((y - yMin) / (yMax - yMin)) * (hi.y - lo.y); };

    const ImU32 frame = ImGui::GetColorU32(ImGuiCol_Border), text = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    const ImU32 grid = ImGui::GetColorU32(ImGuiCol_Border, 0.45f);
    draw->AddRect(lo, hi, frame);
    char buffer[32];
    for (int t = 0; t <= 4; ++t)
    {
        const double x = xMax * t / 4.0, y = yMin + (yMax - yMin) * t / 4.0;
        draw->AddLine(ImVec2(X(x), lo.y), ImVec2(X(x), hi.y), grid);
        draw->AddLine(ImVec2(lo.x, Y(y)), ImVec2(hi.x, Y(y)), grid);
        std::snprintf(buffer, sizeof(buffer), "%.2f", x);
        const ImVec2 size = ImGui::CalcTextSize(buffer);
        draw->AddText(ImVec2(std::clamp(X(x) - size.x * 0.5f, lo.x - size.x * 0.5f, hi.x - size.x), hi.y + 2.0f), text, buffer);
        std::snprintf(buffer, sizeof(buffer), "%.0f", y);
        draw->AddText(ImVec2(lo.x - ImGui::CalcTextSize(buffer).x - 4.0f, Y(y) - ImGui::GetTextLineHeight() * 0.5f), text, buffer);
    }
    const char* xLabel = "mean strain (%)";
    draw->AddText(ImVec2(0.5f * (lo.x + hi.x - ImGui::CalcTextSize(xLabel).x), hi.y + ImGui::GetTextLineHeight() + 2.0f), text, xLabel);

    // The Pareto front, then every match; hover shows it, a click selects it.
    const ImU32 accent = ImGui::GetColorU32(ImGuiCol_CheckMark), muted = ImGui::GetColorU32(ImGuiCol_TextDisabled, 0.75f);
    ImVec2 previous(0, 0);
    bool first = true;
    for (const auto& m : m_matches)
    {
        if (!m.pareto) continue;
        const ImVec2 p(X(m.meanStrain * 100.0), Y(m.atoms));
        if (!first) draw->AddLine(previous, p, accent, 1.5f);
        previous = p;
        first = false;
    }
    const ImVec2 mouse = ImGui::GetIO().MousePos;
    int hovered = -1;
    float nearest = 7.0f * 7.0f;
    for (std::size_t i = 0; i < m_matches.size(); ++i)
    {
        const auto& m = m_matches[i];
        const ImVec2 p(X(m.meanStrain * 100.0), Y(m.atoms));
        draw->AddCircleFilled(p, m.pareto ? 4.0f : 3.0f, m.pareto ? accent : muted);
        const float dx = mouse.x - p.x, dy = mouse.y - p.y;
        if (dx * dx + dy * dy < nearest) { nearest = dx * dx + dy * dy; hovered = static_cast<int>(i); }
    }
    if (m_selected >= 0)
    {
        const auto& m = m_matches[static_cast<std::size_t>(m_selected)];
        draw->AddCircle(ImVec2(X(m.meanStrain * 100.0), Y(m.atoms)), 7.0f, ImGui::GetColorU32(ImGuiCol_Text), 0, 2.0f);
    }
    ImGui::InvisibleButton("##plotarea", ImVec2(std::max(1.0f, available.x), std::max(1.0f, available.y)));
    if (hovered >= 0 && ImGui::IsItemHovered())
    {
        const auto& m = m_matches[static_cast<std::size_t>(hovered)];
        ImGui::SetTooltip("Match %d: %d atoms, strain %.3f%%, area %.1f \xc3\x85\xc2\xb2, rotation %.1f\xc2\xb0", hovered, m.atoms, m.meanStrain * 100.0,
                          m.area, m.thetaDeg);
        if (ImGui::IsItemClicked()) select(hovered);
    }
    ImGui::EndChild();
}

void InterfaceBuilderDialog::drawPreview(float width, float height)
{
    responsive::beginChild("##preview", ImVec2(width, height), ImGuiChildFlags_Borders, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    if (m_preview.atoms.empty() || !m_glReady)
    {
        dialogLayout::emptyState("The selected match is shown here.");
        ImGui::EndChild();
        return;
    }
    ImGui::TextDisabled("%zu atoms; drag to turn, scroll to zoom", m_preview.atoms.size());
    const ImVec2 available = ImGui::GetContentRegionAvail();
    const int w = std::max(1, static_cast<int>(available.x)), h = std::max(1, static_cast<int>(available.y));
    if (w != m_fboWidth || h != m_fboHeight || !m_fbo)
    {
        if (m_fbo) glDeleteFramebuffers(1, &m_fbo);
        if (m_texture) glDeleteTextures(1, &m_texture);
        if (m_depth) glDeleteRenderbuffers(1, &m_depth);
        glGenTextures(1, &m_texture);
        glBindTexture(GL_TEXTURE_2D, m_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
        glGenRenderbuffers(1, &m_depth);
        glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glGenFramebuffers(1, &m_fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texture, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_depth);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        m_fboWidth = w;
        m_fboHeight = h;
    }
    if (m_renderer && m_previewBuffers.atomCount > 0)
    {
        GLint previousFbo = 0, previousViewport[4];
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &previousFbo);
        glGetIntegerv(GL_VIEWPORT, previousViewport);
        Camera camera;
        camera.yaw = m_yaw;
        camera.pitch = m_pitch;
        camera.distance = m_distance;
        FrameView frame;
        frame.framebufferWidth = w;
        frame.framebufferHeight = h;
        buildFrameView(camera, m_previewBuffers, true, frame);
        glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
        glViewport(0, 0, w, h);
        glEnable(GL_DEPTH_TEST);
        const ImVec4& background = ImGui::GetStyle().Colors[ImGuiCol_ChildBg].w > 0.0f ? ImGui::GetStyle().Colors[ImGuiCol_ChildBg] : ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        glClearColor(background.x, background.y, background.z, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        m_renderer->drawBonds(frame.projection, frame.view, frame.lightPosition, frame.cameraPosition, m_previewBuffers.tabCylinderVAO,
                              m_previewBuffers.tabCylinderVertexCount, m_previewBuffers.bondCount);
        m_renderer->drawAtoms(frame.projection, frame.view, frame.lightMVP, frame.lightPosition, frame.cameraPosition, m_shadow,
                              m_previewBuffers.tabSphereVAO, m_previewBuffers.tabSphereIndexCount, m_previewBuffers.atomCount);
        m_renderer->drawBoxLines(frame.projection, frame.view, m_previewBuffers.lineVAO, m_previewBuffers.boxLines.size());
        glBindFramebuffer(GL_FRAMEBUFFER, static_cast<GLuint>(previousFbo));
        glViewport(previousViewport[0], previousViewport[1], previousViewport[2], previousViewport[3]);
    }
    ImGui::Image(static_cast<ImTextureID>(static_cast<intptr_t>(m_texture)), ImVec2(static_cast<float>(w), static_cast<float>(h)), ImVec2(0, 1), ImVec2(1, 0));
    if (ImGui::IsItemHovered())
    {
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
        {
            const ImVec2 delta = ImGui::GetIO().MouseDelta;
            m_yaw -= delta.x * 0.5f;
            m_pitch = std::clamp(m_pitch + delta.y * 0.5f, -89.0f, 89.0f);
        }
        if (const float wheel = ImGui::GetIO().MouseWheel; wheel != 0.0f)
            m_distance = std::clamp(m_distance * (1.0f - 0.1f * wheel), Camera::kMinDistance, Camera::kMaxDistance);
    }
    ImGui::EndChild();
}

void InterfaceBuilderDialog::drawAssembly()
{
    // Two columns of label + field, as one compact block under the results.
    if (!ImGui::BeginTable("##assembly", 4, ImGuiTableFlags_SizingStretchProp)) return;
    ImGui::TableSetupColumn("l1", ImGuiTableColumnFlags_WidthStretch, 0.22f);
    ImGui::TableSetupColumn("f1", ImGuiTableColumnFlags_WidthStretch, 0.28f);
    ImGui::TableSetupColumn("l2", ImGuiTableColumnFlags_WidthStretch, 0.22f);
    ImGui::TableSetupColumn("f2", ImGuiTableColumnFlags_WidthStretch, 0.28f);
    const auto label = [](int column, const char* text) {
        ImGui::TableSetColumnIndex(column);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(text);
        ImGui::TableSetColumnIndex(column + 1);
        ImGui::SetNextItemWidth(-FLT_MIN);
    };
    ImGui::TableNextRow();
    label(0, "Strained crystal");
    int on = m_assembly.strainOn == iface::StrainOn::B ? 0 : m_assembly.strainOn == iface::StrainOn::A ? 1 : 2;
    if (ImGui::Combo("##strainon", &on, "B (A keeps its lattice)\0A (B keeps its lattice)\0Both (to their average)\0"))
    {
        m_assembly.strainOn = on == 0 ? iface::StrainOn::B : on == 1 ? iface::StrainOn::A : iface::StrainOn::Both;
        m_previewDirty = true;
    }
    float gap = static_cast<float>(m_assembly.gap);
    label(2, "Gap (\xc3\x85)");
    if (ImGui::InputFloat("##gap", &gap, 0.1f, 0.5f, "%.2f"))
    {
        m_assembly.gap = std::clamp(gap, 0.0f, 50.0f);
        m_previewDirty = true;
    }
    dialogLayout::hint("Distance between the facing atomic layers of A and B.");
    ImGui::TableNextRow();
    float vacuum = static_cast<float>(m_assembly.vacuum);
    label(0, "Vacuum (\xc3\x85)");
    if (ImGui::InputFloat("##vacuum", &vacuum, 1.0f, 5.0f, "%.1f"))
    {
        m_assembly.vacuum = std::clamp(vacuum, 0.0f, 500.0f);
        m_previewDirty = true;
    }
    dialogLayout::hint("0 keeps the cell periodic along the normal (two interfaces, each at the gap); more adds vacuum above and below.");
    float shift[2] = {static_cast<float>(m_assembly.shift[0]), static_cast<float>(m_assembly.shift[1])};
    label(2, "Shift of B (a, b)");
    if (ImGui::InputFloat2("##shift", shift, "%.3f"))
    {
        m_assembly.shift[0] = shift[0];
        m_assembly.shift[1] = shift[1];
        m_previewDirty = true;
    }
    dialogLayout::hint("In-plane shift of B as fractions of the interface cell vectors: sets the registry of the two crystals.");
    ImGui::TableNextRow();
    int repeats[2] = {m_repeatX, m_repeatY};
    label(0, "Repeat in plane");
    if (ImGui::InputInt2("##repeat", repeats))
    {
        m_repeatX = std::clamp(repeats[0], 1, 50);
        m_repeatY = std::clamp(repeats[1], 1, 50);
        m_previewDirty = true;
    }
    ImGui::EndTable();
}

void InterfaceBuilderDialog::drawDialog(Structure& structure, const std::vector<glm::vec3>& elementColors,
                                        const std::vector<float>& elementRadii, const std::vector<float>& elementShininess,
                                        const std::function<void(Structure&)>& updateBuffers)
{
    (void)elementColors;
    if (!m_pendingDropA.empty()) { loadFile(m_pendingDropA, m_a); m_pendingDropA.clear(); }
    if (!m_pendingDropB.empty()) { loadFile(m_pendingDropB, m_b); m_pendingDropB.clear(); }
    // Editing a step: search with its settings and select its match.
    if (m_stepSearch && !m_a.crystal.atoms.empty() && !m_b.crystal.atoms.empty())
    {
        m_stepSearch = false;
        const int pick = m_stepPick;
        runSearch();
        if (pick >= 0) select(pick);
    }
    if (m_previewDirty) rebuildPreview(elementRadii, elementShininess);

    if (m_openRequested)
    {
        ImGui::OpenPopup("Interface Builder");
        m_openRequested = false;
    }
    m_isOpen = ImGui::IsPopupOpen("Interface Builder");

    // Sized once when it opens; it can then be resized like any dialog.
    responsive::windowSize(ImVec2(1240.0f, 820.0f), ImGuiCond_Appearing);
    responsive::windowConstraints(ImVec2(900.0f, 600.0f), ImVec2(FLT_MAX, FLT_MAX));
    bool dialogOpen = true;
    if (!responsive::beginModal("Interface Builder", &dialogOpen, 0))
    {
        if (responsive::lastModalMinimized()) return;  // hidden, not closed: keeps its state
        if (!dialogOpen || !m_isOpen)
        {
            m_isOpen = false;
            m_step.finish();
        }
        return;
    }
    m_isOpen = true;
    if (m_step.active()) dialogLayout::stepBanner(m_step.step.c_str());

    if (const auto path = m_picker.draw())
        loadFile(*path, m_pickerSlot == 1 ? m_a : m_b);

    const ImGuiStyle& style = ImGui::GetStyle();
    // The body fills exactly what the footer (a rule and one row of buttons)
    // leaves, so the dialog itself never scrolls.
    const float footer = ImGui::GetFrameHeight() + 3.0f * style.ItemSpacing.y + 2.0f;
    const float bodyHeight = std::max(responsive::dp(300), ImGui::GetContentRegionAvail().y - footer);
    responsive::beginChild("##body", ImVec2(ImGui::GetContentRegionAvail().x, bodyHeight), 0, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    if (ImGui::BeginTable("##layout", 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerV, ImVec2(0, bodyHeight)))
    {
        ImGui::TableSetupColumn("settings", ImGuiTableColumnFlags_WidthFixed, responsive::dp(370));
        ImGui::TableSetupColumn("results", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableNextRow();

        // Left: the two crystals and the search.
        ImGui::TableSetColumnIndex(0);
        responsive::beginChild("##settings", ImVec2(0, bodyHeight), 0);
        drawInput("Crystal A (bottom)", m_a, 1, structure);
        ImGui::Spacing();
        drawInput("Crystal B (top)", m_b, 2, structure);
        ImGui::Spacing();
        drawSearchSettings();
        ImGui::Spacing();
        const bool ready = !m_a.crystal.atoms.empty() && !m_b.crystal.atoms.empty();
        ImGui::BeginDisabled(!ready);
        if (dialogLayout::primaryButton("Find Matches", ImVec2(-FLT_MIN, 0)) || (ready && uiScript::take("interface.find"))) runSearch();
        ImGui::EndDisabled();
        if (!ready) ImGui::TextDisabled("Set both crystals first.");
        else if (m_searched && m_stale)
        {
            ImGui::PushStyleColor(ImGuiCol_Text, themeStatusWarn());
            ImGui::TextWrapped("Settings changed: find matches again.");
            ImGui::PopStyleColor();
        }
        ImGui::EndChild();

        // Right: the matches, the plot, the preview and how the interface is put together.
        ImGui::TableSetColumnIndex(1);
        responsive::beginChild("##results", ImVec2(0, bodyHeight), 0, ImGuiWindowFlags_NoScrollbar);
        if (m_matches.empty())
        {
            if (m_searched) dialogLayout::note(dialogLayout::Note::Warning, m_searchMessage.c_str());
            dialogLayout::emptyState(m_searched ? "No matches yet." : "Set both crystals and their planes, then Find Matches.\n"
                                                                     "The matches are listed with the fewest atoms first.");
        }
        else
        {
            ImGui::AlignTextToFramePadding();
            ImGui::Text("%s", m_searchMessage.c_str());
            ImGui::SameLine();
            ImGui::TextDisabled("(fewest atoms first; dimmed atoms: another match has as few and less strain;\xe2\x80\xa2 suggested)");
            const float assemblyHeight = 3.0f * ImGui::GetFrameHeightWithSpacing() + 2.0f * ImGui::GetTextLineHeightWithSpacing() + 3.0f * style.ItemSpacing.y;
            const float rest = ImGui::GetContentRegionAvail().y - assemblyHeight;
            const float tableHeight = std::max(responsive::dp(120), rest * 0.42f);
            drawResultsTable(tableHeight);
            const float lowerHeight = std::max(responsive::dp(160), ImGui::GetContentRegionAvail().y - assemblyHeight);
            const float half = (ImGui::GetContentRegionAvail().x - style.ItemSpacing.x) * 0.5f;
            drawPlot(half * 0.85f, lowerHeight);
            ImGui::SameLine();
            drawPreview(ImGui::GetContentRegionAvail().x, lowerHeight);
            dialogLayout::section("Assembly");
            drawAssembly();
        }
        ImGui::EndChild();
        ImGui::EndTable();
    }
    ImGui::EndChild();

    // Footer: what will be built, then [Close] [Build].
    std::string status;
    if (m_selected >= 0 && m_selected < static_cast<int>(m_matches.size()))
    {
        const auto& m = m_matches[static_cast<std::size_t>(m_selected)];
        char text[200];
        std::snprintf(text, sizeof(text), "Match %d: %zu atoms, mean strain %.3f%% on %s, rotation %.1f\xc2\xb0", m_selected, m_preview.atoms.size(),
                      m.meanStrain * 100.0, m_assembly.strainOn == iface::StrainOn::Both ? "both crystals" : m_assembly.strainOn == iface::StrainOn::A ? "A" : "B",
                      m.thetaDeg);
        status = text;
    }
    ImGui::Separator();
    const bool canBuild = m_step.active() || !m_preview.atoms.empty();
    int pressed = dialogLayout::actionBar(status, dialogLayout::Note::Info, m_step.active() ? "Update Step" : "Build", "Close", canBuild);
    if (pressed == 0 && canBuild && uiScript::take("interface.build")) pressed = 1;
    if (pressed == 1)
    {
        if (m_step.active())
        {
            m_step.commitOptions(stepOptions());
            m_step.finish();
        }
        else
        {
            structure = m_preview;
            updateBuffers(structure);
        }
        m_isOpen = false;
        ImGui::CloseCurrentPopup();
    }
    else if (pressed == 2)
    {
        m_isOpen = false;
        m_step.finish();
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}
