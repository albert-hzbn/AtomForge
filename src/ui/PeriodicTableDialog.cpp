#include "PeriodicTableDialog.h"

#include "imgui.h"

#include <algorithm>
#include <cstdio>

// ---------------------------------------------------------------------------
// Element data
// ---------------------------------------------------------------------------
namespace {

struct ElemDef {
    int         z;
    const char* sym;
    const char* nm;
    int         row; // 1-7 main table, 8 = lanthanides row, 9 = actinides row
    int         col; // 1-18
    int         cat; // see kCatColors
};

// cat: 0=H, 1=Alkali, 2=AlkEarth, 3=Lanthanide, 4=Actinide,
//      5=TransMetal, 6=PostTrans, 7=Metalloid, 8=Nonmetal,
//      9=Halogen, 10=NobleGas

static const ElemDef kElements[] = {
    // Period 1
    {  1,"H",  "Hydrogen",       1,  1,  0},
    {  2,"He", "Helium",         1, 18, 10},
    // Period 2
    {  3,"Li", "Lithium",        2,  1,  1},
    {  4,"Be", "Beryllium",      2,  2,  2},
    {  5,"B",  "Boron",          2, 13,  7},
    {  6,"C",  "Carbon",         2, 14,  8},
    {  7,"N",  "Nitrogen",       2, 15,  8},
    {  8,"O",  "Oxygen",         2, 16,  8},
    {  9,"F",  "Fluorine",       2, 17,  9},
    { 10,"Ne", "Neon",           2, 18, 10},
    // Period 3
    { 11,"Na", "Sodium",         3,  1,  1},
    { 12,"Mg", "Magnesium",      3,  2,  2},
    { 13,"Al", "Aluminum",       3, 13,  6},
    { 14,"Si", "Silicon",        3, 14,  7},
    { 15,"P",  "Phosphorus",     3, 15,  8},
    { 16,"S",  "Sulfur",         3, 16,  8},
    { 17,"Cl", "Chlorine",       3, 17,  9},
    { 18,"Ar", "Argon",          3, 18, 10},
    // Period 4
    { 19,"K",  "Potassium",      4,  1,  1},
    { 20,"Ca", "Calcium",        4,  2,  2},
    { 21,"Sc", "Scandium",       4,  3,  5},
    { 22,"Ti", "Titanium",       4,  4,  5},
    { 23,"V",  "Vanadium",       4,  5,  5},
    { 24,"Cr", "Chromium",       4,  6,  5},
    { 25,"Mn", "Manganese",      4,  7,  5},
    { 26,"Fe", "Iron",           4,  8,  5},
    { 27,"Co", "Cobalt",         4,  9,  5},
    { 28,"Ni", "Nickel",         4, 10,  5},
    { 29,"Cu", "Copper",         4, 11,  5},
    { 30,"Zn", "Zinc",           4, 12,  5},
    { 31,"Ga", "Gallium",        4, 13,  6},
    { 32,"Ge", "Germanium",      4, 14,  7},
    { 33,"As", "Arsenic",        4, 15,  7},
    { 34,"Se", "Selenium",       4, 16,  8},
    { 35,"Br", "Bromine",        4, 17,  9},
    { 36,"Kr", "Krypton",        4, 18, 10},
    // Period 5
    { 37,"Rb", "Rubidium",       5,  1,  1},
    { 38,"Sr", "Strontium",      5,  2,  2},
    { 39,"Y",  "Yttrium",        5,  3,  5},
    { 40,"Zr", "Zirconium",      5,  4,  5},
    { 41,"Nb", "Niobium",        5,  5,  5},
    { 42,"Mo", "Molybdenum",     5,  6,  5},
    { 43,"Tc", "Technetium",     5,  7,  5},
    { 44,"Ru", "Ruthenium",      5,  8,  5},
    { 45,"Rh", "Rhodium",        5,  9,  5},
    { 46,"Pd", "Palladium",      5, 10,  5},
    { 47,"Ag", "Silver",         5, 11,  5},
    { 48,"Cd", "Cadmium",        5, 12,  5},
    { 49,"In", "Indium",         5, 13,  6},
    { 50,"Sn", "Tin",            5, 14,  6},
    { 51,"Sb", "Antimony",       5, 15,  7},
    { 52,"Te", "Tellurium",      5, 16,  7},
    { 53,"I",  "Iodine",         5, 17,  9},
    { 54,"Xe", "Xenon",          5, 18, 10},
    // Period 6
    { 55,"Cs", "Cesium",         6,  1,  1},
    { 56,"Ba", "Barium",         6,  2,  2},
    { 72,"Hf", "Hafnium",        6,  4,  5},
    { 73,"Ta", "Tantalum",       6,  5,  5},
    { 74,"W",  "Tungsten",       6,  6,  5},
    { 75,"Re", "Rhenium",        6,  7,  5},
    { 76,"Os", "Osmium",         6,  8,  5},
    { 77,"Ir", "Iridium",        6,  9,  5},
    { 78,"Pt", "Platinum",       6, 10,  5},
    { 79,"Au", "Gold",           6, 11,  5},
    { 80,"Hg", "Mercury",        6, 12,  5},
    { 81,"Tl", "Thallium",       6, 13,  6},
    { 82,"Pb", "Lead",           6, 14,  6},
    { 83,"Bi", "Bismuth",        6, 15,  6},
    { 84,"Po", "Polonium",       6, 16,  6},
    { 85,"At", "Astatine",       6, 17,  9},
    { 86,"Rn", "Radon",          6, 18, 10},
    // Period 7
    { 87,"Fr", "Francium",       7,  1,  1},
    { 88,"Ra", "Radium",         7,  2,  2},
    {104,"Rf", "Rutherfordium",  7,  4,  5},
    {105,"Db", "Dubnium",        7,  5,  5},
    {106,"Sg", "Seaborgium",     7,  6,  5},
    {107,"Bh", "Bohrium",        7,  7,  5},
    {108,"Hs", "Hassium",        7,  8,  5},
    {109,"Mt", "Meitnerium",     7,  9,  5},
    {110,"Ds", "Darmstadtium",   7, 10,  5},
    {111,"Rg", "Roentgenium",    7, 11,  5},
    {112,"Cn", "Copernicium",    7, 12,  5},
    {113,"Nh", "Nihonium",       7, 13,  6},
    {114,"Fl", "Flerovium",      7, 14,  6},
    {115,"Mc", "Moscovium",      7, 15,  6},
    {116,"Lv", "Livermorium",    7, 16,  6},
    {117,"Ts", "Tennessine",     7, 17,  9},
    {118,"Og", "Oganesson",      7, 18, 10},
    // Lanthanides  (row 8, cols 3-17)
    { 57,"La", "Lanthanum",      8,  3,  3},
    { 58,"Ce", "Cerium",         8,  4,  3},
    { 59,"Pr", "Praseodymium",   8,  5,  3},
    { 60,"Nd", "Neodymium",      8,  6,  3},
    { 61,"Pm", "Promethium",     8,  7,  3},
    { 62,"Sm", "Samarium",       8,  8,  3},
    { 63,"Eu", "Europium",       8,  9,  3},
    { 64,"Gd", "Gadolinium",     8, 10,  3},
    { 65,"Tb", "Terbium",        8, 11,  3},
    { 66,"Dy", "Dysprosium",     8, 12,  3},
    { 67,"Ho", "Holmium",        8, 13,  3},
    { 68,"Er", "Erbium",         8, 14,  3},
    { 69,"Tm", "Thulium",        8, 15,  3},
    { 70,"Yb", "Ytterbium",      8, 16,  3},
    { 71,"Lu", "Lutetium",       8, 17,  3},
    // Actinides (row 9, cols 3-17)
    { 89,"Ac", "Actinium",       9,  3,  4},
    { 90,"Th", "Thorium",        9,  4,  4},
    { 91,"Pa", "Protactinium",   9,  5,  4},
    { 92,"U",  "Uranium",        9,  6,  4},
    { 93,"Np", "Neptunium",      9,  7,  4},
    { 94,"Pu", "Plutonium",      9,  8,  4},
    { 95,"Am", "Americium",      9,  9,  4},
    { 96,"Cm", "Curium",         9, 10,  4},
    { 97,"Bk", "Berkelium",      9, 11,  4},
    { 98,"Cf", "Californium",    9, 12,  4},
    { 99,"Es", "Einsteinium",    9, 13,  4},
    {100,"Fm", "Fermium",        9, 14,  4},
    {101,"Md", "Mendelevium",    9, 15,  4},
    {102,"No", "Nobelium",       9, 16,  4},
    {103,"Lr", "Lawrencium",     9, 17,  4},
};

constexpr int kElemCount = (int)(sizeof(kElements) / sizeof(kElements[0]));

// Category background colors
static const ImVec4 kCatColors[11] = {
    {0.78f, 0.78f, 0.20f, 1.0f}, // 0: H
    {0.92f, 0.40f, 0.25f, 1.0f}, // 1: Alkali metals
    {0.92f, 0.76f, 0.30f, 1.0f}, // 2: Alkaline earth
    {1.00f, 0.72f, 0.77f, 1.0f}, // 3: Lanthanides
    {1.00f, 0.78f, 0.60f, 1.0f}, // 4: Actinides
    {0.72f, 0.72f, 0.95f, 1.0f}, // 5: Transition metals
    {0.55f, 0.85f, 0.80f, 1.0f}, // 6: Post-transition metals
    {0.72f, 0.90f, 0.55f, 1.0f}, // 7: Metalloids
    {0.60f, 0.90f, 0.58f, 1.0f}, // 8: Nonmetals
    {0.55f, 0.90f, 0.70f, 1.0f}, // 9: Halogens
    {0.50f, 0.78f, 0.95f, 1.0f}, // 10: Noble gases
};

static ImU32 toU32(ImVec4 c)
{
    return ImGui::ColorConvertFloat4ToU32(c);
}

static ImVec4 brighter(ImVec4 c, float d)
{
    return { std::min(c.x + d, 1.0f),
             std::min(c.y + d, 1.0f),
             std::min(c.z + d, 1.0f),
             c.w };
}

} // namespace

// ---------------------------------------------------------------------------

void openPeriodicTable()
{
    ImGui::OpenPopup("Periodic Table##picker");
}

bool drawPeriodicTable(std::vector<ElementSelection>& outSelections)
{
    bool applied = false;

    // Static state for tracking single selection
    static int selectedAtomicNumber = -1;
    static bool isPopupOpen = false;

    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(900.0f, 480.0f), ImGuiCond_Always);

    bool open = true;
    if (!ImGui::BeginPopupModal("Periodic Table##picker", &open,
                                ImGuiWindowFlags_NoResize |
                                ImGuiWindowFlags_NoScrollbar))
    {
        // Dialog just closed - reset state
        if (isPopupOpen)
        {
            selectedAtomicNumber = -1;
            isPopupOpen = false;
        }
        return false;
    }

    // Dialog just opened - clear selection
    if (!isPopupOpen)
    {
        selectedAtomicNumber = -1;
        isPopupOpen = true;
    }

    const float cellW = 44.0f;
    const float cellH = 34.0f;
    const float padX  =  6.0f; // left margin
    const float padY  = ImGui::GetCursorPosY() + 4.0f; // start below title bar + small margin
    const float fGap  = 10.0f; // gap between main rows and f-block

    // Returns window-local Y for the top of a given row (1-9).
    auto rowLocalY = [&](int row) -> float {
        if (row <= 7)
            return padY + (row - 1) * cellH;
        return padY + 7.0f * cellH + fGap + (row - 8) * cellH;
    };

    ImDrawList* dl = ImGui::GetWindowDrawList();

    // ---- Render all 118 element buttons ----
    for (int i = 0; i < kElemCount; ++i)
    {
        const ElemDef& e = kElements[i];

        float localX = padX + (e.col - 1) * cellW;
        float localY = rowLocalY(e.row);

        ImGui::SetCursorPos(ImVec2(localX, localY));

        char btnId[16];
        std::snprintf(btnId, sizeof(btnId), "##z%d", e.z);

        bool clicked = ImGui::InvisibleButton(btnId, ImVec2(cellW - 1.0f, cellH - 1.0f));
        bool hovered = ImGui::IsItemHovered();
        
        // Single click to select (replaces previous selection)
        if (clicked)
        {
            selectedAtomicNumber = e.z;
        }

        // Draw background
        ImVec2 pMin = ImGui::GetItemRectMin();
        ImVec2 pMax = ImGui::GetItemRectMax();
        
        bool selected = (selectedAtomicNumber == e.z);
        ImVec4 bg;
        if (selected)
        {
            // Darker, more saturated color for selected element
            bg = brighter(kCatColors[e.cat], 0.35f);
        }
        else if (hovered)
        {
            bg = brighter(kCatColors[e.cat], 0.15f);
        }
        else
        {
            bg = kCatColors[e.cat];
        }
        
        dl->AddRectFilled(pMin, pMax, toU32(bg), 3.0f);
        
        // Thicker border for selected element
        if (selected)
        {
            dl->AddRect(pMin, pMax, IM_COL32(0, 0, 100, 200), 3.0f, 0, 2.0f);
        }
        else
        {
            dl->AddRect(pMin, pMax, IM_COL32(0, 0, 0, 55), 3.0f);
        }

        // Atomic number — small, top-left corner
        char zStr[8];
        std::snprintf(zStr, sizeof(zStr), "%d", e.z);
        dl->AddText(ImVec2(pMin.x + 2.0f, pMin.y + 1.0f),
                    IM_COL32(20, 20, 20, 180), zStr);

        // Symbol — centered
        float symW = ImGui::CalcTextSize(e.sym).x;
        float textH = ImGui::GetTextLineHeight();
        dl->AddText(ImVec2(pMin.x + ((cellW - 1.0f) - symW)  * 0.5f,
                           pMin.y + ((cellH - 1.0f) - textH) * 0.55f),
                    IM_COL32(0, 0, 0, 255), e.sym);

        if (hovered)
        {
            if (selected)
                ImGui::SetTooltip("%d  %s\n%s\n(Selected)", e.z, e.sym, e.nm);
            else
                ImGui::SetTooltip("%d  %s\n%s", e.z, e.sym, e.nm);
        }
    }

    // ---- "*" reference markers at row 6/7, col 3 ----
    static const struct { int row; const char* lbl; } kStars[2] = {
        {6, "57-71"}, {7, "89-103"}
    };
    for (int s = 0; s < 2; ++s)
    {
        float lx = padX + (3 - 1) * cellW;
        float ly = rowLocalY(kStars[s].row);

        ImGui::SetCursorPos(ImVec2(lx, ly));
        // Non-clickable dummy for layout
        ImGui::Dummy(ImVec2(cellW - 1.0f, cellH - 1.0f));

        ImVec2 pMin = ImGui::GetItemRectMin();
        ImVec2 pMax = ImGui::GetItemRectMax();
        dl->AddRectFilled(pMin, pMax, IM_COL32(80, 80, 80, 130), 3.0f);
        dl->AddRect(pMin, pMax, IM_COL32(0, 0, 0, 40), 3.0f);

        float lw   = ImGui::CalcTextSize(kStars[s].lbl).x;
        float lh   = ImGui::GetTextLineHeight();
        dl->AddText(ImVec2(pMin.x + ((cellW - 1.0f) - lw) * 0.5f,
                           pMin.y + ((cellH - 1.0f) - lh) * 0.5f),
                    IM_COL32(210, 210, 210, 255), kStars[s].lbl);
    }

    // ---- Selection summary and buttons ----
    float btnAreaY = rowLocalY(9) + cellH + 10.0f;
    
    // Display selected element
    ImGui::SetCursorPos(ImVec2(padX, btnAreaY));
    
    if (selectedAtomicNumber >= 0)
    {
        // Find the element name
        const char* elemName = "Unknown";
        for (int i = 0; i < kElemCount; ++i)
        {
            if (kElements[i].z == selectedAtomicNumber)
            {
                elemName = kElements[i].sym;
                break;
            }
        }
        ImGui::Text("Selected: Z=%d (%s)", selectedAtomicNumber, elemName);
    }
    else
    {
        ImGui::Text("Selected: None");
    }

    // Buttons
    float btnY = btnAreaY + 25.0f;
    ImGui::SetCursorPos(ImVec2(padX, btnY));
    
    if (ImGui::Button("Apply", ImVec2(80, 0)) && selectedAtomicNumber >= 0)
    {
        // Populate outSelections with the single selection
        outSelections.clear();
        for (int i = 0; i < kElemCount; ++i)
        {
            if (kElements[i].z == selectedAtomicNumber)
            {
                outSelections.push_back({kElements[i].sym, selectedAtomicNumber});
                break;
            }
        }
        applied = true;
        selectedAtomicNumber = -1;
        ImGui::CloseCurrentPopup();
    }

    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(80, 0)))
    {
        selectedAtomicNumber = -1;
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
    return applied;
}

bool drawPeriodicTableInlineSelector(int& selectedAtomicNumber)
{
    bool changed = false;

    const float cellW = 44.0f;
    const float cellH = 34.0f;
    const float padX  =  6.0f;
    const float padY  = 10.0f;
    const float fGap  = 10.0f;

    auto rowLocalY = [&](int row) -> float {
        if (row <= 7)
            return padY + (row - 1) * cellH;
        return padY + 7.0f * cellH + fGap + (row - 8) * cellH;
    };

    ImGui::BeginChild("PeriodicTableInline##selector", ImVec2(900.0f, 390.0f), true);
    ImDrawList* dl = ImGui::GetWindowDrawList();

    for (int i = 0; i < kElemCount; ++i)
    {
        const ElemDef& e = kElements[i];

        float localX = padX + (e.col - 1) * cellW;
        float localY = rowLocalY(e.row);

        ImGui::SetCursorPos(ImVec2(localX, localY));

        char btnId[24];
        std::snprintf(btnId, sizeof(btnId), "##inline_z%d", e.z);

        bool clicked = ImGui::InvisibleButton(btnId, ImVec2(cellW - 1.0f, cellH - 1.0f));
        bool hovered = ImGui::IsItemHovered();
        if (clicked && selectedAtomicNumber != e.z)
        {
            selectedAtomicNumber = e.z;
            changed = true;
        }

        ImVec2 pMin = ImGui::GetItemRectMin();
        ImVec2 pMax = ImGui::GetItemRectMax();
        bool selected = (selectedAtomicNumber == e.z);

        ImVec4 bg;
        if (selected)
            bg = brighter(kCatColors[e.cat], 0.35f);
        else if (hovered)
            bg = brighter(kCatColors[e.cat], 0.15f);
        else
            bg = kCatColors[e.cat];

        dl->AddRectFilled(pMin, pMax, toU32(bg), 3.0f);
        if (selected)
            dl->AddRect(pMin, pMax, IM_COL32(0, 0, 100, 200), 3.0f, 0, 2.0f);
        else
            dl->AddRect(pMin, pMax, IM_COL32(0, 0, 0, 55), 3.0f);

        char zStr[8];
        std::snprintf(zStr, sizeof(zStr), "%d", e.z);
        dl->AddText(ImVec2(pMin.x + 2.0f, pMin.y + 1.0f), IM_COL32(20, 20, 20, 180), zStr);

        float symW = ImGui::CalcTextSize(e.sym).x;
        float textH = ImGui::GetTextLineHeight();
        dl->AddText(ImVec2(pMin.x + ((cellW - 1.0f) - symW)  * 0.5f,
                           pMin.y + ((cellH - 1.0f) - textH) * 0.55f),
                    IM_COL32(0, 0, 0, 255), e.sym);

        if (hovered)
            ImGui::SetTooltip("%d  %s\n%s", e.z, e.sym, e.nm);
    }

    static const struct { int row; const char* lbl; } kStars[2] = {
        {6, "57-71"}, {7, "89-103"}
    };
    for (int s = 0; s < 2; ++s)
    {
        float lx = padX + (3 - 1) * cellW;
        float ly = rowLocalY(kStars[s].row);

        ImGui::SetCursorPos(ImVec2(lx, ly));
        ImGui::Dummy(ImVec2(cellW - 1.0f, cellH - 1.0f));

        ImVec2 pMin = ImGui::GetItemRectMin();
        ImVec2 pMax = ImGui::GetItemRectMax();
        dl->AddRectFilled(pMin, pMax, IM_COL32(80, 80, 80, 130), 3.0f);
        dl->AddRect(pMin, pMax, IM_COL32(0, 0, 0, 40), 3.0f);

        float lw = ImGui::CalcTextSize(kStars[s].lbl).x;
        float lh = ImGui::GetTextLineHeight();
        dl->AddText(ImVec2(pMin.x + ((cellW - 1.0f) - lw) * 0.5f,
                           pMin.y + ((cellH - 1.0f) - lh) * 0.5f),
                    IM_COL32(210, 210, 210, 255), kStars[s].lbl);
    }

    ImGui::EndChild();
    return changed;
}
