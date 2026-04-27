#include "ElementData.h"
#include "io/StructureLoader.h"  // getDefaultElementColor

const char* elementSymbol(int z)
{
    static const char* kSymbols[119] = {
        "",
        "H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
        "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
        "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
        "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
        "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
        "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
        "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
        "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    };

    if (z >= 1 && z <= 118)
        return kSymbols[z];
    return "?";
}

const char* elementName(int z)
{
    static const char* kNames[119] = {
        "",
        "Hydrogen","Helium","Lithium","Beryllium","Boron",
        "Carbon","Nitrogen","Oxygen","Fluorine","Neon",
        "Sodium","Magnesium","Aluminum","Silicon","Phosphorus",
        "Sulfur","Chlorine","Argon","Potassium","Calcium",
        "Scandium","Titanium","Vanadium","Chromium","Manganese",
        "Iron","Cobalt","Nickel","Copper","Zinc",
        "Gallium","Germanium","Arsenic","Selenium","Bromine",
        "Krypton","Rubidium","Strontium","Yttrium","Zirconium",
        "Niobium","Molybdenum","Technetium","Ruthenium","Rhodium",
        "Palladium","Silver","Cadmium","Indium","Tin",
        "Antimony","Tellurium","Iodine","Xenon","Cesium",
        "Barium","Lanthanum","Cerium","Praseodymium","Neodymium",
        "Promethium","Samarium","Europium","Gadolinium","Terbium",
        "Dysprosium","Holmium","Erbium","Thulium","Ytterbium",
        "Lutetium","Hafnium","Tantalum","Tungsten","Rhenium",
        "Osmium","Iridium","Platinum","Gold","Mercury",
        "Thallium","Lead","Bismuth","Polonium","Astatine",
        "Radon","Francium","Radium","Actinium","Thorium",
        "Protactinium","Uranium","Neptunium","Plutonium","Americium",
        "Curium","Berkelium","Californium","Einsteinium","Fermium",
        "Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium",
        "Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium",
        "Roentgenium","Copernicium","Nihonium","Flerovium","Moscovium",
        "Livermorium","Tennessine","Oganesson"
    };

    if (z >= 1 && z <= 118)
        return kNames[z];
    return "Unknown";
}

std::vector<float> makeLiteratureCovalentRadii()
{
    // Covalent radii (Angstrom) from Cordero et al., Dalton Trans. 2008.
    std::vector<float> radii(119, 1.60f);
    radii[0] = 1.0f;

    radii[1]  = 0.31f; radii[2]  = 0.28f; radii[3]  = 1.28f; radii[4]  = 0.96f; radii[5]  = 0.84f;
    radii[6]  = 0.76f; radii[7]  = 0.71f; radii[8]  = 0.66f; radii[9]  = 0.57f; radii[10] = 0.58f;
    radii[11] = 1.66f; radii[12] = 1.41f; radii[13] = 1.21f; radii[14] = 1.11f; radii[15] = 1.07f;
    radii[16] = 1.05f; radii[17] = 1.02f; radii[18] = 1.06f; radii[19] = 2.03f; radii[20] = 1.76f;
    radii[21] = 1.70f; radii[22] = 1.60f; radii[23] = 1.53f; radii[24] = 1.39f; radii[25] = 1.39f;
    radii[26] = 1.32f; radii[27] = 1.26f; radii[28] = 1.24f; radii[29] = 1.32f; radii[30] = 1.22f;
    radii[31] = 1.22f; radii[32] = 1.20f; radii[33] = 1.19f; radii[34] = 1.20f; radii[35] = 1.20f;
    radii[36] = 1.16f; radii[37] = 2.20f; radii[38] = 1.95f; radii[39] = 1.90f; radii[40] = 1.75f;
    radii[41] = 1.64f; radii[42] = 1.54f; radii[43] = 1.47f; radii[44] = 1.46f; radii[45] = 1.42f;
    radii[46] = 1.39f; radii[47] = 1.45f; radii[48] = 1.44f; radii[49] = 1.42f; radii[50] = 1.39f;
    radii[51] = 1.39f; radii[52] = 1.38f; radii[53] = 1.39f; radii[54] = 1.40f; radii[55] = 2.44f;
    radii[56] = 2.15f; radii[57] = 2.07f; radii[58] = 2.04f; radii[59] = 2.03f; radii[60] = 2.01f;
    radii[61] = 1.99f; radii[62] = 1.98f; radii[63] = 1.98f; radii[64] = 1.96f; radii[65] = 1.94f;
    radii[66] = 1.92f; radii[67] = 1.92f; radii[68] = 1.89f; radii[69] = 1.90f; radii[70] = 1.87f;
    radii[71] = 1.87f; radii[72] = 1.75f; radii[73] = 1.70f; radii[74] = 1.62f; radii[75] = 1.51f;
    radii[76] = 1.44f; radii[77] = 1.41f; radii[78] = 1.36f; radii[79] = 1.36f; radii[80] = 1.32f;
    radii[81] = 1.45f; radii[82] = 1.46f; radii[83] = 1.48f; radii[84] = 1.40f; radii[85] = 1.50f;
    radii[86] = 1.50f; radii[87] = 2.60f; radii[88] = 2.21f; radii[89] = 2.15f; radii[90] = 2.06f;
    radii[91] = 2.00f; radii[92] = 1.96f; radii[93] = 1.90f; radii[94] = 1.87f; radii[95] = 1.80f;
    radii[96] = 1.69f;

    return radii;
}

std::vector<glm::vec3> makeDefaultElementColors()
{
    std::vector<glm::vec3> colors(119, glm::vec3(0.7f, 0.7f, 0.7f));
    for (int z = 1; z <= 118; ++z)
    {
        float r, g, b;
        getDefaultElementColor(z, r, g, b);
        colors[z] = glm::vec3(r, g, b);
    }
    return colors;
}

double elementAtomicMass(int z)
{
    // Standard atomic weights (IUPAC 2021), indexed 1–118.
    // Unstable/synthetic elements use mass of most stable isotope.
    static const double kMasses[119] = {
        1.0,          // 0 unused
        1.008,        // 1  H
        4.0026,       // 2  He
        6.941,        // 3  Li
        9.0122,       // 4  Be
        10.811,       // 5  B
        12.011,       // 6  C
        14.007,       // 7  N
        15.999,       // 8  O
        18.998,       // 9  F
        20.180,       // 10 Ne
        22.990,       // 11 Na
        24.305,       // 12 Mg
        26.982,       // 13 Al
        28.086,       // 14 Si
        30.974,       // 15 P
        32.060,       // 16 S
        35.453,       // 17 Cl
        39.948,       // 18 Ar
        39.098,       // 19 K
        40.078,       // 20 Ca
        44.956,       // 21 Sc
        47.867,       // 22 Ti
        50.942,       // 23 V
        51.996,       // 24 Cr
        54.938,       // 25 Mn
        55.845,       // 26 Fe
        58.933,       // 27 Co
        58.693,       // 28 Ni
        63.546,       // 29 Cu
        65.38,        // 30 Zn
        69.723,       // 31 Ga
        72.630,       // 32 Ge
        74.922,       // 33 As
        78.971,       // 34 Se
        79.904,       // 35 Br
        83.798,       // 36 Kr
        85.468,       // 37 Rb
        87.62,        // 38 Sr
        88.906,       // 39 Y
        91.224,       // 40 Zr
        92.906,       // 41 Nb
        95.96,        // 42 Mo
        98.0,         // 43 Tc
        101.07,       // 44 Ru
        102.906,      // 45 Rh
        106.42,       // 46 Pd
        107.868,      // 47 Ag
        112.411,      // 48 Cd
        114.818,      // 49 In
        118.710,      // 50 Sn
        121.760,      // 51 Sb
        127.60,       // 52 Te
        126.904,      // 53 I
        131.293,      // 54 Xe
        132.905,      // 55 Cs
        137.327,      // 56 Ba
        138.905,      // 57 La
        140.116,      // 58 Ce
        140.908,      // 59 Pr
        144.242,      // 60 Nd
        145.0,        // 61 Pm
        150.36,       // 62 Sm
        151.964,      // 63 Eu
        157.25,       // 64 Gd
        158.925,      // 65 Tb
        162.500,      // 66 Dy
        164.930,      // 67 Ho
        167.259,      // 68 Er
        168.934,      // 69 Tm
        173.045,      // 70 Yb
        174.967,      // 71 Lu
        178.49,       // 72 Hf
        180.948,      // 73 Ta
        183.84,       // 74 W
        186.207,      // 75 Re
        190.23,       // 76 Os
        192.217,      // 77 Ir
        195.084,      // 78 Pt
        196.967,      // 79 Au
        200.592,      // 80 Hg
        204.38,       // 81 Tl
        207.2,        // 82 Pb
        208.980,      // 83 Bi
        209.0,        // 84 Po
        210.0,        // 85 At
        222.0,        // 86 Rn
        223.0,        // 87 Fr
        226.0,        // 88 Ra
        227.0,        // 89 Ac
        232.038,      // 90 Th
        231.036,      // 91 Pa
        238.029,      // 92 U
        237.0,        // 93 Np
        244.0,        // 94 Pu
        243.0,        // 95 Am
        247.0,        // 96 Cm
        247.0,        // 97 Bk
        251.0,        // 98 Cf
        252.0,        // 99 Es
        257.0,        // 100 Fm
        258.0,        // 101 Md
        259.0,        // 102 No
        266.0,        // 103 Lr
        267.0,        // 104 Rf
        268.0,        // 105 Db
        269.0,        // 106 Sg
        270.0,        // 107 Bh
        277.0,        // 108 Hs
        278.0,        // 109 Mt
        281.0,        // 110 Ds
        282.0,        // 111 Rg
        285.0,        // 112 Cn
        286.0,        // 113 Nh
        289.0,        // 114 Fl
        290.0,        // 115 Mc
        293.0,        // 116 Lv
        294.0,        // 117 Ts
        294.0         // 118 Og
    };

    if (z >= 1 && z <= 118)
        return kMasses[z];
    return 1.0;
}
