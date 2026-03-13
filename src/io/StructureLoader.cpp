#include "StructureLoader.h"

#include <openbabel3/openbabel/obconversion.h>
#include <openbabel3/openbabel/mol.h>
#include <openbabel3/openbabel/atom.h>
#include <openbabel3/openbabel/elements.h>

#include <iostream>
#include <array>


static void elementColor(int Z,float& r,float& g,float& b)
{
    // CPK coloring scheme for all elements (1..118).
    // Values taken from common molecular visualization defaults.
    static const std::array<std::array<float,3>, 119> colors = {
        std::array<float,3>{0.0f, 0.0f, 0.0f}, // placeholder for index 0
        std::array<float,3>{1.00f, 1.00f, 1.00f}, // 1  H
        std::array<float,3>{0.85f, 1.00f, 1.00f}, // 2  He
        std::array<float,3>{0.80f, 0.50f, 1.00f}, // 3  Li
        std::array<float,3>{0.76f, 1.00f, 0.00f}, // 4  Be
        std::array<float,3>{1.00f, 0.71f, 0.71f}, // 5  B
        std::array<float,3>{0.20f, 0.20f, 0.20f}, // 6  C
        std::array<float,3>{0.00f, 0.00f, 1.00f}, // 7  N
        std::array<float,3>{1.00f, 0.00f, 0.00f}, // 8  O
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 9  F
        std::array<float,3>{0.70f, 0.89f, 0.96f}, // 10 Ne
        std::array<float,3>{0.67f, 0.36f, 0.95f}, // 11 Na
        std::array<float,3>{0.54f, 1.00f, 0.00f}, // 12 Mg
        std::array<float,3>{0.75f, 0.65f, 0.65f}, // 13 Al
        std::array<float,3>{0.94f, 0.78f, 0.62f}, // 14 Si
        std::array<float,3>{1.00f, 0.50f, 0.00f}, // 15 P
        std::array<float,3>{1.00f, 1.00f, 0.00f}, // 16 S
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 17 Cl
        std::array<float,3>{0.50f, 0.82f, 0.89f}, // 18 Ar
        std::array<float,3>{0.56f, 0.00f, 1.00f}, // 19 K
        std::array<float,3>{0.24f, 1.00f, 0.00f}, // 20 Ca
        std::array<float,3>{0.90f, 0.90f, 0.90f}, // 21 Sc
        std::array<float,3>{0.75f, 0.76f, 0.78f}, // 22 Ti
        std::array<float,3>{0.65f, 0.65f, 0.67f}, // 23 V
        std::array<float,3>{0.54f, 0.60f, 0.78f}, // 24 Cr
        std::array<float,3>{0.61f, 0.48f, 0.78f}, // 25 Mn
        std::array<float,3>{0.88f, 0.40f, 0.20f}, // 26 Fe
        std::array<float,3>{0.88f, 0.38f, 0.20f}, // 27 Co
        std::array<float,3>{0.78f, 0.79f, 0.78f}, // 28 Ni
        std::array<float,3>{0.78f, 0.50f, 0.20f}, // 29 Cu
        std::array<float,3>{0.49f, 0.50f, 0.69f}, // 30 Zn
        std::array<float,3>{0.76f, 0.56f, 0.56f}, // 31 Ga
        std::array<float,3>{0.40f, 0.56f, 0.56f}, // 32 Ge
        std::array<float,3>{0.74f, 0.50f, 0.89f}, // 33 As
        std::array<float,3>{1.00f, 0.63f, 0.00f}, // 34 Se
        std::array<float,3>{0.65f, 0.16f, 0.16f}, // 35 Br
        std::array<float,3>{0.36f, 0.72f, 0.82f}, // 36 Kr
        std::array<float,3>{0.44f, 0.18f, 0.69f}, // 37 Rb
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 38 Sr
        std::array<float,3>{0.58f, 1.00f, 1.00f}, // 39 Y
        std::array<float,3>{0.58f, 0.88f, 0.88f}, // 40 Zr
        std::array<float,3>{0.45f, 0.76f, 0.79f}, // 41 Nb
        std::array<float,3>{0.33f, 0.71f, 0.71f}, // 42 Mo
        std::array<float,3>{0.23f, 0.62f, 0.62f}, // 43 Tc
        std::array<float,3>{0.14f, 0.56f, 0.56f}, // 44 Ru
        std::array<float,3>{0.04f, 0.49f, 0.55f}, // 45 Rh
        std::array<float,3>{0.00f, 0.41f, 0.52f}, // 46 Pd
        std::array<float,3>{0.75f, 0.75f, 0.75f}, // 47 Ag
        std::array<float,3>{1.00f, 0.85f, 0.00f}, // 48 Cd
        std::array<float,3>{0.65f, 0.46f, 0.45f}, // 49 In
        std::array<float,3>{0.40f, 0.50f, 0.50f}, // 50 Sn
        std::array<float,3>{0.62f, 0.39f, 0.71f}, // 51 Sb
        std::array<float,3>{0.83f, 0.50f, 0.18f}, // 52 Te
        std::array<float,3>{0.58f, 0.00f, 0.58f}, // 53 I
        std::array<float,3>{0.26f, 0.62f, 0.69f}, // 54 Xe
        std::array<float,3>{0.34f, 0.09f, 0.56f}, // 55 Cs
        std::array<float,3>{0.00f, 0.79f, 0.00f}, // 56 Ba
        std::array<float,3>{0.44f, 0.83f, 1.00f}, // 57 La
        std::array<float,3>{1.00f, 1.00f, 0.78f}, // 58 Ce
        std::array<float,3>{0.85f, 1.00f, 0.78f}, // 59 Pr
        std::array<float,3>{0.78f, 1.00f, 0.78f}, // 60 Nd
        std::array<float,3>{0.64f, 1.00f, 0.78f}, // 61 Pm
        std::array<float,3>{0.56f, 1.00f, 0.78f}, // 62 Sm
        std::array<float,3>{0.38f, 1.00f, 0.78f}, // 63 Eu
        std::array<float,3>{0.27f, 1.00f, 0.78f}, // 64 Gd
        std::array<float,3>{0.19f, 1.00f, 0.78f}, // 65 Tb
        std::array<float,3>{0.12f, 1.00f, 0.78f}, // 66 Dy
        std::array<float,3>{0.10f, 1.00f, 0.61f}, // 67 Ho
        std::array<float,3>{0.10f, 0.86f, 0.56f}, // 68 Er
        std::array<float,3>{0.10f, 0.77f, 0.52f}, // 69 Tm
        std::array<float,3>{0.12f, 0.70f, 0.47f}, // 70 Yb
        std::array<float,3>{0.12f, 0.63f, 0.42f}, // 71 Lu
        std::array<float,3>{0.30f, 0.61f, 0.61f}, // 72 Hf
        std::array<float,3>{0.30f, 0.58f, 0.58f}, // 73 Ta
        std::array<float,3>{0.13f, 0.58f, 0.58f}, // 74 W
        std::array<float,3>{0.14f, 0.54f, 0.54f}, // 75 Re
        std::array<float,3>{0.14f, 0.51f, 0.51f}, // 76 Os
        std::array<float,3>{0.09f, 0.46f, 0.52f}, // 77 Ir
        std::array<float,3>{0.82f, 0.82f, 0.88f}, // 78 Pt
        std::array<float,3>{1.00f, 0.82f, 0.14f}, // 79 Au
        std::array<float,3>{0.72f, 0.72f, 0.82f}, // 80 Hg
        std::array<float,3>{0.65f, 0.33f, 0.30f}, // 81 Tl
        std::array<float,3>{0.34f, 0.35f, 0.38f}, // 82 Pb
        std::array<float,3>{0.62f, 0.31f, 0.71f}, // 83 Bi
        std::array<float,3>{0.67f, 0.36f, 0.95f}, // 84 Po
        std::array<float,3>{0.46f, 0.32f, 0.28f}, // 85 At
        std::array<float,3>{0.26f, 0.42f, 0.47f}, // 86 Rn
        std::array<float,3>{0.26f, 0.00f, 0.40f}, // 87 Fr
        std::array<float,3>{0.00f, 0.42f, 0.00f}, // 88 Ra
        std::array<float,3>{0.00f, 0.70f, 1.00f}, // 89 Ac
        std::array<float,3>{0.00f, 0.73f, 1.00f}, // 90 Th
        std::array<float,3>{0.00f, 0.73f, 1.00f}, // 91 Pa
        std::array<float,3>{0.00f, 0.50f, 1.00f}, // 92 U
        std::array<float,3>{0.00f, 0.47f, 1.00f}, // 93 Np
        std::array<float,3>{0.00f, 0.40f, 1.00f}, // 94 Pu
        std::array<float,3>{0.00f, 0.38f, 1.00f}, // 95 Am
        std::array<float,3>{0.00f, 0.35f, 1.00f}, // 96 Cm
        std::array<float,3>{0.00f, 0.32f, 1.00f}, // 97 Bk
        std::array<float,3>{0.00f, 0.30f, 1.00f}, // 98 Cf
        std::array<float,3>{0.18f, 0.22f, 0.88f}, // 99 Es
        std::array<float,3>{0.20f, 0.20f, 0.90f}, // 100 Fm
        std::array<float,3>{0.20f, 0.20f, 0.88f}, // 101 Md
        std::array<float,3>{0.20f, 0.20f, 0.87f}, // 102 No
        std::array<float,3>{0.20f, 0.20f, 0.86f}, // 103 Lr
        std::array<float,3>{0.18f, 0.20f, 0.83f}, // 104 Rf
        std::array<float,3>{0.16f, 0.19f, 0.80f}, // 105 Db
        std::array<float,3>{0.14f, 0.18f, 0.78f}, // 106 Sg
        std::array<float,3>{0.13f, 0.17f, 0.75f}, // 107 Bh
        std::array<float,3>{0.12f, 0.16f, 0.72f}, // 108 Hs
        std::array<float,3>{0.10f, 0.15f, 0.70f}, // 109 Mt
        std::array<float,3>{0.09f, 0.14f, 0.68f}, // 110 Ds
        std::array<float,3>{0.09f, 0.13f, 0.65f}, // 111 Rg
        std::array<float,3>{0.09f, 0.13f, 0.62f}, // 112 Cn
        std::array<float,3>{0.08f, 0.12f, 0.60f}, // 113 Nh
        std::array<float,3>{0.08f, 0.12f, 0.58f}, // 114 Fl
        std::array<float,3>{0.07f, 0.10f, 0.55f}, // 115 Mc
        std::array<float,3>{0.07f, 0.10f, 0.53f}, // 116 Lv
        std::array<float,3>{0.06f, 0.09f, 0.50f}, // 117 Ts
        std::array<float,3>{0.05f, 0.08f, 0.47f}  // 118 Og
    };

    if (Z >= 1 && Z < (int)colors.size()) {
        r = colors[Z][0];
        g = colors[Z][1];
        b = colors[Z][2];
    } else {
        r = g = b = 0.7f;
    }
}

Structure loadStructure(const std::string& filename)
{
    Structure structure;

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;

    conv.SetInFormat(conv.FormatFromExt(filename.c_str()));
    conv.ReadFile(&mol, filename);

    OpenBabel::OBAtomIterator ai;

    for(OpenBabel::OBAtom* atom = mol.BeginAtom(ai);
        atom;
        atom = mol.NextAtom(ai))
    {
        AtomSite site;

        site.atomicNumber = atom->GetAtomicNum();
        site.symbol = OpenBabel::OBElements::GetSymbol(site.atomicNumber);

        site.x = atom->GetX();
        site.y = atom->GetY();
        site.z = atom->GetZ();

        elementColor(site.atomicNumber, site.r, site.g, site.b);

        structure.atoms.push_back(site);
    }

    return structure;
}
