#include "StructureLoader.h"

#include <openbabel3/openbabel/obconversion.h>
#include <openbabel3/openbabel/mol.h>
#include <openbabel3/openbabel/atom.h>
#include <openbabel3/openbabel/elements.h>
#include <openbabel3/openbabel/generic.h>
#include <openbabel3/openbabel/math/vector3.h>
#include <openbabel3/openbabel/dlhandler.h>
#include <openbabel3/openbabel/plugin.h>

#include <iostream>
#include <array>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <sys/stat.h>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace
{
std::string normalizeSeparators(const std::string& path)
{
    std::string out = path;
    std::replace(out.begin(), out.end(), '\\', '/');
    return out;
}

std::string ensureTrailingSlash(const std::string& path)
{
    if (path.empty())
        return path;
    if (path.back() == '/' || path.back() == '\\')
        return path;
    return path + "/";
}

std::string parentDirectory(const std::string& path)
{
    if (path.empty())
        return std::string();

    std::string out = normalizeSeparators(path);
    while (out.size() > 1 && out.back() == '/')
        out.pop_back();

    std::size_t pos = out.find_last_of('/');
    if (pos == std::string::npos)
        return std::string();
    if (pos == 0)
        return out.substr(0, 1);
    return out.substr(0, pos);
}

std::string joinPath(const std::string& base, const std::string& name)
{
    if (base.empty())
        return name;
    if (base.back() == '/' || base.back() == '\\')
        return base + name;
    return base + "/" + name;
}

bool directoryExists(const std::string& path)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0)
        return false;
#ifdef _WIN32
    return (st.st_mode & _S_IFDIR) != 0;
#else
    return S_ISDIR(st.st_mode);
#endif
}

void setBabelLibDir(const std::string& path)
{
    const std::string normalized = ensureTrailingSlash(normalizeSeparators(path));
#ifdef _WIN32
    _putenv_s("BABEL_LIBDIR", normalized.c_str());
#else
    setenv("BABEL_LIBDIR", normalized.c_str(), 1);
#endif
}

void ensureOpenBabelPlugins()
{
    static bool initialized = false;
    if (initialized)
        return;
    initialized = true;

    const char* currentLibDir = std::getenv("BABEL_LIBDIR");
    if (currentLibDir && *currentLibDir)
    {
        OpenBabel::OBPlugin::LoadAllPlugins();
        return;
    }

    std::vector<std::string> candidates;

    std::string convPath;
    if (DLHandler::getConvDirectory(convPath) && !convPath.empty())
        candidates.push_back(convPath);

#ifdef _WIN32
    HMODULE obModule = nullptr;
    if (GetModuleHandleExA(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCSTR>(&OpenBabel::OBPlugin::LoadAllPlugins),
            &obModule))
    {
        char modulePath[MAX_PATH] = {0};
        DWORD len = GetModuleFileNameA(obModule, modulePath, MAX_PATH);
        if (len > 0 && len < MAX_PATH)
        {
            const std::string binDir = parentDirectory(modulePath);
            const std::string prefixDir = parentDirectory(binDir);
            if (!prefixDir.empty())
                candidates.push_back(joinPath(joinPath(prefixDir, "lib"), "openbabel/3.1.0"));
        }
    }

    candidates.push_back("C:/msys64/ucrt64/lib/openbabel/3.1.0");
    candidates.push_back("C:/msys64/mingw64/lib/openbabel/3.1.0");
#endif

    for (std::size_t i = 0; i < candidates.size(); ++i)
    {
        if (directoryExists(candidates[i]))
        {
            setBabelLibDir(candidates[i]);
            break;
        }
    }

    OpenBabel::OBPlugin::LoadAllPlugins();
}

std::string toLowerCopy(const std::string& value)
{
    std::string out = value;
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = (char)std::tolower((unsigned char)out[i]);
    return out;
}

std::string extractExtension(const std::string& filename)
{
    const std::size_t slashPos = filename.find_last_of("/\\");
    const std::size_t dotPos = filename.find_last_of('.');
    if (dotPos == std::string::npos)
        return std::string();
    if (slashPos != std::string::npos && dotPos < slashPos)
        return std::string();
    return filename.substr(dotPos);
}

bool isSupportedExtension(const std::string& extLower)
{
    static const char* kSupportedExtensions[] = {
        ".cif",
        ".mol",
        ".pdb",
        ".xyz",
        ".sdf",
        ".vasp",
        ".mol2",
        ".pwi",
        ".gjf"
    };

    for (std::size_t i = 0; i < sizeof(kSupportedExtensions) / sizeof(kSupportedExtensions[0]); ++i)
    {
        if (extLower == kSupportedExtensions[i])
            return true;
    }
    return false;
}

std::string supportedExtensionsSummary()
{
    return ".cif, .mol, .pdb, .xyz, .sdf, .vasp, .mol2, .pwi, .gjf";
}
}


void getDefaultElementColor(int Z,float& r,float& g,float& b)
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
        std::array<float,3>{0.39f, 0.13f, 0.67f}, // 104 Rf
        std::array<float,3>{0.46f, 0.13f, 0.64f}, // 105 Db
        std::array<float,3>{0.54f, 0.16f, 0.59f}, // 106 Sg
        std::array<float,3>{0.61f, 0.18f, 0.55f}, // 107 Bh
        std::array<float,3>{0.68f, 0.20f, 0.50f}, // 108 Hs
        std::array<float,3>{0.74f, 0.20f, 0.46f}, // 109 Mt
        std::array<float,3>{0.78f, 0.21f, 0.42f}, // 110 Ds
        std::array<float,3>{0.82f, 0.22f, 0.38f}, // 111 Rg
        std::array<float,3>{0.86f, 0.23f, 0.34f}, // 112 Cn
        std::array<float,3>{0.90f, 0.24f, 0.30f}, // 113 Nh
        std::array<float,3>{0.94f, 0.26f, 0.26f}, // 114 Fl
        std::array<float,3>{0.98f, 0.28f, 0.22f}, // 115 Mc
        std::array<float,3>{1.00f, 0.30f, 0.18f}, // 116 Lv
        std::array<float,3>{1.00f, 0.32f, 0.14f}, // 117 Ts
        std::array<float,3>{1.00f, 0.34f, 0.10f}, // 118 Og
    };

    if (Z < 1) Z = 1;
    if (Z >= (int)colors.size()) Z = (int)colors.size() - 1;

    r = colors[Z][0];
    g = colors[Z][1];
    b = colors[Z][2];
}

bool isSupportedStructureFile(const std::string& filename)
{
    const std::string extLower = toLowerCopy(extractExtension(filename));
    return isSupportedExtension(extLower);
}

bool loadStructureFromFile(const std::string& filename, Structure& structure, std::string& errorMessage)
{
    structure = Structure();
    errorMessage.clear();

    if (filename.empty())
    {
        errorMessage = "No file selected.";
        return false;
    }

    const std::string ext = extractExtension(filename);
    const std::string extLower = toLowerCopy(ext);
    if (!isSupportedExtension(extLower))
    {
        if (extLower.empty())
            errorMessage = "Unsupported file format (missing extension). Supported formats: " + supportedExtensionsSummary();
        else
            errorMessage = "Unsupported file format '" + ext + "'. Supported formats: " + supportedExtensionsSummary();
        return false;
    }

    struct stat fileStat;
    if (stat(filename.c_str(), &fileStat) != 0)
    {
        errorMessage = "File not found: " + filename;
        return false;
    }

    ensureOpenBabelPlugins();

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;

    OpenBabel::OBFormat* inFmt = conv.FormatFromExt(filename.c_str());
    if (!inFmt || !conv.SetInFormat(inFmt))
    {
        errorMessage = "Unsupported file format '" + ext + "'. Supported formats: " + supportedExtensionsSummary();
        return false;
    }

    if (!conv.ReadFile(&mol, filename))
    {
        errorMessage = "Failed to load file. The file may be corrupted or unreadable.";
        return false;
    }

    // If the file provides unit cell / periodic information, make sure we
    // generate the full unit cell (symmetry-equivalent atoms) so that we
    // display the complete crystal structure.
    if (auto *data = mol.GetData(OpenBabel::OBGenericDataType::UnitCell))
    {
        if (auto *cell = dynamic_cast<OpenBabel::OBUnitCell*>(data))
        {
            cell->FillUnitCell(&mol);
            structure.hasUnitCell = true;

            auto vecs = cell->GetCellVectors();
            if (vecs.size() >= 3)
            {
                for (int i = 0; i < 3; ++i)
                {
                    structure.cellVectors[i][0] = vecs[i].GetX();
                    structure.cellVectors[i][1] = vecs[i].GetY();
                    structure.cellVectors[i][2] = vecs[i].GetZ();
                }
            }

            auto off = cell->GetOffset();
            structure.cellOffset[0] = off.GetX();
            structure.cellOffset[1] = off.GetY();
            structure.cellOffset[2] = off.GetZ();
        }
    }

    OpenBabel::OBAtomIterator ai;
    constexpr size_t kWarningAtomCount = 10000;
    constexpr size_t kMaxAtomCount = 500000;

    for(OpenBabel::OBAtom* atom = mol.BeginAtom(ai);
        atom;
        atom = mol.NextAtom(ai))
    {
        if (structure.atoms.size() >= kMaxAtomCount)
        {
            std::cerr << "Warning: Atom count capped at " << kMaxAtomCount 
                      << " to prevent memory exhaustion." << std::endl;
            break;
        }

        AtomSite site;

        site.atomicNumber = atom->GetAtomicNum();
        site.symbol = OpenBabel::OBElements::GetSymbol(site.atomicNumber);

        site.x = atom->GetX();
        site.y = atom->GetY();
        site.z = atom->GetZ();

        getDefaultElementColor(site.atomicNumber, site.r, site.g, site.b);

        structure.atoms.push_back(site);

        if (structure.atoms.size() == kWarningAtomCount)
        {
            std::cerr << "Loading large structure with " << kWarningAtomCount 
                      << " atoms (bond computation: O(n) with spatial hashing)..." << std::endl;
        }
    }

    if (structure.atoms.size() > kWarningAtomCount)
    {
        std::cerr << "Loaded structure with " << structure.atoms.size() 
                  << " atoms." << std::endl;
    }

    if (structure.atoms.empty())
    {
        errorMessage = "File loaded but no atoms were found.";
        return false;
    }

    return true;
}

Structure loadStructure(const std::string& filename)
{
    Structure structure;
    std::string errorMessage;
    if (!filename.empty() && !loadStructureFromFile(filename, structure, errorMessage))
        std::cerr << "Warning: " << errorMessage << std::endl;

    return structure;
}

bool saveStructure(const Structure& structure, const std::string& filename, const std::string& format)
{
    ensureOpenBabelPlugins();

    OpenBabel::OBMol mol;
    mol.BeginModify();

    for (const auto& site : structure.atoms)
    {
        OpenBabel::OBAtom* atom = mol.NewAtom();
        atom->SetAtomicNum(site.atomicNumber);
        atom->SetVector(site.x, site.y, site.z);
    }

    if (structure.hasUnitCell)
    {
        OpenBabel::OBUnitCell* cell = new OpenBabel::OBUnitCell();
        OpenBabel::vector3 va(structure.cellVectors[0][0], structure.cellVectors[0][1], structure.cellVectors[0][2]);
        OpenBabel::vector3 vb(structure.cellVectors[1][0], structure.cellVectors[1][1], structure.cellVectors[1][2]);
        OpenBabel::vector3 vc(structure.cellVectors[2][0], structure.cellVectors[2][1], structure.cellVectors[2][2]);
        cell->SetData(va, vb, vc);
        mol.SetData(cell);
    }

    mol.ConnectTheDots();
    mol.PerceiveBondOrders();
    mol.EndModify();

    OpenBabel::OBConversion conv;
    OpenBabel::OBFormat* outFmt = conv.FindFormat(format.c_str());
    if (!outFmt)
        return false;

    conv.SetOutFormat(outFmt);
    return conv.WriteFile(&mol, filename);
}
