#include "cli/CLIMode.h"

#include "algorithms/AmorphousBuilder.h"
#include "algorithms/BulkCrystalBuilder.h"
#include "algorithms/CSLComputation.h"
#include "algorithms/DislocationBuilder.h"
#include "algorithms/InterfaceBuilder.h"
#include "algorithms/MeshLoader.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "algorithms/PolyCrystalBuilder.h"
#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "io/StructureLoader.h"
#include "util/ElementData.h"
#include "util/PathUtils.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

// ── Argument helpers ─────────────────────────────────────────────────────────

static const char* findArg(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return argv[i + 1];
    return nullptr;
}

static bool hasFlag(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return true;
    return false;
}

// Collect all values after repeated occurrences of `flag`.
static std::vector<std::string> findAllArgs(int argc, char* argv[], const char* flag)
{
    std::vector<std::string> out;
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            out.emplace_back(argv[i + 1]);
    return out;
}

static double argDouble(int argc, char* argv[], const char* flag, double def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v)
        return def;

    try
    {
        std::size_t parsed = 0;
        const std::string text(v);
        const double value = std::stod(text, &parsed);
        if (parsed != text.size() || !std::isfinite(value))
            throw std::invalid_argument("not a finite number");
        return value;
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(std::string("invalid numeric value '") + v + "' for " + flag);
    }
}

static int argInt(int argc, char* argv[], const char* flag, int def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v)
        return def;

    try
    {
        std::size_t parsed = 0;
        const std::string text(v);
        const int value = std::stoi(text, &parsed);
        if (parsed != text.size())
            throw std::invalid_argument("trailing characters");
        return value;
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(std::string("invalid integer value '") + v + "' for " + flag);
    }
}

// Detect Open Babel format string from output filename extension.
static std::string detectFormat(const std::string& path)
{
    const auto dot = path.rfind('.');
    if (dot == std::string::npos)
        return "cif";
    std::string ext = path.substr(dot + 1);
    // lower-case
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

    if (ext == "cif")              return "cif";
    if (ext == "xyz")              return "xyz";
    if (ext == "vasp" || ext == "poscar" || ext == "contcar") return "vasp";
    if (ext == "lmp"  || ext == "lammps") return "lammps";
    if (ext == "pdb")              return "pdb";
    if (ext == "mol2")             return "mol2";
    if (ext == "extxyz")           return "extxyz";
    return ext; // pass-through for anything else
}

// ── Print usage ──────────────────────────────────────────────────────────────

static void printHelpBulk()
{
    std::cout <<
"------------------------------------------------------------------\n"
"BULK CRYSTAL  (--build bulk)\n"
"------------------------------------------------------------------\n"
"  --system     <name>         Crystal system: triclinic | monoclinic |\n"
"                               orthorhombic | tetragonal | trigonal |\n"
"                               hexagonal | cubic   (default: cubic)\n"
"  --spacegroup <N>            Space-group number 1-230  (default: 225)\n"
"  --a   <Ang>                 Lattice parameter a       (default: 4.0)\n"
"  --b   <Ang>                 Lattice parameter b       (default: a)\n"
"  --c   <Ang>                 Lattice parameter c       (default: a)\n"
"  --alpha <deg>               Cell angle alpha          (default: 90)\n"
"  --beta  <deg>               Cell angle beta           (default: 90)\n"
"  --gamma <deg>               Cell angle gamma          (default: 90)\n"
"  --atom  \"SYMBOL fx fy fz\"   Asymmetric unit atom (fractional coords).\n"
"                               May be repeated for multiple atoms.\n"
"                               Example: --atom \"Cu 0 0 0\"\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build bulk --system cubic --spacegroup 225 ^\n"
"            --a 3.61 --atom \"Cu 0 0 0\" --output cu_fcc.cif\n"
<< std::endl;
}

static void printHelpGB()
{
    std::cout <<
"------------------------------------------------------------------\n"
"CSL GRAIN BOUNDARY  (--build gb)\n"
"------------------------------------------------------------------\n"
"  --input  <file>             Reference structure (must have a unit cell)\n"
"  --axis   \"u v w\"            Rotation axis, e.g. \"0 0 1\"   (default: 0 0 1)\n"
"  --sigma  <N>                Exact sigma value to use\n"
"  --sigmamax <N>              Maximum sigma to search when --sigma is omitted\n"
"                               (picks the smallest available sigma, default: 100)\n"
"  --plane  <0|1|2>            GB-plane index within the sigma candidate (default: 0)\n"
"  --uca    <N>                Grain-A unit-cell repetitions along stack dir (default: 1)\n"
"  --ucb    <N>                Grain-B unit-cell repetitions along stack dir (default: 1)\n"
"  --vacuum <Ang>              Vacuum padding on each side  (default: 0)\n"
"  --gap    <Ang>              Interface gap between grains (default: 0)\n"
"  --overlap <Ang>             Overlap removal radius       (default: 0)\n"
"  --conventional              Keep conventional cell (skip primitive reduction)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build gb --input cu_fcc.cif --axis \"0 0 1\" ^\n"
"            --sigma 5 --plane 0 --uca 3 --ucb 3 ^\n"
"            --vacuum 5.0 --output cu_sigma5_gb.cif\n"
<< std::endl;
}

static void printHelpPoly()
{
    std::cout <<
"------------------------------------------------------------------\n"
"POLYCRYSTAL  (--build poly)\n"
"------------------------------------------------------------------\n"
"  --input  <file>             Reference structure (must have a unit cell)\n"
"  --sizex  <Ang>              Box dimension X  (default: 50)\n"
"  --sizey  <Ang>              Box dimension Y  (default: 50)\n"
"  --sizez  <Ang>              Box dimension Z  (default: 50)\n"
"  --grains <N>                Number of Voronoi grains  (default: 8)\n"
"  --seed   <N>                Random seed for reproducibility  (default: 42)\n"
"  --euler  \"phi1 Phi phi2\"    Bunge Euler angles (deg) for grain 0, 1, ...\n"
"                               Specify once per grain in order. If fewer\n"
"                               than --grains values are given, remaining\n"
"                               grains get random orientations.\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build poly --input cu_fcc.cif ^\n"
"            --sizex 100 --sizey 100 --sizez 100 ^\n"
"            --grains 12 --seed 7 --output cu_poly.cif\n"
<< std::endl;
}

static void printHelpNano()
{
    std::cout <<
"------------------------------------------------------------------\n"
"NANOCRYSTAL  (--build nano)\n"
"------------------------------------------------------------------\n"
"  --input  <file>             Reference crystal (must have a unit cell)\n"
"  --shape  <name>             sphere | ellipsoid | box | cylinder |\n"
"                               octahedron | truncated-octahedron |\n"
"                               cuboctahedron   (default: sphere)\n"
"  --radius <Ang>              Sphere/octahedron/cuboctahedron radius (default: 15)\n"
"  --rx <Ang>                  Ellipsoid half-extent X  (default: 15)\n"
"  --ry <Ang>                  Ellipsoid half-extent Y  (default: 12)\n"
"  --rz <Ang>                  Ellipsoid half-extent Z  (default: 10)\n"
"  --hx <Ang>                  Box / truncated-oct half-extent X  (default: 15)\n"
"  --hy <Ang>                  Box / truncated-oct half-extent Y  (default: 15)\n"
"  --hz <Ang>                  Box / truncated-oct half-extent Z  (default: 15)\n"
"  --trunc <Ang>               Truncated-octahedron truncation radius  (default: 12)\n"
"  --cylradius <Ang>           Cylinder radius  (default: 12)\n"
"  --cylheight <Ang>           Cylinder height  (default: 30)\n"
"  --cylaxis   <0|1|2>         Cylinder axis: 0=X 1=Y 2=Z  (default: 2)\n"
"  --vacuum <Ang>              Vacuum padding around particle  (default: 5)\n"
"  --repa <N>                  Manual supercell replication A  (0 = auto)\n"
"  --repb <N>                  Manual supercell replication B  (0 = auto)\n"
"  --repc <N>                  Manual supercell replication C  (0 = auto)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build nano --input cu.cif --shape sphere --radius 20 ^\n"
"            --vacuum 5 --output cu_nano.xyz\n"
<< std::endl;
}

static void printHelpAmorphous()
{
    std::cout <<
"------------------------------------------------------------------\n"
"AMORPHOUS STRUCTURE  (--build amorphous)\n"
"------------------------------------------------------------------\n"
"  --element \"SYMBOL N\"        Element species and atom count.\n"
"                               Repeat for each species in the mixture.\n"
"                               Example: --element \"Si 80\" --element \"O 160\"\n"
"  --density <g/cm3>           Target density for auto box size  (default: 2.0)\n"
"  --boxa <Ang>                Manual box length A (overrides --density)\n"
"  --boxb <Ang>                Manual box length B (overrides --density)\n"
"  --boxc <Ang>                Manual box length C (overrides --density)\n"
"  --scale <factor>            Cell scale factor before packing  (default: 1.0)\n"
"  --mindist \"Z1 Z2 dist\"      Minimum separation for element pair (Ang).\n"
"                               Repeat for multiple pairs.\n"
"  --covtol <frac>             Covalent-radii tolerance fraction  (default: 0.75)\n"
"  --seed   <N>                RNG seed (0 = time-based)  (default: 42)\n"
"  --attempts <N>              Max placement attempts per atom  (default: 1000)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build amorphous --element \"Si 80\" --element \"O 160\" ^\n"
"            --density 2.2 --seed 1 --output sio2.xyz\n"
<< std::endl;
}

[[maybe_unused]] static void printHelpInterface()
{
    std::cout <<
"------------------------------------------------------------------\n"
"HETEROGENEOUS INTERFACE  (--build interface)\n"
"------------------------------------------------------------------\n"
"  --layerA <file>             Structure for layer A (must have a unit cell)\n"
"  --layerB <file>             Structure for layer B (must have a unit cell)\n"
"  --nmax   <N>                Max supercell search index  (default: 4)\n"
"  --maxcells <N>              Max supercell area (unit cells)  (default: 16)\n"
"  --pick   <N>                Index of strain-matched supercell pair to use\n"
"                               (0 = best match, default: 0)\n"
"  --layersA <N>               Z repetitions of layer A  (default: 1)\n"
"  --layersB <N>               Z repetitions of layer B  (default: 1)\n"
"  --gap    <Ang>              Interface gap  (default: 2.0)\n"
"  --vacuum <Ang>              Vacuum above/below  (default: 10.0)\n"
"  --repx   <N>                XY repeat X  (default: 1)\n"
"  --repy   <N>                XY repeat Y  (default: 1)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build interface --layerA cu.cif --layerB ni.cif ^\n"
"            --nmax 4 --layersA 3 --layersB 3 ^\n"
"            --gap 2.0 --vacuum 10.0 --output cu_ni_interface.cif\n"
<< std::endl;
}

static void printHelpCustom()
{
    std::cout <<
"------------------------------------------------------------------\n"
"CUSTOM MESH-FILL  (--build custom)\n"
"------------------------------------------------------------------\n"
"Fill a 3D mesh model with atoms from a reference crystal structure.\n"
"\n"
"  --input   <file>            Reference crystal file (CIF, XYZ, VASP, PDB, ...)\n"
"  --mesh    <file>            3D model file (OBJ or STL)\n"
"  --scale   <factor>          Angstrom per model unit  (default: 1.0)\n"
"  --vacuum  <Ang>             Vacuum padding for output cell  (default: 5.0)\n"
"  --repa <N>                  Manual replication along a  (0 = auto)\n"
"  --repb <N>                  Manual replication along b  (0 = auto)\n"
"  --repc <N>                  Manual replication along c  (0 = auto)\n"
"  --rotx <deg>                Rotate crystal about X before fill  (default: 0)\n"
"  --roty <deg>                Rotate crystal about Y before fill  (default: 0)\n"
"  --rotz <deg>                Rotate crystal about Z before fill  (default: 0)\n"
"  --miller  \"h k l\"          Align crystal direction [hkl] with mesh +Z axis\n"
"                               (overrides --rotx/y/z when given)\n"
"  --output  <file>            Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build custom --input cu.cif --mesh bunny.stl ^\n"
"            --scale 0.1 --vacuum 5 --output cu_bunny.xyz\n"
<< std::endl;
}

static void printHelpSSS()
{
    std::cout <<
"------------------------------------------------------------------\n"
"SUBSTITUTIONAL SOLID SOLUTION  (--build sss)\n"
"------------------------------------------------------------------\n"
"  --input  <file>             Host structure to substitute into\n"
"  --frac   SYM=frac,...       Comma-separated element=fraction pairs.\n"
"                               Fractions must sum to 1.\n"
"                               Example: --frac \"Cu=0.7,Zn=0.3\"\n"
"  --seed   <N>                RNG seed (0 = time-based)  (default: 12345)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build sss --input cu.cif ^\'\n"
"            --frac \"Cu=0.7,Zn=0.3\" --seed 42 --output brass.cif\n"
<< std::endl;
}

static void printHelpDislocation()
{
    std::cout <<
"------------------------------------------------------------------\n"
"DISLOCATION  (--build dislocation)\n"
"------------------------------------------------------------------\n"
"  --input  <file>             Input structure (must have a unit cell)\n"
"  --output <file>             Output file (format from extension)\n"
"  --character <edge|screw|mixed>\n"
"                              Dislocation character (default: edge)\n"
"  --shape <halfplane|cylinder|sphere|ellipsoid|freeform>\n"
"                              Application region shape (default: halfplane)\n"
"  --manual-vectors            Disable automatic FCC/HCP/BCC vector presets\n"
"  --plane   \"h k l\"          Slip-plane Miller indices (default: 1 1 1)\n"
"  --burgers \"u v w\"          Burgers direction indices (default: 1 -1 0)\n"
"  --line    \"u v w\"          Line direction indices (default: 1 1 -2)\n"
"  --line-frac   \"fx fy fz\"   Dislocation line point in fractional coords\n"
"  --line-offset \"x y z\"      Cartesian offset added to line point (A)\n"
"  --bscale <f>                Burgers magnitude scale factor (default: 1.0)\n"
"  --bmag <A>                  Explicit Burgers magnitude in Angstroms (overrides auto-detect)\n"
"  --mixed-angle <deg>         Edge/screw mixing angle for mixed mode\n"
"  --nu <f>                    Poisson ratio (default: 0.33)\n"
"  --core <A>                  Core radius regularization (default: 1.2)\n"
"  --cutoff <A>                Radial attenuation distance, 0=off (default: 20)\n"
"  --line-half <A>             Half-length extent along line (default: 1e6)\n"
"  --cyl-radius <A>            Cylinder radius for shape=cylinder\n"
"  --sphere-radius <A>         Sphere radius for shape=sphere\n"
"  --ellipsoid \"rx ry rz\"     Ellipsoid radii for shape=ellipsoid\n"
"  --poly2d \"x1 y1;x2 y2;...\"  Freeform polygon for shape=freeform\n"
"\n"
"Example:\n"
"  AtomForge --build dislocation --input base.cfg --character edge ^\n"
"            --manual-vectors --line \"0 0 1\" --burgers \"1 0 0\" ^\n"
"            --shape cylinder --cyl-radius 12 --output edge.cfg\n"
<< std::endl;
}

static void printHelp()
{
    std::cout <<
"AtomForge CLI - headless structure builder\n"
"\n"
"Usage:\n"
"  AtomForge --build <mode> [options] --output <file>\n"
"\n"
"Modes:\n"
"  bulk        Build a bulk crystal from a space group and lattice parameters\n"
"  gb          Build a CSL grain boundary bicrystal from an existing structure\n"
"  poly        Build a polycrystalline microstructure from an existing structure\n"
"  nano        Carve a nanocrystal from a bulk reference structure\n"
"  amorphous   Pack an amorphous structure by random sequential addition\n"
"  sss         Build a substitutional solid solution from a host structure\n"
"  dislocation Insert a dislocation displacement field into a structure\n"
"  custom      Fill a 3D mesh model (OBJ/STL) with atoms from a reference crystal\n"
"\n"
"For detailed options per mode run:\n"
"  AtomForge --help bulk\n"
"  AtomForge --help gb\n"
"  AtomForge --help poly\n"
"  AtomForge --help nano\n"
"  AtomForge --help amorphous\n"
"  AtomForge --help sss\n"
"  AtomForge --help dislocation\n"
"  AtomForge --help custom\n"
<< std::endl;
}

static bool parseIvec3(const char* text, glm::ivec3& out)
{
    if (!text)
        return false;
    std::istringstream iss(text);
    int x = 0, y = 0, z = 0;
    if (!(iss >> x >> y >> z))
        return false;
    out = glm::ivec3(x, y, z);
    return true;
}

static bool parseVec3(const char* text, glm::vec3& out)
{
    if (!text)
        return false;
    std::istringstream iss(text);
    float x = 0.0f, y = 0.0f, z = 0.0f;
    if (!(iss >> x >> y >> z))
        return false;
    out = glm::vec3(x, y, z);
    return true;
}

static int runDislocation(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    const char* outPath   = findArg(argc, argv, "--output");

    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build dislocation\n";
        return 1;
    }
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    Structure base;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, base, loadErr))
    {
        std::cerr << "Error loading input structure: " << loadErr << "\n";
        return 1;
    }

    DislocationParams params;
    params.autoDirections = !hasFlag(argc, argv, "--manual-vectors");

    const char* character = findArg(argc, argv, "--character");
    if (character)
    {
        std::string c = character;
        for (char& ch : c) ch = (char)std::tolower((unsigned char)ch);
        if (c == "edge") params.character = DislocationCharacter::Edge;
        else if (c == "screw") params.character = DislocationCharacter::Screw;
        else if (c == "mixed") params.character = DislocationCharacter::Mixed;
        else
        {
            std::cerr << "Error: unsupported --character value '" << character << "'\n";
            return 1;
        }
    }

    const char* shape = findArg(argc, argv, "--shape");
    if (shape)
    {
        std::string s = shape;
        for (char& ch : s) ch = (char)std::tolower((unsigned char)ch);
        if (s == "halfplane") params.shape = DislocationShape::HalfPlane;
        else if (s == "cylinder") params.shape = DislocationShape::Cylinder;
        else if (s == "sphere") params.shape = DislocationShape::Sphere;
        else if (s == "ellipsoid") params.shape = DislocationShape::Ellipsoid;
        else if (s == "freeform") params.shape = DislocationShape::Freeform2D;
        else
        {
            std::cerr << "Error: unsupported --shape value '" << shape << "'\n";
            return 1;
        }
    }

    {
        const char* p = findArg(argc, argv, "--plane");
        if (p && !parseIvec3(p, params.planeHkl))
        {
            std::cerr << "Error: could not parse --plane \"h k l\"\n";
            return 1;
        }
    }
    {
        const char* b = findArg(argc, argv, "--burgers");
        if (b && !parseIvec3(b, params.burgersUvw))
        {
            std::cerr << "Error: could not parse --burgers \"u v w\"\n";
            return 1;
        }
    }
    {
        const char* l = findArg(argc, argv, "--line");
        if (l && !parseIvec3(l, params.lineUvw))
        {
            std::cerr << "Error: could not parse --line \"u v w\"\n";
            return 1;
        }
    }
    {
        const char* lf = findArg(argc, argv, "--line-frac");
        if (lf)
        {
            params.useFractionalLinePoint = true;
            if (!parseVec3(lf, params.linePointFractional))
            {
                std::cerr << "Error: could not parse --line-frac \"fx fy fz\"\n";
                return 1;
            }
        }
    }
    {
        const char* lo = findArg(argc, argv, "--line-offset");
        if (lo && !parseVec3(lo, params.linePointCartesianOffset))
        {
            std::cerr << "Error: could not parse --line-offset \"x y z\"\n";
            return 1;
        }
    }

    params.burgersScale = (float)argDouble(argc, argv, "--bscale", params.burgersScale);
    params.burgersOverrideMagnitude = (float)argDouble(argc, argv, "--bmag", params.burgersOverrideMagnitude);
    params.mixedCharacterAngleDeg = (float)argDouble(argc, argv, "--mixed-angle", params.mixedCharacterAngleDeg);
    params.poissonRatio = (float)argDouble(argc, argv, "--nu", params.poissonRatio);
    params.coreRadius = (float)argDouble(argc, argv, "--core", params.coreRadius);
    params.cutoffRadius = (float)argDouble(argc, argv, "--cutoff", params.cutoffRadius);
    params.lineHalfLength = (float)argDouble(argc, argv, "--line-half", params.lineHalfLength);

    params.cylinderRadius = (float)argDouble(argc, argv, "--cyl-radius", params.cylinderRadius);
    params.sphereRadius = (float)argDouble(argc, argv, "--sphere-radius", params.sphereRadius);
    {
        const char* er = findArg(argc, argv, "--ellipsoid");
        if (er && !parseVec3(er, params.ellipsoidRadii))
        {
            std::cerr << "Error: could not parse --ellipsoid \"rx ry rz\"\n";
            return 1;
        }
    }

    const char* poly = findArg(argc, argv, "--poly2d");
    if (poly)
    {
        params.freeformPoints.clear();
        std::string all = poly;
        std::stringstream ss(all);
        std::string token;
        while (std::getline(ss, token, ';'))
        {
            std::istringstream pss(token);
            float x = 0.0f, y = 0.0f;
            if (!(pss >> x >> y))
            {
                std::cerr << "Error: could not parse polygon point '" << token << "' in --poly2d\n";
                return 1;
            }
            params.freeformPoints.push_back(glm::vec2(x, y));
        }
        if (params.freeformPoints.size() < 3)
        {
            std::cerr << "Error: --poly2d needs at least 3 points\n";
            return 1;
        }
    }

    DislocationResult result = buildDislocation(base, params);
    if (!result.success)
    {
        std::cerr << "Error building dislocation: " << result.message << "\n";
        return 1;
    }

    std::cout << result.message << "\n";
    std::cout << "Shifted atoms: " << result.shiftedAtomCount << "\n";
    std::cout << "Validation: " << result.validation.message << "\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(result.output, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Bulk builder ─────────────────────────────────────────────────────────────

static int runBulk(int argc, char* argv[])
{
    // Crystal system
    CrystalSystem system = CrystalSystem::Cubic;
    const char* sysStr = findArg(argc, argv, "--system");
    if (sysStr)
    {
        std::string s = sysStr;
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if      (s == "triclinic")    system = CrystalSystem::Triclinic;
        else if (s == "monoclinic")   system = CrystalSystem::Monoclinic;
        else if (s == "orthorhombic") system = CrystalSystem::Orthorhombic;
        else if (s == "tetragonal")   system = CrystalSystem::Tetragonal;
        else if (s == "trigonal")     system = CrystalSystem::Trigonal;
        else if (s == "hexagonal")    system = CrystalSystem::Hexagonal;
        else if (s == "cubic")        system = CrystalSystem::Cubic;
        else {
            std::cerr << "Error: unknown crystal system '" << sysStr << "'\n";
            return 1;
        }
    }

    int sgNumber = argInt(argc, argv, "--spacegroup", 225);

    LatticeParameters lp;
    lp.a = argDouble(argc, argv, "--a", 4.0);
    lp.b = argDouble(argc, argv, "--b", lp.a);
    lp.c = argDouble(argc, argv, "--c", lp.a);
    lp.alpha = argDouble(argc, argv, "--alpha", 90.0);
    lp.beta  = argDouble(argc, argv, "--beta",  90.0);
    lp.gamma = argDouble(argc, argv, "--gamma", 90.0);
    applySystemConstraints(system, lp);

    // Validate
    std::string validErr;
    if (!validateParameters(system, lp, validErr))
    {
        std::cerr << "Error: invalid lattice parameters – " << validErr << "\n";
        return 1;
    }

    // Asymmetric atoms: each value is "SYMBOL fx fy fz"
    auto atomStrs = findAllArgs(argc, argv, "--atom");
    auto elementColors = makeDefaultElementColors();

    std::vector<AtomSite> asymAtoms;
    for (const auto& str : atomStrs)
    {
        std::istringstream iss(str);
        std::string sym;
        double fx, fy, fz;
        if (!(iss >> sym >> fx >> fy >> fz))
        {
            std::cerr << "Error: cannot parse --atom value '" << str
                      << "'.  Expected: \"SYMBOL fx fy fz\"\n";
            return 1;
        }
        int z = atomicNumberFromSymbol(sym);
        if (z <= 0)
        {
            std::cerr << "Error: unknown element symbol '" << sym << "'\n";
            return 1;
        }
        AtomSite site;
        site.x = fx; site.y = fy; site.z = fz;
        applyElementToAtom(site, z, elementColors);
        asymAtoms.push_back(site);
    }

    if (asymAtoms.empty())
    {
        std::cerr << "Error: no --atom arguments specified.  "
                     "At least one asymmetric unit atom is required.\n";
        return 1;
    }

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    // Build
    Structure structure;
    BulkBuildResult result = buildBulkCrystal(structure, system, sgNumber,
                                               lp, asymAtoms, elementColors);
    if (!result.success)
    {
        std::cerr << "Error building bulk crystal: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built bulk crystal: " << result.generatedAtoms << " atoms, "
              << "space group " << result.spaceGroupNumber
              << " (" << result.spaceGroupSymbol << ")\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(structure, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── GB builder ───────────────────────────────────────────────────────────────

static int runGB(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build gb\n";
        return 1;
    }

    Structure inputStructure;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, inputStructure, loadErr))
    {
        std::cerr << "Error loading input structure: " << loadErr << "\n";
        return 1;
    }
    if (!inputStructure.hasUnitCell)
    {
        std::cerr << "Error: input structure has no unit cell\n";
        return 1;
    }

    // Rotation axis
    int axis[3] = {0, 0, 1};
    const char* axisStr = findArg(argc, argv, "--axis");
    if (axisStr)
    {
        std::istringstream iss(axisStr);
        if (!(iss >> axis[0] >> axis[1] >> axis[2]))
        {
            std::cerr << "Error: cannot parse --axis value '" << axisStr
                      << "'.  Expected: \"u v w\" (three integers)\n";
            return 1;
        }
    }

    int sigmaTarget = argInt(argc, argv, "--sigma",    0);
    int sigmaMax    = argInt(argc, argv, "--sigmamax", 100);
    int planeIdx    = argInt(argc, argv, "--plane",    0);
    int ucA         = argInt(argc, argv, "--uca",      1);
    int ucB         = argInt(argc, argv, "--ucb",      1);
    float vacuum    = static_cast<float>(argDouble(argc, argv, "--vacuum",  0.0));
    float gap       = static_cast<float>(argDouble(argc, argv, "--gap",     0.0));
    float overlapD  = static_cast<float>(argDouble(argc, argv, "--overlap", 0.0));
    bool conventional = hasFlag(argc, argv, "--conventional");

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    // Compute sigma candidates
    std::vector<SigmaCandidate> candidates = computeGBInfo(axis, sigmaMax);
    if (candidates.empty())
    {
        std::cerr << "Error: no Σ candidates found for axis ["
                  << axis[0] << " " << axis[1] << " " << axis[2]
                  << "] with max Σ=" << sigmaMax << "\n";
        return 1;
    }

    // Select sigma candidate
    int selIdx = 0;
    if (sigmaTarget > 0)
    {
        bool found = false;
        for (int i = 0; i < (int)candidates.size(); ++i)
        {
            if (candidates[i].sigma == sigmaTarget)
            {
                selIdx = i;
                found = true;
                break;
            }
        }
        if (!found)
        {
            std::cerr << "Error: Σ=" << sigmaTarget
                      << " not found for axis ["
                      << axis[0] << " " << axis[1] << " " << axis[2] << "]\n";
            std::cerr << "Available Σ values:";
            for (auto& c : candidates) std::cerr << " " << c.sigma;
            std::cerr << "\n";
            return 1;
        }
    }
    // else selIdx=0 = smallest available sigma

    if (planeIdx < 0 || planeIdx > 2)
    {
        std::cerr << "Error: --plane must be 0, 1, or 2\n";
        return 1;
    }

    const SigmaCandidate& sel = candidates[selIdx];
    const int* plane = sel.plane[planeIdx].data();
    const int direction = planeIdx; // stacking direction = plane index

    std::cout << "Using Σ=" << sel.sigma
              << " (θ=" << sel.thetaDeg << "°)"
              << ", axis [" << axis[0] << " " << axis[1] << " " << axis[2] << "]"
              << ", plane (" << plane[0] << " " << plane[1] << " " << plane[2] << ")\n";

    // Build the bicrystal (mirrors the dialog build logic)
    Grain inputGrain = structureToGrain(inputStructure);

    // CSL^T supercell
    int cslT[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cslT[i][j] = sel.csl[j][i];

    Grain grainA = makeSupercell(inputGrain, cslT);

    if (!isOrthogonal(grainA.cell))
        grainA = setOrthogonalGrain(grainA, direction);

    int scaleB[3] = {1, 1, 1};
    scaleB[direction] = ucB;
    Grain tempA = makeSupercellDiag(grainA, scaleB[0], scaleB[1], scaleB[2]);

    Grain grainB = getBFromA(tempA);
    {
        int ones[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        grainB = makeSupercell(grainB, ones);
    }

    int scaleA[3] = {1, 1, 1};
    scaleA[direction] = ucA;
    grainA = makeSupercellDiag(grainA, scaleA[0], scaleA[1], scaleA[2]);

    Grain gb = stackGrains(grainA, grainB, direction, vacuum, gap);
    int removed = removeOverlaps(gb.atoms, overlapD);

    Structure structure = grainToStructure(gb);

    // PBC boundary tolerance from minimum layer spacing
    {
        double inv[3][3];
        invertCell(gb.cell, inv);
        std::vector<double> layers;
        layers.reserve(structure.atoms.size());
        for (const auto& a : structure.atoms)
        {
            double cart[3] = {a.x, a.y, a.z};
            double frac[3];
            cartToFrac(cart, inv, frac);
            layers.push_back(wrapFrac(frac[direction]));
        }
        std::sort(layers.begin(), layers.end());
        double minSpacing = 1.0;
        for (size_t i = 1; i < layers.size(); i++)
        {
            double d = layers[i] - layers[i - 1];
            if (d > 1e-8 && d < minSpacing) minSpacing = d;
        }
        structure.pbcBoundaryTol = static_cast<float>(minSpacing * 0.5 + 1e-6);
    }

    if (!conventional)
        reduceToPrimitiveGB(structure, direction);

    // Assign grain region IDs
    structure.grainRegionIds.assign(structure.atoms.size(), 0);
    {
        double cell[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                cell[i][j] = structure.cellVectors[i][j];
        double inv[3][3];
        invertCell(cell, inv);
        for (size_t i = 0; i < structure.atoms.size(); ++i)
        {
            double cart[3] = {structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z};
            double frac[3];
            cartToFrac(cart, inv, frac);
            structure.grainRegionIds[i] = (wrapFrac(frac[direction]) >= 0.5) ? 1 : 0;
        }
    }

    std::cout << "Built GB bicrystal: " << (int)structure.atoms.size() << " atoms"
              << " (removed " << removed << " overlapping)\n";
    std::cout << "Boundary type: " << classifyBoundaryType(axis, plane) << "\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(structure, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Polycrystal builder ───────────────────────────────────────────────────────

static int runPoly(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build poly\n";
        return 1;
    }

    Structure reference;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, reference, loadErr))
    {
        std::cerr << "Error loading input structure: " << loadErr << "\n";
        return 1;
    }
    if (!reference.hasUnitCell)
    {
        std::cerr << "Error: input structure has no unit cell\n";
        return 1;
    }

    PolyParams params;
    params.sizeX   = static_cast<float>(argDouble(argc, argv, "--sizex", 50.0));
    params.sizeY   = static_cast<float>(argDouble(argc, argv, "--sizey", 50.0));
    params.sizeZ   = static_cast<float>(argDouble(argc, argv, "--sizez", 50.0));
    params.numGrains = argInt(argc, argv, "--grains", 8);
    params.seed      = argInt(argc, argv, "--seed",   42);

    // Euler angles: each value is "phi1 Phi phi2"
    auto eulerStrs = findAllArgs(argc, argv, "--euler");
    if (!eulerStrs.empty())
    {
        for (int i = 0; i < (int)eulerStrs.size(); ++i)
        {
            std::istringstream iss(eulerStrs[i]);
            float phi1, Phi, phi2;
            if (!(iss >> phi1 >> Phi >> phi2))
            {
                std::cerr << "Error: cannot parse --euler value '" << eulerStrs[i]
                          << "'.  Expected: \"phi1 Phi phi2\"\n";
                return 1;
            }
            GrainOrientation go;
            go.grainIndex = i;
            go.phi1 = phi1;
            go.Phi  = Phi;
            go.phi2 = phi2;
            params.specifiedOrientations.push_back(go);
        }

        if ((int)eulerStrs.size() >= params.numGrains)
            params.orientationMode = GrainOrientationMode::AllSpecified;
        else
            params.orientationMode = GrainOrientationMode::PartialSpecified;
    }

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    auto elementColors = makeDefaultElementColors();
    Structure structure;
    PolyBuildResult result = buildPolycrystal(structure, reference, params, elementColors);
    if (!result.success)
    {
        std::cerr << "Error building polycrystal: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built polycrystal: " << result.outputAtoms << " atoms, "
              << result.numGrains << " grains\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(structure, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Nanocrystal builder ───────────────────────────────────────────────────────

static int runNano(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build nano\n";
        return 1;
    }

    Structure reference;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, reference, loadErr))
    {
        std::cerr << "Error loading input structure: " << loadErr << "\n";
        return 1;
    }
    if (!reference.hasUnitCell)
    {
        std::cerr << "Error: input structure has no unit cell\n";
        return 1;
    }

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    NanoParams params;
    params.generationMode = NanoGenerationMode::Shape;

    const char* shapeStr = findArg(argc, argv, "--shape");
    if (shapeStr)
    {
        std::string s = shapeStr;
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if      (s == "sphere")               params.shape = NanoShape::Sphere;
        else if (s == "ellipsoid")            params.shape = NanoShape::Ellipsoid;
        else if (s == "box")                  params.shape = NanoShape::Box;
        else if (s == "cylinder")             params.shape = NanoShape::Cylinder;
        else if (s == "octahedron")           params.shape = NanoShape::Octahedron;
        else if (s == "truncated-octahedron") params.shape = NanoShape::TruncatedOctahedron;
        else if (s == "cuboctahedron")        params.shape = NanoShape::Cuboctahedron;
        else {
            std::cerr << "Error: unknown shape '" << shapeStr << "'\n";
            return 1;
        }
    }

    float defaultRadius = static_cast<float>(argDouble(argc, argv, "--radius", 15.0));
    params.sphereRadius   = defaultRadius;
    params.octRadius      = defaultRadius;
    params.truncOctRadius = static_cast<float>(argDouble(argc, argv, "--radius", 18.0));
    params.cuboRadius     = defaultRadius;

    params.ellipRx = static_cast<float>(argDouble(argc, argv, "--rx", 15.0));
    params.ellipRy = static_cast<float>(argDouble(argc, argv, "--ry", 12.0));
    params.ellipRz = static_cast<float>(argDouble(argc, argv, "--rz", 10.0));

    params.boxHx = static_cast<float>(argDouble(argc, argv, "--hx", 15.0));
    params.boxHy = static_cast<float>(argDouble(argc, argv, "--hy", 15.0));
    params.boxHz = static_cast<float>(argDouble(argc, argv, "--hz", 15.0));

    params.truncOctTrunc = static_cast<float>(argDouble(argc, argv, "--trunc", 12.0));

    params.cylRadius = static_cast<float>(argDouble(argc, argv, "--cylradius", 12.0));
    params.cylHeight = static_cast<float>(argDouble(argc, argv, "--cylheight", 30.0));
    params.cylAxis   = argInt(argc, argv, "--cylaxis", 2);

    params.vacuumPadding = static_cast<float>(argDouble(argc, argv, "--vacuum", 5.0));
    params.setOutputCell = true;
    params.autoCenterFromAtoms = true;

    int repA = argInt(argc, argv, "--repa", 0);
    int repB = argInt(argc, argv, "--repb", 0);
    int repC = argInt(argc, argv, "--repc", 0);
    if (repA > 0 || repB > 0 || repC > 0)
    {
        params.autoReplicate = false;
        params.repA = (repA > 0) ? repA : 5;
        params.repB = (repB > 0) ? repB : 5;
        params.repC = (repC > 0) ? repC : 5;
    }
    else
    {
        params.autoReplicate = true;
    }

    auto elementColors = makeDefaultElementColors();
    Structure structure;
    NanoBuildResult result = buildNanocrystal(structure, reference, params,
                                               elementColors, {}, {});
    if (!result.success)
    {
        std::cerr << "Error building nanocrystal: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built nanocrystal: " << result.outputAtoms << " atoms"
              << " (shape: " << shapeLabel(params.shape) << ")\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(structure, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Amorphous builder ─────────────────────────────────────────────────────────

static int runAmorphous(int argc, char* argv[])
{
    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    // Element specs: each value is "SYMBOL N"
    auto elemStrs = findAllArgs(argc, argv, "--element");
    if (elemStrs.empty())
    {
        std::cerr << "Error: at least one --element \"SYMBOL N\" is required\n";
        return 1;
    }

    AmorphousParams params;
    for (const auto& str : elemStrs)
    {
        std::istringstream iss(str);
        std::string sym;
        int count;
        if (!(iss >> sym >> count))
        {
            std::cerr << "Error: cannot parse --element value '" << str
                      << "'.  Expected: \"SYMBOL N\"\n";
            return 1;
        }
        int z = atomicNumberFromSymbol(sym);
        if (z <= 0)
        {
            std::cerr << "Error: unknown element symbol '" << sym << "'\n";
            return 1;
        }
        AmorphousElementSpec spec;
        spec.atomicNumber = z;
        spec.count        = count;
        params.elements.push_back(spec);
    }

    // Box mode: manual if any of boxa/boxb/boxc are given, else auto density
    bool hasManualBox = findArg(argc, argv, "--boxa") ||
                        findArg(argc, argv, "--boxb") ||
                        findArg(argc, argv, "--boxc");
    if (hasManualBox)
    {
        params.boxMode = AmorphousBoxMode::Manual;
        params.boxA = static_cast<float>(argDouble(argc, argv, "--boxa", 20.0));
        params.boxB = static_cast<float>(argDouble(argc, argv, "--boxb", 20.0));
        params.boxC = static_cast<float>(argDouble(argc, argv, "--boxc", 20.0));
    }
    else
    {
        params.boxMode        = AmorphousBoxMode::AutoFromDensity;
        params.targetDensity  = static_cast<float>(argDouble(argc, argv, "--density", 2.0));
    }

    params.cellScaleFactor   = static_cast<float>(argDouble(argc, argv, "--scale",    1.0));
    params.covalentTolerance = static_cast<float>(argDouble(argc, argv, "--covtol",   0.75));
    params.seed              = static_cast<unsigned int>(argInt(argc, argv, "--seed",  42));
    params.maxAttempts       = argInt(argc, argv, "--attempts", 1000);

    // Per-pair minimum distances: each value is "Z1 Z2 dist"
    auto pairStrs = findAllArgs(argc, argv, "--mindist");
    for (const auto& str : pairStrs)
    {
        std::istringstream iss(str);
        std::string s1, s2;
        float dist;
        if (!(iss >> s1 >> s2 >> dist))
        {
            std::cerr << "Error: cannot parse --mindist value '" << str
                      << "'.  Expected: \"SYMBOL1 SYMBOL2 dist\"\n";
            return 1;
        }
        int z1 = atomicNumberFromSymbol(s1);
        int z2 = atomicNumberFromSymbol(s2);
        if (z1 <= 0 || z2 <= 0)
        {
            std::cerr << "Error: unknown element symbol in --mindist '" << str << "'\n";
            return 1;
        }
        AmorphousPairDist pd;
        pd.z1 = std::min(z1, z2);
        pd.z2 = std::max(z1, z2);
        pd.minDist = dist;
        params.pairDistances.push_back(pd);
    }

    params.periodic = true;

    auto elementColors  = makeDefaultElementColors();
    auto covalentRadii  = makeLiteratureCovalentRadii();

    AmorphousResult result = buildAmorphousStructure(params, covalentRadii, elementColors);
    if (!result.success)
    {
        std::cerr << "Error building amorphous structure: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built amorphous structure: " << result.placedAtoms << " atoms"
              << " (density: " << result.actualDensity << " g/cm3)\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(result.output, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Interface builder ─────────────────────────────────────────────────────────

[[maybe_unused]] static int runInterface(int argc, char* argv[])
{
    const char* pathA = findArg(argc, argv, "--layerA");
    const char* pathB = findArg(argc, argv, "--layerB");
    if (!pathA || !pathB)
    {
        std::cerr << "Error: --layerA <file> and --layerB <file> are required\n";
        return 1;
    }

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    Structure sA, sB;
    std::string err;
    if (!loadStructureFromFile(pathA, sA, err))
    {
        std::cerr << "Error loading layer A: " << err << "\n";
        return 1;
    }
    if (!loadStructureFromFile(pathB, sB, err))
    {
        std::cerr << "Error loading layer B: " << err << "\n";
        return 1;
    }
    if (!sA.hasUnitCell || !sB.hasUnitCell)
    {
        std::cerr << "Error: both input structures must have a unit cell\n";
        return 1;
    }

    int  nmax      = argInt   (argc, argv, "--nmax",     4);
    int  maxCells  = argInt   (argc, argv, "--maxcells", 16);
    int  pickIdx   = argInt   (argc, argv, "--pick",     0);
    int  layersA   = argInt   (argc, argv, "--layersA",  1);
    int  layersB   = argInt   (argc, argv, "--layersB",  1);
    double gap     = argDouble(argc, argv, "--gap",      2.0);
    double vacuum  = argDouble(argc, argv, "--vacuum",   10.0);
    int  repx      = argInt   (argc, argv, "--repx",     1);
    int  repy      = argInt   (argc, argv, "--repy",     1);

    // Build 2D bases
    double basisA[2][2], basisB[2][2];
    get2DBasis(sA, basisA);
    get2DBasis(sB, basisB);

    // Enumerate supercell candidates for both layers
    auto cellsA = generateUniqueSupercells(basisA, nmax, maxCells);
    auto cellsB = generateUniqueSupercells(basisB, nmax, maxCells);

    if (cellsA.empty() || cellsB.empty())
    {
        std::cerr << "Error: no supercell candidates found (try increasing --nmax or --maxcells)\n";
        return 1;
    }

    // Find best-strain pair
    struct Candidate {
        int iA, iB;
        double strain;
        double vA[2][2], vB[2][2];
    };
    std::vector<Candidate> ranked;
    for (int iA = 0; iA < (int)cellsA.size(); ++iA)
    {
        for (int iB = 0; iB < (int)cellsB.size(); ++iB)
        {
            double exx, eyy, exy;
            if (!strainComponents(cellsA[iA].vecs, cellsB[iB].vecs, exx, eyy, exy))
                continue;
            Candidate c;
            c.iA = iA; c.iB = iB;
            c.strain = meanAbsStrain(exx, eyy, exy);
            std::memcpy(c.vA, cellsA[iA].vecs, sizeof(c.vA));
            std::memcpy(c.vB, cellsB[iB].vecs, sizeof(c.vB));
            ranked.push_back(c);
        }
    }

    if (ranked.empty())
    {
        std::cerr << "Error: no matching supercell pairs found\n";
        return 1;
    }

    std::sort(ranked.begin(), ranked.end(),
              [](const Candidate& a, const Candidate& b){ return a.strain < b.strain; });

    if (pickIdx >= (int)ranked.size())
    {
        std::cerr << "Error: --pick " << pickIdx << " is out of range ("
                  << ranked.size() << " candidates available)\n";
        return 1;
    }

    const Candidate& best = ranked[pickIdx];
    std::cout << "Selected supercell pair " << pickIdx
              << " (mean |strain| = " << best.strain << ")\n";

    // Build layer supercells
    Structure superA = makeSupercell2D(sA, cellsA[best.iA].mat);
    Structure superB = makeSupercell2D(sB, cellsB[best.iB].mat);

    if (layersA > 1) superA = repeatLayersZ(superA, layersA);
    if (layersB > 1) superB = repeatLayersZ(superB, layersB);

    // Strain layer B to match A's 2D cell
    Structure strainedB = applyTransform2D(superB, best.vA);

    // Stack into interface
    Structure iface = assembleInterface(superA, strainedB, gap, vacuum);

    if (repx > 1 || repy > 1)
        iface = repeatInterfaceXY(iface, repx, repy);

    std::cout << "Built interface: " << (int)iface.atoms.size() << " atoms\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(iface, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Custom mesh-fill builder ─────────────────────────────────────────────────

static int runCustom(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    const char* meshPath  = findArg(argc, argv, "--mesh");
    const char* outPath   = findArg(argc, argv, "--output");

    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build custom\n";
        return 1;
    }
    if (!meshPath)
    {
        std::cerr << "Error: --mesh <file> (OBJ or STL) is required for --build custom\n";
        return 1;
    }
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    // Load reference crystal
    Structure reference;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, reference, loadErr))
    {
        std::cerr << "Error loading reference structure: " << loadErr << "\n";
        return 1;
    }

    // Load mesh
    std::vector<glm::vec3>    modelVertices;
    std::vector<unsigned int> modelIndices;
    std::string meshErr;
    {
        std::string mp = meshPath;
        std::string ext;
        auto dot = mp.rfind('.');
        if (dot != std::string::npos)
            ext = mp.substr(dot + 1);
        for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

        bool ok = false;
        if (ext == "obj")
            ok = parseObjMesh(mp, modelVertices, modelIndices, meshErr);
        else if (ext == "stl")
            ok = parseStlMesh(mp, modelVertices, modelIndices, meshErr);
        else
        {
            std::cerr << "Error: unsupported mesh format '." << ext
                      << "'.  Use OBJ or STL.\n";
            return 1;
        }
        if (!ok)
        {
            std::cerr << "Error loading mesh: " << meshErr << "\n";
            return 1;
        }
    }
    if (modelVertices.empty() || modelIndices.empty())
    {
        std::cerr << "Error: mesh file contains no geometry\n";
        return 1;
    }

    // Compute mesh bounding half-extents (needed by buildNanocrystal)
    glm::vec3 mn = modelVertices[0], mx = modelVertices[0];
    for (const auto& v : modelVertices)
    {
        mn = glm::min(mn, v);
        mx = glm::max(mx, v);
    }
    glm::vec3 halfExt = (mx - mn) * 0.5f;

    NanoParams params;
    params.generationMode      = NanoGenerationMode::Shape;
    params.shape               = NanoShape::MeshModel;
    params.modelScale          = static_cast<float>(argDouble(argc, argv, "--scale",  1.0));
    params.modelHx             = halfExt.x;
    params.modelHy             = halfExt.y;
    params.modelHz             = halfExt.z;
    params.vacuumPadding       = static_cast<float>(argDouble(argc, argv, "--vacuum", 5.0));
    params.setOutputCell       = true;
    params.autoCenterFromAtoms = true;

    int repA = argInt(argc, argv, "--repa", 0);
    int repB = argInt(argc, argv, "--repb", 0);
    int repC = argInt(argc, argv, "--repc", 0);
    if (repA > 0 || repB > 0 || repC > 0)
    {
        params.autoReplicate = false;
        params.repA = (repA > 0) ? repA : 5;
        params.repB = (repB > 0) ? repB : 5;
        params.repC = (repC > 0) ? repC : 5;
    }
    else
    {
        params.autoReplicate = true;
    }

    // Orientation
    const char* millerStr = findArg(argc, argv, "--miller");
    if (millerStr)
    {
        params.applyCrystalOrientation = true;
        params.useMillerOrientation    = true;
        int h = 1, k = 0, l = 0;
        std::istringstream iss(millerStr);
        iss >> h >> k >> l;
        params.millerH = h;
        params.millerK = k;
        params.millerL = l;
    }
    else
    {
        float rx = static_cast<float>(argDouble(argc, argv, "--rotx", 0.0));
        float ry = static_cast<float>(argDouble(argc, argv, "--roty", 0.0));
        float rz = static_cast<float>(argDouble(argc, argv, "--rotz", 0.0));
        if (rx != 0.0f || ry != 0.0f || rz != 0.0f)
        {
            params.applyCrystalOrientation = true;
            params.useMillerOrientation    = false;
            params.orientXDeg = rx;
            params.orientYDeg = ry;
            params.orientZDeg = rz;
        }
    }

    auto elementColors = makeDefaultElementColors();
    Structure structure;
    NanoBuildResult result = buildNanocrystal(structure, reference, params,
                                               elementColors,
                                               modelVertices, modelIndices);
    if (!result.success)
    {
        std::cerr << "Error building custom structure: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built custom mesh-fill: " << result.outputAtoms << " atoms\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(structure, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Substitutional solid solution builder ────────────────────────────────────

static int runSSS(int argc, char* argv[])
{
    const char* inputPath = findArg(argc, argv, "--input");
    if (!inputPath)
    {
        std::cerr << "Error: --input <file> is required for --build sss\n";
        return 1;
    }

    const char* outPath = findArg(argc, argv, "--output");
    if (!outPath)
    {
        std::cerr << "Error: --output <file> is required\n";
        return 1;
    }

    Structure base;
    std::string loadErr;
    if (!loadStructureFromFile(inputPath, base, loadErr))
    {
        std::cerr << "Error loading input structure: " << loadErr << "\n";
        return 1;
    }

    // Fractions: single "Cu=0.7,Zn=0.3" string
    const char* fracArg = findArg(argc, argv, "--frac");
    if (!fracArg)
    {
        std::cerr << "Error: --frac \"SYM=frac,...\" is required (e.g. --frac \"Cu=0.7,Zn=0.3\")\n";
        return 1;
    }

    SSSParams params;
    params.seed = static_cast<unsigned int>(argInt(argc, argv, "--seed", 12345));

    {
        std::string fracStr = fracArg;
        std::istringstream stream(fracStr);
        std::string token;
        while (std::getline(stream, token, ','))
        {
            // trim whitespace
            auto b = token.find_first_not_of(" \t");
            auto e = token.find_last_not_of(" \t");
            if (b == std::string::npos) continue;
            token = token.substr(b, e - b + 1);

            auto eq = token.find('=');
            if (eq == std::string::npos)
            {
                std::cerr << "Error: cannot parse --frac token '" << token
                          << "'.  Expected format: SYM=fraction\n";
                return 1;
            }
            std::string sym  = token.substr(0, eq);
            std::string fstr = token.substr(eq + 1);
            float frac;
            try { frac = std::stof(fstr); }
            catch (...)
            {
                std::cerr << "Error: cannot parse fraction '" << fstr
                          << "' in --frac\n";
                return 1;
            }
            int z = atomicNumberFromSymbol(sym);
            if (z <= 0)
            {
                std::cerr << "Error: unknown element symbol '" << sym << "' in --frac\n";
                return 1;
            }
            SSSElementFraction ef;
            ef.atomicNumber = z;
            ef.fraction     = frac;
            params.composition.push_back(ef);
        }
        if (params.composition.empty())
        {
            std::cerr << "Error: --frac string parsed no valid entries\n";
            return 1;
        }
    }

    SSSResult result = buildSubstitutionalSolidSolution(base, params);
    if (!result.success)
    {
        std::cerr << "Error building solid solution: " << result.message << "\n";
        return 1;
    }

    std::cout << "Built solid solution: " << (int)result.output.atoms.size() << " atoms\n";
    for (auto& ec : result.elementCounts)
        std::cout << "  " << elementSymbol(ec.first) << ": " << ec.second << " atoms\n";

    std::string fmt = detectFormat(outPath);
    if (!saveStructure(result.output, outPath, fmt))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

// ── Public interface ──────────────────────────────────────────────────────────

bool isCLIMode(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--build")   == 0) return true;
        if (std::strcmp(argv[i], "--help")    == 0) return true;
        if (std::strcmp(argv[i], "-h")        == 0) return true;
        if (std::strcmp(argv[i], "--version") == 0) return true;
        if (std::strcmp(argv[i], "-v")        == 0) return true;
    }
    return false;
}

int runCLI(int argc, char* argv[])
{
    if (hasFlag(argc, argv, "--version") || hasFlag(argc, argv, "-v"))
    {
        std::cout << "AtomForge " << ATOMFORGE_VERSION << "\n";
        return 0;
    }

    if (hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h"))
    {
        // --help <mode>  →  per-mode detail
        const char* topic = findArg(argc, argv, "--help");
        if (!topic) topic = findArg(argc, argv, "-h");

        if (topic)
        {
            std::string t = topic;
            if      (t == "bulk")      { printHelpBulk();      return 0; }
            else if (t == "gb")        { printHelpGB();        return 0; }
            else if (t == "poly")      { printHelpPoly();      return 0; }
            else if (t == "nano")      { printHelpNano();      return 0; }
            else if (t == "amorphous") { printHelpAmorphous(); return 0; }
            else if (t == "sss")       { printHelpSSS();       return 0; }
            else if (t == "dislocation"){ printHelpDislocation(); return 0; }
            else if (t == "custom")    { printHelpCustom();    return 0; }
            else
            {
                std::cerr << "Unknown help topic '" << t
                          << "'.  Valid topics: bulk | gb | poly | nano | amorphous | sss | dislocation | custom\n";
                return 1;
            }
        }

        // bare --help
        printHelp();
        return 0;
    }

    const char* mode = findArg(argc, argv, "--build");
    if (!mode)
    {
        std::cerr << "Error: --build <mode> is required.  "
                     "Use --help for usage.\n";
        return 1;
    }

    try
    {
        std::string m = mode;
        if      (m == "bulk")        return runBulk       (argc, argv);
        else if (m == "gb")          return runGB         (argc, argv);
        else if (m == "poly")        return runPoly       (argc, argv);
        else if (m == "nano")        return runNano       (argc, argv);
        else if (m == "amorphous")   return runAmorphous  (argc, argv);
        else if (m == "sss")         return runSSS        (argc, argv);
        else if (m == "dislocation") return runDislocation(argc, argv);
        else if (m == "custom")      return runCustom     (argc, argv);

        std::cerr << "Error: unknown build mode '" << m
                  << "'.  Valid modes: bulk | gb | poly | nano | amorphous | sss | dislocation | custom\n";
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
