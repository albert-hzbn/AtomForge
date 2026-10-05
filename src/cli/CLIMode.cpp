#include "cli/CLIMode.h"
#include "cli/AnalysisCLI.h"
#include "cli/RenderCLI.h"
#include "cli/ScienceCLI.h"

#include "algorithms/AmorphousBuilder.h"
#include "algorithms/BulkCrystalBuilder.h"
#include "algorithms/CSLComputation.h"
#include "algorithms/DislocationBuilder.h"
#include "algorithms/InterfaceBuilder.h"
#include "algorithms/StackingFaultBuilder.h"
#include <filesystem>
#include "algorithms/MeshLoader.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "algorithms/NanostructureTools.h"
#include "algorithms/PolyCrystalBuilder.h"
#include "algorithms/SQSBuilder.h"
#include "algorithms/StrainTool.h"
#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "algorithms/SurfaceBuilder.h"
#include "algorithms/VacancyBuilder.h"
#include "io/StructureLoader.h"
#include "util/ElementData.h"
#include "util/PathUtils.h"

#include <algorithm>
#include <array>
#include <string_view>
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
"                               cuboctahedron | wulff (default: sphere)\n"
"  --facet \"h k l energy\"      Wulff facet family; repeat for each family.\n"
"                               --radius sets maximum Wulff plane distance.\n"
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
"  Anisotropic elasticity (Stroh sextic formalism, replaces --nu; elastic\n"
"  constants are given in GPa, in the SAME axes as the input structure's\n"
"  unit cell -- see the BABEL reference in AnisotropicDislocation.h):\n"
"  --anisotropic                Use anisotropic instead of isotropic elasticity\n"
"  --elastic-symmetry <cubic|hexagonal>  (default: cubic)\n"
"  --elastic-c11/--elastic-c12/--elastic-c44 <GPa>   cubic constants\n"
"  --elastic-c13/--elastic-c33 <GPa>   additional hexagonal constants\n"
"                                (hexagonal also needs c11/c12/c44; c-axis\n"
"                                 must be the structure's third cell vector)\n"
"  --elastic-noise <amplitude>  Relative perturbation to break the sextic's\n"
"                                degeneracy at isotropic/high-symmetry\n"
"                                orientations (default: 1e-4)\n"
"\n"
"  Dislocation dipole (adds a second, opposite-Burgers-vector dislocation,\n"
"                       so the pair's net Burgers vector -- and hence its\n"
"                       long-range field -- is periodicity-compatible):\n"
"  --dipole                     Enable dipole mode\n"
"  --dipole-offset \"dx dy\"     Partner offset in the local slip-plane axes,\n"
"                                Angstrom (default: \"15 0\")\n"
"\n"
"Example:\n"
"  AtomForge --build dislocation --input base.cfg --character edge ^\n"
"            --manual-vectors --line \"0 0 1\" --burgers \"1 0 0\" ^\n"
"            --shape cylinder --cyl-radius 12 --output edge.cfg\n"
"  AtomForge --build dislocation --input cu.cfg --character screw ^\n"
"            --anisotropic --elastic-c11 168.4 --elastic-c12 121.4 ^\n"
"            --elastic-c44 75.4 --shape cylinder --cyl-radius 12 --output cu_screw.cfg\n"
<< std::endl;
}

static void printHelp()
{
    std::cout <<
"AtomForge CLI - headless structure builder\n"
"\n"
"Usage:\n"
"  AtomForge --build <mode> [options] --output <file>\n"
"  AtomForge --analyze <cna|rdf|adf|sro|interstitial|sculpt> --input FILE --output FILE\n"
"  AtomForge --convert --input FILE --output FILE [--format FORMAT]\n"
"  AtomForge --render --input FILE --output FILE.png [options]\n"
"  AtomForge --science <tool> --input REQUEST.json --output RESULT.json\n"
"\n"
"Modes:\n"
"  bulk        Build a bulk crystal from a space group and lattice parameters\n"
"  gb          Build a CSL grain boundary bicrystal from an existing structure\n"
"  poly        Build a polycrystalline microstructure from an existing structure\n"
"  nano        Carve a nanocrystal from a bulk reference structure\n"
"  amorphous   Pack an amorphous structure by random sequential addition\n"
"  sss         Build a substitutional solid solution from a host structure\n"
"  dislocation Insert a dislocation displacement field into a structure\n"
"  interface   Match and assemble two periodic layers\n"
"  stacking-fault Generate a sliding stacking-fault sequence\n"
"  custom      Fill a 3D mesh model (OBJ/STL) with atoms from a reference crystal\n"
"  vacancy     Remove atoms to create vacancies at a percentage or count\n"
"  strain      Apply a homogeneous deformation to the cell\n"
"  primitive   Reduce a structure to its standardized primitive cell\n"
"  surface     Cleave a vacuum-padded slab along a Miller plane\n"
"  sqs         Build a special quasirandom alloy on a fixed lattice\n"
"  nanowire    Cut a 1D-periodic wire with a circular or polygonal section\n"
"  core-shell  Relabel a finite particle into core and shell compositions\n"
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
"  AtomForge --help interface\n"
"  AtomForge --help stacking-fault\n"
"  AtomForge --help vacancy | strain | primitive | surface | sqs | nanowire | core-shell\n"
"  AtomForge --analyze cna --help\n"
"  AtomForge --render --help\n"
"  AtomForge --science --help\n"
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

    if (hasFlag(argc, argv, "--dipole"))
    {
        params.dipole = true;
        const char* off = findArg(argc, argv, "--dipole-offset");
        if (off)
        {
            std::istringstream iss(off);
            float dx = 0.0f, dy = 0.0f;
            if (!(iss >> dx >> dy))
            {
                std::cerr << "Error: could not parse --dipole-offset \"dx dy\"\n";
                return 1;
            }
            params.dipoleOffset = glm::vec2(dx, dy);
        }
    }

    if (hasFlag(argc, argv, "--anisotropic"))
    {
        params.anisotropicElasticity = true;
        const char* sym = findArg(argc, argv, "--elastic-symmetry");
        if (sym)
        {
            std::string s = sym;
            for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (s == "cubic") params.elasticSymmetry = DislocationParams::ElasticSymmetry::Cubic;
            else if (s == "hexagonal") params.elasticSymmetry = DislocationParams::ElasticSymmetry::Hexagonal;
            else
            {
                std::cerr << "Error: --elastic-symmetry must be cubic or hexagonal\n";
                return 1;
            }
        }
        params.elasticC11 = argDouble(argc, argv, "--elastic-c11", 0.0);
        params.elasticC12 = argDouble(argc, argv, "--elastic-c12", 0.0);
        params.elasticC44 = argDouble(argc, argv, "--elastic-c44", 0.0);
        params.elasticC13 = argDouble(argc, argv, "--elastic-c13", 0.0);
        params.elasticC33 = argDouble(argc, argv, "--elastic-c33", 0.0);
        params.elasticNoiseAmplitude = argDouble(argc, argv, "--elastic-noise", params.elasticNoiseAmplitude);
        if (params.elasticC11 <= 0.0 || params.elasticC12 < 0.0 || params.elasticC44 <= 0.0)
        {
            std::cerr << "Error: --anisotropic requires positive --elastic-c11/--elastic-c44 "
                         "(and --elastic-c12 >= 0), in GPa\n";
            return 1;
        }
        if (params.elasticSymmetry == DislocationParams::ElasticSymmetry::Hexagonal
            && (params.elasticC13 <= 0.0 || params.elasticC33 <= 0.0))
        {
            std::cerr << "Error: --elastic-symmetry hexagonal also requires positive "
                         "--elastic-c13/--elastic-c33, in GPa\n";
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
        else if (s == "wulff")               params.generationMode = NanoGenerationMode::WulffConstruction;
        else {
            std::cerr << "Error: unknown shape '" << shapeStr << "'\n";
            return 1;
        }
    }

    const auto facets = findAllArgs(argc, argv, "--facet");
    if (params.generationMode == NanoGenerationMode::WulffConstruction)
    {
        if (facets.empty())
            throw std::invalid_argument("Wulff construction requires at least one --facet");
        for (const auto& text : facets)
        {
            WulffPlaneInput plane;
            std::istringstream input(text);
            std::string extra;
            if (!(input >> plane.h >> plane.k >> plane.l >> plane.surfaceEnergy)
                || (input >> extra) || !std::isfinite(plane.surfaceEnergy)
                || plane.surfaceEnergy <= 0 || (plane.h == 0 && plane.k == 0 && plane.l == 0))
                throw std::invalid_argument("Invalid --facet: expected nonzero h k l and positive finite energy");
            params.wulffPlanes.push_back(plane);
        }
        params.wulffMaxRadius = static_cast<float>(argDouble(argc, argv, "--radius", 20.0));
    }
    else if (!facets.empty())
        throw std::invalid_argument("--facet requires --shape wulff");

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
              << " (shape: " << (params.generationMode == NanoGenerationMode::WulffConstruction
                  ? "Wulff" : shapeLabel(params.shape)) << ")\n";

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

    if (nmax < 1 || nmax > 8 || maxCells < 1 || maxCells > 64 || pickIdx < 0 ||
        layersA < 1 || layersB < 1 || layersA > 1000 || layersB > 1000 ||
        repx < 1 || repy < 1 || repx > 100 || repy > 100 || gap < 0 || vacuum < 0)
        throw std::invalid_argument("Invalid interface search bounds, repeats, gap or vacuum");

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
    const auto& u = best.vB;
    const auto& v = best.vA;
    const double determinant = u[0][0]*u[1][1] - u[1][0]*u[0][1];
    if (std::abs(determinant) < 1e-12) throw std::invalid_argument("Singular layer B basis");
    const double transform[2][2] = {
        {(v[0][0]*u[1][1]-v[1][0]*u[0][1])/determinant,
         (-v[0][0]*u[1][0]+v[1][0]*u[0][0])/determinant},
        {(v[0][1]*u[1][1]-v[1][1]*u[0][1])/determinant,
         (-v[0][1]*u[1][0]+v[1][1]*u[0][0])/determinant}};
    Structure strainedB = applyTransform2D(superB, transform);

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
            try {
                std::size_t parsed = 0;
                frac = std::stof(fstr, &parsed);
                if (parsed != fstr.size() || !std::isfinite(frac))
                    throw std::invalid_argument("invalid fraction");
            }
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

static void printHelpStackingFault()
{
    std::cout << "STACKING FAULT (--build stacking-fault)\n"
        "--input FILE --output FILE [--plane 0..6] [--layers 9] [--interval 0.1]\n"
        "[--maximum 2] [--frame 0] [--orthogonal] [--sequence DIRECTORY]\n"
        "Planes: 0 auto, 1 FCC111, 2 HCP basal, 3 HCP prismatic,\n"
        "4 HCP pyramidal, 5 BCC110, 6 BCC112. Sequence writes numbered VASP files.\n";
}

static int runStackingFault(int argc, char* argv[])
{
    const auto input = findArg(argc,argv,"--input");
    const auto output = findArg(argc,argv,"--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input,source,error)) throw std::runtime_error(error);
    StackingFaultParams params;
    const int plane=argInt(argc,argv,"--plane",0);
    if (plane<0 || plane>6) throw std::invalid_argument("Plane must be 0..6");
    params.plane=static_cast<StackingFaultPlane>(plane);
    params.layerCount=argInt(argc,argv,"--layers",9);
    params.interval=argDouble(argc,argv,"--interval",.1);
    params.maxDisplacementFactor=argDouble(argc,argv,"--maximum",2);
    if (hasFlag(argc,argv,"--orthogonal")) params.cellMode=StackingFaultCellMode::OrthogonalCell;
    const auto result=buildStackingFaultSequence(source,params);
    if (!result.success) throw std::runtime_error(result.message);
    const int frame=argInt(argc,argv,"--frame",0);
    if (frame<0 || frame>=(int)result.sequence.size()) throw std::invalid_argument("Frame out of range");
    if (!saveStructure(result.sequence[frame].structure,output,detectFormat(output)))
        throw std::runtime_error("Failed to save stacking fault");
    if (const auto directory=findArg(argc,argv,"--sequence")) {
        std::filesystem::create_directories(directory);
        for (std::size_t i=0;i<result.sequence.size();++i) {
            const auto path=std::filesystem::path(directory)/("frame-"+std::to_string(i)+".vasp");
            if (!saveStructure(result.sequence[i].structure,path.string(),"vasp"))
                throw std::runtime_error("Failed to save sequence frame");
        }
    }
    std::cout << result.message << "\nFrames: " << result.sequence.size() << '\n';
    return 0;
}

// -- Point defects: vacancy generator ----------------------------------------

static void printHelpVacancy()
{
    std::cout << "VACANCY (--build vacancy)\n"
        "--input FILE --output FILE [--element SYMBOL] [--percent P | --count N]\n"
        "[--min-separation D] [--seed S]\n"
        "Randomly removes atoms (optionally restricted to one element) to create\n"
        "point-defect vacancies at a target percentage or exact count.\n";
}

static int runVacancy(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::VacancyParams params;
    if (const auto element = findArg(argc, argv, "--element")) params.element = element;
    params.targetPercentage = argDouble(argc, argv, "--percent", 0.0);
    params.targetCount = argInt(argc, argv, "--count", 0);
    params.minSeparation = argDouble(argc, argv, "--min-separation", 0.0);
    params.seed = (unsigned)argInt(argc, argv, "--seed", 1);
    const auto result = atomforge::buildVacancies(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save vacancy structure");
    std::cout << result.message << '\n';
    return 0;
}

// -- Lattice: homogeneous strain ---------------------------------------------

static void printHelpStrain()
{
    std::cout << "STRAIN (--build strain)\n"
        "--input FILE --output FILE\n"
        "[--exx V] [--eyy V] [--ezz V] [--exy V] [--exz V] [--eyz V]\n"
        "[--matrix \"f11 f12 f13 f21 f22 f23 f31 f32 f33\"]\n"
        "Applies a homogeneous deformation to the cell, holding fractional\n"
        "coordinates fixed. --matrix (a full deformation gradient) overrides the\n"
        "engineering-strain flags when both are given.\n";
}

static int runStrain(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::StrainParams params = atomforge::engineeringStrain(
        argDouble(argc, argv, "--exx", 0.0), argDouble(argc, argv, "--eyy", 0.0), argDouble(argc, argv, "--ezz", 0.0),
        argDouble(argc, argv, "--exy", 0.0), argDouble(argc, argv, "--exz", 0.0), argDouble(argc, argv, "--eyz", 0.0));
    if (const auto matrix = findArg(argc, argv, "--matrix"))
    {
        std::istringstream stream(matrix);
        double f[9];
        for (double& value : f) if (!(stream >> value)) throw std::invalid_argument("--matrix needs 9 numbers");
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) params.f[i][j] = f[i * 3 + j];
    }
    const auto result = atomforge::applyStrain(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save strained structure");
    std::cout << result.message << '\n';
    return 0;
}

// -- Lattice: primitive-cell reduction / symmetrization ----------------------

static void printHelpPrimitive()
{
    std::cout << "PRIMITIVE (--build primitive)\n"
        "--input FILE --output FILE [--symprec 1e-3]\n"
        "Reduces a structure to its symmetry-standardized primitive cell via\n"
        "spglib (requires spglib support at build time).\n";
}

static int runPrimitive(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure structure; std::string error;
    if (!loadStructureFromFile(input, structure, error)) throw std::runtime_error(error);
    const double symprec = argDouble(argc, argv, "--symprec", 1e-3);
    const int before = (int)structure.atoms.size();
    if (!reduceToPrimitive(structure, symprec))
        throw std::runtime_error("Could not reduce to a primitive cell (needs spglib support and a valid unit cell).");
    if (!saveStructure(structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save primitive structure");
    std::cout << "Reduced from " << before << " to " << structure.atoms.size() << " atoms.\n";
    return 0;
}

// -- Surfaces: Miller-plane slab builder --------------------------------------

static void printHelpSurface()
{
    std::cout << "SURFACE (--build surface)\n"
        "--input FILE --output FILE [--h 1] [--k 1] [--l 1]\n"
        "[--layers 4] [--vacuum 15] [--nmax 8] [--no-primitive] [--symprec 1e-3]\n"
        "Cleaves a bulk crystal along Miller plane (h k l) into a 2D-periodic,\n"
        "vacuum-padded slab.\n";
}

static int runSurface(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::SurfaceParams params;
    params.h = argInt(argc, argv, "--h", 1);
    params.k = argInt(argc, argv, "--k", 1);
    params.l = argInt(argc, argv, "--l", 1);
    params.layers = argInt(argc, argv, "--layers", 4);
    params.vacuum = argDouble(argc, argv, "--vacuum", 15.0);
    params.nmax = argInt(argc, argv, "--nmax", 8);
    params.primitiveInput = !hasFlag(argc, argv, "--no-primitive");
    params.primitiveSymprec = argDouble(argc, argv, "--symprec", 1e-3);
    const auto result = atomforge::buildSurface(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save surface");
    std::cout << result.message << '\n';
    return 0;
}

// -- Alloys: SQS-style species optimizer --------------------------------------

static void printHelpSQS()
{
    std::cout << "SQS (--build sqs)\n"
        "--input FILE --output FILE --element \"SYMBOL FRACTION\" (repeatable, >=2)\n"
        "[--shells 2] [--shell-tolerance 0.2] [--steps 3000]\n"
        "[--start-temp 1.0] [--end-temp 0.02] [--seed 1]\n"
        "Simulated-annealing species optimizer: reassigns elements on the input's\n"
        "fixed lattice to a target composition while minimizing Warren-Cowley\n"
        "short-range order toward zero across the first N shells.\n";
}

static int runSQS(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::SQSParams params;
    for (const auto& spec : findAllArgs(argc, argv, "--element"))
    {
        std::istringstream stream(spec);
        std::string symbol; double fraction;
        if (!(stream >> symbol >> fraction))
            throw std::invalid_argument("Cannot parse --element value '" + spec + "'. Expected: \"SYMBOL FRACTION\"");
        params.composition[symbol] = fraction;
    }
    params.shells = argInt(argc, argv, "--shells", 2);
    params.shellTolerance = argDouble(argc, argv, "--shell-tolerance", 0.2);
    params.steps = argInt(argc, argv, "--steps", 3000);
    params.startTemperature = argDouble(argc, argv, "--start-temp", 1.0);
    params.endTemperature = argDouble(argc, argv, "--end-temp", 0.02);
    params.seed = (unsigned)argInt(argc, argv, "--seed", 1);
    const auto result = atomforge::buildSQS(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save SQS structure");
    std::cout << result.message << '\n';
    return 0;
}

// -- Nanostructures: nanowire builder -----------------------------------------

static void printHelpNanowire()
{
    std::cout << "NANOWIRE (--build nanowire)\n"
        "--input FILE --output FILE [--axis 2] [--radius 10]\n"
        "[--sides 0] [--vacuum 10] [--axis-repeats 1]\n"
        "Cuts a 1D-periodic wire from a bulk crystal: periodic along cell vector\n"
        "`axis` (0=a, 1=b, 2=c), bounded by a circular (sides<3) or regular-polygon\n"
        "(sides>=3, radius = apothem) cross-section elsewhere.\n";
}

static int runNanowire(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::NanowireParams params;
    params.axis = argInt(argc, argv, "--axis", 2);
    params.radius = argDouble(argc, argv, "--radius", 10.0);
    params.sides = argInt(argc, argv, "--sides", 0);
    params.vacuum = argDouble(argc, argv, "--vacuum", 10.0);
    params.axisRepeats = argInt(argc, argv, "--axis-repeats", 1);
    const auto result = atomforge::buildNanowire(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save nanowire");
    std::cout << result.message << '\n';
    return 0;
}

// -- Nanostructures: core-shell relabeling ------------------------------------

static void printHelpCoreShell()
{
    std::cout << "CORE-SHELL (--build core-shell)\n"
        "--input FILE --output FILE --core-radius R --core-element SYMBOL\n"
        "--shell-element SYMBOL [--center \"x y z\"]\n"
        "Relabels a finite structure's atoms (e.g. a Nanocrystal Builder output)\n"
        "into a core/shell composition split by radius from its centroid, or from\n"
        "an explicit center.\n";
}

static int runCoreShell(int argc, char* argv[])
{
    const auto input = findArg(argc, argv, "--input");
    const auto output = findArg(argc, argv, "--output");
    if (!input || !output) throw std::invalid_argument("--input and --output are required");
    Structure source; std::string error;
    if (!loadStructureFromFile(input, source, error)) throw std::runtime_error(error);
    atomforge::CoreShellParams params;
    params.coreRadius = argDouble(argc, argv, "--core-radius", 5.0);
    if (const auto e = findArg(argc, argv, "--core-element")) params.coreElement = e;
    if (const auto e = findArg(argc, argv, "--shell-element")) params.shellElement = e;
    if (const auto c = findArg(argc, argv, "--center"))
    {
        std::istringstream stream(c);
        if (!(stream >> params.center[0] >> params.center[1] >> params.center[2]))
            throw std::invalid_argument("Cannot parse --center \"x y z\"");
        params.useCentroid = false;
    }
    const auto result = atomforge::applyCoreShell(source, params);
    if (!result.success) throw std::runtime_error(result.message);
    if (!saveStructure(result.structure, output, detectFormat(output)))
        throw std::runtime_error("Failed to save core-shell structure");
    std::cout << result.message << '\n';
    return 0;
}

namespace
{
struct BuildMode
{
    std::string_view name;
    int (*run)(int argc, char* argv[]);
    void (*printHelp)();
};

// Register a mode once for both execution and topic-specific help.
constexpr std::array<BuildMode, 17> kBuildModes{{
    {"bulk", runBulk, printHelpBulk},
    {"gb", runGB, printHelpGB},
    {"poly", runPoly, printHelpPoly},
    {"nano", runNano, printHelpNano},
    {"amorphous", runAmorphous, printHelpAmorphous},
    {"sss", runSSS, printHelpSSS},
    {"dislocation", runDislocation, printHelpDislocation},
    {"custom", runCustom, printHelpCustom},
    {"interface", runInterface, printHelpInterface},
    {"stacking-fault", runStackingFault, printHelpStackingFault},
    {"vacancy", runVacancy, printHelpVacancy},
    {"strain", runStrain, printHelpStrain},
    {"primitive", runPrimitive, printHelpPrimitive},
    {"surface", runSurface, printHelpSurface},
    {"sqs", runSQS, printHelpSQS},
    {"nanowire", runNanowire, printHelpNanowire},
    {"core-shell", runCoreShell, printHelpCoreShell},
}};

const BuildMode* findBuildMode(std::string_view name)
{
    const auto mode = std::find_if(kBuildModes.begin(), kBuildModes.end(),
                                  [name](const BuildMode& entry) { return entry.name == name; });
    return mode == kBuildModes.end() ? nullptr : &*mode;
}

std::string buildModeNames()
{
    std::string names;
    for (const auto& mode : kBuildModes)
    {
        if (!names.empty()) names += " | ";
        names += mode.name;
    }
    return names;
}
} // namespace

// ── Public interface ──────────────────────────────────────────────────────────

bool isCLIMode(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--build")   == 0) return true;
        if (std::strcmp(argv[i], "--analyze") == 0) return true;
        if (std::strcmp(argv[i], "--convert") == 0) return true;
        if (std::strcmp(argv[i], "--render")  == 0) return true;
        if (std::strcmp(argv[i], "--science") == 0) return true;
        if (std::strcmp(argv[i], "--help")    == 0) return true;
        if (std::strcmp(argv[i], "-h")        == 0) return true;
        if (std::strcmp(argv[i], "--version") == 0) return true;
        if (std::strcmp(argv[i], "-v")        == 0) return true;
    }
    return false;
}

int runCLI(int argc, char* argv[])
{
    if (hasFlag(argc,argv,"--analyze") || hasFlag(argc,argv,"--convert"))
        return runAnalysisCLI(argc,argv);
    if (hasFlag(argc, argv, "--render"))
        return runRenderCLI(argc, argv);
    if (hasFlag(argc, argv, "--science"))
        return runScienceCLI(argc, argv);
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
            if (const auto* mode = findBuildMode(topic))
            {
                mode->printHelp();
                return 0;
            }
            std::cerr << "Unknown help topic '" << topic
                      << "'.  Valid topics: " << buildModeNames() << "\n";
            return 1;
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
        if (const auto* entry = findBuildMode(mode))
            return entry->run(argc, argv);

        std::cerr << "Error: unknown build mode '" << mode
                  << "'.  Valid modes: " << buildModeNames() << "\n";
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
