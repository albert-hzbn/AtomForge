// Interface and defect builders: dislocation, grain boundary, interface and stacking fault.
#include "cli/BuildModes.h"
#include "cli/CliArgs.h"

#include "algorithms/CSLComputation.h"
#include "algorithms/DislocationBuilder.h"
#include "algorithms/InterfaceBuilder.h"
#include "algorithms/InterfaceMatching.h"
#include "algorithms/StackingFaultBuilder.h"
#include "io/StructureLoader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace cli
{
void printHelpDislocation()
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
"  --line-center               Line point at the structure's centre (plus --line-offset)\n"
"                              instead of the fractional point\n"
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
"  AtomForge --build dislocation --input base.xyz --character edge ^\n"
"            --manual-vectors --line \"0 0 1\" --burgers \"1 0 0\" ^\n"
"            --shape cylinder --cyl-radius 12 --output edge.xyz\n"
"  AtomForge --build dislocation --input cu_block.xyz --character screw ^\n"
"            --anisotropic --elastic-c11 168.4 --elastic-c12 121.4 ^\n"
"            --elastic-c44 75.4 --shape cylinder --cyl-radius 12 --output cu_screw.xyz\n"
<< std::endl;
}

int runDislocation(int argc, char* argv[])
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
    if (hasFlag(argc, argv, "--line-center"))
        params.useFractionalLinePoint = false;
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

void printHelpGB()
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

int runGB(int argc, char* argv[])
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
    // An explicit sigma names a coincidence site lattice: an odd number >= 3.
    if (findArg(argc, argv, "--sigma") && (sigmaTarget < 3 || sigmaTarget % 2 == 0))
    {
        std::cerr << "Error: --sigma must be an odd number of at least 3 (omit it to search up to --sigmamax)\n";
        return 1;
    }
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

void printHelpInterface()
{
    std::cout <<
"------------------------------------------------------------------\n"
"HETEROGENEOUS INTERFACE  (--build interface)\n"
"------------------------------------------------------------------\n"
"Matches the surface lattices of two crystals (Stradi et al. 2017, the\n"
"method of QuantumATK's InterfaceBuilder): B is rotated about the interface\n"
"normal, supercells of the two surfaces that coincide within the strain\n"
"limit are listed (fewest atoms first), and one is built with B strained onto\n"
"A (or A onto B, or both onto their average).\n"
"\n"
"  --layerA <file>             Crystal A (bottom), with a unit cell\n"
"  --layerB <file>             Crystal B (top), with a unit cell\n"
"  --hklA \"h k l\"              Plane of A at the interface (default: the file's\n"
"                               a-b plane, for a structure that is already a slab)\n"
"  --hklB \"h k l\"              Plane of B at the interface (default: as for A)\n"
"  --layersA <N>               Surface layers of A  (default: 4)\n"
"  --layersB <N>               Surface layers of B  (default: 4)\n"
"\n"
"  Matching:\n"
"  --max-strain <f>            Largest mean |strain|, (|e11|+|e22|+|e12|)/3  (default: 0.03)\n"
"  --max-area <A2>             Largest interface cell area  (default: 150)\n"
"  --repeats <N>               Largest repetition of B's surface vectors  (default: 6)\n"
"  --nmax <N>                  Older name of --repeats\n"
"  --min-angle <deg>           Smallest angle between the cell vectors  (default: 20)\n"
"  --angle-min <deg>           Rotations of B searched, from  (default: 0)\n"
"  --angle-max <deg>           ... to  (default: 180)\n"
"  --angle-step <deg>          ... in steps of  (default: 1)\n"
"  --list                      Print every match\n"
"  --pick <N>                  Build match N of the list (default: the suggested\n"
"                               match: the Pareto-optimal one with the lowest\n"
"                               strain (%) + area / 25 A^2)\n"
"\n"
"  Assembly:\n"
"  --strain-on <B|A|both>      Side strained to fit the other  (default: B)\n"
"  --gap <Ang>                 Distance between the facing layers  (default: 2.0)\n"
"  --vacuum <Ang>              Vacuum along the normal; 0 makes the cell periodic,\n"
"                               with two interfaces each at the gap  (default: 0)\n"
"  --shift \"fa fb\"             In-plane shift of B, fractions of the cell  (default: 0 0)\n"
"  --repx <N>  --repy <N>      Repeats of the interface cell in plane  (default: 1)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build interface --layerA cu.cif --hklA \"1 1 1\" ^\n"
"            --layerB ni.cif --hklB \"1 1 1\" --layersA 4 --layersB 4 ^\n"
"            --max-strain 0.05 --output cu_ni_interface.cif\n"
<< std::endl;
}

int runInterface(int argc, char* argv[])
{
    namespace iface = atomforge::interfaces;
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

    iface::Side sideA, sideB;
    std::string err;
    if (!loadStructureFromFile(pathA, sideA.crystal, err))
    {
        std::cerr << "Error loading layer A: " << err << "\n";
        return 1;
    }
    if (!loadStructureFromFile(pathB, sideB.crystal, err))
    {
        std::cerr << "Error loading layer B: " << err << "\n";
        return 1;
    }
    const auto plane = [&](const char* flag, iface::Side& side) {
        const char* text = findArg(argc, argv, flag);
        side.cleave = text != nullptr;
        if (!text) return;
        glm::vec3 hkl;
        if (!parseVec3(text, hkl) || hkl != glm::round(hkl))
            throw std::invalid_argument(std::string(flag) + " needs three whole numbers, e.g. \"1 1 1\"");
        side.h = static_cast<int>(hkl.x);
        side.k = static_cast<int>(hkl.y);
        side.l = static_cast<int>(hkl.z);
    };
    plane("--hklA", sideA);
    plane("--hklB", sideB);

    const int layersA = argInt(argc, argv, "--layersA", 4);
    const int layersB = argInt(argc, argv, "--layersB", 4);
    sideA.layers = layersA;
    sideB.layers = layersB;
    iface::SearchParams search;
    search.maxStrain = argDouble(argc, argv, "--max-strain", 0.03);
    search.maxArea = argDouble(argc, argv, "--max-area", 150.0);
    search.maxRepeats = argInt(argc, argv, "--repeats", argInt(argc, argv, "--nmax", 6));
    search.minCellAngle = argDouble(argc, argv, "--min-angle", 20.0);
    search.angleMin = argDouble(argc, argv, "--angle-min", 0.0);
    search.angleMax = argDouble(argc, argv, "--angle-max", 180.0);
    search.angleStep = argDouble(argc, argv, "--angle-step", 1.0);
    iface::Assembly assembly;
    assembly.gap = argDouble(argc, argv, "--gap", 2.0);
    assembly.vacuum = argDouble(argc, argv, "--vacuum", 0.0);
    const std::string strainOn = findArg(argc, argv, "--strain-on") ? findArg(argc, argv, "--strain-on") : "B";
    if (strainOn == "B" || strainOn == "b") assembly.strainOn = iface::StrainOn::B;
    else if (strainOn == "A" || strainOn == "a") assembly.strainOn = iface::StrainOn::A;
    else if (strainOn == "both") assembly.strainOn = iface::StrainOn::Both;
    else throw std::invalid_argument("--strain-on must be B, A or both");
    if (const char* text = findArg(argc, argv, "--shift"))
    {
        std::istringstream in(text);
        if (!(in >> assembly.shift[0] >> assembly.shift[1]))
            throw std::invalid_argument("--shift needs two fractions, e.g. \"0.5 0\"");
    }
    const int repx = argInt(argc, argv, "--repx", 1);
    const int repy = argInt(argc, argv, "--repy", 1);
    if (layersA < 1 || layersB < 1 || layersA > 1000 || layersB > 1000 || search.maxStrain <= 0 || search.maxArea <= 0 ||
        search.maxRepeats < 1 || search.maxRepeats > 30 || search.angleStep <= 0 || search.angleMax < search.angleMin ||
        assembly.gap < 0 || assembly.vacuum < 0 || repx < 1 || repy < 1 || repx > 100 || repy > 100)
        throw std::invalid_argument("Invalid interface layers, search limits, gap, vacuum or repeats");

    Structure slabA, slabB;
    if (!iface::prepareSlab(sideA, slabA, err))
    {
        std::cerr << "Error: layer A " << err << "\n";
        return 1;
    }
    if (!iface::prepareSlab(sideB, slabB, err))
    {
        std::cerr << "Error: layer B " << err << "\n";
        return 1;
    }
    const auto matches = iface::findMatches(slabA, slabB, search);
    if (matches.empty())
    {
        std::cerr << "Error: no matches within " << search.maxStrain * 100 << "% strain and " << search.maxArea
                  << " A^2; raise --max-strain, --max-area or --repeats\n";
        return 1;
    }
    const int suggested = iface::suggestedMatch(matches);
    const int pick = argInt(argc, argv, "--pick", suggested);
    if (pick < 0 || pick >= static_cast<int>(matches.size()))
    {
        std::cerr << "Error: --pick " << pick << " is out of range (" << matches.size() << " matches)\n";
        return 1;
    }
    const bool list = hasFlag(argc, argv, "--list");
    std::cout << matches.size() << " matches (fewest atoms first; * Pareto-optimal, > built):\n"
              << "   #   atoms  strain%   area A2  angle  theta   cells A/B\n";
    for (int i = 0; i < static_cast<int>(matches.size()); ++i)
    {
        if (!list && i >= 12 && i != pick) continue;
        const auto& m = matches[static_cast<std::size_t>(i)];
        char line[160];
        std::snprintf(line, sizeof(line), "%c%c%3d %7d %8.3f %9.2f %6.1f %6.1f   %d/%d\n", i == pick ? '>' : ' ', m.pareto ? '*' : ' ',
                      i, m.atoms, m.meanStrain * 100, m.area, m.cellAngle, m.thetaDeg, m.cellsA, m.cellsB);
        std::cout << line;
    }
    if (!list && matches.size() > 12) std::cout << "  ... (--list prints all)\n";

    Structure result = iface::buildInterface(slabA, slabB, matches[static_cast<std::size_t>(pick)], assembly);
    if (repx > 1 || repy > 1) result = repeatInterfaceXY(result, repx, repy);
    std::cout << "Built interface: " << result.atoms.size() << " atoms\n";
    if (!saveStructure(result, outPath, detectFormat(outPath)))
    {
        std::cerr << "Error: failed to save structure to '" << outPath << "'\n";
        return 1;
    }
    std::cout << "Saved to: " << outPath << "\n";
    return 0;
}

void printHelpStackingFault()
{
    std::cout << "STACKING FAULT (--build stacking-fault)\n"
        "--input FILE --output FILE [--plane 0..6] [--layers 9] [--interval 0.1]\n"
        "[--maximum 2] [--frame 0] [--orthogonal] [--sequence DIRECTORY]\n"
        "Planes: 0 auto, 1 FCC111, 2 HCP basal, 3 HCP prismatic,\n"
        "4 HCP pyramidal, 5 BCC110, 6 BCC112. Sequence writes numbered VASP files.\n";
}

int runStackingFault(int argc, char* argv[])
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
}
