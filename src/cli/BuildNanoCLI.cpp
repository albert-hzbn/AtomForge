// Microstructure and nanostructure builders: polycrystal, nanocrystal, amorphous, nanowire and core-shell.
#include "cli/BuildModes.h"
#include "cli/CliArgs.h"

#include "algorithms/AmorphousBuilder.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "algorithms/NanostructureTools.h"
#include "algorithms/PolyCrystalBuilder.h"
#include "io/StructureLoader.h"
#include "util/ElementData.h"
#include "util/PathUtils.h"

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

namespace cli
{
void printHelpPoly()
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
"  --grain-euler \"N phi1 Phi phi2\"\n"
"                               Euler angles (deg) for grain number N (from 1);\n"
"                               repeat per grain. Grains not given get random\n"
"                               orientations.\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build poly --input cu_fcc.cif ^\n"
"            --sizex 100 --sizey 100 --sizez 100 ^\n"
"            --grains 12 --seed 7 --output cu_poly.cif\n"
<< std::endl;
}

int runPoly(int argc, char* argv[])
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

    // Euler angles of chosen grains: each value is "N phi1 Phi phi2" (N from 1)
    auto grainEulerStrs = findAllArgs(argc, argv, "--grain-euler");
    for (const auto& text : grainEulerStrs)
    {
        std::istringstream iss(text);
        int grain = 0;
        GrainOrientation go;
        if (!(iss >> grain >> go.phi1 >> go.Phi >> go.phi2) || grain < 1)
        {
            std::cerr << "Error: cannot parse --grain-euler value '" << text
                      << "'.  Expected: \"N phi1 Phi phi2\" with N >= 1\n";
            return 1;
        }
        go.grainIndex = grain - 1;
        params.specifiedOrientations.push_back(go);
    }
    if (!grainEulerStrs.empty())
        params.orientationMode = GrainOrientationMode::PartialSpecified;

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

void printHelpNano()
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
"  --no-cell                   No output cell (default: a rectangular cell\n"
"                               around the particle plus --vacuum)\n"
"  --center \"x y z\"            Carving center (Ang)  (default: centre of the atoms)\n"
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

int runNano(int argc, char* argv[])
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
    params.setOutputCell = !hasFlag(argc, argv, "--no-cell");
    params.autoCenterFromAtoms = true;
    if (const char* center = findArg(argc, argv, "--center"))
    {
        glm::vec3 c;
        if (!parseVec3(center, c))
        {
            std::cerr << "Error: cannot parse --center value '" << center
                      << "'.  Expected: \"x y z\"\n";
            return 1;
        }
        params.autoCenterFromAtoms = false;
        params.cx = c.x;
        params.cy = c.y;
        params.cz = c.z;
    }

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

void printHelpAmorphous()
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
"  --no-periodic               Non-periodic box (no unit cell attached)\n"
"  --output <file>             Output file (format from extension)\n"
"\n"
"Example:\n"
"  AtomForge --build amorphous --element \"Si 80\" --element \"O 160\" ^\n"
"            --density 2.2 --seed 1 --output sio2.xyz\n"
<< std::endl;
}

int runAmorphous(int argc, char* argv[])
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

    params.periodic = !hasFlag(argc, argv, "--no-periodic");

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

void printHelpNanowire()
{
    std::cout << "NANOWIRE (--build nanowire)\n"
        "--input FILE --output FILE [--axis 2] [--radius 10]\n"
        "[--sides 0] [--vacuum 10] [--axis-repeats 1]\n"
        "Cuts a 1D-periodic wire from a bulk crystal: periodic along cell vector\n"
        "`axis` (0=a, 1=b, 2=c), bounded by a circular (sides<3) or regular-polygon\n"
        "(sides>=3, radius = apothem) cross-section elsewhere.\n";
}

int runNanowire(int argc, char* argv[])
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

void printHelpCoreShell()
{
    std::cout << "CORE-SHELL (--build core-shell)\n"
        "--input FILE --output FILE --core-radius R --core-element SYMBOL\n"
        "--shell-element SYMBOL [--center \"x y z\"]\n"
        "Relabels a finite structure's atoms (e.g. a Nanocrystal Builder output)\n"
        "into a core/shell composition split by radius from its centroid, or from\n"
        "an explicit center.\n";
}

int runCoreShell(int argc, char* argv[])
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
}
