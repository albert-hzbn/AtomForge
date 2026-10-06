// Crystal builders and editing tools: bulk, custom, solid solution, primitive cell, surface, SQS, vacancy and strain.
#include "cli/BuildModes.h"
#include "cli/CliArgs.h"

#include "algorithms/BulkCrystalBuilder.h"
#include "algorithms/CSLComputation.h"
#include "algorithms/MeshLoader.h"
#include "algorithms/NanoCrystalBuilder.h"
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
void printHelpBulk()
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

int runBulk(int argc, char* argv[])
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

void printHelpCustom()
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

int runCustom(int argc, char* argv[])
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

void printHelpSSS()
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

int runSSS(int argc, char* argv[])
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

void printHelpPrimitive()
{
    std::cout << "PRIMITIVE (--build primitive)\n"
        "--input FILE --output FILE [--symprec 1e-3]\n"
        "Reduces a structure to its symmetry-standardized primitive cell via\n"
        "spglib (requires spglib support at build time).\n";
}

int runPrimitive(int argc, char* argv[])
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

void printHelpSurface()
{
    std::cout << "SURFACE (--build surface)\n"
        "--input FILE --output FILE [--h 1] [--k 1] [--l 1]\n"
        "[--layers 4] [--vacuum 15] [--nmax 8] [--no-primitive] [--symprec 1e-3]\n"
        "Cleaves a bulk crystal along Miller plane (h k l) into a 2D-periodic,\n"
        "vacuum-padded slab.\n";
}

int runSurface(int argc, char* argv[])
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

void printHelpSQS()
{
    std::cout << "SQS (--build sqs)\n"
        "--input FILE --output FILE --element \"SYMBOL FRACTION\" (repeatable, >=2)\n"
        "[--shells 2] [--shell-tolerance 0.2] [--steps 3000]\n"
        "[--start-temp 1.0] [--end-temp 0.02] [--seed 1]\n"
        "Simulated-annealing species optimizer: reassigns elements on the input's\n"
        "fixed lattice to a target composition while minimizing Warren-Cowley\n"
        "short-range order toward zero across the first N shells.\n";
}

int runSQS(int argc, char* argv[])
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

void printHelpVacancy()
{
    std::cout << "VACANCY (--build vacancy)\n"
        "--input FILE --output FILE [--element SYMBOL] [--percent P | --count N]\n"
        "[--min-separation D] [--seed S]\n"
        "Randomly removes atoms (optionally restricted to one element) to create\n"
        "point-defect vacancies at a target percentage or exact count.\n";
}

int runVacancy(int argc, char* argv[])
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

void printHelpStrain()
{
    std::cout << "STRAIN (--build strain)\n"
        "--input FILE --output FILE\n"
        "[--exx V] [--eyy V] [--ezz V] [--exy V] [--exz V] [--eyz V]\n"
        "[--matrix \"f11 f12 f13 f21 f22 f23 f31 f32 f33\"]\n"
        "Applies a homogeneous deformation to the cell, holding fractional\n"
        "coordinates fixed. --matrix (a full deformation gradient) overrides the\n"
        "engineering-strain flags when both are given.\n";
}

int runStrain(int argc, char* argv[])
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
}
