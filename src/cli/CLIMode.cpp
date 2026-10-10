#include "cli/CLIMode.h"
#include "cli/AnalysisCLI.h"
#include "cli/BuildModes.h"
#include "cli/CliArgs.h"
#include "cli/PipelineCLI.h"
#include "cli/RenderCLI.h"
#include "cli/ScienceCLI.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <string>
#include <string_view>

using namespace cli;

// ── Print usage ──────────────────────────────────────────────────────────────

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
"  AtomForge --science-batch BATCH.json --output RESULTS.json [--csv TABLE.csv]\n"
"  AtomForge --pipe \"STEP | STEP\" --input FILE|- --output FILE|-   (see --pipe --help)\n"
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

// ── Bulk builder ─────────────────────────────────────────────────────────────

// ── GB builder ───────────────────────────────────────────────────────────────

// ── Polycrystal builder ───────────────────────────────────────────────────────

// ── Nanocrystal builder ───────────────────────────────────────────────────────

// ── Amorphous builder ─────────────────────────────────────────────────────────

// ── Interface builder ─────────────────────────────────────────────────────────

// ── Custom mesh-fill builder ─────────────────────────────────────────────────

// ── Substitutional solid solution builder ────────────────────────────────────

// -- Point defects: vacancy generator ----------------------------------------

// -- Lattice: homogeneous strain ---------------------------------------------

// -- Lattice: primitive-cell reduction / symmetrization ----------------------

// -- Surfaces: Miller-plane slab builder --------------------------------------

// -- Alloys: SQS-style species optimizer --------------------------------------

// -- Nanostructures: nanowire builder -----------------------------------------

// -- Nanostructures: core-shell relabeling ------------------------------------

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
        if (std::strcmp(argv[i], "--science-batch") == 0) return true;
        if (std::strcmp(argv[i], "--pipe") == 0) return true;
        if (std::strcmp(argv[i], "--pipeline") == 0) return true;
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
    if (hasFlag(argc, argv, "--science") || hasFlag(argc, argv, "--science-batch"))
        return runScienceCLI(argc, argv);
    if (hasFlag(argc, argv, "--pipe") || hasFlag(argc, argv, "--pipeline"))
        return runPipelineCLI(argc, argv);
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
