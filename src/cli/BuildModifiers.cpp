#include "cli/BuildModifiers.h"
#include "cli/AnalysisCLI.h"
#include "cli/BuildModes.h"

#include "io/StructureLoader.h"
#include "pipeline/Options.h"
#include "pipeline/Pipeline.h"

#include <atomic>
#include <chrono>
#include <deque>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using atomforge::pipeline::Json;
using atomforge::pipeline::PipelineData;

namespace
{
// How a mode receives the pipeline's structure.
enum class Input { None, Required, Optional, LayerA };

struct Mode
{
    const char* id;          // --build mode (or "sculpt" for --analyze sculpt)
    const char* title;
    const char* category;    // "Build" (makes a structure) or "Build & edit" (changes it)
    Input input;
    const char* defaults;    // default option flags
    int (*run)(int, char**);
    void (*help)();
};

int runSculpt(int argc, char** argv) { return runAnalysisCLI(argc, argv); }
void helpSculpt()
{
    std::cout << "CELL SCULPTOR (--analyze sculpt)\n--slabs \"h k l lower upper periodic;...\" --nx 1 --ny 1 --nz 1\n"
                 "Keeps the atoms between pairs of (h k l) planes, optionally replicating the cell first (nx ny nz);\n"
                 "the slabs must close a region (three independent planes), which becomes the new cell;\n"
                 "lengths are Angstrom, periodic bounds are integer start/count.\n";
}

const std::vector<Mode>& modes()
{
    using namespace cli;
    static const std::vector<Mode> table = {
        {"bulk", "Bulk crystal", "Build", Input::None, "--system cubic --spacegroup 225 --a 3.615 --atom \"Cu 0 0 0\"", runBulk, printHelpBulk},
        {"custom", "Custom structure", "Build", Input::Optional, "", runCustom, printHelpCustom},
        {"gb", "Grain boundary", "Build", Input::Optional, "", runGB, printHelpGB},
        {"poly", "Polycrystal", "Build", Input::Optional, "", runPoly, printHelpPoly},
        {"nano", "Nanocrystal", "Build", Input::Optional, "", runNano, printHelpNano},
        {"amorphous", "Amorphous structure", "Build", Input::None, "--element \"Si 64\" --density 2.33", runAmorphous, printHelpAmorphous},
        {"interface", "Heterogeneous interface", "Build", Input::LayerA, "", runInterface, printHelpInterface},
        {"nanowire", "Nanowire", "Build", Input::Required, "--axis 2 --radius 10", runNanowire, printHelpNanowire},
        {"core-shell", "Core-shell particle", "Build & edit", Input::Required, "--core-radius 8 --core-element Cu --shell-element Ag", runCoreShell, printHelpCoreShell},
        {"surface", "Surface slab", "Build & edit", Input::Required, "--h 1 --k 1 --l 1 --layers 4 --vacuum 15", runSurface, printHelpSurface},
        {"vacancy", "Vacancies", "Build & edit", Input::Required, "--count 1 --seed 1", runVacancy, printHelpVacancy},
        {"sss", "Substitutional solid solution", "Build & edit", Input::Required, "--frac \"Cu=0.7,Zn=0.3\" --seed 42", runSSS, printHelpSSS},
        {"sqs", "Special quasirandom structure", "Build & edit", Input::Required, "--element \"Cu 0.5\" --element \"Au 0.5\"", runSQS, printHelpSQS},
        {"strain", "Strain the cell", "Build & edit", Input::Required, "--exx 0.01", runStrain, printHelpStrain},
        {"primitive", "Primitive cell", "Build & edit", Input::Required, "", runPrimitive, printHelpPrimitive},
        {"dislocation", "Insert dislocation", "Build & edit", Input::Optional, "", runDislocation, printHelpDislocation},
        {"stacking-fault", "Stacking fault", "Build & edit", Input::Required, "--plane 1 --layers 9", runStackingFault, printHelpStackingFault},
        {"sculpt", "Cell sculptor", "Build & edit", Input::Required, "--slabs \"1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1\"", runSculpt, helpSculpt},
    };
    return table;
}

bool mentions(const std::vector<std::string>& args, const char* flag)
{
    for (const auto& a : args) if (a == flag) return true;
    return false;
}

std::string captureHelp(void (*help)())
{
    std::ostringstream text;
    auto* previous = std::cout.rdbuf(text.rdbuf());
    help();
    std::cout.rdbuf(previous);
    return text.str();
}

// Runs one mode with the pipeline's structure; the result replaces the data.
void runMode(const Mode& mode, const Json& parameters, PipelineData& data)
{
    static std::atomic<unsigned> counter{0};
    const auto folder = std::filesystem::temp_directory_path() / "AtomForge-pipeline" /
        (std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" + std::to_string(counter++));
    std::filesystem::create_directories(folder);
    const auto cleanup = [&] { std::error_code ignored; std::filesystem::remove_all(folder, ignored); };
    try {
        std::vector<std::string> args = {"AtomForge"};
        const bool sculpt = std::string(mode.id) == "sculpt";
        if (sculpt) { args.push_back("--analyze"); args.push_back("sculpt"); }
        else { args.push_back("--build"); args.push_back(mode.id); }
        const auto options = atomforge::pipeline::options::tokens(parameters.contains("options") ? parameters.at("options").string() : mode.defaults);
        args.insert(args.end(), options.begin(), options.end());
        // The pipeline's structure is the input (CIF keeps the cell and the atom order).
        const bool haveAtoms = !data.structure.atoms.empty();
        const bool useInput = parameters.contains("use_input") ? parameters.at("use_input").boolean() : true;
        if ((mode.input == Input::Required || mode.input == Input::LayerA || (mode.input == Input::Optional && useInput)) && haveAtoms) {
            const char* flag = mode.input == Input::LayerA ? "--layerA" : "--input";
            if (!mentions(args, flag)) {
                const auto input = folder / "input.cif";
                if (!saveStructure(data.structure, input.u8string(), "cif")) throw std::runtime_error("Cannot pass the structure to " + std::string(mode.id));
                args.push_back(flag);
                args.push_back(input.u8string());
            }
        } else if (mode.input == Input::Required) throw std::runtime_error(std::string(mode.title) + " needs a structure from earlier steps");
        const auto output = folder / (sculpt ? "output.vasp" : "output.cif");
        args.push_back("--output");
        args.push_back(output.u8string());
        std::vector<char*> argv;
        for (auto& a : args) argv.push_back(a.data());
        // The mode reports on stdout/stderr; keep that out of the pipeline's own output.
        std::ostringstream out, err;
        auto* oldOut = std::cout.rdbuf(out.rdbuf());
        auto* oldErr = std::cerr.rdbuf(err.rdbuf());
        int status = 1;
        try { status = mode.run(static_cast<int>(argv.size()), argv.data()); }
        catch (...) { std::cout.rdbuf(oldOut); std::cerr.rdbuf(oldErr); throw; }
        std::cout.rdbuf(oldOut);
        std::cerr.rdbuf(oldErr);
        if (status != 0 || !std::filesystem::exists(output)) {
            std::string message = err.str().empty() ? out.str() : err.str();
            while (!message.empty() && std::isspace(static_cast<unsigned char>(message.back()))) message.pop_back();
            if (message.rfind("Error: ", 0) == 0) message.erase(0, 7);
            throw std::runtime_error(message.empty() ? std::string(mode.title) + " failed" : message);
        }
        Structure result;
        std::string error;
        if (!loadStructureFromFile(output.u8string(), result, error)) throw std::runtime_error("Cannot read the result of " + std::string(mode.id) + ": " + error);
        data.structure = std::move(result);
        data.clearSelection();
        data.notes.push_back(std::string(mode.title) + ": " + std::to_string(data.structure.atoms.size()) + " atoms");
        cleanup();
    } catch (...) {
        cleanup();
        throw;
    }
}
}

void registerBuildModifiers()
{
    // Registered strings must outlive the registry.
    static std::deque<std::string> texts;
    static bool done = false;
    if (done) return;
    done = true;
    for (const Mode& mode : modes()) {
        texts.push_back(std::string("build-") + mode.id);
        const char* id = texts.back().c_str();
        std::string help = std::string(mode.input == Input::None ? "Makes a new structure (the pipeline's structure is replaced)."
                                       : mode.input == Input::Optional ? "Uses the pipeline's structure as its input when 'use input' is on."
                                       : "Takes the pipeline's structure as its input.") +
                           " Options are the command-line flags of AtomForge --build " + mode.id + ":\n\n" + captureHelp(mode.help);
        texts.push_back(help);
        const char* helpText = texts.back().c_str();
        std::vector<atomforge::pipeline::ModifierParameter> parameters = {{"options", "Options (command-line flags)", "options", mode.defaults, ""}};
        if (mode.input == Input::Optional) parameters.push_back({"use_input", "Use the pipeline's structure as input", "bool", "true", ""});
        const Mode* captured = &mode;
        atomforge::pipeline::registerModifierType({id, mode.title, mode.category, helpText, parameters,
            [captured](const Json& p, PipelineData& data) { runMode(*captured, p, data); }});
    }
}
