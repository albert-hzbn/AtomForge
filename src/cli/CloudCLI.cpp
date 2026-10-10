#include "cli/CloudCLI.h"
#include "cli/CliArgs.h"

#include "cloud/AtomCloud.h"
#include "io/StructureLoader.h"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace cli;
namespace cloud = atomforge::cloud;

namespace
{
void printHelp()
{
    std::cout <<
"AtomForge --cloud - out-of-core atom clouds for very large structures\n"
"\n"
"Atom clouds (.afcloud) hold up to billions of atoms in spatial chunks of\n"
"8 bytes per atom; the desktop opens them with its large-data renderer\n"
"(level of detail and streaming), so their size is limited by the disk, not\n"
"by memory. Building streams the input: memory use does not grow with it.\n"
"\n"
"  AtomForge --cloud build --input FILE --output OUT.afcloud\n"
"      LAMMPS text dumps (x y z, xs ys zs, xu yu zu; element or type columns),\n"
"      XYZ and extended XYZ are streamed; other formats are loaded first.\n"
"      --types \"Cu,Ni\"   element names of LAMMPS types 1, 2, ...\n"
"      --frame N        frame of a multi-frame dump or XYZ file (default 0)\n"
"  AtomForge --cloud generate --lattice fcc|bcc|sc|diamond|hcp --a A [--c C]\n"
"            --cells NX NY NZ --element SYMBOL [--second SYMBOL] --output OUT.afcloud\n"
"      A crystal of NX x NY x NZ cells (e.g. fcc 630^3 cells = 1.0 billion atoms).\n"
"  AtomForge --cloud info FILE.afcloud\n"
"\n"
"Options of build and generate:\n"
"  --chunk N          atoms per chunk (default 524288)\n"
"  --temp DIR         folder for intermediate files (default: next to the output;\n"
"                     needs about 16 bytes per atom of free space)\n"
"  --threads N        worker threads (default: all cores, at most 16)\n"
"  --quiet            no progress output\n"
"\n"
"Example (one billion atoms, about 8 GB):\n"
"  AtomForge --cloud generate --lattice fcc --a 3.615 --cells 630 630 630 --element Cu --output cu1b.afcloud\n"
<< std::endl;
}

std::string bytes(double value)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(value >= 1e9 ? 2 : 1);
    if (value >= 1e9) out << value / 1e9 << " GB";
    else out << value / 1e6 << " MB";
    return out.str();
}

void printInfo(const cloud::CloudInfo& info, const std::filesystem::path& path)
{
    std::error_code ignored;
    std::cout << path.u8string() << "\n"
              << "  atoms    " << info.atoms << "\n"
              << "  chunks   " << info.chunks.size() << "\n"
              << "  species ";
    for (const auto& s : info.species) std::cout << " " << s.name;
    std::cout << "\n  box      [" << info.lower[0] << ", " << info.lower[1] << ", " << info.lower[2] << "] - ["
              << info.upper[0] << ", " << info.upper[1] << ", " << info.upper[2] << "] Angstrom\n"
              << "  size     " << bytes(static_cast<double>(std::filesystem::file_size(path, ignored))) << "\n";
    if (!info.source.empty()) std::cout << "  source   " << info.source << "\n";
}
}

int runCloudCLI(int argc, char* argv[])
{
    const char* command = findArg(argc, argv, "--cloud");
    if (!command || hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h")) { printHelp(); return command ? 0 : 1; }
    try {
        cloud::BuildOptions options;
        if (const char* chunk = findArg(argc, argv, "--chunk")) options.chunkAtoms = std::stoull(chunk);
        if (const char* temp = findArg(argc, argv, "--temp")) options.temporary = std::filesystem::u8path(temp);
        if (const char* threads = findArg(argc, argv, "--threads")) options.threads = static_cast<unsigned>(std::stoul(threads));
        const bool quiet = hasFlag(argc, argv, "--quiet");
        const auto start = std::chrono::steady_clock::now();
        int lastPercent = -1;
        const cloud::Progress progress = [&](double fraction, const std::string& message) {
            if (quiet || fraction < 0) return true;
            const int percent = static_cast<int>(fraction * 100);
            if (percent != lastPercent) {
                lastPercent = percent;
                std::cerr << "\r  " << std::setw(3) << percent << "%  " << message << "          " << std::flush;
            }
            return true;
        };
        const auto finished = [&](const cloud::CloudInfo& info, const std::filesystem::path& output) {
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (!quiet) std::cerr << "\r  100%  done in " << std::fixed << std::setprecision(1) << seconds << " s ("
                                  << std::setprecision(1) << static_cast<double>(info.atoms) / std::max(seconds, 1e-9) / 1e6 << " M atoms/s)          \n";
            printInfo(info, output);
        };
        const std::string mode = command;
        if (mode == "info") {
            // The file is the argument after "info", or --input.
            const char* file = findArg(argc, argv, "--input");
            for (int i = 1; !file && i + 1 < argc; ++i)
                if (std::strcmp(argv[i], "info") == 0) file = argv[i + 1];
            if (!file) throw std::invalid_argument("--cloud info FILE.afcloud");
            const auto path = std::filesystem::u8path(file);
            cloud::CloudFile cloudFile(path);
            printInfo(cloudFile.info(), path);
            return 0;
        }
        const char* output = findArg(argc, argv, "--output");
        if (!output) throw std::invalid_argument("--output OUT.afcloud is required");
        const auto out = std::filesystem::u8path(output);
        if (mode == "generate") {
            const char* lattice = findArg(argc, argv, "--lattice");
            const char* a = findArg(argc, argv, "--a");
            const char* element = findArg(argc, argv, "--element");
            if (!lattice || !a || !element) throw std::invalid_argument("--lattice, --a and --element are required");
            std::array<std::uint64_t, 3> cells{1, 1, 1};
            for (int i = 1; i + 3 < argc; ++i)
                if (std::strcmp(argv[i], "--cells") == 0)
                    for (int k = 0; k < 3; ++k) cells[static_cast<std::size_t>(k)] = std::stoull(argv[i + 1 + k]);
            const char* c = findArg(argc, argv, "--c");
            const char* second = findArg(argc, argv, "--second");
            finished(cloud::generateCrystal(out, lattice, std::stod(a), c ? std::stod(c) : 0.0, cells, element, second ? second : "", progress, options), out);
            return 0;
        }
        if (mode == "build") {
            const char* input = findArg(argc, argv, "--input");
            if (!input) throw std::invalid_argument("--input FILE is required");
            const auto in = std::filesystem::u8path(input);
            std::vector<std::string> types;
            if (const char* names = findArg(argc, argv, "--types")) {
                std::stringstream list(names);
                for (std::string name; std::getline(list, name, ',');) types.push_back(name);
            }
            const std::size_t frame = findArg(argc, argv, "--frame") ? std::stoull(findArg(argc, argv, "--frame")) : 0;
            // Streamed formats by content, others through the structure loaders.
            std::ifstream probe(in);
            std::string first;
            std::getline(probe, first);
            probe.close();
            std::string extension = in.extension().string();
            for (char& ch : extension) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
            if (first.rfind("ITEM: TIMESTEP", 0) == 0) finished(cloud::buildFromLammpsDump(in, out, types, frame, progress, options), out);
            else if (extension == ".xyz" || extension == ".extxyz") finished(cloud::buildFromXyz(in, out, frame, progress, options), out);
            else {
                Structure structure;
                std::string error;
                if (!loadStructureFromFile(in.u8string(), structure, error)) throw std::runtime_error("Cannot read " + in.u8string() + ": " + error);
                std::array<double, 3> lower{1e300, 1e300, 1e300}, upper{-1e300, -1e300, -1e300};
                for (const auto& atom : structure.atoms) {
                    const double p[3] = {atom.x, atom.y, atom.z};
                    for (std::size_t k = 0; k < 3; ++k) { lower[k] = std::min(lower[k], p[k]); upper[k] = std::max(upper[k], p[k]); }
                }
                if (structure.atoms.empty()) throw std::runtime_error(in.u8string() + " has no atoms");
                cloud::CloudBuilder builder(out, lower, upper, structure.atoms.size(), options);
                builder.setSource(in.filename().u8string());
                if (structure.hasUnitCell) builder.setCell(structure.cellVectors, structure.cellOffset);
                for (const auto& atom : structure.atoms) builder.add({atom.x, atom.y, atom.z, builder.species(atom.symbol)});
                finished(builder.finish(progress), out);
            }
            return 0;
        }
        throw std::invalid_argument("Unknown --cloud command " + mode + " (build, generate or info)");
    } catch (const std::exception& error) {
        std::cerr << "\nError: " << error.what() << "\n";
        return 1;
    }
}
