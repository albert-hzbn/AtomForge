#include "cli/ScienceCLI.h"
#include "io/StructureLoader.h"
#include "science/ScienceCatalog.h"
#include "science/Batch.h"
#include "science/ScienceTools.h"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
using atomforge::science::Json;

void printScienceHelp()
{
    std::cout <<
"AtomForge scientific tools (native, no Python required)\n"
"\n"
"Usage:\n"
"  AtomForge --science <tool> --input REQUEST.json --output RESULT.json\n"
"            [--report REPORT.txt] [--structures FRAMES.extxyz] [--files DIR] [--overwrite]\n"
"  AtomForge --science --catalog      Print tools, parameters and defaults as JSON\n"
"  AtomForge --science-batch BATCH.json --output RESULTS.json [--csv TABLE.csv]\n"
"            BATCH.json: {\"tool\": ..., \"base\": {...}, \"sweep\": {\"param\": [values]},\n"
"                         \"files\": {\"param\": \"data/*.xyz\"}, \"collect\": [\"result.key\"]}\n"
"            Runs the Cartesian product of sweeps and file matches; exit 2 if a run failed.\n"
"\n"
"REQUEST.json is an object of tool parameters. Data values are inline JSON\n"
"arrays or {\"file\": PATH, \"field\": NAME, \"column\": N} references to JSON,\n"
"NPY, CSV/TXT/DAT, EIGENVAL, XYZ/extXYZ, POSCAR/XDATCAR, LAMMPS dump or other\n"
"structure files; paths are relative to REQUEST.json.\n"
"Simulation tools take an interatomic potential such as {\"potential\": \"EMT\"}\n"
"or {\"potential\": \"LennardJones\", \"epsilon\": 0.0104, \"sigma\": 3.4, \"cutoff\": 8.5}\n"
"or {\"potential\": \"EAM\", \"file\": \"Cu.eam.alloy\"} (setfl, eam/fs or funcfl tables).\n"
"\n"
"Tools:\n";
    for (const auto& tool : scienceToolCatalog())
        std::cout << "  " << tool.id << std::string(26 - std::min<std::size_t>(25, std::strlen(tool.id)), ' ') << tool.title << '\n';
    std::cout << std::endl;
}

Json catalogJson()
{
    Json tools = Json::array();
    int priority = 1;
    for (const auto& tool : scienceToolCatalog()) {
        Json entry = Json::object();
        entry["priority"] = priority++;
        entry["id"] = tool.id;
        entry["name"] = tool.title;
        entry["category"] = tool.category;
        entry["description"] = tool.help;
        Json parameters = Json::array();
        for (const auto& parameter : tool.parameters) {
            Json item = Json::object();
            item["name"] = parameter.name;
            item["label"] = parameter.label;
            item["kind"] = parameter.kind;
            item["required"] = parameter.required;
            item["default"] = parameter.value[0] ? Json::parse(parameter.value) : Json();
            parameters.push(item);
        }
        entry["parameters"] = parameters;
        tools.push(entry);
    }
    return tools;
}

void writeText(const std::filesystem::path& path, const std::string& text)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    // Write beside the destination first so a failure never truncates an existing file.
    const auto temporary = path.parent_path() / (".atomforge-" + path.filename().u8string() + ".tmp");
    {
        std::ofstream output(temporary, std::ios::binary);
        output << text;
        if (!output) throw std::runtime_error("Cannot write " + path.u8string());
    }
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    std::filesystem::rename(temporary, path);
}
}

int runScienceCLI(int argc, char* argv[])
{
    std::map<std::string, std::string> args;
    std::string tool, batchFile;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--science") {
            if (i + 1 < argc && argv[i + 1][0] != '-') tool = argv[++i];
        } else if (key == "--science-batch") {
            if (i + 1 >= argc) { std::cerr << "Error: --science-batch requires a batch file\n"; return 1; }
            batchFile = argv[++i];
        } else if (key == "--csv") {
            if (i + 1 >= argc) { std::cerr << "Error: --csv requires a value\n"; return 1; }
            args[key] = argv[++i];
        } else if (key == "--catalog" || key == "--overwrite" || key == "--help" || key == "-h") args[key] = "";
        else if (key == "--input" || key == "--output" || key == "--report" || key == "--structures" || key == "--files") {
            if (i + 1 >= argc) { std::cerr << "Error: " << key << " requires a value\n"; return 1; }
            args[key] = argv[++i];
        } else { std::cerr << "Error: unknown science option " << key << "\n"; return 1; }
    }
    try {
        if (args.count("--catalog")) { std::cout << catalogJson().dump(2) << std::endl; return 0; }
        if (!batchFile.empty()) {
            if (!args.count("--output")) throw std::invalid_argument("--output is required for a batch");
            const auto input = std::filesystem::u8path(batchFile);
            std::ifstream stream(input, std::ios::binary);
            if (!stream) throw std::runtime_error("Cannot open " + input.u8string());
            std::ostringstream text;
            text << stream.rdbuf();
            const auto reader = [](const std::filesystem::path& path) {
                Structure structure;
                std::string error;
                if (!loadStructureFromFile(path.u8string(), structure, error)) throw std::runtime_error(error);
                return structure;
            };
            const auto output = std::filesystem::u8path(args.at("--output"));
            if (std::filesystem::exists(output) && !args.count("--overwrite")) throw std::runtime_error(output.u8string() + " exists; use --overwrite");
            const Json result = atomforge::science::runBatch(Json::parse(text.str()), std::filesystem::absolute(input).parent_path(), reader);
            writeText(output, result.dump(2) + "\n");
            if (args.count("--csv")) writeText(std::filesystem::u8path(args.at("--csv")), atomforge::science::batchCsv(result));
            std::cout << "Batch: " << result.at("run_count").number() << " runs, " << result.at("failures").number() << " failed; saved " << output.u8string() << std::endl;
            return result.at("failures").number() > 0 ? 2 : 0;
        }
        if (args.count("--help") || args.count("-h") || tool.empty()) { printScienceHelp(); return tool.empty() && !args.count("--help") && !args.count("-h") ? 1 : 0; }
        if (!args.count("--input") || !args.count("--output")) throw std::invalid_argument("--input and --output are required");
        const auto input = std::filesystem::u8path(args.at("--input"));
        const auto output = std::filesystem::u8path(args.at("--output"));
        const bool overwrite = args.count("--overwrite") > 0;
        std::vector<std::filesystem::path> destinations = {output};
        if (args.count("--report")) destinations.push_back(std::filesystem::u8path(args.at("--report")));
        if (args.count("--structures")) destinations.push_back(std::filesystem::u8path(args.at("--structures")));
        for (std::size_t i = 0; i < destinations.size(); ++i) {
            if (std::filesystem::exists(destinations[i]) && !overwrite)
                throw std::runtime_error(destinations[i].u8string() + " exists; choose another file or use --overwrite");
            for (std::size_t j = 0; j < i; ++j)
                if (std::filesystem::absolute(destinations[i]).lexically_normal() == std::filesystem::absolute(destinations[j]).lexically_normal())
                    throw std::runtime_error("Result, report and structure files must have distinct paths");
            if (std::filesystem::absolute(destinations[i]).lexically_normal() == std::filesystem::absolute(input).lexically_normal())
                throw std::runtime_error("Input and output files must differ");
        }
        std::ifstream stream(input, std::ios::binary);
        if (!stream) throw std::runtime_error("Cannot open " + input.u8string());
        std::ostringstream text;
        text << stream.rdbuf();
        const Json request = Json::parse(text.str());
        const auto reader = [](const std::filesystem::path& path) {
            Structure structure;
            std::string error;
            if (!loadStructureFromFile(path.u8string(), structure, error)) throw std::runtime_error(error);
            return structure;
        };
        const auto result = atomforge::science::runTool(tool, request, std::filesystem::absolute(input).parent_path(), reader);
        if (args.count("--structures") && !result.frames.empty())
            atomforge::science::writeExtxyz(std::filesystem::u8path(args.at("--structures")), result.frames, result.velocities, result.times);
        writeText(output, atomforge::science::resultDocument(tool, request, result).dump(2) + "\n");
        if (args.count("--report")) writeText(std::filesystem::u8path(args.at("--report")), atomforge::science::resultReport(tool, result.result));
        // Generated input files (e.g. KPOINTS, POSCAR) go to the --files directory.
        if (args.count("--files"))
            if (const Json* files = result.result.find("files"); files && files->isObject()) {
                const auto directory = std::filesystem::u8path(args.at("--files"));
                for (const auto& [name, content] : files->members()) {
                    const auto target = directory / std::filesystem::u8path(name);
                    if (std::filesystem::exists(target) && !overwrite) throw std::runtime_error(target.u8string() + " exists; use --overwrite");
                    writeText(target, content.string());
                }
            }
        std::cout << "Saved " << output.u8string() << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AtomForge science: " << error.what() << std::endl;
        return 1;
    }
}
