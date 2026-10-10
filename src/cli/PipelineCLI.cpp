#include "cli/PipelineCLI.h"
#include "cli/CliArgs.h"
#include "ui/MenuParity.h"

#include "io/StructureLoader.h"
#include "pipeline/Pipeline.h"
#include "science/ScienceData.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace cli;
using atomforge::pipeline::Pipeline;
using atomforge::pipeline::Json;

namespace
{
void printHelp()
{
    std::cout <<
"AtomForge --pipe - non-destructive structure editing pipeline\n"
"\n"
"Usage:\n"
"  AtomForge --pipe \"STEP | STEP | ...\" --input FILE --output FILE\n"
"  AtomForge --pipeline STEPS.json|STEPS.txt --input FILE --output FILE\n"
"  AtomForge --pipe --list                 List modifiers and their parameters\n"
"\n"
"  --input FILE | -       Structure file, or extended XYZ on stdin\n"
"  --output FILE | -      Output file (format from extension), or extended XYZ on\n"
"                         stdout (with selection and per-atom property columns)\n"
"  --save-pipeline FILE   Also save the pipeline as JSON (open it in the desktop)\n"
"  --report FILE          Write the per-step results as JSON\n"
"  --list-json            The modifier descriptions as JSON\n"
"  --quiet                No per-step report on stderr\n"
"\n"
"Each STEP is a modifier id followed by its parameters, positionally or as\n"
"name=value; vectors take three numbers (or a,b,c), an expression or element\n"
"list takes the rest of its step. Modifiers act on the selection made by the\n"
"selection modifiers above them. The input is never modified.\n"
"\n"
"Examples:\n"
"  AtomForge --pipe \"replicate 4 4 4 | select-random 0.05 | delete-selected\" ^\n"
"            --input cu.cif --output cu_vacancies.xyz\n"
"  AtomForge --pipe \"add-vacuum c 15 | select-expression fz > 0.6 | assign-element Ni\" ^\n"
"            --input slab.vasp --output - | AtomForge --pipe \"wrap\" --input - --output out.cif\n"
<< std::endl;
}

void listModifiers()
{
    std::string category;
    for (const auto& type : atomforge::pipeline::modifierTypes()) {
        if (category != type.category) { category = type.category; std::cout << "\n" << category << "\n"; }
        std::cout << "  " << std::left << std::setw(20) << type.id << type.title << "\n";
        for (const auto& p : type.parameters) {
            std::cout << "      " << std::setw(16) << p.name << std::setw(12) << p.kind << "default " << p.value;
            if (*p.options) std::cout << "  (" << p.options << ")";
            std::cout << "\n";
        }
    }
    std::cout << std::endl;
}

Structure readInput(const std::string& path)
{
    if (path == "-") {
        auto frames = atomforge::science::readExtxyzFrames(std::cin, 1);
        return frames.front().structure;
    }
    Structure structure;
    std::string error;
    if (!loadStructureFromFile(path, structure, error)) throw std::runtime_error("Cannot read " + path + ": " + error);
    return structure;
}

// Extended XYZ with the pipeline's selection and per-atom property.
void writeExtxyz(std::ostream& out, const atomforge::pipeline::PipelineData& data)
{
    const Structure& s = data.structure;
    const bool property = s.atomProperty.size() == s.atoms.size() && !s.atoms.empty();
    out << std::setprecision(12) << s.atoms.size() << "\n";
    if (s.hasUnitCell) {
        out << "Lattice=\"";
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) out << (r || k ? " " : "") << s.cellVectors[r][k];
        out << "\" pbc=\"T T T\" ";
    }
    out << "Properties=species:S:1:pos:R:3:selection:I:1" << (property ? ":property:R:1" : "") << "\n";
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const auto& a = s.atoms[i];
        out << a.symbol << ' ' << a.x << ' ' << a.y << ' ' << a.z << ' ' << static_cast<int>(data.selected[i]);
        if (property) out << ' ' << (std::isfinite(s.atomProperty[i]) ? s.atomProperty[i] : 0.0);
        out << "\n";
    }
}
}

int runPipelineCLI(int argc, char* argv[])
{
    if (hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h")) { printHelp(); return 0; }
    if (hasFlag(argc, argv, "--list")) { listModifiers(); return 0; }
    if (hasFlag(argc, argv, "--list-json")) {
        // Modifier descriptions for tools (the Python package uses this).
        Json types = Json::array();
        for (const auto& type : atomforge::pipeline::modifierTypes()) {
            Json entry = Json::object();
            entry["id"] = type.id;
            entry["title"] = type.title;
            entry["category"] = type.category;
            entry["help"] = type.help;
            entry["menu"] = menuPathOf(type.id);  // the same feature in the desktop menus
            const auto defaults = atomforge::pipeline::makeModifier(type.id).parameters;
            Json parameters = Json::array();
            for (const auto& parameter : type.parameters) {
                Json item = Json::object();
                item["name"] = parameter.name;
                item["label"] = parameter.label;
                item["kind"] = parameter.kind;
                item["default"] = defaults.at(parameter.name);
                item["options"] = parameter.options;
                parameters.push(item);
            }
            entry["parameters"] = parameters;
            types.push(entry);
        }
        std::cout << types.dump() << std::endl;
        return 0;
    }
    try {
        Pipeline pipeline;
        if (const char* text = findArg(argc, argv, "--pipe"); text && std::strcmp(text, "--list") != 0)
            pipeline = Pipeline::parse(text);
        else if (const char* file = findArg(argc, argv, "--pipeline")) {
            std::ifstream in(std::filesystem::u8path(file));
            if (!in) throw std::runtime_error(std::string("Cannot read ") + file);
            std::stringstream content;
            content << in.rdbuf();
            const std::string body = content.str();
            const auto start = body.find_first_not_of(" \t\r\n");
            pipeline = start != std::string::npos && body[start] == '{' ? Pipeline::fromJson(atomforge::science::Json::parse(body)) : Pipeline::parse(body);
        } else throw std::invalid_argument("--pipe \"STEPS\" or --pipeline FILE is required (see --pipe --help)");
        const char* input = findArg(argc, argv, "--input");
        const char* output = findArg(argc, argv, "--output");
        if (!input || !output) throw std::invalid_argument("--input FILE|- and --output FILE|- are required");
        const bool quiet = hasFlag(argc, argv, "--quiet");
        if (const char* save = findArg(argc, argv, "--save-pipeline")) {
            std::ofstream out(std::filesystem::u8path(save));
            out << pipeline.toJson().dump(2) << "\n";
            if (!out) throw std::runtime_error(std::string("Cannot write ") + save);
        }
        const Structure source = readInput(input);
        std::vector<atomforge::pipeline::StageResult> stages;
        const auto data = pipeline.evaluate(source, &stages);
        if (const char* report = findArg(argc, argv, "--report")) {
            // Per-step results as JSON (atoms, selection, time, error, notes).
            Json list = Json::array();
            for (std::size_t i = 0; i < stages.size(); ++i) {
                Json stage = Json::object();
                stage["type"] = pipeline.modifiers[i].type;
                stage["atoms"] = stages[i].atoms;
                stage["selected"] = stages[i].selected;
                stage["milliseconds"] = stages[i].milliseconds;
                stage["skipped"] = stages[i].skipped;
                stage["error"] = stages[i].error;
                Json notes = Json::array();
                for (const auto& note : stages[i].notes) notes.push(note);
                stage["notes"] = notes;
                list.push(stage);
            }
            Json document = Json::object();
            document["stages"] = list;
            document["property_name"] = data.structure.atomPropertyName;
            std::ofstream out(std::filesystem::u8path(report));
            out << document.dump(2) << "\n";
        }
        // The report goes to stderr so that stdout can carry the structure.
        for (std::size_t i = 0; i < stages.size(); ++i) {
            const auto& stage = stages[i];
            if (!stage.error.empty()) {
                std::cerr << "Error: step " << i + 1 << " (" << pipeline.modifiers[i].type << "): " << stage.error << "\n";
                return 1;
            }
            if (quiet) continue;
            std::cerr << "  " << i + 1 << ' ' << std::left << std::setw(20) << pipeline.modifiers[i].type;
            if (stage.skipped) std::cerr << "skipped (disabled)\n";
            else std::cerr << stage.atoms << " atoms, " << stage.selected << " selected (" << std::fixed << std::setprecision(1) << stage.milliseconds << " ms)\n";
        }
        if (std::strcmp(output, "-") == 0) {
            writeExtxyz(std::cout, data);
            std::cout.flush();
        } else {
            if (data.structure.atoms.empty()) throw std::runtime_error("The pipeline removed every atom; nothing to write");
            if (!saveStructure(data.structure, output, detectFormat(output))) throw std::runtime_error(std::string("Cannot write ") + output);
            if (!quiet) std::cerr << "Saved " << data.structure.atoms.size() << " atoms to " << output << "\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
