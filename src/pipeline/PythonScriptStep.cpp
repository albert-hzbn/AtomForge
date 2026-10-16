// The python-script step: runs a user's Python script on the structure, as
// OVITO's Python script modifier.
#include "pipeline/ModifierSupport.h"
#include "pipeline/PipelineSteps.h"
#include "science/ScienceData.h"
#include "science/ScienceTools.h"
#include "util/ElementData.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace atomforge::pipeline
{
namespace
{
using namespace support;

// Runs a Python script on the structure (OVITO's Python script modifier). The
// script is called as  python SCRIPT INPUT.json OUTPUT.json  and writes the
// modified data: {"structure": {"symbols", "positions", "cell"}, "pbc",
// "selected", "property", "property_name"}; atomforge.pipeline.script_step
// does the reading and writing.
void pythonScript(const Json& p, PipelineData& data)
{
    const std::string script = text(p, "python-script", "script");
    if (script.empty() || !std::filesystem::is_regular_file(std::filesystem::u8path(script))) throw std::runtime_error("Choose an existing Python script file");
    const auto directory = std::filesystem::temp_directory_path() / ("AtomForge-python-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    const auto input = directory / "input.json", output = directory / "output.json", log = directory / "log.txt";
    Structure& s = data.structure;
    Json document = Json::object();
    document["structure"] = science::structureJson(s);
    document["pbc"] = Json::array({science::periodicAxes(s)[0], science::periodicAxes(s)[1], science::periodicAxes(s)[2]});
    Json selected = Json::array();
    for (char c : data.selected) selected.push(c ? 1.0 : 0.0);
    document["selected"] = selected;
    if (s.atomProperty.size() == s.atoms.size() && !s.atoms.empty()) {
        document["property"] = science::toJson(s.atomProperty);
        document["property_name"] = s.atomPropertyName;
    }
    document["arguments"] = text(p, "python-script", "arguments");
    { std::ofstream(input) << document.dump(); }
    std::string interpreter = text(p, "python-script", "python");
    if (interpreter.empty()) {
        const char* configured = std::getenv("ATOMFORGE_PYTHON");
#ifdef _WIN32
        interpreter = configured && *configured ? configured : "python";
#else
        interpreter = configured && *configured ? configured : "python3";
#endif
    }
    const auto quoted = [](const std::string& word) { return "\"" + word + "\""; };
    std::string command = quoted(interpreter) + " " + quoted(script) + " " + quoted(input.u8string()) + " " + quoted(output.u8string()) + " > " +
                          quoted(log.u8string()) + " 2>&1";
#ifdef _WIN32
    command = "\"" + command + "\"";  // cmd /c strips one pair of quotes
#endif
    const int status = std::system(command.c_str());
    std::string messages;
    {
        std::ifstream in(log);
        std::stringstream buffer;
        buffer << in.rdbuf();
        messages = buffer.str();
    }
    std::error_code ignored;  // temporary files; a leftover is harmless
    if (status != 0 || !std::filesystem::exists(output)) {
        std::filesystem::remove_all(directory, ignored);
        const std::string tail = messages.size() > 1500 ? messages.substr(messages.size() - 1500) : messages;
        throw std::runtime_error("The Python script failed (" + interpreter + "):\n" + tail);
    }
    std::string text;
    {
        std::ifstream in(output);
        std::stringstream buffer;
        buffer << in.rdbuf();
        text = buffer.str();
    }
    const Json result = Json::parse(text);
    std::filesystem::remove_all(directory, ignored);
    if (const Json* structure = result.find("structure")) {
        Structure changed = science::structureFromJson(*structure);
        changed.pbc = s.pbc;
        // Atoms keep their colours where they still exist.
        for (std::size_t i = 0; i < changed.atoms.size(); ++i) {
            if (i < s.atoms.size() && s.atoms[i].symbol == changed.atoms[i].symbol) {
                changed.atoms[i].r = s.atoms[i].r; changed.atoms[i].g = s.atoms[i].g; changed.atoms[i].b = s.atoms[i].b;
                changed.atoms[i].fixed = s.atoms[i].fixed;
            } else getDefaultElementColor(changed.atoms[i].atomicNumber, changed.atoms[i].r, changed.atoms[i].g, changed.atoms[i].b);
        }
        if (changed.atoms.size() == s.atoms.size()) {
            changed.grainColors = s.grainColors;
            changed.grainRegionIds = s.grainRegionIds;
            changed.atomProperty = s.atomProperty;
            changed.atomPropertyName = s.atomPropertyName;
            changed.frozenProperty = s.frozenProperty;
            changed.frozenPropertyName = s.frozenPropertyName;
        } else data.selected.assign(changed.atoms.size(), 0);
        changed.cellOffset = s.cellOffset;
        s = std::move(changed);
    }
    if (const Json* pbc = result.find("pbc"); pbc && pbc->isArray() && pbc->size() == 3)
        for (std::size_t k = 0; k < 3; ++k) s.pbc[k] = pbc->items()[k].boolean();
    if (const Json* mask = result.find("selected"); mask && mask->isArray()) {
        if (mask->size() != s.atoms.size()) throw std::runtime_error("The script's selection has " + std::to_string(mask->size()) + " entries for " + std::to_string(s.atoms.size()) + " atoms");
        data.selected.clear();
        for (const auto& v : mask->items()) data.selected.push_back(v.isBool() ? v.boolean() : v.number() != 0);
    }
    if (const Json* values = result.find("property"); values && values->isArray()) {
        if (values->size() != s.atoms.size()) throw std::runtime_error("The script's property has " + std::to_string(values->size()) + " values for " + std::to_string(s.atoms.size()) + " atoms");
        s.atomProperty.clear();
        for (const auto& v : values->items()) s.atomProperty.push_back(v.isNumber() ? v.number() : std::numeric_limits<double>::quiet_NaN());
        if (const Json* name = result.find("property_name"); name && name->isString()) s.atomPropertyName = name->string();
    }
    std::istringstream lines(messages);
    for (std::string line; std::getline(lines, line);)
        if (!line.empty()) data.notes.push_back(line);
}
}

ModifierType pythonScriptStepType()
{
    return
        {"python-script", "Python script", "Modification",
         "Runs a Python script on the structure: python SCRIPT INPUT.json OUTPUT.json. With atomforge.pipeline.script_step the script is a function of the\nstructure (and its selection and property) returning the changed data; what it prints becomes the step's notes.",
         {{"script", "Python script file", "string", "", ""}, {"arguments", "Arguments for the script (text)", "string", "", ""},
          {"python", "Python interpreter (empty: ATOMFORGE_PYTHON or python)", "string", "", ""}}, pythonScript};
}
}
