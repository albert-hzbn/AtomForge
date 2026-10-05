#include "science/ScienceTools.h"
#include "science/Analysis.h"
#include "science/LammpsExport.h"
#include "science/VaspElectronic.h"
#include "science/Diffraction.h"
#include "science/Clusters.h"
#include "science/DislocationLines.h"
#include "science/Phonons.h"
#include "science/ReciprocalPath.h"
#include "science/ScienceCatalog.h"
#include "science/Simulation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <set>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const ScienceToolDef& definition(const std::string& tool)
{
    for (const auto& candidate : scienceToolCatalog())
        if (tool == candidate.id) return candidate;
    throw std::runtime_error("Unknown scientific tool: " + tool);
}

void validateKeys(const ScienceToolDef& tool, const Json& request)
{
    if (!request.isObject()) throw std::runtime_error("parameters must be an object");
    std::set<std::string> known;
    for (const auto& parameter : tool.parameters) known.insert(parameter.name);
    for (const auto& member : request.members())
        if (!known.count(member.first)) throw std::runtime_error(std::string(tool.id) + " has no parameter '" + member.first + "'");
    for (const auto& parameter : tool.parameters)
        if (parameter.required && !request.contains(parameter.name))
            throw std::runtime_error(std::string(tool.id) + " requires parameter '" + parameter.name + "'");
}

DynamicsOptions dynamicsOptions(const Parameters& p, bool npt)
{
    DynamicsOptions options;
    options.steps = integer(p.number("steps", 1000), "steps");
    options.timestepFs = positive(p.number("timestep_fs", 1.0), "timestep_fs");
    options.temperatureK = positive(p.number("temperature_K", 300.0), "temperature_K");
    options.thermostatFs = positive(p.number("thermostat_fs", 100.0), "thermostat_fs");
    options.seed = static_cast<unsigned long long>(integer(p.number("seed", 0), "seed", 0));
    options.sampleInterval = integer(p.number("sample_interval", 10), "sample_interval");
    options.productionSteps = integer(p.number("production_steps", 0), "production_steps", 0);
    options.equilibrationFs = p.number("equilibration_fs", 0.0);
    if (options.equilibrationFs < 0) throw std::runtime_error("equilibration_fs must be nonnegative");
    if (npt) {
        options.pressureGPa = p.number("pressure_GPa", 0.0);
        options.barostatFs = positive(p.number("barostat_fs", 1000.0), "barostat_fs");
    }
    return options;
}

ToolOutput numeric(Json result)
{
    ToolOutput output;
    output.result = std::move(result);
    return output;
}

bool isStructure(const Json& value)
{
    return value.isObject() && value.contains("symbols") && value.contains("positions");
}

bool isNumeric(const Json& value)
{
    if (value.isArray()) {
        for (const auto& item : value.items())
            if (!isNumeric(item)) return false;
        return true;
    }
    return value.isNumber() || value.isNull() || value.isBool();
}

void collect(const Json& value, std::vector<double>& out)
{
    if (value.isArray()) for (const auto& item : value.items()) collect(item, out);
    else if (value.isNumber()) out.push_back(value.number());
    else if (value.isBool()) out.push_back(value.boolean() ? 1 : 0);
    else out.push_back(NAN);
}

std::string number(double value)
{
    if (!std::isfinite(value)) return "nan";
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.7g", value);
    return buffer;
}

std::string compact(const Json& value)
{
    if (value.isArray()) {
        std::string out = "[";
        for (std::size_t i = 0; i < value.size(); ++i) out += (i ? " " : "") + compact(value.items()[i]);
        return out + "]";
    }
    if (value.isBool()) return value.boolean() ? "True" : "False";
    if (value.isNumber()) return number(value.number());
    return "nan";
}

std::string shape(const Json& value)
{
    std::string out = "(";
    const Json* current = &value;
    bool first = true;
    while (current->isArray()) {
        out += (first ? "" : ", ") + std::to_string(current->size());
        first = false;
        if (!current->size()) break;
        current = &current->items()[0];
    }
    if (out == "(" ) return "()";
    return out + (out.find(',') == std::string::npos ? ",)" : ")");
}

void describe(const std::string& key, const Json& value, std::vector<std::string>& lines)
{
    if (isStructure(value)) {
        lines.push_back(key + ": " + std::to_string(value.at("symbols").size()) + " atoms");
    } else if (value.isObject()) {
        for (const auto& [child, item] : value.members()) describe(key.empty() ? child : key + "." + child, item, lines);
    } else if (value.isArray() && value.size() && isStructure(value.items()[0])) {
        lines.push_back(key + ": " + std::to_string(value.size()) + " structural frames");
    } else if (value.isArray() && value.size() && isNumeric(value)) {
        std::vector<double> data;
        collect(value, data);
        lines.push_back(key + " (shape " + shape(value) + "):");
        if (data.size() <= 36) lines.push_back(compact(value));
        else {
            double low = HUGE_VAL, high = -HUGE_VAL, sum = 0;
            std::size_t valid = 0;
            for (double v : data)
                if (std::isfinite(v)) { low = std::min(low, v); high = std::max(high, v); sum += v; ++valid; }
            if (valid)
                lines.push_back("  min=" + number(low) + ", max=" + number(high) + ", mean=" + number(sum / static_cast<double>(valid)) +
                                "; " + std::to_string(valid) + " valid / " + std::to_string(data.size()) + " values");
            else lines.push_back("  No valid values");
        }
    } else if (value.isArray()) {
        std::string text = value.dump();
        if (text.size() > 1000) text = text.substr(0, 1000);
        lines.push_back(key + ": " + text);
    } else if (value.isString()) {
        const std::string text = value.string();
        lines.push_back(key + ": " + (text.find('\n') != std::string::npos ? "(" + std::to_string(std::count(text.begin(), text.end(), '\n')) + "-line file; use Save generated files)" : text));
    } else if (value.isBool()) {
        lines.push_back(key + ": " + std::string(value.boolean() ? "True" : "False"));
    } else if (value.isNumber()) {
        char buffer[40];
        std::snprintf(buffer, sizeof(buffer), "%.12g", value.number());
        lines.push_back(key + ": " + buffer);
    } else {
        lines.push_back(key + ": None");
    }
}
}

ToolOutput runTool(const std::string& tool, const Json& request, const std::filesystem::path& base, const StructureReader& reader)
{
    const ScienceToolDef& def = definition(tool);
    validateKeys(def, request);
    const Parameters p(tool, request, base, reader);
    if (tool == "msd") return numeric(meanSquareDisplacement(p));
    if (tool == "diffusion") return numeric(diffusionCoefficient(p));
    if (tool == "vacf") return numeric(velocityAutocorrelation(p));
    if (tool == "vibrational-spectrum") return numeric(vibrationalSpectrum(p));
    if (tool == "local-strain") return numeric(localStrain(p));
    if (tool == "centrosymmetry") return numeric(centrosymmetry(p));
    if (tool == "structure-type") return numeric(structureType(p));
    if (tool == "cluster-analysis") return numeric(clusterAnalysis(p));
    if (tool == "void-analysis") return numeric(voidAnalysis(p));
    if (tool == "dislocation-lines") return dislocationLines(p.structure("reference"), p.structure("structure"), p);
    if (tool == "bond-order") return numeric(bondOrder(p));
    if (tool == "wigner-seitz") return numeric(wignerSeitz(p));
    if (tool == "structure-factor") return numeric(staticStructureFactor(p));
    if (tool == "band-gap") return numeric(bandGap(p));
    if (tool == "effective-mass") return numeric(effectiveMass(p));
    if (tool == "work-function") return numeric(workFunction(p));
    if (tool == "equation-of-state") return numeric(equationOfState(p));
    if (tool == "elastic-tensor") return numeric(elasticTensor(p));
    if (tool == "phonon-dos") return numeric(phononDos(p));
    if (tool == "harmonic-thermodynamics") return numeric(harmonicThermodynamics(p));
    if (tool == "neb") {
        NebOptions options;
        if (p.has("restart_images"))
            for (const auto& frame : p.frames("restart_images")) options.restart.push_back(configurationFrom(frame.structure, frame.pbc));
        else if (!p.has("initial") || !p.has("final"))
            throw std::runtime_error("NEB needs initial and final structures, or restart_images");
        const Configuration initial = options.restart.empty() ? configurationFrom(p.structure("initial").structure, p.structure("initial").pbc) : options.restart.front();
        const Configuration final = options.restart.empty() ? configurationFrom(p.structure("final").structure, p.structure("final").pbc) : options.restart.back();
        options.images = integer(p.number("images", 7), "images", 0);
        options.imageSpacing = p.number("image_spacing_A", 0.5);
        options.relaxEndpoints = p.boolean("relax_endpoints", false);
        options.endpointFmax = positive(p.number("endpoint_fmax", 0.01), "endpoint_fmax");
        options.fmax = positive(p.number("fmax", 0.03), "fmax");
        options.steps = integer(p.number("steps", 300), "steps");
        options.spring = positive(p.number("spring_eV_per_A2", 0.1), "spring_eV_per_A2");
        options.climb = p.boolean("climb", true);
        options.mic = p.boolean("mic", false);
        return migrationPath(initial, final, potentialFactory(p.json("calculator_factory"), base), options);
    }
    if (tool == "relax") {
        const auto& structure = p.structure("structure");
        const auto potential = potentialFactory(p.json("calculator"), base)();
        RelaxOptions options;
        options.fmax = positive(p.number("fmax", 0.01), "fmax");
        options.steps = integer(p.number("steps", 500), "steps", 0);
        options.relaxCell = p.boolean("relax_cell", false);
        options.pressureGPa = p.number("pressure_GPa", 0.0);
        return relaxStructure(configurationFrom(structure.structure, structure.pbc), *potential, options);
    }
    if (tool == "phonons") {
        const auto potential = potentialFactory(p.json("calculator"), base)();
        return phononCalculation(p.structure("structure"), *potential, p);
    }
    if (tool == "nvt" || tool == "npt") {
        const auto& structure = p.structure("structure");
        const auto potential = potentialFactory(p.json("calculator"), base)();
        const auto configuration = configurationFrom(structure.structure, structure.pbc);
        return tool == "nvt" ? nvtDynamics(configuration, *potential, dynamicsOptions(p, false))
                             : nptDynamics(configuration, *potential, dynamicsOptions(p, true));
    }
    if (tool == "vasp-electronic") return vaspElectronic(p);
    if (tool == "lammps-export") return lammpsExport(p);
    if (tool == "powder-xrd") return powderXrd(p.structure("structure"), p);
    if (tool == "electron-diffraction") return electronDiffraction(p.structure("structure"), p);
    if (tool == "dft-inputs")
        return dftInputs(p.structure("structure"), static_cast<int>(integer(p.number("points_per_segment", 40), "points_per_segment", 2)),
                         p.number("symprec_A", 1e-5), p.boolean("time_reversal", true));
    if (tool == "reciprocal-path")
        return reciprocalPath(p.structure("structure"), p.number("spacing_inv_A", 0.025), p.number("symprec_A", 1e-5),
                              p.boolean("time_reversal", true));
    throw std::runtime_error("Unknown scientific tool: " + tool);
}

std::string resultReport(const std::string& tool, const Json& result)
{
    std::vector<std::string> lines = {definition(tool).title, "", "Full numerical data are available through Save results.", ""};
    describe("", result, lines);
    std::string text;
    for (const auto& line : lines) text += line + "\n";
    return text;
}

Json resultDocument(const std::string& tool, const Json& request, const ToolOutput& output)
{
    Json document = Json::object();
    document["tool"] = tool;
    document["parameters"] = request;
    document["result"] = output.result;
    return document;
}
}
