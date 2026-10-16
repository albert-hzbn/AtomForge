#include "science/ScienceTools.h"
#include "util/Suggest.h"
#include "science/Analysis.h"
#include "science/LammpsExport.h"
#include "science/Symmetry.h"
#include "science/Workflows.h"
#include "science/TrajectoryStructure.h"
#include "science/VaspElectronic.h"
#include "science/Diffraction.h"
#include "science/Clusters.h"
#include "science/LocalEnvironment.h"
#include "science/PropertyAnalysis.h"
#include "science/TemplateMatching.h"
#include "science/VoronoiAnalysis.h"
#include "science/TrajectoryTools.h"
#include "science/SurfaceMesh.h"
#include "science/BondsAndRegions.h"
#include "science/DislocationLines.h"
#include "science/Phonons.h"
#include "science/ReciprocalPath.h"
#include "science/ScienceCatalog.h"
#include "science/Simulation.h"

#include <map>
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
    if (const ScienceToolDef* found = findScienceTool(tool)) return *found;
    std::vector<std::string> ids;
    for (const auto& known : scienceToolCatalog()) ids.push_back(known.id);
    throw std::runtime_error("Unknown scientific tool: " + tool + suggest::didYouMean(tool, ids) +
                             "; AtomForge --science --help lists them");
}

void validateKeys(const ScienceToolDef& tool, const Json& request)
{
    if (!request.isObject()) throw std::runtime_error("parameters must be an object");
    std::set<std::string> known;
    for (const auto& parameter : tool.parameters) known.insert(parameter.name);
    for (const auto& member : request.members())
        if (!known.count(member.first))
            throw std::runtime_error(std::string(tool.id) + " has no parameter '" + member.first + "'" +
                                     suggest::didYouMean(member.first, known) + "; it takes " + suggest::list(known));
    for (const auto& parameter : tool.parameters)
        if (parameter.required && !request.contains(parameter.name))
            throw std::runtime_error(std::string(tool.id) + " requires parameter '" + parameter.name + "'" +
                                     (parameter.value && *parameter.value
                                          ? std::string(", for example \"") + parameter.name + "\": " + parameter.value
                                          : std::string()));
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

// One entry per catalog tool; ScienceTools registry tests check the two agree.
using ToolRunner = ToolOutput (*)(const Parameters&, const std::filesystem::path& base);

// Runners of tools that only need their parameters: numerical results, or
// results with structures (frames).
template <Json (*Analyse)(const Parameters&)>
ToolOutput numericTool(const Parameters& p, const std::filesystem::path&) { return numeric(Analyse(p)); }
template <ToolOutput (*Run)(const Parameters&)>
ToolOutput outputTool(const Parameters& p, const std::filesystem::path&) { return Run(p); }

ToolOutput runNeb(const Parameters& p, const std::filesystem::path& base)
{
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

ToolOutput runRelax(const Parameters& p, const std::filesystem::path& base)
{
    const auto& structure = p.structure("structure");
    const auto potential = potentialFactory(p.json("calculator"), base)();
    RelaxOptions options;
    options.fmax = positive(p.number("fmax", 0.01), "fmax");
    options.steps = integer(p.number("steps", 500), "steps", 0);
    options.relaxCell = p.boolean("relax_cell", false);
    options.pressureGPa = p.number("pressure_GPa", 0.0);
    return relaxStructure(configurationFrom(structure.structure, structure.pbc), *potential, options);
}

ToolOutput runPhonons(const Parameters& p, const std::filesystem::path& base)
{
    const auto potential = potentialFactory(p.json("calculator"), base)();
    return phononCalculation(p.structure("structure"), *potential, p);
}

ToolOutput runDynamics(const Parameters& p, const std::filesystem::path& base, bool npt)
{
    const auto& structure = p.structure("structure");
    const auto potential = potentialFactory(p.json("calculator"), base)();
    const auto configuration = configurationFrom(structure.structure, structure.pbc);
    return !npt ? nvtDynamics(configuration, *potential, dynamicsOptions(p, false))
                         : nptDynamics(configuration, *potential, dynamicsOptions(p, true));
}

const std::map<std::string, ToolRunner>& toolRunners()
{
    static const std::map<std::string, ToolRunner> runners = {
        {"msd", numericTool<meanSquareDisplacement>},
        {"diffusion", numericTool<diffusionCoefficient>},
        {"vacf", numericTool<velocityAutocorrelation>},
        {"vibrational-spectrum", numericTool<vibrationalSpectrum>},
        {"local-strain", numericTool<localStrain>},
        {"centrosymmetry", numericTool<centrosymmetry>},
        {"structure-type", numericTool<structureType>},
        {"diamond-structure", numericTool<diamondStructure>},
        {"cluster-analysis", numericTool<clusterAnalysis>},
        {"void-analysis", numericTool<voidAnalysis>},
        {"dislocation-lines", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return dislocationLines(p.structure("reference"), p.structure("structure"), p); }},
        {"bond-order", numericTool<bondOrder>},
        {"wigner-seitz", numericTool<wignerSeitz>},
        {"trajectory-average", outputTool<trajectoryAverage>},
        {"unwrap-trajectory", outputTool<unwrapTrajectory>},
        {"frame-difference", numericTool<frameDifference>},
        {"time-series", numericTool<timeSeries>},
        {"surface-mesh", numericTool<surfaceMesh>},
        {"voronoi", numericTool<voronoiAnalysis>},
        {"vorotop", numericTool<voroTop>},
        {"cna", numericTool<adaptiveCna>},
        {"ptm", numericTool<polyhedralTemplateMatching>},
        {"grain-segmentation", numericTool<grainSegmentation>},
        {"elastic-strain", numericTool<elasticStrain>},
        {"fcc-planar-faults", numericTool<fccPlanarFaults>},
        {"chill-plus", numericTool<chillPlus>},
        {"local-entropy", numericTool<localEntropy>},
        {"bond-distributions", numericTool<bondDistributions>},
        {"bond-orders", numericTool<bondOrders>},
        {"lammps-regions", numericTool<lammpsRegions>},
        {"find-rings", numericTool<findRings>},
        {"displacements", numericTool<displacementVectors>},
        {"spatial-binning", numericTool<spatialBinning>},
        {"property-histogram", numericTool<propertyHistogram>},
        {"scatter-plot", numericTool<scatterPlot>},
        {"ambient-occlusion", numericTool<ambientOcclusion>},
        {"spatial-correlation", numericTool<spatialCorrelation>},
        {"structure-factor", numericTool<staticStructureFactor>},
        {"band-gap", numericTool<bandGap>},
        {"effective-mass", numericTool<effectiveMass>},
        {"work-function", numericTool<workFunction>},
        {"trajectory-structure", numericTool<trajectoryStructure>},
        {"elastic-constants", [](const Parameters& p, const std::filesystem::path& base) -> ToolOutput {
            return elasticConstants(p.structure("structure"), *potentialFactory(p.json("calculator"), base)(), p); }},
        {"eos-scan", [](const Parameters& p, const std::filesystem::path& base) -> ToolOutput {
            return equationOfStateScan(p.structure("structure"), *potentialFactory(p.json("calculator"), base)(), p); }},
        {"formation-energy", [](const Parameters& p, const std::filesystem::path& base) -> ToolOutput {
            return formationEnergy(p.structure("structure"), p.structure("bulk"), *potentialFactory(p.json("calculator"), base)(), p); }},
        {"planar-defect-energy", [](const Parameters& p, const std::filesystem::path& base) -> ToolOutput {
            return planarDefectEnergy(p.structure("structure"), p.structure("bulk"), *potentialFactory(p.json("calculator"), base)(), p); }},
        {"equation-of-state", numericTool<equationOfState>},
        {"elastic-tensor", numericTool<elasticTensor>},
        {"phonon-dos", numericTool<phononDos>},
        {"harmonic-thermodynamics", numericTool<harmonicThermodynamics>},
        {"neb", runNeb},
        {"relax", runRelax},
        {"phonons", runPhonons},
        {"nvt", [](const Parameters& p, const std::filesystem::path& base) { return runDynamics(p, base, false); }},
        {"npt", [](const Parameters& p, const std::filesystem::path& base) { return runDynamics(p, base, true); }},
        {"vasp-electronic", outputTool<vaspElectronic>},
        {"symmetry", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return symmetryAnalysis(p.structure("structure"), p); }},
        {"lammps-export", outputTool<lammpsExport>},
        {"powder-xrd", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return powderXrd(p.structure("structure"), p); }},
        {"electron-diffraction", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return electronDiffraction(p.structure("structure"), p); }},
        {"dft-inputs", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return dftInputs(p.structure("structure"), static_cast<int>(integer(p.number("points_per_segment", 40), "points_per_segment", 2)), p.number("symprec_A", 1e-5), p.boolean("time_reversal", true)); }},
        {"reciprocal-path", [](const Parameters& p, const std::filesystem::path&) -> ToolOutput { return reciprocalPath(p.structure("structure"), p.number("spacing_inv_A", 0.025), p.number("symprec_A", 1e-5), p.boolean("time_reversal", true)); }},
    };
    return runners;
}
}

ToolOutput runTool(const std::string& tool, const Json& request, const std::filesystem::path& base, const StructureReader& reader)
{
    const ScienceToolDef& def = definition(tool);
    validateKeys(def, request);
    const Parameters p(tool, request, base, reader);
    const auto runner = toolRunners().find(tool);
    if (runner == toolRunners().end()) throw std::runtime_error("Unknown scientific tool: " + tool);
    return runner->second(p, base);
}

Pbc periodicAxes(const Structure& structure)
{
    return {structure.hasUnitCell && structure.pbc[0], structure.hasUnitCell && structure.pbc[1], structure.hasUnitCell && structure.pbc[2]};
}

Json atomsRequest(const Structure& structure, const Pbc& pbc)
{
    Json positions = Json::array();
    for (const auto& atom : structure.atoms) positions.push(Json::array({atom.x, atom.y, atom.z}));
    Json request = Json::object();
    request["positions"] = positions;
    if (structure.hasUnitCell) {
        Json cell = Json::array();
        for (const auto& row : structure.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
        request["cell"] = cell;
        request["pbc"] = Json::array({pbc[0], pbc[1], pbc[2]});
    }
    return request;
}

Json perAtomRequest(const std::string& tool, const Structure& structure, const Pbc& pbc, double cutoffA)
{
    const ScienceToolDef& def = definition(tool);
    const auto parameter = [&](const char* name) -> const ScienceParameterDef* {
        for (const auto& candidate : def.parameters)
            if (std::strcmp(candidate.name, name) == 0) return &candidate;
        return nullptr;
    };
    Json request = atomsRequest(structure, pbc);
    // The cutoff goes to tools that require one. Others keep their default or
    // automatic range, except that the automatic range of the Ackland-Jones and
    // diamond structure types comes from the cell, so open structures need it.
    if (const ScienceParameterDef* cutoff = parameter("cutoff_A")) {
        const bool automaticNeedsCell = tool == "structure-type" || tool == "diamond-structure";
        if (cutoff->required || (automaticNeedsCell && !structure.hasUnitCell)) request["cutoff_A"] = cutoffA;
    }
    if (parameter("types")) {
        Json species = Json::array();
        for (const auto& atom : structure.atoms) species.push(static_cast<double>(atom.atomicNumber));
        request["types"] = species;
    }
    // The global density needs a fully periodic cell.
    if (parameter("local_density") && !(pbc[0] && pbc[1] && pbc[2])) request["local_density"] = true;
    return request;
}

AtomProperty analysePerAtom(const Structure& structure, const std::string& tool, double cutoffA, const std::string& property)
{
    Json request = perAtomRequest(tool, structure, periodicAxes(structure), cutoffA);
    if (tool == "bond-order") request["degrees"] = Json::array({6});  // the q6 property
    const auto properties = perAtomProperties(tool, runTool(tool, request).result);
    for (const auto& candidate : properties)
        if (property.empty() || candidate.name == property) return candidate;
    throw std::runtime_error(tool + (property.empty() ? " returned no per-atom values" : " has no per-atom property " + property));
}

std::vector<std::string> runnableTools()
{
    std::vector<std::string> tools;
    for (const auto& entry : toolRunners()) tools.push_back(entry.first);
    return tools;
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
