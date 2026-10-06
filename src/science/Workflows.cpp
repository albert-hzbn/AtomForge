#include "science/Workflows.h"
#include "science/Analysis.h"
#include "science/Simulation.h"
#include "util/TaskControl.h"

#include <cmath>
#include <map>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kEvPerA3ToGPa = 160.21766208;

Configuration periodicConfiguration(const StructureInput& input, const char* tool)
{
    if (!input.structure.hasUnitCell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || input.structure.atoms.empty())
        throw std::runtime_error(std::string(tool) + " needs a nonempty fully periodic structure");
    return configurationFrom(input.structure, input.pbc);
}

// Homogeneous deformation: positions and cell rows times (I + strain).
Configuration deformed(const Configuration& c, const Mat3& deformation)
{
    Configuration d = c;
    for (auto& x : d.positions) x = rowTimes(x, deformation);
    for (int r = 0; r < 3; ++r) d.cell[r] = rowTimes(c.cell[r], deformation);
    return d;
}

// Ions relaxed at a fixed cell (when requested), returning the final energy and stress.
PotentialResult evaluate(Configuration& c, const Potential& potential, bool relaxIons, double fmax, long long steps, bool& converged)
{
    taskCheckpoint();
    if (!relaxIons) {
        converged = true;
        return potential.compute(c, true);
    }
    RelaxOptions options;
    options.fmax = fmax;
    options.steps = steps;
    const RelaxResult relaxed = relaxConfiguration(c, potential, options);
    c = relaxed.configuration;
    converged = relaxed.converged;
    return relaxed.forces.hasStress ? relaxed.forces : potential.compute(c, true);
}

Json array(const std::vector<double>& values) { return toJson(values); }
}

ToolOutput elasticConstants(const StructureInput& input, const Potential& potential, const Parameters& p)
{
    const Configuration reference = periodicConfiguration(input, "Elastic constants");
    const double maxStrain = positive(p.number("max_strain", 0.01), "max_strain");
    if (maxStrain > 0.05) throw std::runtime_error("max_strain should be at most 0.05 for the linear (small-strain) regime");
    const long long amplitudes = integer(p.number("strain_steps", 2), "strain_steps", 1);
    const bool relaxIons = p.boolean("relax_ions", true);
    const double fmax = positive(p.number("fmax", 1e-4), "fmax");
    const long long steps = integer(p.number("steps", 1000), "steps");

    // Engineering-Voigt strains xx, yy, zz, yz, xz, xy; shear components are 2 eps_ij.
    Json strains = Json::array(), stresses = Json::array();
    long long unconverged = 0, sample = 0;
    const long long total = 6 * 2 * amplitudes + 1;
    auto addSample = [&](const std::array<double, 6>& voigt) {
        taskProgress(static_cast<double>(sample++) / static_cast<double>(total));
        Mat3 deformation = identity();
        const int pairs[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
        for (int k = 0; k < 6; ++k) {
            const auto [a, b] = pairs[k];
            const double e = k < 3 ? voigt[static_cast<std::size_t>(k)] : 0.5 * voigt[static_cast<std::size_t>(k)];
            deformation[a][b] += e;
            if (a != b) deformation[b][a] += e;
        }
        Configuration c = deformed(reference, deformation);
        bool converged = true;
        const PotentialResult result = evaluate(c, potential, relaxIons, fmax, steps, converged);
        unconverged += !converged;
        const Mat3& s = result.stress;
        strains.push(Json::array({voigt[0], voigt[1], voigt[2], voigt[3], voigt[4], voigt[5]}));
        stresses.push(Json::array({s[0][0] * kEvPerA3ToGPa, s[1][1] * kEvPerA3ToGPa, s[2][2] * kEvPerA3ToGPa,
                                   s[1][2] * kEvPerA3ToGPa, s[0][2] * kEvPerA3ToGPa, s[0][1] * kEvPerA3ToGPa}));
    };
    addSample({0, 0, 0, 0, 0, 0});
    for (int component = 0; component < 6; ++component)
        for (long long n = 1; n <= amplitudes; ++n)
            for (int sign : {-1, 1}) {
                std::array<double, 6> voigt{};
                voigt[static_cast<std::size_t>(component)] = sign * maxStrain * static_cast<double>(n) / static_cast<double>(amplitudes);
                addSample(voigt);
            }
    Json request = Json::object();
    request["strains"] = strains;
    request["stresses_GPa"] = stresses;
    Json result = elasticTensor(Parameters("elastic-tensor", request, "."));
    result["strains"] = strains;
    result["stresses_GPa"] = stresses;
    result["max_strain"] = maxStrain;
    result["relaxed_ions"] = relaxIons;
    result["unconverged_relaxations"] = unconverged;
    result["potential"] = potential.description();
    ToolOutput output;
    output.result = result;
    return output;
}

ToolOutput equationOfStateScan(const StructureInput& input, const Potential& potential, const Parameters& p)
{
    const Configuration reference = periodicConfiguration(input, "The equation-of-state scan");
    const double range = positive(p.number("volume_range", 0.1), "volume_range");
    if (range >= 0.5) throw std::runtime_error("volume_range must be below 0.5");
    const long long points = integer(p.number("points", 11), "points", 5);
    const bool relaxIons = p.boolean("relax_ions", true);
    const double fmax = positive(p.number("fmax", 1e-3), "fmax");
    const long long steps = integer(p.number("steps", 500), "steps");
    const double volume0 = cellVolume(reference.cell);
    std::vector<double> volumes, energies, pressures;
    long long unconverged = 0;
    for (long long n = 0; n < points; ++n) {
        taskProgress(static_cast<double>(n) / static_cast<double>(points));
        const double fraction = 1 - range + 2 * range * static_cast<double>(n) / static_cast<double>(points - 1);
        Configuration c = deformed(reference, [&] { Mat3 m = identity(); for (int k = 0; k < 3; ++k) m[k][k] = std::cbrt(fraction); return m; }());
        bool converged = true;
        const PotentialResult result = evaluate(c, potential, relaxIons, fmax, steps, converged);
        unconverged += !converged;
        volumes.push_back(volume0 * fraction);
        energies.push_back(result.energy);
        pressures.push_back(-(result.stress[0][0] + result.stress[1][1] + result.stress[2][2]) / 3 * kEvPerA3ToGPa);
    }
    Json request = Json::object();
    request["volumes_A3"] = array(volumes);
    request["energies_eV"] = array(energies);
    Json result = equationOfState(Parameters("equation-of-state", request, "."));
    result["volumes_A3"] = array(volumes);
    result["energies_eV"] = array(energies);
    result["pressures_GPa"] = array(pressures);
    result["atom_count"] = static_cast<long long>(reference.size());
    result["relaxed_ions"] = relaxIons;
    result["unconverged_relaxations"] = unconverged;
    result["potential"] = potential.description();
    ToolOutput output;
    output.result = result;
    return output;
}
}

namespace atomforge::science
{
namespace
{
const double kEvPerA2ToJPerM2 = 16.02176634;

struct Relaxed
{
    Configuration configuration;
    double energy = 0;
    bool converged = true;
};

Relaxed relaxed(const Configuration& start, const Potential& potential, const Parameters& p, bool relaxCell)
{
    Relaxed out;
    if (!p.boolean("relax", true) && !relaxCell) {
        out.configuration = start;
        out.energy = potential.compute(start, false).energy;
        return out;
    }
    RelaxOptions options;
    options.fmax = positive(p.number("fmax", 1e-3), "fmax");
    options.steps = integer(p.number("steps", 2000), "steps");
    options.relaxCell = relaxCell;
    const RelaxResult result = relaxConfiguration(start, potential, options);
    out.configuration = result.configuration;
    out.energy = result.forces.energy;
    out.converged = result.converged;
    return out;
}

std::map<std::string, long long> composition(const Configuration& c)
{
    std::map<std::string, long long> counts;
    for (const auto& symbol : c.symbols) ++counts[symbol];
    return counts;
}
}

ToolOutput formationEnergy(const StructureInput& defectInput, const StructureInput& bulkInput, const Potential& potential, const Parameters& p)
{
    const Configuration defect = periodicConfiguration(defectInput, "The defect cell");
    const Configuration bulk = periodicConfiguration(bulkInput, "The bulk reference");
    taskProgress(0.0);
    const Relaxed bulkRelaxed = relaxed(bulk, potential, p, false);
    taskProgress(0.3);
    const Relaxed defectRelaxed = relaxed(defect, potential, p, p.boolean("relax_cell", false));
    const auto nBulk = composition(bulk), nDefect = composition(defect);
    // Chemical potentials: given, else the bulk energy per atom of a one-element reference.
    std::map<std::string, double> mu;
    if (p.has("chemical_potentials_eV"))
        for (const auto& [element, value] : p.json("chemical_potentials_eV").members()) mu[element] = value.number();
    std::map<std::string, long long> elements = nBulk;
    for (const auto& entry : nDefect) elements.emplace(entry.first, 0);
    double formation = defectRelaxed.energy - bulkRelaxed.energy;
    Json change = Json::object(), used = Json::object();
    for (const auto& entry : elements) {
        const std::string& element = entry.first;
        const long long dn = (nDefect.count(element) ? nDefect.at(element) : 0) - (nBulk.count(element) ? nBulk.at(element) : 0);
        change[element] = dn;
        if (dn == 0) continue;
        if (!mu.count(element)) {
            if (nBulk.size() == 1 && nBulk.count(element)) mu[element] = bulkRelaxed.energy / static_cast<double>(bulk.size());
            else throw std::runtime_error("Give chemical_potentials_eV for " + element + " (the bulk reference does not fix it)");
        }
        used[element] = mu.at(element);
        formation -= static_cast<double>(dn) * mu.at(element);
    }
    Json result = Json::object();
    result["formation_energy_eV"] = formation;
    result["defect_energy_eV"] = defectRelaxed.energy;
    result["bulk_energy_eV"] = bulkRelaxed.energy;
    result["bulk_energy_per_atom_eV"] = bulkRelaxed.energy / static_cast<double>(bulk.size());
    result["atom_change"] = change;
    result["chemical_potentials_eV"] = used;
    result["defect_converged"] = defectRelaxed.converged;
    result["bulk_converged"] = bulkRelaxed.converged;
    result["relaxation_volume_A3"] = cellVolume(defectRelaxed.configuration.cell) - cellVolume(defect.cell);
    result["potential"] = potential.description();
    ToolOutput output;
    output.result = result;
    output.frames.push_back(structureFrom(defectRelaxed.configuration));
    return output;
}

ToolOutput planarDefectEnergy(const StructureInput& cellInput, const StructureInput& bulkInput, const Potential& potential, const Parameters& p)
{
    if (!cellInput.structure.hasUnitCell || cellInput.structure.atoms.empty()) throw std::runtime_error("The defect cell needs atoms and lattice vectors");
    const Configuration cell = configurationFrom(cellInput.structure, cellInput.pbc);
    const Configuration bulk = periodicConfiguration(bulkInput, "The bulk reference");
    const long long interfaces = integer(p.number("interfaces", 2), "interfaces");
    const int axis = static_cast<int>(integer(p.number("normal_axis", 2), "normal_axis", 0));
    if (axis > 2) throw std::runtime_error("normal_axis must be 0 (a), 1 (b) or 2 (c)");
    // The cell must hold a whole multiple of the bulk composition.
    const auto nCell = composition(cell), nBulk = composition(bulk);
    const double ratio = static_cast<double>(cell.size()) / static_cast<double>(bulk.size());
    bool stoichiometric = nCell.size() == nBulk.size();
    for (const auto& [element, count] : nBulk) {
        const auto found = nCell.find(element);
        stoichiometric = stoichiometric && found != nCell.end() && std::abs(static_cast<double>(found->second) - ratio * static_cast<double>(count)) < 1e-6;
    }
    if (!stoichiometric) throw std::runtime_error("The defect cell is not stoichiometric with the bulk reference; planar energies need the same composition");
    taskProgress(0.0);
    const Relaxed bulkRelaxed = relaxed(bulk, potential, p, false);
    taskProgress(0.3);
    const Relaxed cellRelaxed = relaxed(cell, potential, p, false);
    const Vec3 u = cell.cell[static_cast<std::size_t>((axis + 1) % 3)], v = cell.cell[static_cast<std::size_t>((axis + 2) % 3)];
    const double area = norm(cross(u, v));
    const double excess = cellRelaxed.energy - ratio * bulkRelaxed.energy;
    const double gamma = excess / (static_cast<double>(interfaces) * area);
    Json result = Json::object();
    result["energy_eV_per_A2"] = gamma;
    result["energy_J_per_m2"] = gamma * kEvPerA2ToJPerM2;
    result["energy_mJ_per_m2"] = gamma * kEvPerA2ToJPerM2 * 1000;
    result["excess_energy_eV"] = excess;
    result["area_A2"] = area;
    result["interfaces"] = interfaces;
    result["cell_energy_eV"] = cellRelaxed.energy;
    result["bulk_energy_per_atom_eV"] = bulkRelaxed.energy / static_cast<double>(bulk.size());
    result["cell_converged"] = cellRelaxed.converged;
    result["bulk_converged"] = bulkRelaxed.converged;
    result["potential"] = potential.description();
    ToolOutput output;
    output.result = result;
    output.frames.push_back(structureFrom(cellRelaxed.configuration));
    return output;
}
}
