#include "science/Simulation.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kBoltzmann = 8.617333262145e-5;      // eV/K
// 1 eV/(Angstrom amu) expressed in Angstrom/fs^2.
const double kAcceleration = 1.602176634e-19 / (1e-10 * 1.66053906660e-27) * 1e10 * 1e-30;
const double kEvPerA3ToGPa = 160.2176634;

double maxAtomForce(const std::vector<Vec3>& forces)
{
    double result = 0;
    for (const auto& force : forces) result = std::max(result, norm(force));
    return result;
}

bool sameCell(const Configuration& a, const Configuration& b)
{
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            if (std::abs(a.cell[r][c] - b.cell[r][c]) > 1e-10) return false;
    return true;
}

// FIRE minimiser (Bitzek et al. 2006) with ASE's default parameters; returns
// the displacement for the supplied generalised forces.
class Fire
{
public:
    std::vector<Vec3> step(const std::vector<Vec3>& forces)
    {
        if (m_velocity.size() != forces.size()) m_velocity.assign(forces.size(), {0, 0, 0});
        else {
            double vf = 0, ff = 0, vv = 0;
            for (std::size_t k = 0; k < forces.size(); ++k) {
                vf += dot(forces[k], m_velocity[k]); ff += dot(forces[k], forces[k]); vv += dot(m_velocity[k], m_velocity[k]);
            }
            if (vf > 0) {
                for (std::size_t k = 0; k < forces.size(); ++k)
                    m_velocity[k] = add(scale(m_velocity[k], 1 - m_alpha), scale(forces[k], m_alpha / std::sqrt(ff) * std::sqrt(vv)));
                if (m_positiveSteps > 5) { m_dt = std::min(m_dt * 1.1, 1.0); m_alpha *= 0.99; }
                ++m_positiveSteps;
            } else {
                for (auto& v : m_velocity) v = {0, 0, 0};
                m_alpha = 0.1; m_dt *= 0.5; m_positiveSteps = 0;
            }
        }
        std::vector<Vec3> displacement(forces.size());
        double length = 0;
        for (std::size_t k = 0; k < forces.size(); ++k) {
            m_velocity[k] = add(m_velocity[k], scale(forces[k], m_dt));
            displacement[k] = scale(m_velocity[k], m_dt);
            length += dot(displacement[k], displacement[k]);
        }
        length = std::sqrt(length);
        if (length > m_maxStep) for (auto& d : displacement) d = scale(d, m_maxStep / length);
        return displacement;
    }
private:
    std::vector<Vec3> m_velocity;
    double m_dt = 0.1, m_alpha = 0.1, m_maxStep = 0.2;
    int m_positiveSteps = 0;
};

// Reproducible across platforms: SplitMix64 with Box-Muller normals.
class Random
{
public:
    explicit Random(std::uint64_t seed) : m_state(seed) {}
    double uniform()
    {
        std::uint64_t z = (m_state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z ^= z >> 31;
        return (static_cast<double>(z >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    }
    double normal()
    {
        if (m_hasSpare) { m_hasSpare = false; return m_spare; }
        const double radius = std::sqrt(-2.0 * std::log(uniform()));
        const double angle = 2.0 * std::acos(-1.0) * uniform();
        m_spare = radius * std::sin(angle);
        m_hasSpare = true;
        return radius * std::cos(angle);
    }
private:
    std::uint64_t m_state;
    double m_spare = 0;
    bool m_hasSpare = false;
};

double sinhc(double x)
{
    return std::abs(x) < 1e-4 ? 1 + x * x / 6 + x * x * x * x / 120 : std::sinh(x) / x;
}

struct Dynamics
{
    Configuration atoms;
    std::vector<double> mass;  // eV fs^2 / Angstrom^2
    std::vector<Vec3> velocity;  // Angstrom/fs
    PotentialResult forces;

    double kinetic() const
    {
        double result = 0;
        for (std::size_t i = 0; i < velocity.size(); ++i) result += 0.5 * mass[i] * dot(velocity[i], velocity[i]);
        return result;
    }
    double temperature() const { return 2 * kinetic() / (3.0 * static_cast<double>(atoms.size()) * kBoltzmann); }
};

void validateDynamics(const Configuration& start, const DynamicsOptions& options)
{
    integer(static_cast<double>(options.steps), "steps");
    integer(static_cast<double>(options.sampleInterval), "sample_interval");
    positive(options.timestepFs, "timestep_fs");
    positive(options.temperatureK, "temperature_K");
    positive(options.thermostatFs, "thermostat_fs");
    if (start.size() < 2) throw std::runtime_error("At least two atoms are required");
}

Dynamics initialise(const Configuration& start, const Potential& potential, double temperature, Random& random, bool stress)
{
    Dynamics state;
    state.atoms = start;
    for (double m : start.masses) state.mass.push_back(m / kAcceleration);
    const double kT = kBoltzmann * temperature;
    for (double m : state.mass) {
        const double width = std::sqrt(kT / m);
        // Draw components separately to keep the stream order explicit.
        const double x = random.normal(), y = random.normal(), z = random.normal();
        state.velocity.push_back({x * width, y * width, z * width});
    }
    state.forces = potential.compute(state.atoms, stress);
    return state;
}

struct Recorder
{
    ToolOutput output;
    std::vector<double> temperatures, energies, pressures;
    Json volumes = Json::array();

    void record(const Dynamics& state, double time, bool npt)
    {
        const double total = state.forces.energy + state.kinetic();
        bool finite = std::isfinite(total);
        for (const auto& p : state.atoms.positions) finite = finite && std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
        if (!finite) throw std::runtime_error("Dynamics became nonfinite; inspect timestep and potential");
        output.frames.push_back(structureFrom(state.atoms));
        output.times.push_back(time);
        output.velocities.push_back(state.velocity);
        temperatures.push_back(state.temperature());
        energies.push_back(total);
        if (std::abs(determinant(state.atoms.cell)) > 1e-12) volumes.push(cellVolume(state.atoms.cell));
        else volumes.push(Json());
        if (npt) {
            const Mat3& s = state.forces.stress;
            const double volume = cellVolume(state.atoms.cell);
            const double pressure = -(s[0][0] + s[1][1] + s[2][2]) / 3 + 2 * state.kinetic() / (3 * volume);
            pressures.push_back(pressure * kEvPerA3ToGPa);
        }
    }

    ToolOutput finish(const char* ensemble, double target)
    {
        Json result = Json::object();
        Json frames = Json::array();
        for (const auto& frame : output.frames) frames.push(structureJson(frame));
        Json velocities = Json::array();
        for (const auto& frame : output.velocities) {
            Json rows = Json::array();
            for (const auto& v : frame) rows.push(toJson(v));
            velocities.push(rows);
        }
        result["frames"] = frames;
        result["time_fs"] = toJson(output.times);
        result["temperature_K"] = toJson(temperatures);
        result["total_energy_eV"] = toJson(energies);
        result["volume_A3"] = volumes;
        result["velocities_A_per_fs"] = velocities;
        result["pressure_GPa"] = toJson(pressures);
        result["ensemble"] = ensemble;
        result["target_temperature_K"] = target;
        output.result = result;
        return std::move(output);
    }
};
}

ToolOutput migrationPath(const Configuration& initial, const Configuration& final,
                         const PotentialFactory& factory, const NebOptions& options)
{
    integer(static_cast<double>(options.images), "images", 3);
    integer(static_cast<double>(options.steps), "steps");
    positive(options.fmax, "fmax");
    positive(options.spring, "spring_eV_per_A2");
    if (!initial.size() || initial.symbols != final.symbols)
        throw std::runtime_error("NEB requires identical ordered species at both endpoints");
    if (initial.pbc != final.pbc || !sameCell(initial, final))
        throw std::runtime_error("Variable-cell NEB is not supported; endpoint cells must match");
    const std::size_t count = static_cast<std::size_t>(options.images), atoms = initial.size();
    std::vector<Configuration> chain(count, initial);
    chain.back() = final;
    std::vector<std::unique_ptr<Potential>> potentials;
    for (std::size_t i = 0; i < count; ++i) {
        potentials.push_back(factory());
        if (!potentials.back()) throw std::runtime_error("The potential factory must return a potential for each image");
    }
    const MinimumImage mic = options.mic ? MinimumImage(initial.cell, initial.pbc) : MinimumImage();
    for (std::size_t a = 0; a < atoms; ++a) {
        const Vec3 delta = scale(mic(sub(final.positions[a], initial.positions[a])), 1.0 / static_cast<double>(count - 1));
        for (std::size_t i = 1; i + 1 < count; ++i)
            chain[i].positions[a] = add(initial.positions[a], scale(delta, static_cast<double>(i)));
    }
    std::vector<PotentialResult> results(count);
    results.front() = potentials.front()->compute(chain.front(), false);
    results.back() = potentials.back()->compute(chain.back(), false);
    auto tangentVector = [&](std::size_t from, std::size_t to) {
        std::vector<Vec3> t(atoms);
        for (std::size_t a = 0; a < atoms; ++a) t[a] = mic(sub(chain[to].positions[a], chain[from].positions[a]));
        return t;
    };
    auto vectorNorm = [](const std::vector<Vec3>& v) { double s = 0; for (const auto& x : v) s += dot(x, x); return std::sqrt(s); };
    // NEB forces on the movable images, flattened image-by-image.
    auto nebForces = [&]() {
        for (std::size_t i = 1; i + 1 < count; ++i) results[i] = potentials[i]->compute(chain[i], false);
        std::size_t highest = 1;
        for (std::size_t i = 1; i + 1 < count; ++i)
            if (results[i].energy > results[highest].energy) highest = i;
        std::vector<Vec3> total;
        for (std::size_t i = 1; i + 1 < count; ++i) {
            const auto minus = tangentVector(i - 1, i), plus = tangentVector(i, i + 1);
            const double eMinus = results[i - 1].energy, e = results[i].energy, ePlus = results[i + 1].energy;
            std::vector<Vec3> tangent(atoms);
            if (ePlus > e && e > eMinus) tangent = plus;
            else if (ePlus < e && e < eMinus) tangent = minus;
            else {
                const double high = std::max(std::abs(ePlus - e), std::abs(eMinus - e));
                const double low = std::min(std::abs(ePlus - e), std::abs(eMinus - e));
                for (std::size_t a = 0; a < atoms; ++a)
                    tangent[a] = ePlus > eMinus ? add(scale(plus[a], high), scale(minus[a], low))
                                                : add(scale(plus[a], low), scale(minus[a], high));
            }
            const double length = vectorNorm(tangent);
            if (length > 0) for (auto& t : tangent) t = scale(t, 1 / length);
            double projection = 0;
            for (std::size_t a = 0; a < atoms; ++a) projection += dot(results[i].forces[a], tangent[a]);
            const bool climbing = options.climb && i == highest;
            const double springForce = options.spring * (vectorNorm(plus) - vectorNorm(minus));
            for (std::size_t a = 0; a < atoms; ++a) {
                Vec3 f = sub(results[i].forces[a], scale(tangent[a], (climbing ? 2 : 1) * projection));
                if (!climbing) f = add(f, scale(tangent[a], springForce));
                total.push_back(f);
            }
        }
        return total;
    };
    Fire fire;
    long long steps = 0;
    auto forces = nebForces();
    bool converged = maxAtomForce(forces) < options.fmax;
    while (!converged && steps < options.steps) {
        taskProgress(static_cast<double>(steps) / static_cast<double>(options.steps));
        const auto step = fire.step(forces);
        for (std::size_t i = 1, k = 0; i + 1 < count; ++i)
            for (std::size_t a = 0; a < atoms; ++a, ++k) chain[i].positions[a] = add(chain[i].positions[a], step[k]);
        ++steps;
        forces = nebForces();
        converged = maxAtomForce(forces) < options.fmax;
    }
    std::vector<double> energies;
    for (const auto& r : results) energies.push_back(r.energy);
    const double peak = *std::max_element(energies.begin(), energies.end());
    ToolOutput output;
    Json images = Json::array();
    for (const auto& image : chain) {
        output.frames.push_back(structureFrom(image));
        images.push(structureJson(output.frames.back()));
    }
    Json result = Json::object();
    result["images"] = images;
    result["energies_eV"] = toJson(energies);
    std::vector<double> coordinate = {0.0};
    for (std::size_t i = 1; i < count; ++i) coordinate.push_back(coordinate.back() + vectorNorm(tangentVector(i - 1, i)));
    result["reaction_coordinate_A"] = toJson(coordinate);
    result["forward_barrier_eV"] = peak - energies.front();
    result["reverse_barrier_eV"] = peak - energies.back();
    result["reaction_energy_eV"] = energies.back() - energies.front();
    result["converged"] = converged;
    result["steps"] = steps;
    result["endpoint_max_forces_eV_per_A"] = Json::array({maxAtomForce(results.front().forces), maxAtomForce(results.back().forces)});
    result["potential"] = potentials.front()->description();
    output.result = result;
    return output;
}

ToolOutput nvtDynamics(const Configuration& start, const Potential& potential, const DynamicsOptions& options)
{
    validateDynamics(start, options);
    Random random(options.seed);
    Dynamics state = initialise(start, potential, options.temperatureK, random, false);
    const double dt = options.timestepFs, gamma = 1 / options.thermostatFs;
    const double c1 = std::exp(-gamma * dt);
    const double kT = kBoltzmann * options.temperatureK;
    Recorder recorder;
    recorder.record(state, 0, false);
    const std::size_t n = state.atoms.size();
    for (long long step = 1; step <= options.steps; ++step) {
        if (step % 16 == 0) taskProgress(static_cast<double>(step) / static_cast<double>(options.steps));
        for (std::size_t i = 0; i < n; ++i) {
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
            state.atoms.positions[i] = add(state.atoms.positions[i], scale(state.velocity[i], 0.5 * dt));
            const double width = std::sqrt((1 - c1 * c1) * kT / state.mass[i]);
            const double x = random.normal(), y = random.normal(), z = random.normal();
            state.velocity[i] = add(scale(state.velocity[i], c1), {x * width, y * width, z * width});
            state.atoms.positions[i] = add(state.atoms.positions[i], scale(state.velocity[i], 0.5 * dt));
        }
        state.forces = potential.compute(state.atoms, false);
        for (std::size_t i = 0; i < n; ++i)
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
        if (step % options.sampleInterval == 0) recorder.record(state, static_cast<double>(step) * dt, false);
    }
    if (options.steps % options.sampleInterval) recorder.record(state, static_cast<double>(options.steps) * dt, false);
    ToolOutput output = recorder.finish("nvt", options.temperatureK);
    output.result["potential"] = potential.description();
    return output;
}

ToolOutput nptDynamics(const Configuration& start, const Potential& potential, const DynamicsOptions& options)
{
    validateDynamics(start, options);
    positive(options.barostatFs, "barostat_fs");
    if (!std::isfinite(options.pressureGPa)) throw std::runtime_error("pressure_GPa must be finite");
    if (!start.fullyPeriodic() || std::abs(determinant(start.cell)) < 1e-12) throw std::runtime_error("NPT requires a full periodic cell");
    Random random(options.seed);
    Dynamics state = initialise(start, potential, options.temperatureK, random, true);
    if (!state.forces.hasStress) throw std::runtime_error("NPT requires a potential that provides stress");
    const std::size_t n = state.atoms.size();
    const double dt = options.timestepFs, kT = kBoltzmann * options.temperatureK;
    const double dof = 3.0 * static_cast<double>(n);
    const double alpha = 1 + 3 / dof;
    const double externalPressure = options.pressureGPa / kEvPerA3ToGPa;
    const double barostatMass = (dof + 3) * kT * options.barostatFs * options.barostatFs;
    const int chain = 3;
    std::vector<double> particleMass(chain, kT * options.thermostatFs * options.thermostatFs);
    particleMass[0] = dof * kT * options.thermostatFs * options.thermostatFs;
    std::vector<double> barostatChainMass(chain, kT * options.barostatFs * options.barostatFs);
    std::vector<double> particleChain(chain, 0.0), barostatChain(chain, 0.0);
    double barostatMomentum = 0.0;
    // Nose-Hoover chain propagation over dt/2 (Martyna, Tuckerman, Klein 1996).
    auto chainHalf = [&](std::vector<double>& p, const std::vector<double>& q, double twoK, double g) {
        auto force = [&](int j, double kinetic) { return j == 0 ? kinetic - g * kT : p[j - 1] * p[j - 1] / q[j - 1] - kT; };
        const int m = static_cast<int>(p.size());
        p[m - 1] += 0.25 * dt * force(m - 1, twoK);
        for (int j = m - 2; j >= 0; --j) {
            const double damping = std::exp(-0.125 * dt * p[j + 1] / q[j + 1]);
            p[j] = p[j] * damping * damping + 0.25 * dt * force(j, twoK) * damping;
        }
        const double scaling = std::exp(-0.5 * dt * p[0] / q[0]);
        twoK *= scaling * scaling;
        for (int j = 0; j < m - 1; ++j) {
            const double damping = std::exp(-0.125 * dt * p[j + 1] / q[j + 1]);
            p[j] = p[j] * damping * damping + 0.25 * dt * force(j, twoK) * damping;
        }
        p[m - 1] += 0.25 * dt * force(m - 1, twoK);
        return scaling;
    };
    auto thermostats = [&]() {
        const double barostatScale = chainHalf(barostatChain, barostatChainMass, barostatMomentum * barostatMomentum / barostatMass, 1.0);
        barostatMomentum *= barostatScale;
        const double particleScale = chainHalf(particleChain, particleMass, 2 * state.kinetic(), dof);
        for (auto& v : state.velocity) v = scale(v, particleScale);
    };
    auto barostatForce = [&]() {
        const Mat3& s = state.forces.stress;
        const double volume = cellVolume(state.atoms.cell);
        const double virialPressure = -(s[0][0] + s[1][1] + s[2][2]) / 3;
        return alpha * 2 * state.kinetic() + 3 * volume * (virialPressure - externalPressure);
    };
    auto velocityHalf = [&]() {
        const double rate = alpha * barostatMomentum / barostatMass;
        const double decay = std::exp(-rate * 0.5 * dt), mid = std::exp(-rate * 0.25 * dt) * sinhc(rate * 0.25 * dt);
        for (std::size_t i = 0; i < n; ++i)
            state.velocity[i] = add(scale(state.velocity[i], decay), scale(state.forces.forces[i], 0.5 * dt / state.mass[i] * mid));
    };
    Recorder recorder;
    recorder.record(state, 0, true);
    for (long long step = 1; step <= options.steps; ++step) {
        if (step % 4 == 0) taskProgress(static_cast<double>(step) / static_cast<double>(options.steps));
        thermostats();
        barostatMomentum += 0.5 * dt * barostatForce();
        velocityHalf();
        const double strainRate = barostatMomentum / barostatMass;
        const double expansion = std::exp(strainRate * dt);
        const double drift = dt * std::exp(0.5 * strainRate * dt) * sinhc(0.5 * strainRate * dt);
        for (std::size_t i = 0; i < n; ++i)
            state.atoms.positions[i] = add(scale(state.atoms.positions[i], expansion), scale(state.velocity[i], drift));
        for (auto& row : state.atoms.cell) row = scale(row, expansion);
        state.forces = potential.compute(state.atoms, true);
        velocityHalf();
        barostatMomentum += 0.5 * dt * barostatForce();
        thermostats();
        if (step % options.sampleInterval == 0) recorder.record(state, static_cast<double>(step) * dt, true);
    }
    if (options.steps % options.sampleInterval) recorder.record(state, static_cast<double>(options.steps) * dt, true);
    ToolOutput output = recorder.finish("npt", options.temperatureK);
    output.result["potential"] = potential.description();
    return output;
}
}

namespace atomforge::science
{
RelaxResult relaxConfiguration(const Configuration& start, const Potential& potential, const RelaxOptions& options)
{
    positive(options.fmax, "fmax");
    integer(static_cast<double>(options.steps), "steps", 0);
    if (!std::isfinite(options.pressureGPa)) throw std::runtime_error("pressure_GPa must be finite");
    if (start.size() == 0) throw std::runtime_error("The structure has no atoms");
    const bool periodic = start.fullyPeriodic() && std::abs(determinant(start.cell)) > 1e-12;
    if (options.relaxCell && !periodic) throw std::runtime_error("Cell relaxation requires a full periodic cell");
    const std::size_t n = start.size();
    const double pressure = options.pressureGPa / kEvPerA3ToGPa;
    const double cellFactor = static_cast<double>(n);
    // Generalised coordinates: undeformed positions x = r F^-T and cellFactor * F.
    std::vector<Vec3> undeformed = start.positions;
    Mat3 deformation = identity();
    RelaxResult result;
    result.configuration = start;
    auto evaluate = [&]() {
        const Mat3 ft = transpose(deformation);
        for (std::size_t i = 0; i < n; ++i) result.configuration.positions[i] = rowTimes(undeformed[i], ft);
        for (int r = 0; r < 3; ++r) result.configuration.cell[r] = rowTimes(start.cell[r], ft);
        result.forces = potential.compute(result.configuration, periodic);
        std::vector<Vec3> generalised(n);
        for (std::size_t i = 0; i < n; ++i) generalised[i] = rowTimes(result.forces.forces[i], deformation);
        double enthalpy = result.forces.energy;
        if (options.relaxCell) {
            const double volume = cellVolume(result.configuration.cell);
            enthalpy += pressure * volume;
            Mat3 virial{};
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) virial[a][b] = -volume * (result.forces.stress[a][b] + (a == b ? pressure : 0.0));
            // Pull the virial back to the reference frame: virial F^-T.
            const Mat3 pulled = multiply(virial, transpose(inverse(deformation)));
            for (int r = 0; r < 3; ++r) generalised.push_back(scale(pulled[r], 1.0 / cellFactor));
        }
        result.enthalpies.push_back(enthalpy);
        result.maxForces.push_back(maxAtomForce(generalised));
        return generalised;
    };
    Fire fire;
    auto forces = evaluate();
    result.converged = result.maxForces.back() < options.fmax;
    while (!result.converged && result.steps < options.steps) {
        taskProgress(static_cast<double>(result.steps) / static_cast<double>(std::max<long long>(1, options.steps)));
        const auto step = fire.step(forces);
        for (std::size_t i = 0; i < n; ++i) undeformed[i] = add(undeformed[i], step[i]);
        if (options.relaxCell)
            for (int r = 0; r < 3; ++r) deformation[r] = add(deformation[r], scale(step[n + r], 1.0 / cellFactor));
        ++result.steps;
        forces = evaluate();
        if (!std::isfinite(result.enthalpies.back())) throw std::runtime_error("Relaxation became nonfinite; inspect the structure and potential");
        result.converged = result.maxForces.back() < options.fmax;
    }
    return result;
}

ToolOutput relaxStructure(const Configuration& start, const Potential& potential, const RelaxOptions& options)
{
    const RelaxResult relaxed = relaxConfiguration(start, potential, options);
    const PotentialResult initial = potential.compute(start, false);
    ToolOutput output;
    output.frames = {structureFrom(start), structureFrom(relaxed.configuration)};
    Json result = Json::object();
    result["converged"] = relaxed.converged;
    result["steps"] = relaxed.steps;
    result["initial_energy_eV"] = initial.energy;
    result["energy_eV"] = relaxed.forces.energy;
    result["energy_change_eV"] = relaxed.forces.energy - initial.energy;
    result["max_force_eV_per_A"] = maxAtomForce(relaxed.forces.forces);
    result["max_generalised_force"] = relaxed.maxForces.back();
    if (relaxed.forces.hasStress) {
        Mat3 stress{};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) stress[a][b] = relaxed.forces.stress[a][b] * kEvPerA3ToGPa;
        result["stress_GPa"] = toJson(stress);
        result["pressure_GPa"] = -(stress[0][0] + stress[1][1] + stress[2][2]) / 3;
        result["volume_A3"] = cellVolume(relaxed.configuration.cell);
        result["initial_volume_A3"] = cellVolume(start.cell);
        result["enthalpy_eV"] = relaxed.enthalpies.back();
        result["cell_A"] = toJson(relaxed.configuration.cell);
    }
    result["enthalpy_history_eV"] = toJson(relaxed.enthalpies);
    result["max_force_history"] = toJson(relaxed.maxForces);
    result["relaxed_structure"] = structureJson(output.frames.back());
    result["potential"] = potential.description();
    output.result = result;
    return output;
}
}
