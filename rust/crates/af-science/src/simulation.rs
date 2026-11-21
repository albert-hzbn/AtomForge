//! Simulations with the native potentials (the C++ science/Simulation):
//! climbing-image NEB, FIRE relaxation of positions and cell, Langevin NVT and
//! Martyna-Tobias-Klein NPT dynamics with optional NVE production.

use af_core::json::Json;

use crate::analysis::jv;
use crate::core::*;
use crate::data::{structure_json, Parameters};
use crate::potentials::{configuration_from, potential_factory, structure_from, Configuration, Potential, PotentialResult};
use crate::tools::ToolOutput;

const BOLTZMANN: f64 = 8.617333262145e-5; // eV/K
/// 1 eV/(Angstrom amu) expressed in Angstrom/fs^2.
const ACCELERATION: f64 = 1.602176634e-19 / (1e-10 * 1.66053906660e-27) * 1e10 * 1e-30;
pub const EV_PER_A3_TO_GPA: f64 = 160.2176634;

fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

fn min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

pub fn max_atom_force(forces: &[Vec3]) -> f64 {
    let mut result = 0.0;
    for force in forces {
        result = max(result, norm(*force));
    }
    result
}

fn same_cell(a: &Configuration, b: &Configuration) -> bool {
    (0..3).all(|r| (0..3).all(|c| (a.cell[r][c] - b.cell[r][c]).abs() <= 1e-10))
}

/// FIRE minimiser (Bitzek et al. 2006) with ASE's default parameters; returns
/// the displacement for the supplied generalised forces.
struct Fire {
    velocity: Vec<Vec3>,
    dt: f64,
    alpha: f64,
    max_step: f64,
    positive_steps: i32,
}

impl Fire {
    fn new() -> Self {
        Self { velocity: Vec::new(), dt: 0.1, alpha: 0.1, max_step: 0.2, positive_steps: 0 }
    }

    fn step(&mut self, forces: &[Vec3]) -> Vec<Vec3> {
        if self.velocity.len() != forces.len() {
            self.velocity = vec![[0.0; 3]; forces.len()];
        } else {
            let (mut vf, mut ff, mut vv) = (0.0, 0.0, 0.0);
            for k in 0..forces.len() {
                vf += dot(forces[k], self.velocity[k]);
                ff += dot(forces[k], forces[k]);
                vv += dot(self.velocity[k], self.velocity[k]);
            }
            if vf > 0.0 {
                for k in 0..forces.len() {
                    self.velocity[k] = add(scale(self.velocity[k], 1.0 - self.alpha), scale(forces[k], self.alpha / ff.sqrt() * vv.sqrt()));
                }
                if self.positive_steps > 5 {
                    self.dt = min(self.dt * 1.1, 1.0);
                    self.alpha *= 0.99;
                }
                self.positive_steps += 1;
            } else {
                for v in self.velocity.iter_mut() {
                    *v = [0.0; 3];
                }
                self.alpha = 0.1;
                self.dt *= 0.5;
                self.positive_steps = 0;
            }
        }
        let mut displacement = vec![[0.0; 3]; forces.len()];
        let mut length = 0.0;
        for k in 0..forces.len() {
            self.velocity[k] = add(self.velocity[k], scale(forces[k], self.dt));
            displacement[k] = scale(self.velocity[k], self.dt);
            length += dot(displacement[k], displacement[k]);
        }
        length = f64::sqrt(length);
        if length > self.max_step {
            for d in displacement.iter_mut() {
                *d = scale(*d, self.max_step / length);
            }
        }
        displacement
    }
}

/// Reproducible across platforms: SplitMix64 with Box-Muller normals.
struct Random {
    state: u64,
    spare: f64,
    has_spare: bool,
}

impl Random {
    fn new(seed: u64) -> Self {
        Self { state: seed, spare: 0.0, has_spare: false }
    }

    fn uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) * (1.0 / 9007199254740992.0)
    }

    fn normal(&mut self) -> f64 {
        if self.has_spare {
            self.has_spare = false;
            return self.spare;
        }
        let radius = (-2.0 * self.uniform().ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * self.uniform();
        self.spare = radius * angle.sin();
        self.has_spare = true;
        radius * angle.cos()
    }
}

fn sinhc(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        1.0 + x * x / 6.0 + x * x * x * x / 120.0
    } else {
        x.sinh() / x
    }
}

struct Dynamics {
    atoms: Configuration,
    /// eV fs^2 / Angstrom^2.
    mass: Vec<f64>,
    /// Angstrom/fs.
    velocity: Vec<Vec3>,
    forces: PotentialResult,
}

impl Dynamics {
    fn kinetic(&self) -> f64 {
        let mut result = 0.0;
        for i in 0..self.velocity.len() {
            result += 0.5 * self.mass[i] * dot(self.velocity[i], self.velocity[i]);
        }
        result
    }

    fn temperature(&self) -> f64 {
        2.0 * self.kinetic() / (3.0 * self.atoms.len() as f64 * BOLTZMANN)
    }
}

#[derive(Clone, Debug)]
pub struct DynamicsOptions {
    pub steps: i64,
    pub timestep_fs: f64,
    pub temperature_k: f64,
    pub thermostat_fs: f64,
    pub pressure_gpa: f64,
    pub barostat_fs: f64,
    pub seed: u64,
    pub sample_interval: i64,
    /// NVE (velocity Verlet) steps after the thermostatted run, at the final cell.
    pub production_steps: i64,
    /// Samples before this time are excluded from the statistics.
    pub equilibration_fs: f64,
}

impl Default for DynamicsOptions {
    fn default() -> Self {
        Self {
            steps: 1000,
            timestep_fs: 1.0,
            temperature_k: 300.0,
            thermostat_fs: 100.0,
            pressure_gpa: 0.0,
            barostat_fs: 1000.0,
            seed: 0,
            sample_interval: 10,
            production_steps: 0,
            equilibration_fs: 0.0,
        }
    }
}

fn validate_dynamics(start: &Configuration, options: &DynamicsOptions) -> Result<(), String> {
    integer(options.steps as f64, "steps", 1)?;
    integer(options.sample_interval as f64, "sample_interval", 1)?;
    positive(options.timestep_fs, "timestep_fs")?;
    positive(options.temperature_k, "temperature_K")?;
    positive(options.thermostat_fs, "thermostat_fs")?;
    if start.len() < 2 {
        return Err("At least two atoms are required".into());
    }
    Ok(())
}

fn initialise(start: &Configuration, potential: &dyn Potential, temperature: f64, random: &mut Random, stress: bool) -> Result<Dynamics, String> {
    let mass: Vec<f64> = start.masses.iter().map(|m| m / ACCELERATION).collect();
    let k_t = BOLTZMANN * temperature;
    let mut velocity = Vec::new();
    for m in &mass {
        let width = (k_t / m).sqrt();
        // Draw components separately to keep the stream order explicit.
        let x = random.normal();
        let y = random.normal();
        let z = random.normal();
        velocity.push([x * width, y * width, z * width]);
    }
    let forces = potential.compute(start, stress)?;
    Ok(Dynamics { atoms: start.clone(), mass, velocity, forces })
}

#[derive(Default)]
struct Recorder {
    output: ToolOutput,
    temperatures: Vec<f64>,
    energies: Vec<f64>,
    pressures: Vec<f64>,
    volumes: Vec<Json>,
    production_start: f64,
    production_index: usize,
}

impl Recorder {
    fn new() -> Self {
        Self { production_start: -1.0, ..Default::default() }
    }

    fn record(&mut self, state: &Dynamics, time: f64, npt: bool) -> Result<(), String> {
        let total = state.forces.energy + state.kinetic();
        let mut finite = total.is_finite();
        for p in &state.atoms.positions {
            finite = finite && p[0].is_finite() && p[1].is_finite() && p[2].is_finite();
        }
        if !finite {
            return Err("Dynamics became nonfinite; inspect timestep and potential".into());
        }
        self.output.frames.push(structure_from(&state.atoms));
        self.output.times.push(time);
        self.output.velocities.push(state.velocity.clone());
        self.temperatures.push(state.temperature());
        self.energies.push(total);
        if determinant(&state.atoms.cell).abs() > 1e-12 {
            self.volumes.push(Json::from(cell_volume(&state.atoms.cell)));
        } else {
            self.volumes.push(Json::Null);
        }
        if npt {
            let s = &state.forces.stress;
            let volume = cell_volume(&state.atoms.cell);
            let pressure = -(s[0][0] + s[1][1] + s[2][2]) / 3.0 + 2.0 * state.kinetic() / (3.0 * volume);
            self.pressures.push(pressure * EV_PER_A3_TO_GPA);
        }
        Ok(())
    }

    fn finish(&mut self, ensemble: &str, target: f64) -> ToolOutput {
        let mut result = Json::object();
        result.set("frames", Json::Array(self.output.frames.iter().map(structure_json).collect()));
        result.set("time_fs", jv(&self.output.times));
        result.set("temperature_K", jv(&self.temperatures));
        result.set("total_energy_eV", jv(&self.energies));
        result.set("volume_A3", Json::Array(self.volumes.clone()));
        result.set(
            "velocities_A_per_fs",
            Json::Array(self.output.velocities.iter().map(|frame| Json::Array(frame.iter().map(|v| jv(v)).collect())).collect()),
        );
        result.set("pressure_GPa", jv(&self.pressures));
        result.set("ensemble", ensemble);
        result.set("target_temperature_K", target);
        let mut output = self.output.clone();
        output.result = result;
        output
    }
}

/// Mean, standard deviation and the standard error of the mean from block
/// averages (contiguous blocks), which accounts for correlated samples.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlockStatistics {
    pub mean: f64,
    pub standard_deviation: f64,
    pub standard_error: f64,
    pub samples: usize,
    pub blocks: usize,
}

pub fn block_statistics(values: &[f64], blocks: usize) -> BlockStatistics {
    let mut result = BlockStatistics::default();
    let finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    result.samples = finite.len();
    if finite.is_empty() {
        result.mean = f64::NAN;
        result.standard_deviation = f64::NAN;
        result.standard_error = f64::NAN;
        return result;
    }
    let mut sum = 0.0;
    for v in &finite {
        sum += v;
    }
    result.mean = sum / finite.len() as f64;
    let mut variance = 0.0;
    for v in &finite {
        variance += (v - result.mean) * (v - result.mean);
    }
    result.standard_deviation = if finite.len() > 1 { (variance / (finite.len() - 1) as f64).sqrt() } else { 0.0 };
    if blocks < 2 || finite.len() < 2 * blocks {
        result.standard_error = f64::NAN;
        return result;
    }
    let size = finite.len() / blocks;
    let mut means = Vec::new();
    for b in 0..blocks {
        let mut total = 0.0;
        for v in &finite[b * size..(b + 1) * size] {
            total += v;
        }
        means.push(total / size as f64);
    }
    let mut mean_of_means = 0.0;
    for m in &means {
        mean_of_means += m;
    }
    mean_of_means /= blocks as f64;
    let mut spread = 0.0;
    for m in &means {
        spread += (m - mean_of_means) * (m - mean_of_means);
    }
    result.standard_error = (spread / (blocks - 1) as f64 / blocks as f64).sqrt();
    result.blocks = blocks;
    result
}

/// Continues with energy-conserving velocity Verlet at fixed cell.
fn run_production(state: &mut Dynamics, potential: &dyn Potential, options: &DynamicsOptions, recorder: &mut Recorder, npt: bool) -> Result<(), String> {
    if options.production_steps <= 0 {
        return Ok(());
    }
    integer(options.production_steps as f64, "production_steps", 1)?;
    let dt = options.timestep_fs;
    let start = options.steps as f64 * dt;
    recorder.production_start = start;
    recorder.production_index = recorder.output.times.len() - 1;
    let n = state.atoms.len();
    for step in 1..=options.production_steps {
        for i in 0..n {
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
            state.atoms.positions[i] = add(state.atoms.positions[i], scale(state.velocity[i], dt));
        }
        state.forces = potential.compute(&state.atoms, npt)?;
        for i in 0..n {
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
        }
        if step % options.sample_interval == 0 || step == options.production_steps {
            recorder.record(state, start + step as f64 * dt, npt)?;
        }
    }
    Ok(())
}

fn statistics_json(s: &BlockStatistics) -> Json {
    let mut result = Json::object();
    result.set("mean", s.mean);
    result.set("standard_deviation", s.standard_deviation);
    result.set("standard_error", if s.standard_error.is_finite() { Json::from(s.standard_error) } else { Json::Null });
    result.set("samples", s.samples);
    result.set("blocks", s.blocks);
    result
}

fn add_statistics(output: &mut ToolOutput, recorder: &Recorder, options: &DynamicsOptions, npt: bool) {
    // Statistics use production samples when present, else samples after equilibration.
    let production = recorder.production_start >= 0.0;
    let selected: Vec<usize> = (0..recorder.output.times.len())
        .filter(|&i| if production { i >= recorder.production_index } else { recorder.output.times[i] >= options.equilibration_fs })
        .collect();
    let pick = |values: &[f64]| -> Vec<f64> { selected.iter().filter(|&&i| i < values.len()).map(|&i| values[i]).collect() };
    let mut statistics = Json::object();
    statistics.set("temperature_K", statistics_json(&block_statistics(&pick(&recorder.temperatures), 5)));
    statistics.set("total_energy_eV", statistics_json(&block_statistics(&pick(&recorder.energies), 5)));
    if npt {
        let volumes: Vec<f64> = recorder.volumes.iter().map(|v| if v.is_number() { v.number().unwrap() } else { f64::NAN }).collect();
        statistics.set("volume_A3", statistics_json(&block_statistics(&pick(&volumes), 5)));
        statistics.set("pressure_GPa", statistics_json(&block_statistics(&pick(&recorder.pressures), 5)));
    }
    output.result.set("statistics", statistics);
    output.result.set("equilibration_fs", options.equilibration_fs);
    if production {
        output.result.set("production_start_fs", recorder.production_start);
        let energies = pick(&recorder.energies);
        let span = recorder.output.times.last().unwrap() - recorder.production_start;
        let atoms = recorder.output.frames.last().unwrap().atoms.len() as f64;
        if energies.len() > 1 && span > 0.0 {
            output.result.set("nve_energy_drift_eV_per_atom_per_ps", (energies.last().unwrap() - energies[0]) / atoms / (span / 1000.0));
        }
    }
}

/// Langevin NVT (BAOAB splitting) with Maxwell-Boltzmann initial velocities.
pub fn nvt_dynamics(start: &Configuration, potential: &dyn Potential, options: &DynamicsOptions) -> Result<ToolOutput, String> {
    validate_dynamics(start, options)?;
    let mut random = Random::new(options.seed);
    let mut state = initialise(start, potential, options.temperature_k, &mut random, false)?;
    let dt = options.timestep_fs;
    let gamma = 1.0 / options.thermostat_fs;
    let c1 = (-gamma * dt).exp();
    let k_t = BOLTZMANN * options.temperature_k;
    let mut recorder = Recorder::new();
    recorder.record(&state, 0.0, false)?;
    let n = state.atoms.len();
    for step in 1..=options.steps {
        for i in 0..n {
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
            state.atoms.positions[i] = add(state.atoms.positions[i], scale(state.velocity[i], 0.5 * dt));
            let width = ((1.0 - c1 * c1) * k_t / state.mass[i]).sqrt();
            let x = random.normal();
            let y = random.normal();
            let z = random.normal();
            state.velocity[i] = add(scale(state.velocity[i], c1), [x * width, y * width, z * width]);
            state.atoms.positions[i] = add(state.atoms.positions[i], scale(state.velocity[i], 0.5 * dt));
        }
        state.forces = potential.compute(&state.atoms, false)?;
        for i in 0..n {
            state.velocity[i] = add(state.velocity[i], scale(state.forces.forces[i], 0.5 * dt / state.mass[i]));
        }
        if step % options.sample_interval == 0 {
            recorder.record(&state, step as f64 * dt, false)?;
        }
    }
    if options.steps % options.sample_interval != 0 {
        recorder.record(&state, options.steps as f64 * dt, false)?;
    }
    run_production(&mut state, potential, options, &mut recorder, false)?;
    let mut output = recorder.finish("nvt", options.temperature_k);
    add_statistics(&mut output, &recorder, options, false);
    output.result.set("potential", potential.description());
    Ok(output)
}

/// Nose-Hoover chain propagation over dt/2 (Martyna, Tuckerman, Klein 1996).
fn chain_half(p: &mut [f64], q: &[f64], mut two_k: f64, g: f64, dt: f64, k_t: f64) -> f64 {
    let force = |p: &[f64], j: usize, kinetic: f64| if j == 0 { kinetic - g * k_t } else { p[j - 1] * p[j - 1] / q[j - 1] - k_t };
    let m = p.len();
    p[m - 1] += 0.25 * dt * force(p, m - 1, two_k);
    for j in (0..m - 1).rev() {
        let damping = (-0.125 * dt * p[j + 1] / q[j + 1]).exp();
        p[j] = p[j] * damping * damping + 0.25 * dt * force(p, j, two_k) * damping;
    }
    let scaling = (-0.5 * dt * p[0] / q[0]).exp();
    two_k *= scaling * scaling;
    for j in 0..m - 1 {
        let damping = (-0.125 * dt * p[j + 1] / q[j + 1]).exp();
        p[j] = p[j] * damping * damping + 0.25 * dt * force(p, j, two_k) * damping;
    }
    p[m - 1] += 0.25 * dt * force(p, m - 1, two_k);
    scaling
}

/// Isotropic Martyna-Tobias-Klein NPT with Nose-Hoover chains (Tuckerman et al. 2006 splitting).
pub fn npt_dynamics(start: &Configuration, potential: &dyn Potential, options: &DynamicsOptions) -> Result<ToolOutput, String> {
    validate_dynamics(start, options)?;
    positive(options.barostat_fs, "barostat_fs")?;
    if !options.pressure_gpa.is_finite() {
        return Err("pressure_GPa must be finite".into());
    }
    if !start.fully_periodic() || determinant(&start.cell).abs() < 1e-12 {
        return Err("NPT requires a full periodic cell".into());
    }
    let mut random = Random::new(options.seed);
    let mut state = initialise(start, potential, options.temperature_k, &mut random, true)?;
    if !state.forces.has_stress {
        return Err("NPT requires a potential that provides stress".into());
    }
    let n = state.atoms.len();
    let dt = options.timestep_fs;
    let k_t = BOLTZMANN * options.temperature_k;
    let dof = 3.0 * n as f64;
    let alpha = 1.0 + 3.0 / dof;
    let external_pressure = options.pressure_gpa / EV_PER_A3_TO_GPA;
    let barostat_mass = (dof + 3.0) * k_t * options.barostat_fs * options.barostat_fs;
    let chain = 3;
    let mut particle_mass = vec![k_t * options.thermostat_fs * options.thermostat_fs; chain];
    particle_mass[0] = dof * k_t * options.thermostat_fs * options.thermostat_fs;
    let barostat_chain_mass = vec![k_t * options.barostat_fs * options.barostat_fs; chain];
    let mut particle_chain = vec![0.0; chain];
    let mut barostat_chain = vec![0.0; chain];
    let mut barostat_momentum = 0.0;
    let thermostats = |state: &mut Dynamics, barostat_momentum: &mut f64, particle_chain: &mut Vec<f64>, barostat_chain: &mut Vec<f64>| {
        let barostat_scale = chain_half(barostat_chain, &barostat_chain_mass, *barostat_momentum * *barostat_momentum / barostat_mass, 1.0, dt, k_t);
        *barostat_momentum *= barostat_scale;
        let particle_scale = chain_half(particle_chain, &particle_mass, 2.0 * state.kinetic(), dof, dt, k_t);
        for v in state.velocity.iter_mut() {
            *v = scale(*v, particle_scale);
        }
    };
    let barostat_force = |state: &Dynamics| {
        let s = &state.forces.stress;
        let volume = cell_volume(&state.atoms.cell);
        let virial_pressure = -(s[0][0] + s[1][1] + s[2][2]) / 3.0;
        alpha * 2.0 * state.kinetic() + 3.0 * volume * (virial_pressure - external_pressure)
    };
    let velocity_half = |state: &mut Dynamics, barostat_momentum: f64| {
        let rate = alpha * barostat_momentum / barostat_mass;
        let decay = (-rate * 0.5 * dt).exp();
        let mid = (-rate * 0.25 * dt).exp() * sinhc(rate * 0.25 * dt);
        for i in 0..n {
            state.velocity[i] = add(scale(state.velocity[i], decay), scale(state.forces.forces[i], 0.5 * dt / state.mass[i] * mid));
        }
    };
    let mut recorder = Recorder::new();
    recorder.record(&state, 0.0, true)?;
    for step in 1..=options.steps {
        thermostats(&mut state, &mut barostat_momentum, &mut particle_chain, &mut barostat_chain);
        barostat_momentum += 0.5 * dt * barostat_force(&state);
        velocity_half(&mut state, barostat_momentum);
        let strain_rate = barostat_momentum / barostat_mass;
        let expansion = (strain_rate * dt).exp();
        let drift = dt * (0.5 * strain_rate * dt).exp() * sinhc(0.5 * strain_rate * dt);
        for i in 0..n {
            state.atoms.positions[i] = add(scale(state.atoms.positions[i], expansion), scale(state.velocity[i], drift));
        }
        for row in state.atoms.cell.iter_mut() {
            *row = scale(*row, expansion);
        }
        state.forces = potential.compute(&state.atoms, true)?;
        velocity_half(&mut state, barostat_momentum);
        barostat_momentum += 0.5 * dt * barostat_force(&state);
        thermostats(&mut state, &mut barostat_momentum, &mut particle_chain, &mut barostat_chain);
        if step % options.sample_interval == 0 {
            recorder.record(&state, step as f64 * dt, true)?;
        }
    }
    if options.steps % options.sample_interval != 0 {
        recorder.record(&state, options.steps as f64 * dt, true)?;
    }
    run_production(&mut state, potential, options, &mut recorder, true)?;
    let mut output = recorder.finish("npt", options.temperature_k);
    add_statistics(&mut output, &recorder, options, true);
    output.result.set("potential", potential.description());
    Ok(output)
}

#[derive(Clone, Debug)]
pub struct RelaxOptions {
    /// eV/Angstrom, also applied to the scaled cell forces.
    pub fmax: f64,
    pub steps: i64,
    /// Full 3x3 cell (shape and volume).
    pub relax_cell: bool,
    /// Target hydrostatic pressure for cell relaxation.
    pub pressure_gpa: f64,
}

impl Default for RelaxOptions {
    fn default() -> Self {
        Self { fmax: 0.01, steps: 500, relax_cell: false, pressure_gpa: 0.0 }
    }
}

pub struct RelaxResult {
    pub configuration: Configuration,
    /// At the final configuration (stress when periodic).
    pub forces: PotentialResult,
    pub converged: bool,
    pub steps: i64,
    /// E + P V per evaluated step (eV).
    pub enthalpies: Vec<f64>,
    /// Largest generalised force per step.
    pub max_forces: Vec<f64>,
}

/// FIRE minimisation of atomic positions and, optionally, the cell through a
/// deformation gradient with ASE UnitCellFilter scaling (cell factor = atom count).
pub fn relax_configuration(start: &Configuration, potential: &dyn Potential, options: &RelaxOptions) -> Result<RelaxResult, String> {
    positive(options.fmax, "fmax")?;
    integer(options.steps as f64, "steps", 0)?;
    if !options.pressure_gpa.is_finite() {
        return Err("pressure_GPa must be finite".into());
    }
    if start.is_empty() {
        return Err("The structure has no atoms".into());
    }
    let periodic = start.fully_periodic() && determinant(&start.cell).abs() > 1e-12;
    if options.relax_cell && !periodic {
        return Err("Cell relaxation requires a full periodic cell".into());
    }
    let n = start.len();
    let pressure = options.pressure_gpa / EV_PER_A3_TO_GPA;
    let cell_factor = n as f64;
    // Generalised coordinates: undeformed positions x = r F^-T and cell_factor * F.
    let mut undeformed = start.positions.clone();
    let mut deformation = identity();
    let mut result = RelaxResult {
        configuration: start.clone(),
        forces: PotentialResult::default(),
        converged: false,
        steps: 0,
        enthalpies: Vec::new(),
        max_forces: Vec::new(),
    };
    let evaluate = |result: &mut RelaxResult, undeformed: &[Vec3], deformation: &Mat3| -> Result<Vec<Vec3>, String> {
        let ft = transpose(deformation);
        for i in 0..n {
            result.configuration.positions[i] = row_times(undeformed[i], &ft);
        }
        for r in 0..3 {
            result.configuration.cell[r] = row_times(start.cell[r], &ft);
        }
        result.forces = potential.compute(&result.configuration, periodic)?;
        let mut generalised: Vec<Vec3> = (0..n).map(|i| row_times(result.forces.forces[i], deformation)).collect();
        let mut enthalpy = result.forces.energy;
        if options.relax_cell {
            let volume = cell_volume(&result.configuration.cell);
            enthalpy += pressure * volume;
            let mut virial = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    virial[a][b] = -volume * (result.forces.stress[a][b] + if a == b { pressure } else { 0.0 });
                }
            }
            // Pull the virial back to the reference frame: virial F^-T.
            let pulled = multiply(&virial, &transpose(&inverse(deformation)?));
            for r in 0..3 {
                generalised.push(scale(pulled[r], 1.0 / cell_factor));
            }
        }
        result.enthalpies.push(enthalpy);
        result.max_forces.push(max_atom_force(&generalised));
        Ok(generalised)
    };
    let mut fire = Fire::new();
    let mut forces = evaluate(&mut result, &undeformed, &deformation)?;
    result.converged = *result.max_forces.last().unwrap() < options.fmax;
    while !result.converged && result.steps < options.steps {
        let step = fire.step(&forces);
        for i in 0..n {
            undeformed[i] = add(undeformed[i], step[i]);
        }
        if options.relax_cell {
            for r in 0..3 {
                deformation[r] = add(deformation[r], scale(step[n + r], 1.0 / cell_factor));
            }
        }
        result.steps += 1;
        forces = evaluate(&mut result, &undeformed, &deformation)?;
        if !result.enthalpies.last().unwrap().is_finite() {
            return Err("Relaxation became nonfinite; inspect the structure and potential".into());
        }
        result.converged = *result.max_forces.last().unwrap() < options.fmax;
    }
    Ok(result)
}

pub fn relax_structure(start: &Configuration, potential: &dyn Potential, options: &RelaxOptions) -> Result<ToolOutput, String> {
    let relaxed = relax_configuration(start, potential, options)?;
    let initial = potential.compute(start, false)?;
    let mut output = ToolOutput { frames: vec![structure_from(start), structure_from(&relaxed.configuration)], ..Default::default() };
    let mut result = Json::object();
    result.set("converged", relaxed.converged);
    result.set("steps", relaxed.steps);
    result.set("initial_energy_eV", initial.energy);
    result.set("energy_eV", relaxed.forces.energy);
    result.set("energy_change_eV", relaxed.forces.energy - initial.energy);
    result.set("max_force_eV_per_A", max_atom_force(&relaxed.forces.forces));
    result.set("max_generalised_force", *relaxed.max_forces.last().unwrap());
    if relaxed.forces.has_stress {
        let stress: Mat3 = std::array::from_fn(|a| std::array::from_fn(|b| relaxed.forces.stress[a][b] * EV_PER_A3_TO_GPA));
        result.set("stress_GPa", mat3_json(&stress));
        result.set("pressure_GPa", -(stress[0][0] + stress[1][1] + stress[2][2]) / 3.0);
        result.set("volume_A3", cell_volume(&relaxed.configuration.cell));
        result.set("initial_volume_A3", cell_volume(&start.cell));
        result.set("enthalpy_eV", *relaxed.enthalpies.last().unwrap());
        result.set("cell_A", mat3_json(&relaxed.configuration.cell));
    }
    result.set("enthalpy_history_eV", jv(&relaxed.enthalpies));
    result.set("max_force_history", jv(&relaxed.max_forces));
    result.set("relaxed_structure", structure_json(output.frames.last().unwrap()));
    result.set("potential", potential.description());
    output.result = result;
    Ok(output)
}

#[derive(Clone, Debug)]
pub struct NebOptions {
    /// Including both endpoints; 0 chooses the count from the endpoint separation.
    pub images: i64,
    pub fmax: f64,
    pub steps: i64,
    pub spring: f64,
    pub climb: bool,
    pub mic: bool,
    /// Angstrom of total displacement per image (images == 0).
    pub image_spacing: f64,
    /// Fixed-cell FIRE relaxation of both endpoints first.
    pub relax_endpoints: bool,
    pub endpoint_fmax: f64,
    /// Optional restart chain (all images, endpoints included); replaces interpolation.
    pub restart: Vec<Configuration>,
}

impl Default for NebOptions {
    fn default() -> Self {
        Self {
            images: 7,
            fmax: 0.03,
            steps: 300,
            spring: 0.1,
            climb: true,
            mic: false,
            image_spacing: 0.5,
            relax_endpoints: false,
            endpoint_fmax: 0.01,
            restart: Vec::new(),
        }
    }
}

/// Fixed-cell climbing-image NEB with the improved tangent (Henkelman and
/// Jonsson 2000) and FIRE optimisation, with ASE's default FIRE parameters.
pub fn migration_path(initial: &Configuration, last_input: &Configuration, potential: &dyn Potential, options: &NebOptions) -> Result<ToolOutput, String> {
    integer(options.images as f64, "images", 0)?;
    integer(options.steps as f64, "steps", 1)?;
    positive(options.fmax, "fmax")?;
    positive(options.spring, "spring_eV_per_A2")?;
    let restarting = !options.restart.is_empty();
    let mut first = if restarting { options.restart[0].clone() } else { initial.clone() };
    let mut last = if restarting { options.restart.last().unwrap().clone() } else { last_input.clone() };
    if first.is_empty() || first.symbols != last.symbols {
        return Err("NEB requires identical ordered species at both endpoints".into());
    }
    if first.pbc != last.pbc || !same_cell(&first, &last) {
        return Err("Variable-cell NEB is not supported; endpoint cells must match".into());
    }
    for image in &options.restart {
        if image.symbols != first.symbols || !same_cell(image, &first) {
            return Err("Restart images must share the endpoints' species, order and cell".into());
        }
    }
    if restarting && options.restart.len() < 3 {
        return Err("A restart needs at least three images".into());
    }
    let mut endpoint_relaxation = Json::Null;
    if options.relax_endpoints {
        let relax = RelaxOptions { fmax: options.endpoint_fmax, steps: 2000, ..Default::default() };
        let a = relax_configuration(&first, potential, &relax)?;
        let b = relax_configuration(&last, potential, &relax)?;
        first = a.configuration;
        last = b.configuration;
        let mut e = Json::object();
        e.set("initial_converged", a.converged);
        e.set("final_converged", b.converged);
        e.set("initial_steps", a.steps);
        e.set("final_steps", b.steps);
        endpoint_relaxation = e;
    }
    let mic = if options.mic { MinimumImage::new(&first.cell, first.pbc)? } else { MinimumImage::default() };
    let atoms = first.len();
    let mut count = if restarting { options.restart.len() } else { options.images as usize };
    if !restarting && count == 0 {
        positive(options.image_spacing, "image_spacing_A")?;
        let mut squared = 0.0;
        for a in 0..atoms {
            let d = mic.apply(sub(last.positions[a], first.positions[a]));
            squared += dot(d, d);
        }
        count = (((squared.sqrt() / options.image_spacing).ceil() as i64) + 1).clamp(5, 31) as usize;
    }
    if count < 3 {
        return Err("images must be an integer >= 3 (or 0 for automatic)".into());
    }
    let mut chain: Vec<Configuration>;
    if restarting {
        chain = options.restart.clone();
        chain[0] = first.clone();
        *chain.last_mut().unwrap() = last.clone();
    } else {
        chain = vec![first.clone(); count];
        chain[count - 1] = last.clone();
        for a in 0..atoms {
            let delta = scale(mic.apply(sub(last.positions[a], first.positions[a])), 1.0 / (count - 1) as f64);
            for i in 1..count - 1 {
                chain[i].positions[a] = add(first.positions[a], scale(delta, i as f64));
            }
        }
    }
    let mut results = vec![PotentialResult::default(); count];
    results[0] = potential.compute(&chain[0], false)?;
    results[count - 1] = potential.compute(&chain[count - 1], false)?;
    let tangent_vector = |chain: &[Configuration], from: usize, to: usize| -> Vec<Vec3> {
        (0..atoms).map(|a| mic.apply(sub(chain[to].positions[a], chain[from].positions[a]))).collect()
    };
    let vector_norm = |v: &[Vec3]| {
        let mut s = 0.0;
        for x in v {
            s += dot(*x, *x);
        }
        s.sqrt()
    };
    // NEB forces on the movable images, flattened image-by-image.
    let neb_forces = |chain: &[Configuration], results: &mut Vec<PotentialResult>| -> Result<Vec<Vec3>, String> {
        for i in 1..count - 1 {
            results[i] = potential.compute(&chain[i], false)?;
        }
        let mut highest = 1;
        for i in 1..count - 1 {
            if results[i].energy > results[highest].energy {
                highest = i;
            }
        }
        let mut total = Vec::new();
        for i in 1..count - 1 {
            let minus = tangent_vector(chain, i - 1, i);
            let plus = tangent_vector(chain, i, i + 1);
            let (e_minus, e, e_plus) = (results[i - 1].energy, results[i].energy, results[i + 1].energy);
            let mut tangent = vec![[0.0; 3]; atoms];
            if e_plus > e && e > e_minus {
                tangent = plus.clone();
            } else if e_plus < e && e < e_minus {
                tangent = minus.clone();
            } else {
                let high = max((e_plus - e).abs(), (e_minus - e).abs());
                let low = min((e_plus - e).abs(), (e_minus - e).abs());
                for a in 0..atoms {
                    tangent[a] = if e_plus > e_minus { add(scale(plus[a], high), scale(minus[a], low)) } else { add(scale(plus[a], low), scale(minus[a], high)) };
                }
            }
            let length = vector_norm(&tangent);
            if length > 0.0 {
                for t in tangent.iter_mut() {
                    *t = scale(*t, 1.0 / length);
                }
            }
            let mut projection = 0.0;
            for a in 0..atoms {
                projection += dot(results[i].forces[a], tangent[a]);
            }
            let climbing = options.climb && i == highest;
            let spring_force = options.spring * (vector_norm(&plus) - vector_norm(&minus));
            for a in 0..atoms {
                let mut f = sub(results[i].forces[a], scale(tangent[a], if climbing { 2.0 } else { 1.0 } * projection));
                if !climbing {
                    f = add(f, scale(tangent[a], spring_force));
                }
                total.push(f);
            }
        }
        Ok(total)
    };
    let mut fire = Fire::new();
    let mut steps: i64 = 0;
    let mut forces = neb_forces(&chain, &mut results)?;
    let mut converged = max_atom_force(&forces) < options.fmax;
    while !converged && steps < options.steps {
        let step = fire.step(&forces);
        let mut k = 0;
        for image in chain.iter_mut().take(count - 1).skip(1) {
            for a in 0..atoms {
                image.positions[a] = add(image.positions[a], step[k]);
                k += 1;
            }
        }
        steps += 1;
        forces = neb_forces(&chain, &mut results)?;
        converged = max_atom_force(&forces) < options.fmax;
    }
    let energies: Vec<f64> = results.iter().map(|r| r.energy).collect();
    // std::max_element: the first maximum.
    let mut peak = energies[0];
    for &e in &energies[1..] {
        if peak < e {
            peak = e;
        }
    }
    let mut output = ToolOutput::default();
    let mut images = Vec::new();
    for image in &chain {
        output.frames.push(structure_from(image));
        images.push(structure_json(output.frames.last().unwrap()));
    }
    let mut result = Json::object();
    result.set("images", Json::Array(images));
    result.set("energies_eV", jv(&energies));
    let mut coordinate = vec![0.0];
    for i in 1..count {
        coordinate.push(coordinate.last().unwrap() + vector_norm(&tangent_vector(&chain, i - 1, i)));
    }
    result.set("reaction_coordinate_A", jv(&coordinate));
    result.set("forward_barrier_eV", peak - energies[0]);
    result.set("reverse_barrier_eV", peak - energies[count - 1]);
    result.set("reaction_energy_eV", energies[count - 1] - energies[0]);
    result.set("converged", converged);
    result.set("steps", steps);
    result.set("endpoint_max_forces_eV_per_A", jv(&[max_atom_force(&results[0].forces), max_atom_force(&results[count - 1].forces)]));
    result.set("potential", potential.description());
    result.set("image_count", count);
    result.set("restarted", restarting);
    if !endpoint_relaxation.is_null() {
        result.set("endpoint_relaxation", endpoint_relaxation);
    }
    output.result = result;
    Ok(output)
}

// ------------------------------------------------------- tool runners

pub fn dynamics_options(p: &Parameters, npt: bool) -> Result<DynamicsOptions, String> {
    let mut options = DynamicsOptions {
        steps: integer(p.number_or("steps", 1000.0)?, "steps", 1)?,
        timestep_fs: positive(p.number_or("timestep_fs", 1.0)?, "timestep_fs")?,
        temperature_k: positive(p.number_or("temperature_K", 300.0)?, "temperature_K")?,
        thermostat_fs: positive(p.number_or("thermostat_fs", 100.0)?, "thermostat_fs")?,
        seed: integer(p.number_or("seed", 0.0)?, "seed", 0)? as u64,
        sample_interval: integer(p.number_or("sample_interval", 10.0)?, "sample_interval", 1)?,
        production_steps: integer(p.number_or("production_steps", 0.0)?, "production_steps", 0)?,
        equilibration_fs: p.number_or("equilibration_fs", 0.0)?,
        ..Default::default()
    };
    if options.equilibration_fs < 0.0 {
        return Err("equilibration_fs must be nonnegative".into());
    }
    if npt {
        options.pressure_gpa = p.number_or("pressure_GPa", 0.0)?;
        options.barostat_fs = positive(p.number_or("barostat_fs", 1000.0)?, "barostat_fs")?;
    }
    Ok(options)
}

pub fn run_neb(p: &Parameters, base: &std::path::Path) -> Result<ToolOutput, String> {
    let mut options = NebOptions::default();
    if p.has("restart_images") {
        for frame in p.frames("restart_images")? {
            options.restart.push(configuration_from(&frame.structure, &frame.pbc)?);
        }
    } else if !p.has("initial") || !p.has("final") {
        return Err("NEB needs initial and final structures, or restart_images".into());
    }
    let initial = if options.restart.is_empty() {
        let s = p.structure("initial")?;
        configuration_from(&s.structure, &s.pbc)?
    } else {
        options.restart[0].clone()
    };
    let last = if options.restart.is_empty() {
        let s = p.structure("final")?;
        configuration_from(&s.structure, &s.pbc)?
    } else {
        options.restart.last().unwrap().clone()
    };
    options.images = integer(p.number_or("images", 7.0)?, "images", 0)?;
    options.image_spacing = p.number_or("image_spacing_A", 0.5)?;
    options.relax_endpoints = p.boolean("relax_endpoints", false)?;
    options.endpoint_fmax = positive(p.number_or("endpoint_fmax", 0.01)?, "endpoint_fmax")?;
    options.fmax = positive(p.number_or("fmax", 0.03)?, "fmax")?;
    options.steps = integer(p.number_or("steps", 300.0)?, "steps", 1)?;
    options.spring = positive(p.number_or("spring_eV_per_A2", 0.1)?, "spring_eV_per_A2")?;
    options.climb = p.boolean("climb", true)?;
    options.mic = p.boolean("mic", false)?;
    let potential = potential_factory(p.json("calculator_factory")?, base)?;
    migration_path(&initial, &last, potential.as_ref(), &options)
}

pub fn run_relax(p: &Parameters, base: &std::path::Path) -> Result<ToolOutput, String> {
    let structure = p.structure("structure")?;
    let potential = potential_factory(p.json("calculator")?, base)?;
    let options = RelaxOptions {
        fmax: positive(p.number_or("fmax", 0.01)?, "fmax")?,
        steps: integer(p.number_or("steps", 500.0)?, "steps", 0)?,
        relax_cell: p.boolean("relax_cell", false)?,
        pressure_gpa: p.number_or("pressure_GPa", 0.0)?,
    };
    relax_structure(&configuration_from(&structure.structure, &structure.pbc)?, potential.as_ref(), &options)
}

pub fn run_dynamics(p: &Parameters, base: &std::path::Path, npt: bool) -> Result<ToolOutput, String> {
    let structure = p.structure("structure")?;
    let potential = potential_factory(p.json("calculator")?, base)?;
    let configuration = configuration_from(&structure.structure, &structure.pbc)?;
    let options = dynamics_options(p, npt)?;
    if npt {
        npt_dynamics(&configuration, potential.as_ref(), &options)
    } else {
        nvt_dynamics(&configuration, potential.as_ref(), &options)
    }
}
