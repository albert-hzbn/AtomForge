//! Harmonic phonons by finite displacements (the C++ science/Phonons): force
//! constants in a supercell, frequencies along the HPKOT path or given
//! q-points, a mesh, the DOS and harmonic thermodynamics.
//!
//! The dynamical matrix is diagonalised with nalgebra's Hermitian solver, not
//! Eigen's SelfAdjointEigenSolver<MatrixXcd>, so frequencies agree with the
//! C++ to rounding; near-zero acoustic modes (sqrt of rounding noise) agree
//! only absolutely (docs/rust-port-plan.md 6a).

use std::f64::consts::PI;
use std::path::Path;

use af_core::json::Json;
use nalgebra::{Complex, DMatrix};

use crate::analysis::{harmonic_thermodynamics, jv, phonon_dos};
use crate::core::*;
use crate::data::{Parameters, StructureInput};
use crate::potentials::{configuration_from, structure_from, Configuration, Potential};
use crate::reciprocal_path::reciprocal_path;
use crate::tools::ToolOutput;

/// eV/(A^2 amu) -> 1/fs^2.
const ACCELERATION: f64 = 1.602176634e-19 / (1e-10 * 1.66053906660e-27) * 1e10 * 1e-30;
/// h * 1 THz in eV.
pub const THZ_TO_EV: f64 = 4.135667696e-3;

fn reciprocal(cell: &Mat3) -> Result<Mat3, String> {
    let inv = inverse(cell)?;
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| 2.0 * PI * inv[j][i])))
}

fn frequency_from_eigenvalue(lambda: f64) -> f64 {
    let omega = (lambda.abs() * ACCELERATION).sqrt(); // rad/fs
    (if lambda < 0.0 { -1.0 } else { 1.0 }) * omega / (2.0 * PI) * 1000.0 // THz
}

/// Repetitions giving each supercell plane spacing at least `minimum_length`.
pub fn supercell_for(unit: &Configuration, minimum_length: f64) -> Result<[i32; 3], String> {
    positive(minimum_length, "supercell_min_A")?;
    let inv = inverse(&unit.cell)?;
    Ok(std::array::from_fn(|k| {
        let spacing = 1.0 / norm([inv[0][k], inv[1][k], inv[2][k]]);
        1.max((minimum_length / spacing - 1e-9).ceil() as i32)
    }))
}

pub struct ForceConstants {
    pub unit: Configuration,
    pub supercell: [i32; 3],
    pub sup: Configuration,
    /// phi[a][b]: d^2 E / d u_a d u_b for unit-cell atom a and supercell atom b.
    pub phi: Vec<Vec<Mat3>>,
}

pub fn force_constants(unit: &Configuration, potential: &dyn Potential, supercell: [i32; 3], displacement: f64) -> Result<ForceConstants, String> {
    if !unit.fully_periodic() || determinant(&unit.cell).abs() < 1e-12 {
        return Err("Phonons require a full periodic cell".into());
    }
    positive(displacement, "displacement_A")?;
    for n in supercell {
        if !(1..=50).contains(&n) {
            return Err("Supercell repetitions must be between 1 and 50".into());
        }
    }
    let mut sup = Configuration { pbc: [true; 3], ..Default::default() };
    for r in 0..3 {
        sup.cell[r] = scale(unit.cell[r], supercell[r] as f64);
    }
    for i in 0..supercell[0] {
        for j in 0..supercell[1] {
            for k in 0..supercell[2] {
                let shift = row_times([i as f64, j as f64, k as f64], &unit.cell);
                for a in 0..unit.len() {
                    sup.symbols.push(unit.symbols[a].clone());
                    sup.numbers.push(unit.numbers[a]);
                    sup.masses.push(unit.masses[a]);
                    sup.positions.push(add(unit.positions[a], shift));
                }
            }
        }
    }
    let (n, total) = (unit.len(), sup.len());
    let mut phi = vec![vec![[[0.0; 3]; 3]; total]; n];
    for a in 0..n {
        for alpha in 0..3 {
            let mut plus = sup.clone();
            let mut minus = sup.clone();
            plus.positions[a][alpha] += displacement;
            minus.positions[a][alpha] -= displacement;
            let f_plus = potential.compute(&plus, false)?.forces;
            let f_minus = potential.compute(&minus, false)?.forces;
            for b in 0..total {
                for beta in 0..3 {
                    phi[a][b][alpha][beta] = -(f_plus[b][beta] - f_minus[b][beta]) / (2.0 * displacement);
                }
            }
        }
    }
    // Acoustic sum rule: a rigid translation produces no force.
    for a in 0..n {
        for alpha in 0..3 {
            for beta in 0..3 {
                let mut sum = 0.0;
                for b in 0..total {
                    sum += phi[a][b][alpha][beta];
                }
                phi[a][a][alpha][beta] -= sum;
            }
        }
    }
    Ok(ForceConstants { unit: unit.clone(), supercell, sup, phi })
}

/// Frequencies (THz, imaginary as negative) at the Cartesian wave vector q (1/Angstrom, 2 pi included).
pub fn phonon_frequencies(constants: &ForceConstants, q: Vec3) -> Vec<f64> {
    let (n, total) = (constants.unit.len(), constants.sup.len());
    let super_cell = &constants.sup.cell;
    let mut dynamical = DMatrix::<Complex<f64>>::zeros(3 * n, 3 * n);
    for a in 0..n {
        for bs in 0..total {
            let b = bs % n;
            let base = sub(constants.sup.positions[bs], constants.unit.positions[a]);
            // Minimum-image vectors (shared equally between equidistant images).
            let mut images: Vec<Vec3> = Vec::new();
            let mut best = f64::INFINITY;
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let d = add(base, row_times([i as f64, j as f64, k as f64], super_cell));
                        let length = norm(d);
                        if length < best - 1e-5 {
                            best = length;
                            images = vec![d];
                        } else if (length - best).abs() <= 1e-5 {
                            images.push(d);
                        }
                    }
                }
            }
            let (mut re, mut im) = (0.0, 0.0);
            for d in &images {
                let angle = dot(q, *d);
                re += angle.cos();
                im += angle.sin();
            }
            let count = images.len() as f64;
            let phase = Complex::new(re / count, im / count);
            let mass_factor = 1.0 / (constants.unit.masses[a] * constants.unit.masses[b]).sqrt();
            let f = &constants.phi[a][bs];
            for alpha in 0..3 {
                for beta in 0..3 {
                    let w = f[alpha][beta] * mass_factor;
                    dynamical[(3 * a + alpha, 3 * b + beta)] += Complex::new(w * phase.re, w * phase.im);
                }
            }
        }
    }
    let hermitian = (&dynamical + dynamical.adjoint()) / Complex::new(2.0, 0.0);
    let mut eigenvalues: Vec<f64> = hermitian.symmetric_eigenvalues().iter().copied().collect();
    eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    eigenvalues.into_iter().map(frequency_from_eigenvalue).collect()
}

pub fn phonon_calculation(structure: &StructureInput, potential: &dyn Potential, p: &Parameters) -> Result<ToolOutput, String> {
    let mut result = Json::object();
    let mut unit = configuration_from(&structure.structure, &structure.pbc)?;
    let mut path_fractional: Vec<Vec3> = Vec::new();
    let mut path_distance: Vec<f64> = Vec::new();
    let mut labels = Json::Array(Vec::new());
    if p.boolean("symmetry_path", true)? {
        // Phonons are computed in the standardized primitive cell so the
        // HPKOT k-points refer to the same reciprocal basis.
        let symprec = p.number_or("symprec_A", 1e-5)?;
        let spacing = p.number_or("spacing_inv_A", 0.05)?;
        let path = reciprocal_path(structure, spacing, symprec, true)?;
        unit = configuration_from(&path.frames[0], &[true; 3])?;
        let k = to_array(path.result.at("kpoints_fractional")?, "kpoints")?;
        for i in 0..k.shape[0] {
            path_fractional.push([k.at2(i, 0), k.at2(i, 1), k.at2(i, 2)]);
        }
        path_distance = to_array(path.result.at("distance_inv_A")?, "distance")?.values;
        labels = path.result.at("labels")?.clone();
        result.set("spacegroup_symbol", path.result.at("spacegroup_symbol")?.clone());
        result.set("path", path.result.at("path")?.clone());
        result.set("primitive_structure", path.result.at("primitive_structure")?.clone());
    } else if p.has("q_points_fractional") {
        for q in points(p.array("q_points_fractional")?, "q_points_fractional")? {
            if let Some(&previous) = path_fractional.last() {
                let b = reciprocal(&unit.cell)?;
                let next = path_distance.last().unwrap() + norm(row_times(sub(q, previous), &b));
                path_distance.push(next);
            } else {
                path_distance.push(0.0);
            }
            path_fractional.push(q);
            labels.push("");
        }
    }
    let supercell = if p.has("supercell") {
        let values = finite_array(p.array("supercell")?, "supercell", 1)?.values;
        if values.len() != 3 {
            return Err("supercell needs three repetitions".into());
        }
        let mut s = [0i32; 3];
        for k in 0..3 {
            s[k] = integer(values[k], "supercell", 1)? as i32;
        }
        s
    } else {
        supercell_for(&unit, p.number_or("supercell_min_A", 12.0)?)?
    };
    let constants = force_constants(&unit, potential, supercell, p.number_or("displacement_A", 0.01)?)?;
    let b = reciprocal(&unit.cell)?;

    let bands: Vec<Json> = path_fractional.iter().map(|q| jv(&phonon_frequencies(&constants, row_times(*q, &b)))).collect();
    let mut mesh = [8i32, 8, 8];
    if p.has("mesh") {
        let values = finite_array(p.array("mesh")?, "mesh", 1)?.values;
        if values.len() != 3 {
            return Err("mesh needs three sizes".into());
        }
        for k in 0..3 {
            mesh[k] = integer(values[k], "mesh", 1)? as i32;
        }
    }
    let mut mesh_energies = Vec::new();
    let (mut lowest, mut highest) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut imaginary = 0i32;
    for i in 0..mesh[0] {
        for j in 0..mesh[1] {
            for k in 0..mesh[2] {
                let q = [i as f64 / mesh[0] as f64, j as f64 / mesh[1] as f64, k as f64 / mesh[2] as f64];
                let mut row = Vec::new();
                for f in phonon_frequencies(&constants, row_times(q, &b)) {
                    if f < lowest {
                        lowest = f;
                    }
                    if highest < f {
                        highest = f;
                    }
                    imaginary += (f < -0.05) as i32;
                    row.push(f * THZ_TO_EV);
                }
                mesh_energies.push(jv(&row));
            }
        }
    }
    let mesh_energies = Json::Array(mesh_energies);
    result.set("supercell", Json::Array(supercell.iter().map(|&v| Json::from(v)).collect()));
    result.set("unit_cell_atoms", unit.len());
    result.set("displacement_A", p.number_or("displacement_A", 0.01)?);
    if !path_fractional.is_empty() {
        result.set("q_points_fractional", Json::Array(path_fractional.iter().map(|k| jv(k)).collect()));
        result.set("distance_inv_A", jv(&path_distance));
        result.set("labels", labels);
        result.set("frequencies_THz", Json::Array(bands));
    }
    result.set("mesh", Json::Array(mesh.iter().map(|&v| Json::from(v)).collect()));
    result.set("mesh_energies_eV", mesh_energies.clone());
    result.set("min_frequency_THz", lowest);
    result.set("max_frequency_THz", highest);
    result.set("imaginary_mesh_modes", imaginary);

    // Density of states and thermodynamics reuse the tabulated-mode tools.
    let sigma = positive(p.number_or("sigma_eV", 0.001)?, "sigma_eV")?;
    let from = (if lowest * THZ_TO_EV < 0.0 { lowest * THZ_TO_EV } else { 0.0 }) - 5.0 * sigma;
    let to = highest * THZ_TO_EV + 5.0 * sigma;
    let grid: Vec<Json> = (0..=1000).map(|i| Json::from(from + (to - from) * i as f64 / 1000.0)).collect();
    let mut dos_request = Json::object();
    dos_request.set("energies_eV", mesh_energies.clone());
    dos_request.set("energy_grid_eV", Json::Array(grid));
    dos_request.set("sigma_eV", sigma);
    result.set("dos", phonon_dos(&Parameters::new("phonon-dos", &dos_request, Path::new("."), None)?)?);
    if imaginary == 0 {
        let mut thermo_request = Json::object();
        thermo_request.set("energies_eV", mesh_energies);
        if p.has("temperatures_K") {
            thermo_request.set("temperatures_K", array_json(&p.array("temperatures_K")?));
        } else {
            thermo_request.set("temperatures_K", Json::Array((0..=1000).step_by(25).map(|t: i32| Json::from(t)).collect()));
        }
        // Acoustic modes at Gamma are zero only to finite-difference accuracy.
        thermo_request.set("zero_tolerance_eV", 2e-4);
        result.set("thermodynamics", harmonic_thermodynamics(&Parameters::new("harmonic-thermodynamics", &thermo_request, Path::new("."), None)?)?);
    } else {
        result.set("thermodynamics", Json::Null);
        result.set(
            "thermodynamics_note",
            "Imaginary modes on the mesh: the structure is dynamically unstable with this potential; harmonic thermodynamics is not defined.",
        );
    }
    result.set("potential", potential.description());
    Ok(ToolOutput { result, frames: vec![structure_from(&unit)], ..Default::default() })
}
