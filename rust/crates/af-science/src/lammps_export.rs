//! LAMMPS input decks (the C++ science/LammpsExport): data file in the LAMMPS
//! frame, pair style for the calculator, and minimize, NVT, NPT or NEB input.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use af_core::format::format_g;
use af_core::json::Json;

use crate::core::*;
use crate::data::{atomic_mass, Parameters};
use crate::tools::ToolOutput;

/// Box parameters of a cell in the LAMMPS frame (a along x, b in the xy plane).
fn lammps_box(c: &Mat3) -> Mat3 {
    let lx = norm(c[0]);
    let ax = scale(c[0], 1.0 / lx);
    let xy = dot(c[1], ax);
    let ly = norm(cross(ax, c[1]));
    let xz = dot(c[2], ax);
    let yz = (dot(c[1], c[2]) - xy * xz) / ly;
    let squared = dot(c[2], c[2]) - xz * xz - yz * yz;
    let lz = (if 0.0 < squared { squared } else { 0.0 }).sqrt();
    [[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]]
}

/// The same lattice with tilts reduced by lattice translations (original frame).
fn reduce_tilts(mut cell: Mat3) -> Mat3 {
    for _ in 0..16 {
        let mut changed = false;
        let mut b = lammps_box(&cell);
        if b[2][1].abs() > 0.5 * b[1][1] + 1e-9 {
            cell[2] = sub(cell[2], scale(cell[1], (b[2][1] / b[1][1]).round()));
            changed = true;
            b = lammps_box(&cell);
        }
        if b[2][0].abs() > 0.5 * b[0][0] + 1e-9 {
            cell[2] = sub(cell[2], scale(cell[0], (b[2][0] / b[0][0]).round()));
            changed = true;
            b = lammps_box(&cell);
        }
        if b[1][0].abs() > 0.5 * b[0][0] + 1e-9 {
            cell[1] = sub(cell[1], scale(cell[0], (b[1][0] / b[0][0]).round()));
            changed = true;
        }
        if !changed {
            break;
        }
    }
    cell
}

/// The LAMMPS box (rows: (lx,0,0), (xy,ly,0), (xz,yz,lz)) of a right-handed cell.
pub fn lammps_cell(input: &Mat3) -> Result<Mat3, String> {
    if determinant(input) <= 0.0 {
        return Err("LAMMPS needs a right-handed cell; reorder or negate a cell vector".into());
    }
    Ok(lammps_box(&reduce_tilts(*input)))
}

fn file_name(path: &str) -> String {
    Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
}

fn g6(x: f64) -> String {
    format_g(x, 6)
}

pub fn lammps_export(p: &Parameters) -> Result<ToolOutput, String> {
    let input = p.structure("structure")?;
    let structure = &input.structure;
    if structure.atoms.is_empty() {
        return Err("The structure has no atoms".into());
    }
    let periodic = structure.has_unit_cell && input.pbc[0] && input.pbc[1] && input.pbc[2];
    let mut species: Vec<String> = Vec::new();
    for atom in &structure.atoms {
        if !species.contains(&atom.symbol) {
            species.push(atom.symbol.clone());
        }
    }
    // Box and coordinates in the LAMMPS frame.
    let mut positions: Vec<Vec3> = Vec::new();
    let mut reduced = [[0.0; 3]; 3]; // lattice basis of the LAMMPS box, in the original frame
    let box_: Mat3;
    if periodic {
        let cell = structure.cell_vectors;
        box_ = lammps_cell(&cell)?;
        reduced = reduce_tilts(cell);
        for atom in &structure.atoms {
            let mut f = fractional([atom.x, atom.y, atom.z], &reduced)?;
            for v in f.iter_mut() {
                *v -= v.floor();
            }
            positions.push(row_times(f, &box_));
        }
    } else {
        let first = &structure.atoms[0];
        let mut high = [first.x, first.y, first.z];
        let mut low = high;
        for atom in &structure.atoms {
            let x = [atom.x, atom.y, atom.z];
            for k in 0..3 {
                if x[k] < low[k] {
                    low[k] = x[k];
                }
                if high[k] < x[k] {
                    high[k] = x[k];
                }
            }
        }
        let pad = 10.0;
        for k in 0..3 {
            low[k] -= pad;
            high[k] += pad;
        }
        box_ = [[high[0] - low[0], 0.0, 0.0], [0.0, high[1] - low[1], 0.0], [0.0, 0.0, high[2] - low[2]]];
        for atom in &structure.atoms {
            positions.push(sub([atom.x, atom.y, atom.z], low));
        }
    }
    // Potential.
    let element_list: String = species.iter().map(|s| format!(" {s}")).collect();
    let mut pair = String::new();
    let mut notes: Vec<String> = Vec::new();
    let mut charges: BTreeMap<String, f64> = BTreeMap::new(); // atom_style charge when non-empty
    let mut extra_files: BTreeMap<String, String> = BTreeMap::new();
    let potential = if p.has("calculator") { p.json("calculator")?.clone() } else { Json::parse(r#"{"potential": "EMT"}"#)? };
    let name = if potential.contains("potential") { potential.at("potential")?.string()?.to_string() } else { "EMT".to_string() };
    let key: String = name.chars().filter(|c| c.is_ascii_alphanumeric()).map(|c| c.to_ascii_lowercase()).collect();
    let number_or = |k: &str, fallback: f64| -> Result<f64, String> { if potential.contains(k) { Ok(potential.at(k)?.number()?) } else { Ok(fallback) } };
    match key.as_str() {
        "eam" | "eamalloy" | "eamfs" | "finnissinclair" => {
            let file = file_name(potential.at("file")?.string()?);
            let mut format = if potential.contains("format") { potential.at("format")?.string()?.to_string() } else { "auto".to_string() };
            if format == "auto" {
                let lower = file.to_ascii_lowercase();
                format = if lower.contains(".fs") {
                    "fs"
                } else if lower.contains("alloy") || lower.contains("setfl") {
                    "setfl"
                } else {
                    "funcfl"
                }
                .to_string();
            }
            if format == "funcfl" {
                let _ = writeln!(pair, "pair_style eam\npair_coeff 1 1 {file}");
            } else {
                let _ = writeln!(pair, "pair_style eam/{}\npair_coeff * * {file}{element_list}", if format == "fs" { "fs" } else { "alloy" });
            }
            notes.push(format!("Copy {file} beside in.lammps."));
        }
        "lennardjones" | "lj" => {
            let sigma = number_or("sigma", 1.0)?;
            let epsilon = number_or("epsilon", 1.0)?;
            let cutoff = number_or("cutoff", 3.0 * sigma)?;
            let _ = writeln!(pair, "pair_style lj/cut {}\npair_coeff * * {} {}\npair_modify shift yes", g6(cutoff), g6(epsilon), g6(sigma));
        }
        "tersoff" | "stillingerweber" | "sw" => {
            let tersoff = key == "tersoff";
            let mut file = if potential.contains("file") && !potential.at("file")?.string()?.is_empty() {
                file_name(potential.at("file")?.string()?)
            } else {
                String::new()
            };
            if file.is_empty() {
                // The built-in Si parameters, in LAMMPS's own file format.
                file = if tersoff { "Si.tersoff" } else { "Si.sw" }.to_string();
                extra_files.insert(
                    file.clone(),
                    if tersoff {
                        "# Si, Tersoff PRB 38 9902 (1988)\nSi Si Si 3.0 1.0 0.0 1.0039e5 16.217 -0.59825 0.78734 1.1e-6 1.7322 471.18 2.85 0.15 2.4799 1830.8\n"
                    } else {
                        "# Si, Stillinger and Weber PRB 31 5262 (1985)\nSi Si Si 2.1683 2.0951 1.80 21.0 1.20 -0.333333333333 7.049556277 0.6022245584 4.0 0.0 0.0\n"
                    }
                    .to_string(),
                );
            } else {
                notes.push(format!("Copy {file} beside in.lammps."));
            }
            let _ = writeln!(pair, "pair_style {}\npair_coeff * * {file}{element_list}", if tersoff { "tersoff" } else { "sw" });
        }
        "buckingham" | "buck" => {
            let cutoff = number_or("cutoff", 10.0)?;
            if let Some(table) = potential.find("charges").filter(|t| t.is_object()) {
                for (element, q) in table.members()? {
                    charges.insert(element.clone(), q.number()?);
                }
            }
            let _ = writeln!(pair, "{}{}", if charges.is_empty() { "pair_style buck " } else { "pair_style buck/coul/long " }, g6(cutoff));
            // Unlisted pairs interact only through their charges.
            for a in 0..species.len() {
                for b in a..species.len() {
                    let (mut big_a, mut rho, mut c) = (0.0, 1.0, 0.0);
                    if let Some(list) = potential.find("pairs").filter(|l| l.is_array()) {
                        for entry in list.items()? {
                            let e = entry.at("elements")?.items()?;
                            let (e0, e1) = (e[0].string()?, e[1].string()?);
                            if (e0 == species[a] && e1 == species[b]) || (e0 == species[b] && e1 == species[a]) {
                                big_a = entry.at("A")?.number()?;
                                rho = entry.at("rho")?.number()?;
                                c = if entry.contains("C") { entry.at("C")?.number()? } else { 0.0 };
                            }
                        }
                    }
                    let _ = writeln!(pair, "pair_coeff {} {} {} {} {}  # {}-{}", a + 1, b + 1, g6(big_a), g6(rho), g6(c), species[a], species[b]);
                }
            }
            if !charges.is_empty() {
                if !periodic {
                    return Err("LAMMPS export of Coulomb interactions needs a periodic structure".into());
                }
                let _ = writeln!(pair, "kspace_style ewald {}", g6(number_or("ewald_accuracy", 1e-6)?));
            }
        }
        _ => {
            let _ = writeln!(pair, "# EMT has no LAMMPS pair style: choose a potential for{element_list}, e.g.\n# pair_style eam/alloy\n# pair_coeff * * FILE.eam.alloy{element_list}");
            notes.push("EMT is not available in LAMMPS; the pair style is left as a commented placeholder.".into());
        }
    }
    let g12 = |x: f64| format_g(x, 12);
    let triclinic = periodic && (box_[1][0].abs() > 1e-12 || box_[2][0].abs() > 1e-12 || box_[2][1].abs() > 1e-12);
    let mut data = String::new();
    let _ = write!(data, "LAMMPS data file written by AtomForge\n\n{} atoms\n{} atom types\n\n", structure.atoms.len(), species.len());
    let _ = write!(data, "0.0 {} xlo xhi\n0.0 {} ylo yhi\n0.0 {} zlo zhi\n", g12(box_[0][0]), g12(box_[1][1]), g12(box_[2][2]));
    if triclinic {
        let _ = writeln!(data, "{} {} {} xy xz yz", g12(box_[1][0]), g12(box_[2][0]), g12(box_[2][1]));
    }
    data.push_str("\nMasses\n\n");
    for (t, s) in species.iter().enumerate() {
        let _ = writeln!(data, "{} {}  # {s}", t + 1, g12(atomic_mass(s)?));
    }
    data.push_str(if charges.is_empty() { "\nAtoms  # atomic\n\n" } else { "\nAtoms  # charge\n\n" });
    for (i, atom) in structure.atoms.iter().enumerate() {
        let type_ = species.iter().position(|s| s == &atom.symbol).unwrap() + 1;
        let _ = write!(data, "{} {type_} ", i + 1);
        if !charges.is_empty() {
            let Some(charge) = charges.get(&atom.symbol) else {
                return Err(format!("No charge given for {}", atom.symbol));
            };
            let _ = write!(data, "{} ", g12(*charge));
        }
        let _ = writeln!(data, "{} {} {}", g12(positions[i][0]), g12(positions[i][1]), g12(positions[i][2]));
    }

    let task = if p.has("task") { p.json("task")?.string()?.to_string() } else { "minimize".to_string() };
    let temperature = p.number_or("temperature_K", 300.0)?;
    let pressure_bar = p.number_or("pressure_GPa", 0.0)? * 10000.0;
    let timestep_ps = positive(p.number_or("timestep_fs", 1.0)?, "timestep_fs")? / 1000.0;
    let steps = integer(p.number_or("steps", 10000.0)?, "steps", 1)?;
    let interval = integer(p.number_or("sample_interval", 100.0)?, "sample_interval", 1)?;
    let thermostat = positive(p.number_or("thermostat_fs", 100.0)?, "thermostat_fs")? / 1000.0;
    let barostat = positive(p.number_or("barostat_fs", 1000.0)?, "barostat_fs")? / 1000.0;
    let seed = integer(p.number_or("seed", 12345.0)?, "seed", 1)?;
    let g10 = |x: f64| format_g(x, 10);
    let mut inp = String::new();
    let _ = write!(
        inp,
        "# LAMMPS input written by AtomForge ({task})\nunits metal\natom_style {}\nboundary {}\n",
        if charges.is_empty() { "atomic" } else { "charge" },
        if periodic { "p p p" } else { "f f f" }
    );
    if task == "neb" {
        inp.push_str("atom_modify map array\n");
    }
    let _ = write!(inp, "read_data data.lammps\n\n{pair}\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nthermo {interval}\n");
    inp.push_str("thermo_style custom step temp pe etotal press vol\n");
    // Frames with element names are read back by AtomForge's trajectory playback.
    let dump = format!("dump traj all custom {interval} dump.lammpstrj id element x y z\ndump_modify traj element{element_list} sort id\n");
    let mut files = Json::object();
    match task.as_str() {
        "minimize" => {
            inp.push_str(&dump);
            inp.push_str("min_style cg\n");
            if p.boolean("relax_cell", false)? {
                if !periodic {
                    return Err("Cell relaxation needs a periodic structure".into());
                }
                let _ = writeln!(inp, "fix relax all box/relax tri {} vmax 0.001", g10(pressure_bar));
            }
            let _ = write!(inp, "minimize 0.0 {} {steps} {}\nwrite_data relaxed.data\n", g10(p.number_or("fmax", 0.01)?), 10 * steps);
        }
        "nvt" | "npt" => {
            let _ = write!(inp, "velocity all create {} {seed} dist gaussian\ntimestep {}\n{dump}", g10(temperature), g10(timestep_ps));
            if task == "nvt" {
                let _ = writeln!(inp, "fix md all nvt temp {} {} {}", g10(temperature), g10(temperature), g10(thermostat));
            } else {
                if !periodic {
                    return Err("NPT needs a periodic structure".into());
                }
                let _ = writeln!(
                    inp,
                    "fix md all npt temp {} {} {} iso {} {} {}",
                    g10(temperature),
                    g10(temperature),
                    g10(thermostat),
                    g10(pressure_bar),
                    g10(pressure_bar),
                    g10(barostat)
                );
            }
            let _ = writeln!(inp, "run {steps}");
        }
        "neb" => {
            if !p.has("final") {
                return Err("NEB export needs the final structure".into());
            }
            let last = &p.structure("final")?.structure;
            if last.atoms.len() != structure.atoms.len() {
                return Err("NEB endpoints need the same atoms".into());
            }
            let mut neb = format!("{}\n", last.atoms.len());
            // Final coordinates in the same (LAMMPS) frame: apply the initial structure's mapping.
            for (i, (a, s0)) in last.atoms.iter().zip(&structure.atoms).enumerate() {
                // Displacements rotate with the cell into the LAMMPS frame.
                let shift = sub([a.x, a.y, a.z], [s0.x, s0.y, s0.z]);
                let x = if periodic { add(positions[i], row_times(row_times(shift, &inverse(&reduced)?), &box_)) } else { add(positions[i], shift) };
                let _ = writeln!(neb, "{} {} {} {}", i + 1, g12(x[0]), g12(x[1]), g12(x[2]));
            }
            files.set("final.neb", neb);
            let images = integer(p.number_or("images", 7.0)?, "images", 3)?;
            let _ = write!(
                inp,
                "variable u uloop {images}\nfix neb all neb {}\nmin_style quickmin\ntimestep 0.01\n",
                g10(p.number_or("spring_eV_per_A2", 0.1)?)
            );
            let _ = writeln!(inp, "neb 0.0 {} {steps} {steps} {interval} final final.neb", g10(p.number_or("fmax", 0.03)?));
            notes.push(format!("Run NEB with one partition per image, e.g. mpirun -np {images} lmp -partition {images}x1 -in in.lammps"));
        }
        _ => return Err("task must be minimize, nvt, npt or neb".into()),
    }
    for (file, text) in extra_files {
        files.set(&file, text);
    }
    files.set("data.lammps", data);
    files.set("in.lammps", inp);
    let mut result = Json::object();
    result.set("task", task);
    result.set("atom_types", Json::Array(species.iter().map(|s| Json::from(s.as_str())).collect()));
    result.set("box", mat3_json(&box_));
    result.set("triclinic", triclinic);
    result.set("notes", Json::Array(notes.into_iter().map(Json::from).collect()));
    result.set("files", files);
    Ok(ToolOutput { result, ..Default::default() })
}
