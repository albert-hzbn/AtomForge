//! Opening and saving structure files (the C++ loadStructureFromFile and
//! saveStructure): format by extension (or a POSCAR/CONTCAR name), the
//! cell wrap and orientation colours after reading, the sidecar after writing.

use std::path::Path;

use af_core::Structure;

use crate::{cif, ipf, lammps, molecular, occupancy, poscar, simulation, writers, xyz};

const SUPPORTED: [&str; 15] =
    [".cif", ".mol", ".pdb", ".xyz", ".sdf", ".vasp", ".mol2", ".pwi", ".gjf", ".extxyz", ".data", ".lmp", ".xsf", ".gen", ".cell"];

pub fn supported_extensions_summary() -> &'static str {
    ".cif, .xyz, .extxyz, .vasp (or POSCAR/CONTCAR), .data/.lmp (LAMMPS data), .pdb, .mol, .mol2, .sdf, .pwi, .gjf, \
     .xsf, .gen (DFTB+), .cell (CASTEP)"
}

/// The extension with its dot ("" when none, or when the last dot is in a directory name).
pub fn extension(filename: &str) -> &str {
    let slash = filename.rfind(['/', '\\']);
    match filename.rfind('.') {
        Some(dot) if slash.map_or(true, |s| dot >= s) => &filename[dot..],
        _ => "",
    }
}

pub fn is_poscar_name(filename: &str) -> bool {
    let base = &filename[filename.rfind(['/', '\\']).map_or(0, |s| s + 1)..];
    let upper = base.to_ascii_uppercase();
    upper.contains("POSCAR") || upper.contains("CONTCAR")
}

pub fn is_supported_structure_file(filename: &str) -> bool {
    SUPPORTED.contains(&extension(filename).to_ascii_lowercase().as_str()) || is_poscar_name(filename)
}

pub fn is_large_data_file(filename: &str) -> bool {
    matches!(extension(filename).to_ascii_lowercase().as_str(), ".afcloud" | ".dump" | ".lammpstrj")
}

fn read_text(filename: &str) -> Result<String, String> {
    let bytes = std::fs::read(filename).map_err(|_| format!("File not found: {filename}"))?;
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// Reads a structure file.
pub fn load_structure_from_file(filename: &str) -> Result<Structure, String> {
    if filename.is_empty() {
        return Err("No file selected.".into());
    }
    let ext = extension(filename);
    let lower = ext.to_ascii_lowercase();
    let supported = SUPPORTED.contains(&lower.as_str());
    let detected_poscar = !supported && is_poscar_name(filename);
    if !supported && !detected_poscar {
        return Err(if lower.is_empty() {
            format!("Unsupported file format (missing extension). Supported formats: {}", supported_extensions_summary())
        } else {
            format!("Unsupported file format '{ext}'. Supported formats: {}", supported_extensions_summary())
        });
    }
    if !Path::new(filename).exists() {
        return Err(format!("File not found: {filename}"));
    }
    let text = read_text(filename)?;
    let mut s = Structure::default();
    match lower.as_str() {
        ".xyz" | ".extxyz" => xyz::read(&text, &mut s)?,
        ".vasp" | ".poscar" => poscar::read(&text, &mut s)?,
        _ if detected_poscar => poscar::read(&text, &mut s)?,
        ".cif" => {
            if !cif::read(&text, &mut s)? {
                return Err("This CIF lists no atom sites (_atom_site_fract_x or _atom_site_Cartn_x).".into());
            }
        }
        ".data" | ".lmp" => lammps::read(&text, &mut s)?,
        ".xsf" | ".gen" | ".cell" | ".pwi" => s = simulation::read_simulation_format(&text, &lower[1..])?,
        ".pdb" => molecular::read_pdb(&text, &mut s)?,
        ".mol" | ".sdf" => molecular::read_mol(&text, &mut s)?,
        ".mol2" => molecular::read_mol2(&text, &mut s)?,
        ".gjf" => molecular::read_gjf(&text, &mut s)?,
        _ => unreachable!(),
    }
    finish_loaded_structure(filename, &mut s)?;
    Ok(s)
}

fn finish_loaded_structure(filename: &str, s: &mut Structure) -> Result<(), String> {
    if s.has_unit_cell {
        occupancy::wrap_atoms_into_primary_cell(s);
    }
    let restored = ipf::restore_sidecar(filename, s);
    let recovered = !restored && ipf::recover_from_geometry(s);
    s.ipf_load_status = if restored {
        "IPF loaded from saved metadata.".into()
    } else if recovered {
        "IPF reconstructed from geometry fallback.".into()
    } else {
        "IPF not found; using non-orientation colors.".into()
    };
    if s.atoms.is_empty() {
        return Err("File loaded but no atoms were found.".into());
    }
    Ok(())
}

/// The text a structure is saved as, in a format (cif, xyz/extxyz, vasp/poscar,
/// data/lmp/lammps/lmpdat, pwi, cell, gen, xsf, cp2k, pdb, mol, sdf, mol2, gjf).
pub fn structure_text(s: &Structure, format: &str) -> Result<String, String> {
    let mut f = format.to_ascii_lowercase();
    if !s.has_unit_cell && (f == "cif" || f == "vasp" || f == "poscar") {
        return Err(format!(
            "{} files need a unit cell, and this structure has none. Save it as XYZ, PDB or MOL, \
             or give it a cell with Edit > Structure Operations > Add Vacuum.",
            if f == "cif" { "CIF" } else { "VASP" }
        ));
    }
    if f == "extxyz" {
        f = "xyz".into();
    }
    if simulation::is_simulation_format(&f) {
        return simulation::write_simulation_format(s, &f);
    }
    match f.as_str() {
        "cif" => writers::write_cif(s),
        "data" | "lmp" | "lammps" | "lmpdat" => writers::write_lammps_data(s),
        "xyz" => Ok(writers::write_xyz(s)),
        "vasp" | "poscar" => writers::write_poscar(s),
        "pdb" | "ent" => Ok(molecular::write_pdb(s)),
        "mol" | "mdl" => Ok(molecular::write_mol(s, false)),
        "sdf" | "sd" => Ok(molecular::write_mol(s, true)),
        "mol2" | "ml2" => Ok(molecular::write_mol2(s)),
        "gjf" | "com" | "gau" => Ok(molecular::write_gjf(s)),
        _ => Err(format!("The {format} format is not available.")),
    }
}

/// Saves a structure (and its orientation colours beside it).
pub fn save_structure(s: &Structure, filename: &str, format: &str) -> Result<(), String> {
    let text = structure_text(s, format)?;
    std::fs::write(filename, text).map_err(|_| format!("Could not write {filename}; check that the folder exists and is writable."))?;
    let f = format.to_ascii_lowercase();
    // The C++ writes no sidecar for the simulation-code formats.
    if !simulation::is_simulation_format(&f) {
        ipf::save_sidecar(s, filename);
    }
    Ok(())
}
