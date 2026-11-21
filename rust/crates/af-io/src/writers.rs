//! The native writers of the C++ saveStructure: CIF (P1, fractional
//! coordinates to twelve decimals), LAMMPS data, (extended) XYZ and POSCAR.

use af_core::format::format_g;
use af_core::{elements, Structure};
use glam::{DMat3, DVec3};

fn clean(v: f64) -> f64 {
    if v.abs() < 5e-11 {
        0.0
    } else {
        v
    }
}

fn cell_of(s: &Structure) -> DMat3 {
    let [a, b, c] = s.cell_vectors.map(DVec3::from);
    DMat3::from_cols(a, b, c)
}

fn symbol_of(a: &af_core::AtomSite) -> String {
    if a.symbol.is_empty() {
        elements::symbol(a.atomic_number).to_string()
    } else {
        a.symbol.clone()
    }
}

pub fn write_cif(s: &Structure) -> Result<String, String> {
    let cell = cell_of(s);
    if cell.determinant().abs() < 1e-9 {
        return Err("The unit cell has no volume, so it cannot be written as CIF.".into());
    }
    let to_fractional = cell.inverse();
    let (la, lb, lc) = (cell.col(0).length(), cell.col(1).length(), cell.col(2).length());
    let angle = |u: DVec3, v: DVec3| (u.dot(v) / (u.length() * v.length())).clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI;
    let mut out = String::from("# CIF file written by AtomForge\ndata_AtomForge\n");
    out.push_str(&format!(
        "_cell_length_a {}\n_cell_length_b {}\n_cell_length_c {}\n_cell_angle_alpha {}\n_cell_angle_beta {}\n_cell_angle_gamma {}\n",
        format_g(la, 10),
        format_g(lb, 10),
        format_g(lc, 10),
        format_g(angle(cell.col(1), cell.col(2)), 10),
        format_g(angle(cell.col(0), cell.col(2)), 10),
        format_g(angle(cell.col(0), cell.col(1)), 10)
    ));
    out.push_str(
        "_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\nloop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n\
         loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n",
    );
    let mut serial: std::collections::BTreeMap<String, i32> = Default::default();
    for site in &s.atoms {
        let symbol = symbol_of(site);
        let f = to_fractional * DVec3::new(site.x, site.y, site.z);
        let n = serial.entry(symbol.clone()).or_default();
        *n += 1;
        let label = format!("{symbol}{n}");
        out.push_str(&format!("{label:<7} {symbol:<3} {:>16.12} {:>16.12} {:>16.12} 1\n", clean(f.x), clean(f.y), clean(f.z)));
    }
    Ok(out)
}

pub fn write_lammps_data(s: &Structure) -> Result<String, String> {
    let (cell, origin) = if s.has_unit_cell {
        (cell_of(s), DVec3::from(s.cell_offset))
    } else {
        let mut low = DVec3::ZERO;
        let mut high = DVec3::ZERO;
        for (i, a) in s.atoms.iter().enumerate() {
            let q = DVec3::new(a.x, a.y, a.z);
            low = if i > 0 { low.min(q) } else { q };
            high = if i > 0 { high.max(q) } else { q };
        }
        let size = high - low + DVec3::splat(20.0);
        (DMat3::from_diagonal(size), low - DVec3::splat(10.0))
    };
    let det = cell.determinant();
    if det.abs() < 1e-9 {
        return Err("The unit cell has no volume, so it cannot be written as a LAMMPS box.".into());
    }
    let (a, b, c) = (cell.col(0), cell.col(1), cell.col(2));
    let ax = a.length();
    let ahat = a / ax;
    let bx = b.dot(ahat);
    let by = ahat.cross(b).length();
    let cx = c.dot(ahat);
    let cy = (b.dot(c) - bx * cx) / by;
    let cz2 = c.dot(c) - cx * cx - cy * cy;
    let cz = cz2.max(0.0).sqrt() * if det < 0.0 { -1.0 } else { 1.0 };
    let lammps = DMat3::from_cols(DVec3::new(ax, 0.0, 0.0), DVec3::new(bx, by, 0.0), DVec3::new(cx, cy, cz));
    let to_fractional = cell.inverse();
    if cz < 0.0 {
        return Err("The cell is left-handed; LAMMPS needs a right-handed box (swap two cell vectors first).".into());
    }
    let mut kinds: Vec<i32> = Vec::new();
    for site in &s.atoms {
        if !kinds.contains(&site.atomic_number) {
            kinds.push(site.atomic_number);
        }
    }
    kinds.sort_unstable();
    let mut out = format!("LAMMPS data file written by AtomForge\n\n{} atoms\n{} atom types\n\n", s.atoms.len(), kinds.len());
    out.push_str(&format!(
        "0.0 {} xlo xhi\n0.0 {} ylo yhi\n0.0 {} zlo zhi\n",
        format_g(clean(ax), 10),
        format_g(clean(by), 10),
        format_g(clean(cz), 10)
    ));
    if bx.abs() > 1e-10 || cx.abs() > 1e-10 || cy.abs() > 1e-10 {
        out.push_str(&format!("{} {} {} xy xz yz\n", format_g(clean(bx), 10), format_g(clean(cx), 10), format_g(clean(cy), 10)));
    }
    out.push_str("\nMasses\n\n");
    for (t, &z) in kinds.iter().enumerate() {
        out.push_str(&format!("{} {}  # {}\n", t + 1, format_g(elements::atomic_mass(z), 6), elements::symbol(z)));
    }
    let charges = s.atoms.iter().any(|a| a.charge != 0.0);
    let velocities = s.atoms.iter().any(|a| a.velocity.iter().any(|&v| v != 0.0));
    out.push_str(if charges { "\nAtoms  # charge\n\n" } else { "\nAtoms  # atomic\n\n" });
    for (i, site) in s.atoms.iter().enumerate() {
        let f = to_fractional * (DVec3::new(site.x, site.y, site.z) - origin);
        let q = lammps * f;
        let t = kinds.iter().position(|&z| z == site.atomic_number).unwrap() + 1;
        if charges {
            out.push_str(&format!("{} {} {} {:.10} {:.10} {:.10}\n", i + 1, t, format_g(site.charge, 10), clean(q.x), clean(q.y), clean(q.z)));
        } else {
            out.push_str(&format!("{} {} {:.10} {:.10} {:.10}\n", i + 1, t, clean(q.x), clean(q.y), clean(q.z)));
        }
    }
    if velocities {
        out.push_str("\nVelocities\n\n");
        let rotate = lammps * to_fractional;
        for (i, site) in s.atoms.iter().enumerate() {
            let v = rotate * DVec3::from(site.velocity) * 1000.0;
            out.push_str(&format!("{} {} {} {}\n", i + 1, format_g(clean(v.x), 10), format_g(clean(v.y), 10), format_g(clean(v.z), 10)));
        }
    }
    Ok(out)
}

pub fn write_xyz(s: &Structure) -> String {
    let moments = s.atoms.iter().any(|a| a.magmom != 0.0);
    let charges = s.atoms.iter().any(|a| a.charge != 0.0);
    let velocities = s.atoms.iter().any(|a| a.velocity.iter().any(|&v| v != 0.0));
    let mut properties = String::from("species:S:1:pos:R:3");
    if moments {
        properties.push_str(":initial_magmoms:R:1");
    }
    if charges {
        properties.push_str(":initial_charges:R:1");
    }
    if velocities {
        properties.push_str(":velo:R:3");
    }
    let flag = |p: bool| if p { 'T' } else { 'F' };
    let mut out = format!("{}\n", s.atoms.len());
    if s.has_unit_cell {
        out.push_str("Lattice=\"");
        for r in 0..3 {
            for c in 0..3 {
                if r > 0 || c > 0 {
                    out.push(' ');
                }
                out.push_str(&format_g(clean(s.cell_vectors[r][c]), 10));
            }
        }
        out.push_str(&format!("\" Properties={properties} pbc=\"{} {} {}\"\n", flag(s.pbc[0]), flag(s.pbc[1]), flag(s.pbc[2])));
    } else if moments || charges || velocities {
        out.push_str(&format!("Properties={properties} pbc=\"F F F\"\n"));
    } else {
        out.push_str("AtomForge structure\n");
    }
    for site in &s.atoms {
        out.push_str(&format!("{:<3} {:>16.10} {:>16.10} {:>16.10}", symbol_of(site), clean(site.x), clean(site.y), clean(site.z)));
        if moments {
            out.push_str(&format!(" {}", format_g(site.magmom, 10)));
        }
        if charges {
            out.push_str(&format!(" {}", format_g(site.charge, 10)));
        }
        if velocities {
            out.push_str(&format!(" {} {} {}", format_g(site.velocity[0], 10), format_g(site.velocity[1], 10), format_g(site.velocity[2], 10)));
        }
        out.push('\n');
    }
    out
}

/// POSCAR (VASP 5 with the species also on the title line), atoms grouped by
/// element in increasing atomic number, Direct coordinates.
pub fn write_poscar(s: &Structure) -> Result<String, String> {
    let cell = cell_of(s);
    if cell.determinant().abs() < 1e-9 {
        return Err("The unit cell has no volume, so it cannot be written as POSCAR.".into());
    }
    let to_fractional = cell.inverse();
    let mut kinds: Vec<i32> = Vec::new();
    for site in &s.atoms {
        if !kinds.contains(&site.atomic_number) {
            kinds.push(site.atomic_number);
        }
    }
    kinds.sort_unstable();
    let symbols: Vec<&str> = kinds.iter().map(|&z| elements::symbol(z)).collect();
    let counts: Vec<String> = kinds.iter().map(|&z| s.atoms.iter().filter(|a| a.atomic_number == z).count().to_string()).collect();
    let mut out = format!("{}\n1.0\n", symbols.join(" "));
    for r in 0..3 {
        out.push_str(&format!(
            "  {:>20.12} {:>20.12} {:>20.12}\n",
            clean(s.cell_vectors[r][0]),
            clean(s.cell_vectors[r][1]),
            clean(s.cell_vectors[r][2])
        ));
    }
    let selective = s.atoms.iter().any(|a| a.fixed != 0);
    out.push_str(&format!("{}\n{}\n{}Direct\n", symbols.join(" "), counts.join(" "), if selective { "Selective dynamics\n" } else { "" }));
    for &z in &kinds {
        for site in s.atoms.iter().filter(|a| a.atomic_number == z) {
            let f = to_fractional * DVec3::new(site.x, site.y, site.z);
            out.push_str(&format!("  {:>16.12} {:>16.12} {:>16.12}", clean(f.x), clean(f.y), clean(f.z)));
            if selective {
                for k in 0..3 {
                    out.push_str(if (site.fixed >> k) & 1 != 0 { " F" } else { " T" });
                }
            }
            out.push('\n');
        }
    }
    Ok(out)
}
