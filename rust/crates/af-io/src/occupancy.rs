//! Wrapping atoms into the cell and resolving partially occupied sites (the
//! C++ wrapAtomsIntoPrimaryCell and resolvePartialOccupancy).

use std::collections::{BTreeMap, BTreeSet};

use af_core::cell::cell_matrices;
use af_core::Structure;
use glam::{DMat3, DVec3, IVec3};

/// Moves atoms outside the cell by whole cell vectors into it; atoms inside
/// keep their coordinates exactly.
pub fn wrap_atoms_into_primary_cell(structure: &mut Structure) {
    if !structure.has_unit_cell || structure.atoms.is_empty() {
        return;
    }
    let origin = DVec3::from(structure.cell_offset);
    let [a, b, c] = structure.cell_vectors.map(DVec3::from);
    let cell = DMat3::from_cols(a, b, c);
    if cell.determinant().abs() <= 1e-8 {
        return;
    }
    let inverse = cell.inverse();
    const TOLERANCE: f64 = 1e-5;
    let vectors = [a, b, c];
    for atom in structure.atoms.iter_mut() {
        let position = DVec3::new(atom.x, atom.y, atom.z);
        let frac = inverse * (position - origin);
        let inside = (0..3).all(|k| frac[k] >= -1e-12 && frac[k] < 1.0 - TOLERANCE);
        if inside {
            continue;
        }
        let mut cells = DVec3::ZERO;
        for k in 0..3 {
            cells[k] = frac[k].floor();
            if frac[k] - cells[k] >= 1.0 - TOLERANCE {
                cells[k] += 1.0;
            }
        }
        let mut wrapped = position - (cells.x * a + cells.y * b + cells.z * c);
        let moved = inverse * (wrapped - origin);
        for k in 0..3 {
            if cells[k] != 0.0 && moved[k] < 0.0 && moved[k] > -TOLERANCE {
                wrapped -= moved[k] * vectors[k];
            }
        }
        atom.x = wrapped.x;
        atom.y = wrapped.y;
        atom.z = wrapped.z;
    }
}

/// One atom per site of a disordered structure (the majority species; one
/// of a split position), with a note.
pub fn resolve_partial_occupancy(structure: &mut Structure, occupancies: &[f64]) {
    if occupancies.len() != structure.atoms.len() {
        return;
    }
    if !occupancies.iter().any(|&o| o < 1.0 - 1e-6) {
        return;
    }
    if structure.has_unit_cell {
        wrap_atoms_into_primary_cell(structure);
    }
    let matrices = cell_matrices(structure);
    let periodic = matrices.is_some();
    let (cell, inverse) = matrices.unwrap_or((DMat3::IDENTITY, DMat3::IDENTITY));
    let atoms = &structure.atoms;
    let separation = |i: usize, j: usize| {
        let mut d = DVec3::new(atoms[j].x - atoms[i].x, atoms[j].y - atoms[i].y, atoms[j].z - atoms[i].z);
        if periodic {
            let mut f = inverse * d;
            f -= DVec3::new(f.x.round(), f.y.round(), f.z.round());
            d = cell * f;
        }
        d.length()
    };
    let partials: Vec<usize> = (0..atoms.len()).filter(|&i| occupancies[i] < 1.0 - 1e-6).collect();
    let mut low = DVec3::ZERO;
    let mut extent = DVec3::ONE;
    let mut bins = IVec3::ONE;
    if periodic {
        let volume = cell.determinant().abs();
        for k in 0..3 {
            let height = volume / cell.col((k + 1) % 3).cross(cell.col((k + 2) % 3)).length();
            bins[k] = ((height / 0.5) as i32).clamp(1, 1000).max(1);
        }
    } else if !partials.is_empty() {
        let first = &atoms[partials[0]];
        let mut high = DVec3::new(first.x, first.y, first.z);
        low = high;
        for &i in &partials {
            let p = DVec3::new(atoms[i].x, atoms[i].y, atoms[i].z);
            low = low.min(p);
            high = high.max(p);
        }
        extent = (high - low).max(DVec3::splat(1e-9));
        for k in 0..3 {
            bins[k] = ((extent[k] / 0.5) as i32).clamp(1, 1000).max(1);
        }
    }
    let coordinate = |i: usize| {
        let p = DVec3::new(atoms[i].x, atoms[i].y, atoms[i].z);
        if periodic {
            let f = inverse * p;
            f - DVec3::new(f.x.floor(), f.y.floor(), f.z.floor())
        } else {
            (p - low) / extent
        }
    };
    let bin_of = |i: usize| {
        let f = coordinate(i);
        let mut b = IVec3::ZERO;
        for k in 0..3 {
            b[k] = (bins[k] - 1).min(0.max((f[k] * bins[k] as f64) as i32));
        }
        b
    };
    let mut binned: BTreeMap<(i32, i32, i32), Vec<usize>> = BTreeMap::new();
    for &i in &partials {
        let b = bin_of(i);
        binned.entry((b.x, b.y, b.z)).or_default().push(i);
    }
    let candidates = |i: usize| {
        let b = bin_of(i);
        let mut found = Vec::new();
        let mut seen = BTreeSet::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let mut c = b + IVec3::new(dx, dy, dz);
                    let mut inside = true;
                    for k in 0..3 {
                        if periodic {
                            c[k] = ((c[k] % bins[k]) + bins[k]) % bins[k];
                        } else if c[k] < 0 || c[k] >= bins[k] {
                            inside = false;
                        }
                    }
                    if !inside || !seen.insert((c.x, c.y, c.z)) {
                        continue;
                    }
                    if let Some(list) = binned.get(&(c.x, c.y, c.z)) {
                        found.extend_from_slice(list);
                    }
                }
            }
        }
        found.sort_unstable();
        found
    };
    let mut site_of = vec![-1i32; atoms.len()];
    let mut sites: Vec<Vec<usize>> = Vec::new();
    for &i in &partials {
        if site_of[i] >= 0 {
            continue;
        }
        site_of[i] = sites.len() as i32;
        sites.push(vec![i]);
        let mut total = occupancies[i];
        for j in candidates(i) {
            if site_of[j] >= 0 {
                continue;
            }
            let distance = separation(i, j);
            if distance < 0.05 || (distance < 0.5 && total + occupancies[j] <= 1.0 + 0.02) {
                site_of[j] = site_of[i];
                sites.last_mut().unwrap().push(j);
                total += occupancies[j];
            }
        }
    }
    let mut keep = vec![true; atoms.len()];
    let (mut shared, mut split, mut partial) = (0usize, 0usize, 0usize);
    for members in &sites {
        if members.len() > 1 {
            let one_species = members.iter().all(|&m| atoms[m].atomic_number == atoms[members[0]].atomic_number);
            if one_species {
                split += 1;
            } else {
                shared += 1;
            }
            let mut best = members[0];
            for &m in members {
                if occupancies[m] > occupancies[best] + 1e-9 {
                    best = m;
                }
            }
            for &m in members {
                keep[m] = m == best;
            }
        } else if occupancies[members[0]] < 1.0 - 1e-6 {
            partial += 1;
        }
    }
    if shared > 0 || split > 0 {
        let kept: Vec<_> = structure.atoms.iter().zip(&keep).filter(|(_, &k)| k).map(|(a, _)| a.clone()).collect();
        structure.atoms = kept;
    }
    let mut note = String::new();
    if shared > 0 {
        note = format!(
            "{shared}{} more than one species (partial occupancy); each keeps its majority species \
             (Build > Substitutional Solid Solution makes a random mixture).",
            if shared == 1 { " site holds" } else { " sites hold" }
        );
    }
    if split > 0 {
        if !note.is_empty() {
            note.push(' ');
        }
        note.push_str(&format!(
            "{split}{} split over nearby positions (positional disorder); one position is kept.",
            if split == 1 { " atom is" } else { " atoms are" }
        ));
    }
    if partial > 0 {
        if !note.is_empty() {
            note.push(' ');
        }
        note.push_str(&format!(
            "{partial}{}",
            if partial == 1 { " site is partially occupied; its atom is kept." } else { " sites are partially occupied; their atoms are kept." }
        ));
    }
    structure.load_note = note;
}
