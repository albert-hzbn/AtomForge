//! The atomic structure (the C++ model/Structure.h): atoms with per-atom
//! metadata, an optional unit cell with periodicity, and optional per-atom
//! arrays from builders and analyses.

use serde::{Deserialize, Serialize};

use crate::elements;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AtomSite {
    pub symbol: String,
    pub atomic_number: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub r: f32,
    pub g: f32,
    pub b: f32,
    /// Coordinates held fixed in a relaxation: bit 0 x, bit 1 y, bit 2 z.
    pub fixed: u8,
    /// Initial magnetic moment (Bohr magnetons).
    pub magmom: f64,
    /// Charge (e).
    pub charge: f64,
    /// Velocity (Angstrom/fs).
    pub velocity: [f64; 3],
}

impl Default for AtomSite {
    fn default() -> Self {
        Self {
            symbol: String::new(),
            atomic_number: 0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            r: 1.0,
            g: 1.0,
            b: 1.0,
            fixed: 0,
            magmom: 0.0,
            charge: 0.0,
            velocity: [0.0; 3],
        }
    }
}

impl AtomSite {
    /// An atom of an element at a position, in the element's default colour.
    pub fn new(symbol: &str, position: [f64; 3]) -> Self {
        let z = elements::atomic_number(symbol);
        let [r, g, b] = elements::default_colour(z);
        Self { symbol: symbol.to_string(), atomic_number: z, x: position[0], y: position[1], z: position[2], r, g, b, ..Default::default() }
    }

    pub fn position(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    pub fn set_position(&mut self, p: [f64; 3]) {
        self.x = p[0];
        self.y = p[1];
        self.z = p[2];
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Structure {
    pub atoms: Vec<AtomSite>,
    pub has_unit_cell: bool,
    /// Rows are the cell vectors a, b and c.
    pub cell_vectors: [[f64; 3]; 3],
    /// Periodic along a, b and c (with a cell).
    pub pbc: [bool; 3],
    pub cell_offset: [f64; 3],
    /// When > 0, the tolerance for periodic-boundary image atoms.
    pub pbc_boundary_tol: f32,
    /// Per-atom grain orientation colours (IPF-Z), same length as atoms or empty.
    pub grain_colors: Vec<[f32; 3]>,
    /// Per-atom grain or region ids, same length as atoms or empty.
    pub grain_region_ids: Vec<i32>,
    /// Per-atom scalar from an analysis (NaN = invalid), same length as atoms or empty.
    pub atom_property: Vec<f64>,
    pub atom_property_name: String,
    pub frozen_property: Vec<f64>,
    pub frozen_property_name: String,
    pub ipf_load_status: String,
    pub load_note: String,
    /// Dislocation loop outline from the dislocation builder.
    pub dislocation_loop_points: Vec<[f32; 3]>,
    pub dislocation_detection_done: bool,
}

impl Default for Structure {
    fn default() -> Self {
        Self {
            atoms: Vec::new(),
            has_unit_cell: false,
            cell_vectors: [[0.0; 3]; 3],
            pbc: [true, true, true],
            cell_offset: [0.0; 3],
            pbc_boundary_tol: 0.0,
            grain_colors: Vec::new(),
            grain_region_ids: Vec::new(),
            atom_property: Vec::new(),
            atom_property_name: String::new(),
            frozen_property: Vec::new(),
            frozen_property_name: String::new(),
            ipf_load_status: String::new(),
            load_note: String::new(),
            dislocation_loop_points: Vec::new(),
            dislocation_detection_done: false,
        }
    }
}

impl Structure {
    /// Removes an atom, keeping the aligned per-atom arrays aligned (arrays of
    /// another length are cleared), as the C++ eraseAtom.
    pub fn erase_atom(&mut self, index: usize) {
        if index >= self.atoms.len() {
            return;
        }
        let n = self.atoms.len();
        if self.grain_colors.len() == n {
            self.grain_colors.remove(index);
        } else {
            self.grain_colors.clear();
        }
        if self.grain_region_ids.len() == n {
            self.grain_region_ids.remove(index);
        } else {
            self.grain_region_ids.clear();
        }
        if self.atom_property.len() == n {
            self.atom_property.remove(index);
        } else {
            self.atom_property.clear();
        }
        if self.frozen_property.len() == n {
            self.frozen_property.remove(index);
        } else {
            self.frozen_property.clear();
        }
        self.atoms.remove(index);
        self.dislocation_loop_points.clear();
        self.dislocation_detection_done = false;
    }

    /// Periodic along at least one cell vector (with a cell).
    pub fn is_periodic(&self) -> bool {
        self.has_unit_cell && self.pbc.iter().any(|&p| p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erase_keeps_arrays_aligned() {
        let mut s = Structure::default();
        for i in 0..4 {
            s.atoms.push(AtomSite::new("Cu", [i as f64, 0.0, 0.0]));
        }
        s.atom_property = vec![0.0, 1.0, 2.0, 3.0];
        s.grain_region_ids = vec![1, 2];  // misaligned: cleared
        s.erase_atom(1);
        assert_eq!(s.atoms.len(), 3);
        assert_eq!(s.atom_property, vec![0.0, 2.0, 3.0]);
        assert!(s.grain_region_ids.is_empty());
        assert_eq!(s.atoms[1].x, 2.0);
        assert_eq!(s.atoms[0].atomic_number, 29);
    }
}
