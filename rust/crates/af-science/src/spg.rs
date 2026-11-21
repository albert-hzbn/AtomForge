//! Space-group detection and standardisation through spglib 2.7.0 (the
//! library the C++ AtomForge links, vendored by af-spglib), so results match
//! the C++ bit for bit. Lattices are rows here; af-spglib transposes.

pub use af_spglib::{delaunay_reduce, niggli_reduce, spacegroup_type, symmetry_from_database, Dataset};

#[derive(Clone, Debug, PartialEq)]
pub struct SpgDataset {
    pub spacegroup_number: i32,
    pub hall_number: i32,
    /// Short international symbol of the standard setting ("Fm-3m").
    pub international_symbol: String,
    /// Symmetry operations of the input cell: rotation (rows) and translation.
    pub operations: Vec<([[i32; 3]; 3], [f64; 3])>,
}

/// A cell: lattice vectors as rows, fractional positions and species.
pub type SpgCell = af_spglib::Cell;

/// The idealised primitive standardized cell: spg_standardize_cell with
/// to_primitive = 1 and no_idealize = 0.
pub fn standardize_primitive(lattice: &[[f64; 3]; 3], positions: &[[f64; 3]], types: &[i32], symprec: f64) -> Option<SpgCell> {
    af_spglib::standardize_cell(lattice, positions, types, true, false, symprec)
}

/// The full spg_get_dataset result (lattice vectors as rows).
pub fn dataset(lattice: &[[f64; 3]; 3], positions: &[[f64; 3]], types: &[i32], symprec: f64) -> Option<Dataset> {
    af_spglib::get_dataset(lattice, positions, types, symprec)
}

/// The symmetry of a cell (lattice vectors as rows, fractional positions,
/// species numbers) within `symprec` (A), or None.
pub fn get_dataset(lattice: &[[f64; 3]; 3], positions: &[[f64; 3]], types: &[i32], symprec: f64) -> Option<SpgDataset> {
    let d = af_spglib::get_dataset(lattice, positions, types, symprec)?;
    Some(SpgDataset {
        spacegroup_number: d.spacegroup_number,
        hall_number: d.hall_number,
        international_symbol: d.international_symbol,
        operations: d.rotations.into_iter().zip(d.translations).collect(),
    })
}
