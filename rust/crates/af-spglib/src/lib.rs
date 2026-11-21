//! spglib 2.7.0, the library the C++ AtomForge links, so symmetry results
//! match it bit for bit. Safe wrappers take lattice vectors as ROWS (the
//! AtomForge convention) and hand spglib columns.

use std::ffi::{c_char, c_int, CStr};

pub mod ffi {
    use std::ffi::{c_char, c_double, c_int};

    #[repr(C)]
    pub struct SpglibDataset {
        pub spacegroup_number: c_int,
        pub hall_number: c_int,
        pub international_symbol: [c_char; 11],
        pub hall_symbol: [c_char; 17],
        pub choice: [c_char; 6],
        pub transformation_matrix: [[c_double; 3]; 3],
        pub origin_shift: [c_double; 3],
        pub n_operations: c_int,
        pub rotations: *mut [[c_int; 3]; 3],
        pub translations: *mut [c_double; 3],
        pub n_atoms: c_int,
        pub wyckoffs: *mut c_int,
        pub site_symmetry_symbols: *mut [c_char; 7],
        pub equivalent_atoms: *mut c_int,
        pub crystallographic_orbits: *mut c_int,
        pub primitive_lattice: [[c_double; 3]; 3],
        pub mapping_to_primitive: *mut c_int,
        pub n_std_atoms: c_int,
        pub std_lattice: [[c_double; 3]; 3],
        pub std_types: *mut c_int,
        pub std_positions: *mut [c_double; 3],
        pub std_rotation_matrix: [[c_double; 3]; 3],
        pub std_mapping_to_primitive: *mut c_int,
        pub pointgroup_symbol: [c_char; 6],
    }

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub struct SpglibSpacegroupType {
        pub number: c_int,
        pub international_short: [c_char; 11],
        pub international_full: [c_char; 20],
        pub international: [c_char; 32],
        pub schoenflies: [c_char; 7],
        pub hall_number: c_int,
        pub hall_symbol: [c_char; 17],
        pub choice: [c_char; 6],
        pub pointgroup_international: [c_char; 6],
        pub pointgroup_schoenflies: [c_char; 4],
        pub arithmetic_crystal_class_number: c_int,
        pub arithmetic_crystal_class_symbol: [c_char; 7],
    }

    extern "C" {
        pub fn spg_get_version() -> *const c_char;
        pub fn spg_get_error_code() -> c_int;
        pub fn spg_get_dataset(lattice: *const [c_double; 3], position: *const [c_double; 3], types: *const c_int, num_atom: c_int, symprec: c_double) -> *mut SpglibDataset;
        pub fn spgat_get_dataset(
            lattice: *const [c_double; 3],
            position: *const [c_double; 3],
            types: *const c_int,
            num_atom: c_int,
            symprec: c_double,
            angle_tolerance: c_double,
        ) -> *mut SpglibDataset;
        pub fn spg_free_dataset(dataset: *mut SpglibDataset);
        pub fn spg_standardize_cell(
            lattice: *mut [c_double; 3],
            position: *mut [c_double; 3],
            types: *mut c_int,
            num_atom: c_int,
            to_primitive: c_int,
            no_idealize: c_int,
            symprec: c_double,
        ) -> c_int;
        pub fn spg_niggli_reduce(lattice: *mut [c_double; 3], symprec: c_double) -> c_int;
        pub fn spg_delaunay_reduce(lattice: *mut [c_double; 3], symprec: c_double) -> c_int;
        pub fn spg_get_symmetry_from_database(rotations: *mut [[c_int; 3]; 3], translations: *mut [c_double; 3], hall_number: c_int) -> c_int;
        pub fn spg_get_spacegroup_type(hall_number: c_int) -> SpglibSpacegroupType;
    }
}

pub type Mat3 = [[f64; 3]; 3];

fn transpose(m: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i]))
}

fn text(chars: &[c_char]) -> String {
    let bytes: Vec<u8> = chars.iter().take_while(|&&c| c != 0).map(|&c| c as u8).collect();
    String::from_utf8_lossy(&bytes).trim_end().to_string()
}

pub fn version() -> String {
    unsafe { CStr::from_ptr(ffi::spg_get_version()) }.to_string_lossy().into_owned()
}

/// spglib's error code after the last call (0 = success).
pub fn error_code() -> i32 {
    unsafe { ffi::spg_get_error_code() }
}

/// An owned copy of spg_get_dataset's result; lattices as rows.
#[derive(Clone, Debug, PartialEq)]
pub struct Dataset {
    pub spacegroup_number: i32,
    pub hall_number: i32,
    pub international_symbol: String,
    pub hall_symbol: String,
    pub choice: String,
    pub transformation_matrix: Mat3,
    pub origin_shift: [f64; 3],
    pub rotations: Vec<[[i32; 3]; 3]>,
    pub translations: Vec<[f64; 3]>,
    pub wyckoffs: Vec<i32>,
    pub site_symmetry_symbols: Vec<String>,
    pub equivalent_atoms: Vec<i32>,
    pub crystallographic_orbits: Vec<i32>,
    /// Rows (spglib's columns transposed).
    pub primitive_lattice: Mat3,
    pub mapping_to_primitive: Vec<i32>,
    /// Rows.
    pub std_lattice: Mat3,
    pub std_types: Vec<i32>,
    pub std_positions: Vec<[f64; 3]>,
    pub std_rotation_matrix: Mat3,
    pub std_mapping_to_primitive: Vec<i32>,
    pub pointgroup_symbol: String,
}

unsafe fn slice<'a, T>(p: *const T, n: c_int) -> &'a [T] {
    if p.is_null() || n <= 0 {
        &[]
    } else {
        std::slice::from_raw_parts(p, n as usize)
    }
}

/// spg_get_dataset; `lattice` rows are the cell vectors, positions fractional.
pub fn get_dataset(lattice: &Mat3, positions: &[[f64; 3]], types: &[i32], symprec: f64) -> Option<Dataset> {
    get_dataset_with_angle(lattice, positions, types, symprec, -1.0)
}

/// spgat_get_dataset: as get_dataset with an angle tolerance in degrees
/// (negative: spglib's default).
pub fn get_dataset_with_angle(lattice: &Mat3, positions: &[[f64; 3]], types: &[i32], symprec: f64, angle_tolerance: f64) -> Option<Dataset> {
    if positions.len() != types.len() || positions.is_empty() {
        return None;
    }
    let columns = transpose(lattice);
    unsafe {
        let d = ffi::spgat_get_dataset(columns.as_ptr(), positions.as_ptr(), types.as_ptr(), positions.len() as c_int, symprec, angle_tolerance);
        if d.is_null() {
            return None;
        }
        let r = &*d;
        let out = Dataset {
            spacegroup_number: r.spacegroup_number,
            hall_number: r.hall_number,
            international_symbol: text(&r.international_symbol),
            hall_symbol: text(&r.hall_symbol),
            choice: text(&r.choice),
            transformation_matrix: r.transformation_matrix,
            origin_shift: r.origin_shift,
            rotations: slice(r.rotations, r.n_operations).to_vec(),
            translations: slice(r.translations, r.n_operations).to_vec(),
            wyckoffs: slice(r.wyckoffs, r.n_atoms).to_vec(),
            site_symmetry_symbols: slice(r.site_symmetry_symbols, r.n_atoms).iter().map(|s| text(s)).collect(),
            equivalent_atoms: slice(r.equivalent_atoms, r.n_atoms).to_vec(),
            crystallographic_orbits: slice(r.crystallographic_orbits, r.n_atoms).to_vec(),
            primitive_lattice: transpose(&r.primitive_lattice),
            mapping_to_primitive: slice(r.mapping_to_primitive, r.n_atoms).to_vec(),
            std_lattice: transpose(&r.std_lattice),
            std_types: slice(r.std_types, r.n_std_atoms).to_vec(),
            std_positions: slice(r.std_positions, r.n_std_atoms).to_vec(),
            std_rotation_matrix: r.std_rotation_matrix,
            std_mapping_to_primitive: slice(r.std_mapping_to_primitive, r.n_std_atoms).to_vec(),
            pointgroup_symbol: text(&r.pointgroup_symbol),
        };
        ffi::spg_free_dataset(d);
        Some(out)
    }
}

/// A cell: lattice vectors as rows, fractional positions and species.
#[derive(Clone, Debug, PartialEq)]
pub struct Cell {
    pub lattice: Mat3,
    pub positions: Vec<[f64; 3]>,
    pub types: Vec<i32>,
}

/// spg_standardize_cell; None when spglib returns 0 atoms.
pub fn standardize_cell(lattice: &Mat3, positions: &[[f64; 3]], types: &[i32], to_primitive: bool, no_idealize: bool, symprec: f64) -> Option<Cell> {
    let n = positions.len();
    if n == 0 || n != types.len() {
        return None;
    }
    // A conventional cell can hold up to four times the atoms of a primitive input.
    let capacity = n * 4;
    let mut columns = transpose(lattice);
    let mut pos = positions.to_vec();
    pos.resize(capacity, [0.0; 3]);
    let mut ty = types.to_vec();
    ty.resize(capacity, 0);
    let count = unsafe {
        ffi::spg_standardize_cell(columns.as_mut_ptr(), pos.as_mut_ptr(), ty.as_mut_ptr(), n as c_int, to_primitive as c_int, no_idealize as c_int, symprec)
    };
    if count <= 0 {
        return None;
    }
    pos.truncate(count as usize);
    ty.truncate(count as usize);
    Some(Cell { lattice: transpose(&columns), positions: pos, types: ty })
}

/// spg_niggli_reduce of row vectors; None on failure.
pub fn niggli_reduce(lattice: &Mat3, symprec: f64) -> Option<Mat3> {
    let mut columns = transpose(lattice);
    (unsafe { ffi::spg_niggli_reduce(columns.as_mut_ptr(), symprec) } != 0).then(|| transpose(&columns))
}

/// spg_delaunay_reduce of row vectors; None on failure.
pub fn delaunay_reduce(lattice: &Mat3, symprec: f64) -> Option<Mat3> {
    let mut columns = transpose(lattice);
    (unsafe { ffi::spg_delaunay_reduce(columns.as_mut_ptr(), symprec) } != 0).then(|| transpose(&columns))
}

/// The symmetry operations of a Hall setting (1..=530).
pub fn symmetry_from_database(hall_number: i32) -> Vec<([[i32; 3]; 3], [f64; 3])> {
    let mut rotations = vec![[[0i32; 3]; 3]; 192];
    let mut translations = vec![[0.0f64; 3]; 192];
    let n = unsafe { ffi::spg_get_symmetry_from_database(rotations.as_mut_ptr(), translations.as_mut_ptr(), hall_number) };
    (0..n.max(0) as usize).map(|i| (rotations[i], translations[i])).collect()
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpacegroupType {
    pub number: i32,
    pub international_short: String,
    pub international_full: String,
    pub international: String,
    pub schoenflies: String,
    pub hall_number: i32,
    pub hall_symbol: String,
    pub choice: String,
    pub pointgroup_international: String,
    pub pointgroup_schoenflies: String,
    pub arithmetic_crystal_class_number: i32,
    pub arithmetic_crystal_class_symbol: String,
}

pub fn spacegroup_type(hall_number: i32) -> SpacegroupType {
    let t = unsafe { ffi::spg_get_spacegroup_type(hall_number) };
    SpacegroupType {
        number: t.number,
        international_short: text(&t.international_short),
        international_full: text(&t.international_full),
        international: text(&t.international),
        schoenflies: text(&t.schoenflies),
        hall_number: t.hall_number,
        hall_symbol: text(&t.hall_symbol),
        choice: text(&t.choice),
        pointgroup_international: text(&t.pointgroup_international),
        pointgroup_schoenflies: text(&t.pointgroup_schoenflies),
        arithmetic_crystal_class_number: t.arithmetic_crystal_class_number,
        arithmetic_crystal_class_symbol: text(&t.arithmetic_crystal_class_symbol),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_spglib_2_7_0() {
        // The vendored build on Windows; elsewhere the system's spglib 2.x.
        if cfg!(all(windows, target_env = "msvc")) {
            assert_eq!(version(), "2.7.0");
        } else {
            assert!(version().starts_with("2."), "spglib {}", version());
        }
    }

    #[test]
    fn fcc_dataset_and_primitive() {
        let a = 3.6;
        let lattice = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];
        let types = [29; 4];
        let d = get_dataset(&lattice, &positions, &types, 1e-3).unwrap();
        assert_eq!(d.spacegroup_number, 225);
        assert_eq!(d.international_symbol, "Fm-3m");
        assert_eq!(d.rotations.len(), 192);
        let p = standardize_cell(&lattice, &positions, &types, true, false, 1e-3).unwrap();
        assert_eq!(p.positions.len(), 1);
        assert_eq!(spacegroup_type(d.hall_number).number, 225);
        assert_eq!(symmetry_from_database(d.hall_number).len(), 192);
    }
}
