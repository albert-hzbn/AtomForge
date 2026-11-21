//! AtomForge file input and output: structure formats (native readers and
//! writers for every format, including those the C++ read through Open
//! Babel), orientation-colour sidecars and simulation-code inputs. Port of
//! the C++ src/io.

pub mod cif;
pub mod ipf;
pub mod lammps;
pub mod loader;
pub mod molecular;
pub mod occupancy;
pub mod poscar;
pub mod simulation;
pub mod spacegroups;
pub mod text;
pub mod trajectory;
pub mod writers;
pub mod xyz;

pub use loader::{is_large_data_file, is_supported_structure_file, load_structure_from_file, save_structure, structure_text};
