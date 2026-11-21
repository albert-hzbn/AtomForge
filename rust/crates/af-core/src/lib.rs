//! AtomForge core: the structure model, element data, cell maths, neighbour
//! search, display settings and the C-compatible number formatting shared by
//! every other crate. Port of the C++ src/model, src/math and src/util.

pub mod cell;
pub mod crmath;
pub mod cxx;
pub mod cxx_random;
pub mod display;
pub mod dump;
pub mod eigen;
mod element_tables;
pub mod elements;
pub mod format;
pub mod glm;
pub mod json;
pub mod neighbours;
pub mod occlusion;
pub mod picking;
pub mod structure;
pub mod suggest;
pub mod supercell;
pub mod task;

pub use structure::{AtomSite, Structure};
