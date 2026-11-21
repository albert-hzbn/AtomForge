//! Integer supercells (the C++ graphics/StructureInstanceBuilder
//! buildSupercell): the rows of the transform give the new cell vectors as
//! combinations of a, b and c; atoms fill the new cell, moved by whole cell
//! vectors in double precision.

use glam::{DMat3, DVec3};

use crate::glm::{determinant_d, inverse_d, mul_d};
use crate::Structure;

/// The supercell of `structure` for the integer transform (rows: new
/// vectors). Unchanged without a cell or for a singular transform.
pub fn build_supercell(structure: &Structure, transform: &[[i32; 3]; 3]) -> Structure {
    if !structure.has_unit_cell {
        return structure.clone();
    }
    let eps = 1e-5f32 as f64;
    let origin = DVec3::from(structure.cell_offset);
    let [a, b, c] = structure.cell_vectors.map(DVec3::from);
    let cell = DMat3::from_cols(a, b, c);
    if determinant_d(&cell).abs() <= 1e-12 {
        return structure.clone();
    }
    let inverse_cell = inverse_d(&cell);
    // transform[col][row] = T[row][col] in glm terms.
    let mut m = [[0.0f64; 3]; 3];
    for row in 0..3 {
        for col in 0..3 {
            m[col][row] = transform[row][col] as f64;
        }
    }
    let t = DMat3::from_cols_array_2d(&m);
    if determinant_d(&t).abs() <= eps {
        return structure.clone();
    }
    let inverse_t = inverse_d(&t);
    let (mut n_min, mut n_max) = ([0i32; 3], [0i32; 3]);
    for row in 0..3 {
        let (mut low, mut high) = (0, 0);
        for col in 0..3 {
            let v = transform[row][col];
            if v < 0 {
                low += v;
            }
            if v > 0 {
                high += v;
            }
        }
        n_min[row] = low - 1;
        n_max[row] = high + 1;
    }
    let mut result = Structure { has_unit_cell: true, cell_offset: structure.cell_offset, ..Default::default() };
    let n = structure.atoms.len();
    let has_colours = structure.grain_colors.len() == n;
    let has_regions = structure.grain_region_ids.len() == n;
    let row = |r: usize| transform[r][0] as f64 * a + transform[r][1] as f64 * b + transform[r][2] as f64 * c;
    result.cell_vectors = [row(0).to_array(), row(1).to_array(), row(2).to_array()];
    for (index, atom) in structure.atoms.iter().enumerate() {
        let world = DVec3::new(atom.x, atom.y, atom.z);
        let raw = mul_d(&inverse_cell, world - origin);
        let cell_of = raw.floor();
        let frac = raw - cell_of;
        for ix in n_min[0]..=n_max[0] {
            for iy in n_min[1]..=n_max[1] {
                for iz in n_min[2]..=n_max[2] {
                    let shift = DVec3::new(ix as f64, iy as f64, iz as f64);
                    let mapped = mul_d(&inverse_t, frac + shift);
                    if mapped.x < -eps || mapped.y < -eps || mapped.z < -eps {
                        continue;
                    }
                    if mapped.x >= 1.0 - eps || mapped.y >= 1.0 - eps || mapped.z >= 1.0 - eps {
                        continue;
                    }
                    // Moved by whole cells: (ix, iy, iz) less the cell it was in.
                    let cells = shift - cell_of;
                    let moved = world + cells.x * a + cells.y * b + cells.z * c;
                    let mut site = atom.clone();
                    site.x = moved.x;
                    site.y = moved.y;
                    site.z = moved.z;
                    result.atoms.push(site);
                    if has_colours {
                        result.grain_colors.push(structure.grain_colors[index]);
                    }
                    if has_regions {
                        result.grain_region_ids.push(structure.grain_region_ids[index]);
                    }
                }
            }
        }
    }
    result
}
