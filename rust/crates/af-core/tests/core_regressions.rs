//! Port of tests/core_regressions.cpp (the picking, metadata and cell parts;
//! the sculptor, SSS and amorphous checks move to af-algorithms).

use af_core::cell::cell_matrices;
use af_core::picking::{pick_atom, pick_ray_dir, pick_ray_origin};
use af_core::{AtomSite, Structure};
use glam::{Mat4, Vec3};

#[test]
fn picking() {
    let origin = Vec3::ZERO;
    let dir = Vec3::Z;
    let positions = [Vec3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 6.0)];
    assert_eq!(pick_atom(origin, dir, &positions, &[0.5, 2.0], 1.0), Some(1), "the nearest surface, not centre");
    assert_eq!(pick_atom(origin, 3.0 * dir, &positions, &[0.5, 2.0], 1.0), Some(1), "a non-unit direction");
    assert_eq!(pick_atom(origin, origin, &positions, &[], 1.0), None, "a zero ray misses");
    assert_eq!(pick_atom(origin, dir, &[Vec3::new(0.0, 0.0, -0.5)], &[1.0], 1.0), Some(0), "a ray starting inside a sphere hits");
    assert_eq!(pick_atom(origin, dir, &[Vec3::new(0.0, 0.0, -5.0)], &[1.0], 1.0), None, "a sphere behind the ray misses");
}

#[test]
fn cursor_rays_go_through_the_scene() {
    let view = Mat4::look_at_rh(Vec3::new(0.0, 0.0, 10.0), Vec3::ZERO, Vec3::Y);
    let projection = Mat4::perspective_rh(40f32.to_radians(), 1.5, 0.1, 100.0);
    // The window centre looks straight down -z.
    let d = pick_ray_dir(300.0, 200.0, 600, 400, &projection, &view);
    assert!((d - Vec3::NEG_Z).length() < 1e-5);
    let o = pick_ray_origin(300.0, 200.0, 600, 400, &projection, &view);
    assert!((o - Vec3::new(0.0, 0.0, 9.9)).length() < 1e-4, "{o}");
    // A cursor right of centre points to +x.
    assert!(pick_ray_dir(500.0, 200.0, 600, 400, &projection, &view).x > 0.1);
}

#[test]
fn metadata_and_cells() {
    let mut base = Structure::default();
    base.atoms = vec![AtomSite::default(); 3];
    base.grain_colors = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    base.grain_region_ids = vec![10, 20, 30];
    base.erase_atom(1);
    assert!(base.atoms.len() == 2 && base.grain_region_ids == vec![10, 30] && base.grain_colors[1][2] == 1.0);
    base.has_unit_cell = true;
    base.cell_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(cell_matrices(&base).is_some(), "a valid cell is accepted");
    base.cell_vectors[0][0] = f64::NAN;
    assert!(cell_matrices(&base).is_none(), "a non-finite cell is rejected");
}
