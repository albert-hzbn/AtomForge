//! Grab mode (G, the C++ app/InteractionHandlers Blender-style grab): the
//! selected atoms follow the mouse in the view plane through their centroid,
//! optionally along one axis; periodic structures wrap the moved atoms into
//! the cell and keep the boundary images. Confirm keeps the move, cancel
//! restores every position.

use std::collections::{HashMap, HashSet};

use af_core::glm;
use af_core::Structure;
use glam::{Mat3, Vec3, Vec4};

use crate::scene::{FrameView, Scene};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Axis {
    #[default]
    Free,
    X,
    Y,
    Z,
}

impl Axis {
    /// X, Y or Z again releases the axis.
    pub fn toggled(self, axis: Axis) -> Axis {
        if self == axis {
            Axis::Free
        } else {
            axis
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Axis::Free => "",
            Axis::X => " [X axis]",
            Axis::Y => " [Y axis]",
            Axis::Z => " [Z axis]",
        }
    }
}

const PBC_TOL: f32 = 1e-4;

/// A fractional coordinate in [0, 1); within the tolerance of 0 or 1 it is 0.
pub fn wrap_frac(v: f32) -> f32 {
    let v = v - v.floor();
    if v.abs() <= PBC_TOL || (1.0 - v).abs() <= PBC_TOL {
        0.0
    } else {
        v
    }
}

#[derive(Clone, Debug)]
struct Cell {
    a: Vec3,
    b: Vec3,
    c: Vec3,
    matrix: Mat3,
    inverse: Mat3,
    origin: Vec3,
    tolerance: f32,
}

#[derive(Clone, Debug)]
pub struct Grab {
    pub axis: Axis,
    start_mouse: [f32; 2],
    /// The selected instances' positions when the grab started (selection order).
    pub original_positions: Vec<Vec3>,
    /// Every instance of each grabbed atom, with its first position (instance order).
    siblings: HashMap<i32, Vec<(usize, Vec3)>>,
    original_base: HashMap<i32, Vec3>,
    cell: Option<Cell>,
}

impl Grab {
    /// startGrab: the selected atoms and all their periodic images.
    pub fn start(scene: &Scene, s: &Structure, selection: &[usize], mouse: [f32; 2]) -> Grab {
        let positions = &scene.instances.positions;
        let indices = &scene.instances.atom_indices;
        let original_positions = selection.iter().map(|&i| positions.get(i).copied().unwrap_or(Vec3::ZERO)).collect();
        let cell = if s.has_unit_cell {
            let v = |k: usize| Vec3::new(s.cell_vectors[k][0] as f32, s.cell_vectors[k][1] as f32, s.cell_vectors[k][2] as f32);
            let (a, b, c) = (v(0), v(1), v(2));
            let matrix = Mat3::from_cols(a, b, c);
            (glm::determinant(&matrix).abs() > 1e-8).then(|| Cell {
                a,
                b,
                c,
                matrix,
                inverse: glm::inverse(&matrix),
                origin: Vec3::new(s.cell_offset[0] as f32, s.cell_offset[1] as f32, s.cell_offset[2] as f32),
                tolerance: if s.pbc_boundary_tol > 0.0 { s.pbc_boundary_tol } else { PBC_TOL },
            })
        } else {
            None
        };
        let grabbed: HashSet<i32> = selection.iter().filter_map(|&i| indices.get(i).copied()).collect();
        let mut original_base = HashMap::new();
        for &base in &grabbed {
            if base >= 0 {
                if let Some(a) = s.atoms.get(base as usize) {
                    original_base.insert(base, Vec3::new(a.x as f32, a.y as f32, a.z as f32));
                }
            }
        }
        let mut siblings: HashMap<i32, Vec<(usize, Vec3)>> = HashMap::new();
        for (i, &base) in indices.iter().enumerate() {
            if grabbed.contains(&base) {
                siblings.entry(base).or_default().push((i, positions[i]));
            }
        }
        Grab { axis: Axis::Free, start_mouse: mouse, original_positions, siblings, original_base, cell }
    }

    /// The selected atoms' centroid when the grab started.
    fn centroid(&self) -> Vec3 {
        if self.original_positions.is_empty() {
            return Vec3::ZERO;
        }
        self.original_positions.iter().fold(Vec3::ZERO, |a, &p| a + p) / self.original_positions.len() as f32
    }

    /// The world displacement for the mouse at `mouse` (pixels), axis applied.
    pub fn delta(&self, frame: &FrameView, mouse: [f32; 2]) -> Vec3 {
        let d = screen_to_world_delta([mouse[0] - self.start_mouse[0], mouse[1] - self.start_mouse[1]], self.centroid(), frame);
        match self.axis {
            Axis::Free => d,
            Axis::X => Vec3::new(d.x, 0.0, 0.0),
            Axis::Y => Vec3::new(0.0, d.y, 0.0),
            Axis::Z => Vec3::new(0.0, 0.0, d.z),
        }
    }

    /// Moves the selected atoms (and their images) by `delta` from where they started.
    pub fn apply(&self, scene: &mut Scene, s: &mut Structure, selection: &[usize], delta: Vec3) {
        let mut done = HashSet::new();
        for &instance in selection {
            let Some(&base) = scene.instances.atom_indices.get(instance) else { continue };
            if !done.insert(base) {
                continue;
            }
            let (Some(siblings), Some(&original)) = (self.siblings.get(&base), self.original_base.get(&base)) else { continue };
            let moved = original + delta;
            let Some(cell) = &self.cell else {
                for &(i, p) in siblings {
                    scene.instances.positions[i] = p + delta;
                }
                set_position(s, base, moved);
                continue;
            };
            let f = glm::mul(&cell.inverse, moved - cell.origin);
            let f = Vec3::new(wrap_frac(f.x), wrap_frac(f.y), wrap_frac(f.z));
            let wrapped = cell.origin + glm::mul(&cell.matrix, f);
            set_position(s, base, wrapped);
            // The cell and its boundary images: the first instances show them,
            // the others wait at the wrapped position.
            let shifts = |v: f32| if v.abs() <= cell.tolerance { vec![0.0, 1.0] } else { vec![0.0] };
            let mut needed = Vec::new();
            for sx in shifts(f.x) {
                for sy in shifts(f.y) {
                    for sz in shifts(f.z) {
                        needed.push(wrapped + cell.a * sx + cell.b * sy + cell.c * sz);
                    }
                }
            }
            for (j, &(i, _)) in siblings.iter().enumerate() {
                scene.instances.positions[i] = needed.get(j).copied().unwrap_or(wrapped);
            }
        }
    }

    /// Esc or a right click: every instance and atom back where it was.
    pub fn cancel(&self, scene: &mut Scene, s: &mut Structure) {
        for siblings in self.siblings.values() {
            for &(i, p) in siblings {
                scene.instances.positions[i] = p;
            }
        }
        for (&base, &p) in &self.original_base {
            set_position(s, base, p);
        }
    }

    /// A left click or Enter: the moved atoms are wrapped into the cell (the scene is rebuilt after).
    pub fn confirm(&self, s: &mut Structure) {
        let Some(cell) = &self.cell else { return };
        for &base in self.original_base.keys() {
            let Some(a) = (base >= 0).then(|| s.atoms.get(base as usize)).flatten() else { continue };
            let f = glm::mul(&cell.inverse, Vec3::new(a.x as f32, a.y as f32, a.z as f32) - cell.origin);
            let f = Vec3::new(wrap_frac(f.x), wrap_frac(f.y), wrap_frac(f.z));
            set_position(s, base, cell.origin + glm::mul(&cell.matrix, f));
        }
    }
}

fn set_position(s: &mut Structure, base: i32, p: Vec3) {
    if let Some(a) = (base >= 0).then(|| s.atoms.get_mut(base as usize)).flatten() {
        a.x = p.x as f64;
        a.y = p.y as f64;
        a.z = p.z as f64;
    }
}

/// A screen displacement (pixels, y down) as a world displacement in the
/// plane through `reference` facing the camera (screenToWorldDelta).
pub fn screen_to_world_delta(screen: [f32; 2], reference: Vec3, frame: &FrameView) -> Vec3 {
    let vp = frame.projection * frame.view;
    let clip = vp * Vec4::new(reference.x, reference.y, reference.z, 1.0);
    if clip.w <= 1e-6 {
        return Vec3::ZERO;
    }
    let ndc = clip.truncate() / clip.w;
    let (dx, dy) = (screen[0] / frame.width as f32 * 2.0, -(screen[1] / frame.height as f32) * 2.0);
    let inverse = vp.inverse();
    let a = inverse * Vec4::new(ndc.x, ndc.y, ndc.z, 1.0);
    let b = inverse * Vec4::new(ndc.x + dx, ndc.y + dy, ndc.z, 1.0);
    if a.w.abs() < 1e-8 || b.w.abs() < 1e-8 {
        return Vec3::ZERO;
    }
    (b / b.w).truncate() - (a / a.w).truncate()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{build_frame_view, build_instance_data, BondStyle, OrbitCamera};
    use af_core::display::{default_look, DisplaySettings};
    use af_core::structure::AtomSite;

    fn scene_of(s: &Structure) -> Scene {
        let look = default_look();
        Scene::new(build_instance_data(s, false, &[[1, 0, 0], [0, 1, 0], [0, 0, 1]], &look.radii, &vec![32.0; 119]), &BondStyle::default())
    }

    fn front_view(scene: &Scene) -> FrameView {
        let camera = OrbitCamera { yaw: 0.0, pitch: 0.0, roll: 0.0, distance: 20.0, ..Default::default() };
        build_frame_view(&camera, scene, true, 400, 300, &DisplaySettings::default())
    }

    #[test]
    fn wrapping() {
        assert_eq!(wrap_frac(1.25), 0.25);
        assert_eq!(wrap_frac(-0.25), 0.75);
        assert_eq!(wrap_frac(0.99999), 0.0);
        assert_eq!(wrap_frac(-0.00001), 0.0);
        assert_eq!(Axis::X.toggled(Axis::X), Axis::Free);
        assert_eq!(Axis::Free.toggled(Axis::Y), Axis::Y);
    }

    #[test]
    fn moves_in_the_view_plane_and_cancels() {
        let mut s = Structure::default();
        for x in [0.0, 3.0] {
            s.atoms.push(AtomSite { symbol: "C".into(), atomic_number: 6, x, ..Default::default() });
        }
        let mut scene = scene_of(&s);
        let frame = front_view(&scene);
        // Seen from +z: x to the right, y up on screen.
        let grab = Grab::start(&scene, &s, &[1], [100.0, 100.0]);
        let right = grab.delta(&frame, [140.0, 100.0]);
        assert!(right.x > 0.0 && right.y.abs() < 1e-4 && right.z.abs() < 1e-4, "{right}");
        let up = grab.delta(&frame, [100.0, 60.0]);
        assert!(up.y > 0.0 && up.x.abs() < 1e-4, "{up}");
        // 40 pixels both ways move the same distance.
        assert!((right.x - up.y).abs() < 1e-4);
        let original = scene.instances.positions[1];
        grab.apply(&mut scene, &mut s, &[1], right);
        assert_eq!(scene.instances.positions[1], original + right);
        assert_eq!(s.atoms[1].x, (3.0 + right.x) as f64);
        assert_eq!(s.atoms[0].x, 0.0);
        grab.cancel(&mut scene, &mut s);
        assert_eq!(scene.instances.positions[1], original);
        assert_eq!(s.atoms[1].x, 3.0);
        // An axis keeps only that component.
        let mut z = Grab::start(&scene, &s, &[1], [0.0, 0.0]);
        z.axis = Axis::Z;
        assert_eq!(z.delta(&frame, [50.0, 50.0]), Vec3::ZERO);
    }

    #[test]
    fn periodic_atoms_wrap_and_keep_their_images() {
        let a = 4.0;
        let mut s = Structure { has_unit_cell: true, cell_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], ..Default::default() };
        s.atoms.push(AtomSite { symbol: "Cu".into(), atomic_number: 29, x: 1.0, y: 1.0, z: 1.0, ..Default::default() });
        s.atoms.push(AtomSite { symbol: "Cu".into(), atomic_number: 29, x: 2.0, y: 2.0, z: 2.0, ..Default::default() });
        let mut scene = scene_of(&s);
        let grab = Grab::start(&scene, &s, &[0], [0.0, 0.0]);
        // Moved past the cell edge: wrapped back in.
        grab.apply(&mut scene, &mut s, &[0], Vec3::new(3.5, 0.0, 0.0));
        assert!((s.atoms[0].x - 0.5).abs() < 1e-5, "{}", s.atoms[0].x);
        // Onto the face x = 0: the image at x = a appears too.
        grab.apply(&mut scene, &mut s, &[0], Vec3::new(-1.0, 0.0, 0.0));
        assert_eq!(s.atoms[0].x, 0.0);
        let shown: Vec<Vec3> = scene.instances.atom_indices.iter().zip(&scene.instances.positions).filter(|(&b, _)| b == 0).map(|(_, &p)| p).collect();
        // The atom had one instance (inside the cell); it shows the cell position.
        assert_eq!(shown, vec![Vec3::new(0.0, 1.0, 1.0)]);
        grab.confirm(&mut s);
        assert_eq!(s.atoms[0].x, 0.0);
        assert_eq!((s.atoms[1].x, s.atoms[1].y), (2.0, 2.0));
    }
}
