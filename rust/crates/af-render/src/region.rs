//! Box and lasso selection (the C++ app/InteractionHandlers
//! handleBoxSelection and handleLassoSelection): the instances whose centres
//! project inside a screen rectangle or polygon (pixels, y down), in front of
//! the far plane and behind the near one.

use glam::Vec4;

use crate::scene::{FrameView, Scene};

/// A drag smaller than 4 pixels in either direction is a click.
pub fn drag_large_enough(start: [f32; 2], end: [f32; 2]) -> bool {
    (end[0] - start[0]).abs() >= 4.0 && (end[1] - start[1]).abs() >= 4.0
}

pub fn inside_rect(p: [f32; 2], a: [f32; 2], b: [f32; 2]) -> bool {
    p[0] >= a[0].min(b[0]) && p[0] <= a[0].max(b[0]) && p[1] >= a[1].min(b[1]) && p[1] <= a[1].max(b[1])
}

/// Ray casting (even-odd) point in polygon; fewer than 3 points: outside.
pub fn inside_lasso(p: [f32; 2], polygon: &[[f32; 2]]) -> bool {
    if polygon.len() < 3 {
        return false;
    }
    let mut inside = false;
    let n = polygon.len();
    let mut j = n - 1;
    for i in 0..n {
        let ([xi, yi], [xj, yj]) = (polygon[i], polygon[j]);
        if (yi > p[1]) != (yj > p[1]) && p[0] < (xj - xi) * (p[1] - yi) / (yj - yi) + xi {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Screen positions (pixels) of the instances in the view volume.
fn screen_positions<'a>(scene: &'a Scene, frame: &FrameView) -> impl Iterator<Item = (usize, [f32; 2])> + 'a {
    let m = frame.projection * frame.view;
    let (w, h) = (frame.width as f32, frame.height as f32);
    scene.instances.positions.iter().enumerate().filter_map(move |(i, p)| {
        let clip = m * Vec4::new(p.x, p.y, p.z, 1.0);
        if clip.w <= 1e-6 {
            return None;
        }
        let ndc = clip.truncate() / clip.w;
        if ndc.z < -1.0 || ndc.z > 1.0 {
            return None;
        }
        Some((i, [(ndc.x * 0.5 + 0.5) * w, (1.0 - (ndc.y * 0.5 + 0.5)) * h]))
    })
}

/// Adds the instances inside the rectangle to `selection` (in index order,
/// without duplicates); without `additive` the selection starts empty.
pub fn select_box(selection: &mut Vec<usize>, scene: &Scene, frame: &FrameView, start: [f32; 2], end: [f32; 2], additive: bool) {
    if !drag_large_enough(start, end) {
        return;
    }
    if !additive {
        selection.clear();
    }
    for (i, p) in screen_positions(scene, frame) {
        if inside_rect(p, start, end) && !selection.contains(&i) {
            selection.push(i);
        }
    }
}

/// As select_box for a lasso polygon (at least 3 points).
pub fn select_lasso(selection: &mut Vec<usize>, scene: &Scene, frame: &FrameView, polygon: &[[f32; 2]], additive: bool) {
    if polygon.len() < 3 {
        return;
    }
    if !additive {
        selection.clear();
    }
    for (i, p) in screen_positions(scene, frame) {
        if inside_lasso(p, polygon) && !selection.contains(&i) {
            selection.push(i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::{build_frame_view, build_instance_data, BondStyle, OrbitCamera};
    use af_core::display::{default_look, DisplaySettings};
    use af_core::structure::{AtomSite, Structure};

    #[test]
    fn polygon_and_rectangle_tests() {
        let square = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]];
        assert!(inside_lasso([5.0, 5.0], &square));
        assert!(!inside_lasso([15.0, 5.0], &square));
        assert!(!inside_lasso([5.0, 5.0], &square[..2]));
        // A concave L shape.
        let l = [[0.0, 0.0], [10.0, 0.0], [10.0, 4.0], [4.0, 4.0], [4.0, 10.0], [0.0, 10.0]];
        assert!(inside_lasso([2.0, 8.0], &l));
        assert!(!inside_lasso([8.0, 8.0], &l));
        assert!(inside_rect([3.0, 3.0], [5.0, 5.0], [1.0, 1.0]));
        assert!(!drag_large_enough([0.0, 0.0], [3.0, 30.0]));
    }

    #[test]
    fn selects_the_atoms_seen_inside() {
        // A row of atoms along x, seen from +z.
        let mut s = Structure::default();
        for i in 0..5 {
            s.atoms.push(AtomSite { symbol: "C".into(), atomic_number: 6, x: i as f64 * 2.0, ..Default::default() });
        }
        let look = default_look();
        let scene = Scene::new(build_instance_data(&s, false, &[[1, 0, 0], [0, 1, 0], [0, 0, 1]], &look.radii, &vec![32.0; 119]), &BondStyle::default());
        let camera = OrbitCamera { yaw: 0.0, pitch: 0.0, roll: 0.0, distance: 20.0, ..Default::default() };
        let frame = build_frame_view(&camera, &scene, true, 400, 300, &DisplaySettings::default());
        let xs: Vec<f32> = screen_positions(&scene, &frame).map(|(_, p)| p[0]).collect();
        assert_eq!(xs.len(), 5);
        assert!(xs.windows(2).all(|w| w[1] > w[0]), "{xs:?}");
        // A box around the middle three.
        let mut selection = vec![4];
        select_box(&mut selection, &scene, &frame, [(xs[0] + xs[1]) / 2.0, 0.0], [(xs[3] + xs[4]) / 2.0, 300.0], false);
        assert_eq!(selection, vec![1, 2, 3]);
        // Additive: the last one joins, nothing twice.
        select_box(&mut selection, &scene, &frame, [(xs[3] + xs[4]) / 2.0, 0.0], [400.0, 300.0], true);
        assert_eq!(selection, vec![1, 2, 3, 4]);
        // A triangle around the first two.
        let tip = |k: usize| (xs[k] + xs[k + 1]) / 2.0;
        select_lasso(&mut selection, &scene, &frame, &[[0.0, 0.0], [tip(0), 150.0], [0.0, 300.0]], false);
        assert_eq!(selection, vec![0]);
        select_lasso(&mut selection, &scene, &frame, &[[0.0, -300.0], [tip(1), 150.0], [0.0, 600.0]], false);
        assert_eq!(selection, vec![0, 1]);
    }
}
