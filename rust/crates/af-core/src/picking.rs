//! Picking atoms under the cursor (the C++ graphics/Picking). The projection
//! uses the wgpu depth range (0 at the near plane, 1 at the far plane), as
//! af-render's cameras do; the C++ used OpenGL's -1 to 1.

use glam::{Mat4, Vec3, Vec4};

fn ray_endpoints(mx: f64, my: f64, w: i32, h: i32, projection: &Mat4, view: &Mat4) -> (Vec3, Vec3) {
    let ndc_x = 2.0 * mx as f32 / w as f32 - 1.0;
    let ndc_y = -2.0 * my as f32 / h as f32 + 1.0;
    let inverse = (*projection * *view).inverse();
    let near = inverse * Vec4::new(ndc_x, ndc_y, 0.0, 1.0);
    let far = inverse * Vec4::new(ndc_x, ndc_y, 1.0, 1.0);
    (near.truncate() / near.w, far.truncate() / far.w)
}

/// The cursor's ray on the near plane, in world space.
pub fn pick_ray_origin(mx: f64, my: f64, w: i32, h: i32, projection: &Mat4, view: &Mat4) -> Vec3 {
    ray_endpoints(mx, my, w, h, projection, view).0
}

/// The cursor's ray direction (unit), in world space.
pub fn pick_ray_dir(mx: f64, my: f64, w: i32, h: i32, projection: &Mat4, view: &Mat4) -> Vec3 {
    let (near, far) = ray_endpoints(mx, my, w, h, projection, view);
    (far - near).normalize()
}

/// The atom whose sphere surface the ray meets first, or None. Radii come
/// from `radii`, or `fallback_radius` where missing or not positive.
pub fn pick_atom(origin: Vec3, dir: Vec3, positions: &[Vec3], radii: &[f32], fallback_radius: f32) -> Option<usize> {
    let length = dir.length();
    if !length.is_finite() || length <= 0.0 {
        return None;
    }
    let ray = dir / length;
    let mut best = None;
    let mut best_t = 1e30f32;
    for (i, &p) in positions.iter().enumerate() {
        let oc = p - origin;
        let t = oc.dot(ray);
        let radius = match radii.get(i) {
            Some(&r) if r > 0.0 => r,
            _ => fallback_radius,
        };
        let d2 = oc.dot(oc) - t * t;
        if radius <= 0.0 || d2 > radius * radius {
            continue;
        }
        let half_chord = (radius * radius - d2).max(0.0).sqrt();
        if t + half_chord < 0.0 {
            continue;
        }
        // Compare visible surfaces, since larger atoms may cover nearer centres.
        let hit = (t - half_chord).max(0.0);
        if hit < best_t {
            best = Some(i);
            best_t = hit;
        }
    }
    best
}
