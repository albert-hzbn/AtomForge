// The structure renderer's shaders: the C++ graphics/Renderer GLSL (atoms as
// instanced sphere meshes, bonds as two-coloured cylinders, a PCF shadow map,
// six-light shading, depth cueing and outlines, and the cell box as wide
// lines). Matrices are glm's (clip depth -1..1); to_clip maps them to wgpu's
// 0..1 depth.

struct Frame {
    projection: mat4x4<f32>,
    view: mat4x4<f32>,
    light_mvp: mat4x4<f32>,
    light_pos: vec4<f32>,
    view_pos: vec4<f32>,
    // ambient, saturation, contrast, specular intensity
    material: vec4<f32>,
    // shadow strength, shininess scale, shininess floor, unused
    material2: vec4<f32>,
    // near, far (camera distance), strength (0: off), unused
    fog: vec4<f32>,
    fog_color: vec4<f32>,
    // width (0: off)
    outline: vec4<f32>,
    outline_color: vec4<f32>,
    // rgb, width in pixels
    line: vec4<f32>,
    // framebuffer width, height
    viewport: vec4<f32>,
    // selection halo colour (rgb)
    selection: vec4<f32>,
};

@group(0) @binding(0) var<uniform> frame: Frame;
@group(0) @binding(1) var shadow_map: texture_depth_2d;
@group(0) @binding(2) var shadow_sampler: sampler_comparison;

fn to_clip(p: vec4<f32>) -> vec4<f32> {
    return vec4<f32>(p.x, p.y, 0.5 * p.z + 0.5 * p.w, p.w);
}

fn finish_colour(lit_in: vec3<f32>, world_pos: vec3<f32>, norm: vec3<f32>) -> vec3<f32> {
    var lit = lit_in;
    if (frame.outline.x > 0.0) {
        let facing = abs(dot(norm, normalize(frame.view_pos.xyz - world_pos)));
        lit = mix(lit, frame.outline_color.xyz, 1.0 - smoothstep(frame.outline.x, frame.outline.x + 0.08, facing));
    }
    if (frame.fog.z > 0.0) {
        let t = clamp((length(frame.view_pos.xyz - world_pos) - frame.fog.x) / max(frame.fog.y - frame.fog.x, 1e-4), 0.0, 1.0);
        lit = mix(lit, frame.fog_color.xyz, t * frame.fog.z);
    }
    return lit;
}

fn compute_shadow(pos: vec4<f32>, cos_theta: f32) -> f32 {
    var proj = pos.xyz / pos.w;
    proj = proj * 0.5 + 0.5;
    if (proj.x < 0.0 || proj.x > 1.0 || proj.y < 0.0 || proj.y > 1.0 || proj.z < 0.0 || proj.z > 1.0) {
        return 0.0;
    }
    // Slope-scale bias: larger at grazing angles to prevent self-shadowing acne.
    let bias = max(0.012 * (1.0 - cos_theta), 0.002);
    let texel = vec2<f32>(1.0 / 2048.0);
    // wgpu textures start at the top row; GL's at the bottom.
    let uv = vec2<f32>(proj.x, 1.0 - proj.y);
    var shadow = 0.0;
    for (var x = -2; x <= 2; x++) {
        for (var y = -2; y <= 2; y++) {
            shadow += 1.0 - textureSampleCompareLevel(shadow_map, shadow_sampler, uv + vec2<f32>(f32(x), f32(y)) * texel, proj.z - bias);
        }
    }
    return shadow / 25.0;
}

struct Lighting {
    diff: f32,
    spec: f32,
};

// Six lights at cube-octant positions: the key light (+x,+y,+z) and five fills.
fn six_lights(norm: vec3<f32>, world_pos: vec3<f32>, spec_power: f32) -> Lighting {
    let lp = frame.light_pos.xyz;
    let l0 = normalize(lp - world_pos);
    let l1 = normalize(vec3<f32>(-lp.x, lp.y, lp.z) - world_pos);
    let l2 = normalize(vec3<f32>(lp.x, -lp.y, lp.z) - world_pos);
    let l3 = normalize(vec3<f32>(-lp.x, -lp.y, lp.z) - world_pos);
    let l4 = normalize(vec3<f32>(lp.x, lp.y, -lp.z) - world_pos);
    let l5 = normalize(vec3<f32>(-lp.x, lp.y, -lp.z) - world_pos);
    let v = normalize(frame.view_pos.xyz - world_pos);
    let diff = 0.22 * max(dot(norm, l0), 0.0) + 0.17 * max(dot(norm, l1), 0.0) + 0.17 * max(dot(norm, l2), 0.0)
        + 0.15 * max(dot(norm, l3), 0.0) + 0.15 * max(dot(norm, l4), 0.0) + 0.14 * max(dot(norm, l5), 0.0);
    let s0 = pow(max(dot(norm, normalize(l0 + v)), 0.0), spec_power);
    let s1 = pow(max(dot(norm, normalize(l1 + v)), 0.0), spec_power);
    let s2 = pow(max(dot(norm, normalize(l2 + v)), 0.0), spec_power);
    let s3 = pow(max(dot(norm, normalize(l3 + v)), 0.0), spec_power);
    let s4 = pow(max(dot(norm, normalize(l4 + v)), 0.0), spec_power);
    let s5 = pow(max(dot(norm, normalize(l5 + v)), 0.0), spec_power);
    let spec = 0.22 * s0 + 0.17 * s1 + 0.17 * s2 + 0.15 * s3 + 0.15 * s4 + 0.14 * s5;
    return Lighting(diff, spec);
}

fn graded(colour: vec3<f32>) -> vec3<f32> {
    let luma = dot(colour, vec3<f32>(0.2126, 0.7152, 0.0722));
    let saturated = mix(vec3<f32>(luma), colour, frame.material.y);
    return clamp((saturated - 0.5) * frame.material.z + 0.5, vec3<f32>(0.0), vec3<f32>(1.0));
}

// ---- atoms

struct AtomIn {
    @location(0) position: vec3<f32>,
    @location(1) centre: vec3<f32>,
    @location(2) scale: f32,
    @location(3) colour: vec3<f32>,
    @location(4) shininess: f32,
};

struct AtomOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) colour: vec3<f32>,
    @location(1) world_pos: vec3<f32>,
    @location(2) normal: vec3<f32>,
    @location(3) shininess: f32,
    @location(4) light_pos: vec4<f32>,
};

@vertex
fn atom_vs(v: AtomIn) -> AtomOut {
    let world = v.position * v.scale + v.centre;
    var o: AtomOut;
    o.clip = to_clip(frame.projection * frame.view * vec4<f32>(world, 1.0));
    o.light_pos = frame.light_mvp * vec4<f32>(world, 1.0);
    o.colour = v.colour;
    o.world_pos = world;
    o.normal = normalize(v.position);
    o.shininess = v.shininess;
    return o;
}

@fragment
fn atom_fs(f: AtomOut) -> @location(0) vec4<f32> {
    let norm = normalize(f.normal);
    let l0 = normalize(frame.light_pos.xyz - f.world_pos);
    let cos_theta = clamp(dot(norm, l0), 0.0, 1.0);
    let shadow = compute_shadow(f.light_pos, cos_theta);
    let light = six_lights(norm, f.world_pos, max(f.shininess * frame.material2.y, frame.material2.z));
    let base = graded(f.colour);
    let lit_factor = frame.material.x + (1.0 - frame.material.x) * light.diff * (1.0 - shadow * frame.material2.x);
    let diffuse = base * lit_factor;
    let specular = vec3<f32>(frame.material.w * light.spec * (1.0 - shadow * frame.material2.x));
    return vec4<f32>(finish_colour(diffuse + specular, f.world_pos, norm), 1.0);
}

@vertex
fn atom_shadow_vs(v: AtomIn) -> @builtin(position) vec4<f32> {
    return to_clip(frame.light_mvp * vec4<f32>(v.position * v.scale + v.centre, 1.0));
}

// ---- bonds

struct BondIn {
    @location(0) position: vec3<f32>,
    @location(1) start: vec3<f32>,
    @location(2) radius: f32,
    @location(3) end: vec3<f32>,
    @location(4) shininess_a: f32,
    @location(5) colour_a: vec3<f32>,
    @location(6) shininess_b: f32,
    @location(7) colour_b: vec3<f32>,
};

struct BondOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) colour_a: vec3<f32>,
    @location(1) colour_b: vec3<f32>,
    @location(2) axis: f32,
    @location(3) world_pos: vec3<f32>,
    @location(4) normal: vec3<f32>,
    @location(5) shininess_a: f32,
    @location(6) shininess_b: f32,
};

fn bond_world(v: BondIn) -> array<vec3<f32>, 2> {
    let axis = v.end - v.start;
    let length_axis = length(axis);
    var dir = vec3<f32>(0.0, 0.0, 1.0);
    if (length_axis > 1e-5) {
        dir = axis / length_axis;
    }
    var helper = vec3<f32>(1.0, 0.0, 0.0);
    if (abs(dir.y) < 0.999) {
        helper = vec3<f32>(0.0, 1.0, 0.0);
    }
    let tangent = normalize(cross(helper, dir));
    let bitangent = cross(dir, tangent);
    let local = vec3<f32>(v.position.xy * v.radius, v.position.z * length_axis);
    let centre = 0.5 * (v.start + v.end);
    let world = centre + tangent * local.x + bitangent * local.y + dir * local.z;
    let n = normalize(vec3<f32>(v.position.xy, 0.0));
    let normal = normalize(tangent * n.x + bitangent * n.y + dir * n.z);
    return array<vec3<f32>, 2>(world, normal);
}

@vertex
fn bond_vs(v: BondIn) -> BondOut {
    let w = bond_world(v);
    var o: BondOut;
    o.clip = to_clip(frame.projection * frame.view * vec4<f32>(w[0], 1.0));
    o.colour_a = v.colour_a;
    o.colour_b = v.colour_b;
    o.axis = v.position.z;
    o.world_pos = w[0];
    o.normal = w[1];
    o.shininess_a = v.shininess_a;
    o.shininess_b = v.shininess_b;
    return o;
}

@fragment
fn bond_fs(f: BondOut) -> @location(0) vec4<f32> {
    let norm = normalize(f.normal);
    var colour = f.colour_b;
    var shininess = f.shininess_b;
    if (f.axis < 0.0) {
        colour = f.colour_a;
        shininess = f.shininess_a;
    }
    let light = six_lights(norm, f.world_pos, max(shininess * frame.material2.y, frame.material2.z));
    let lit_factor = frame.material.x + (1.0 - frame.material.x) * light.diff;
    let diffuse = graded(colour) * lit_factor;
    let specular = vec3<f32>(frame.material.w * light.spec);
    return vec4<f32>(finish_colour(diffuse + specular, f.world_pos, norm), 1.0);
}

@vertex
fn bond_shadow_vs(v: BondIn) -> @builtin(position) vec4<f32> {
    return to_clip(frame.light_mvp * vec4<f32>(bond_world(v)[0], 1.0));
}

// ---- cell box: GL wide lines (non-antialiased: widened across the minor axis)

struct LineIn {
    @location(0) a: vec3<f32>,
    @location(1) b: vec3<f32>,
};

@vertex
fn line_vs(@builtin(vertex_index) index: u32, v: LineIn) -> @builtin(position) vec4<f32> {
    let ca = to_clip(frame.projection * frame.view * vec4<f32>(v.a, 1.0));
    let cb = to_clip(frame.projection * frame.view * vec4<f32>(v.b, 1.0));
    let size = frame.viewport.xy;
    let pa = ca.xy / ca.w * size;
    let pb = cb.xy / cb.w * size;
    var offset = vec2<f32>(0.0, frame.line.w / size.y);
    if (abs(pb.y - pa.y) > abs(pb.x - pa.x)) {
        offset = vec2<f32>(frame.line.w / size.x, 0.0);
    }
    // Two triangles: (a-, a+, b+), (a-, b+, b-).
    let corner = array<u32, 6>(0u, 1u, 3u, 0u, 3u, 2u)[index];
    var c = ca;
    if (corner >= 2u) {
        c = cb;
    }
    var side = -1.0;
    if (corner % 2u == 1u) {
        side = 1.0;
    }
    return vec4<f32>(c.xy + side * offset * c.w, c.z, c.w);
}

@fragment
fn line_fs() -> @location(0) vec4<f32> {
    return vec4<f32>(frame.line.xyz, 1.0);
}

// Selection halo (the C++ kSelHaloVS / kSelHaloFS): a camera-facing quad per
// selected atom, scaled to the sphere's projected silhouette (larger than the
// radius under perspective) so the ring hugs the outline at any zoom; a crisp
// accent ring just outside the outline, a thin dark edge for contrast and a
// soft glow that fades out. The atom itself is untouched.
const HALO_EXTENT: f32 = 1.6;

struct HaloOut {
    @builtin(position) pos: vec4<f32>,
    @location(0) local: vec2<f32>,
};

@vertex
fn halo_vs(@builtin(vertex_index) vi: u32, @location(0) centre: vec3<f32>, @location(1) radius: f32) -> HaloOut {
    var corners = array<vec2<f32>, 6>(vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, -1.0), vec2<f32>(1.0, 1.0),
                                      vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, 1.0), vec2<f32>(-1.0, 1.0));
    let corner = corners[vi];
    let v = frame.view;
    let right = vec3<f32>(v[0][0], v[1][0], v[2][0]);
    let up = vec3<f32>(v[0][1], v[1][1], v[2][1]);
    var silhouette = 1.0;
    if (frame.projection[3][3] < 0.5) {
        let d = length((v * vec4<f32>(centre, 1.0)).xyz);
        silhouette = d / sqrt(max(d * d - radius * radius, 1e-6 * d * d));
    }
    let world = centre + (right * corner.x + up * corner.y) * radius * silhouette * HALO_EXTENT;
    // At the sphere's front surface, so atoms behind the selected one never cover its ring.
    var eye = v * vec4<f32>(world, 1.0);
    eye.z = eye.z + radius;
    var out: HaloOut;
    out.pos = to_clip(frame.projection * eye);
    out.local = corner * HALO_EXTENT;
    return out;
}

fn halo_colour(local: vec2<f32>, opacity: f32) -> vec4<f32> {
    let d = length(local);
    if (d < 0.985 || d > HALO_EXTENT) {
        discard;
    }
    let pixel = fwidth(d);
    let width = max(0.12, 2.0 * pixel);
    let inner = 1.0;
    let outer = 1.0 + width;
    let ring = smoothstep(inner - pixel, inner + pixel, d) * (1.0 - smoothstep(outer - pixel, outer + pixel, d));
    let edge = smoothstep(outer - pixel, outer, d) * (1.0 - smoothstep(outer + pixel, outer + 2.0 * pixel, d));
    let glow = 0.38 * pow(1.0 - smoothstep(outer, HALO_EXTENT, d), 2.0) * step(outer, d);
    let colour = mix(frame.selection.rgb, vec3<f32>(0.06), edge * 0.75);
    let alpha = max(max(ring, edge * 0.75), glow);
    if (alpha < 0.01) {
        discard;
    }
    return vec4<f32>(colour, alpha * opacity);
}

@fragment
fn halo_fs(in: HaloOut) -> @location(0) vec4<f32> {
    return halo_colour(in.local, 1.0);
}

// Where the halo is hidden behind other atoms: faint, so selected atoms stay findable.
@fragment
fn halo_hidden_fs(in: HaloOut) -> @location(0) vec4<f32> {
    return halo_colour(in.local, 0.3);
}
