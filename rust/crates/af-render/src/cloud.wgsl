// The atom-cloud renderer's shaders (the C++ graphics/CloudRenderer GLSL):
// each atom is a camera-facing quad (wgpu has no point sizes) drawn as a lit
// sphere impostor that writes the depth of its surface. Atoms are 8 bytes
// (x, y, z quantized in their chunk's box, species). Matrices are glm's
// (clip depth -1..1); to_clip maps them to wgpu's 0..1 depth.

struct Frame {
    view: mat4x4<f32>,
    projection: mat4x4<f32>,
    // light direction (view space, normalised), ambient
    light: vec4<f32>,
    // pixels per unit (projection[1][1] * height / 2), orthographic, radius scale, colour axis
    params: vec4<f32>,
    // colour range along the axis (lower, upper)
    colour_range: vec4<f32>,
    // box line colour
    line_colour: vec4<f32>,
};

struct Chunk {
    // xyz: the chunk box's lower corner; w: the point enlargement of a subsample
    lower: vec4<f32>,
    // xyz: (upper - lower) / 65535
    scale: vec4<f32>,
};

@group(0) @binding(0) var<uniform> frame: Frame;
// 256 species per row: rgb colour, a radius
@group(0) @binding(1) var species: texture_2d<f32>;
@group(1) @binding(0) var<uniform> chunk: Chunk;

fn to_clip(p: vec4<f32>) -> vec4<f32> {
    return vec4<f32>(p.x, p.y, 0.5 * p.z + 0.5 * p.w, p.w);
}

fn ramp(t_in: f32) -> vec3<f32> {
    let t = clamp(t_in, 0.0, 1.0);
    return clamp(vec3<f32>(1.5 - abs(4.0 * t - 3.0), 1.5 - abs(4.0 * t - 2.0), 1.5 - abs(4.0 * t - 1.0)), vec3<f32>(0.0), vec3<f32>(1.0));
}

struct PointOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) colour: vec3<f32>,
    @location(1) centre: vec3<f32>,
    @location(2) radius: f32,
    @location(3) corner: vec2<f32>,
};

@vertex
fn point_vs(@builtin(vertex_index) index: u32, @location(0) packed: vec4<u32>) -> PointOut {
    var o: PointOut;
    let p = chunk.lower.xyz + vec3<f32>(packed.xyz) * chunk.scale.xyz;
    let s = i32(packed.w);
    let sp = textureLoad(species, vec2<i32>(s & 255, s >> 8u), 0);
    o.colour = sp.rgb;
    let axis = i32(frame.params.w);
    if (axis > 0) {
        o.colour = ramp((p[axis - 1] - frame.colour_range.x) / max(frame.colour_range.y - frame.colour_range.x, 1e-6));
    }
    o.radius = sp.a * frame.params.z;
    let view_position = frame.view * vec4<f32>(p, 1.0);
    o.centre = view_position.xyz;
    // The quad: the GL point sprite of side 2 * pixels * fill (at least a pixel).
    var units_per_pixel = 1.0 / frame.params.x;
    if (frame.params.y < 0.5) {
        units_per_pixel = max(-view_position.z, 1e-3) / frame.params.x;
    }
    let half = max(o.radius * chunk.lower.w, 0.5 * units_per_pixel);
    var corner = vec2<f32>(-1.0, -1.0);
    if ((index & 1u) != 0u) {
        corner.x = 1.0;
    }
    if ((index & 2u) != 0u) {
        corner.y = 1.0;
    }
    o.corner = corner;
    var clip = frame.projection * vec4<f32>(view_position.xy + corner * half, view_position.z, 1.0);
    // Rasterise at the sphere's front, so that the depth written later is
    // never nearer than this (as the C++).
    let front = frame.projection * vec4<f32>(view_position.xy, view_position.z + o.radius, 1.0);
    clip.z = front.z / front.w * clip.w;
    o.clip = to_clip(clip);
    return o;
}

struct PointFragment {
    @location(0) colour: vec4<f32>,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn point_fs(i: PointOut) -> PointFragment {
    let c = i.corner;
    let d2 = dot(c, c);
    if (d2 > 1.0) {
        discard;
    }
    let n = vec3<f32>(c, sqrt(1.0 - d2));
    // Depth of the sphere surface, so that impostors intersect correctly.
    let clip = frame.projection * vec4<f32>(i.centre + n * i.radius, 1.0);
    var o: PointFragment;
    o.depth = clamp(clip.z / clip.w * 0.5 + 0.5, 0.0, 1.0);
    let light = frame.light.xyz;
    let ambient = frame.light.w;
    let diffuse = max(dot(n, light), 0.0);
    let specular = pow(max(dot(reflect(-light, n), vec3<f32>(0.0, 0.0, 1.0)), 0.0), 40.0) * 0.35;
    var colour = i.colour * (ambient + (1.0 - ambient) * diffuse) + vec3<f32>(specular);
    // Darken the rim for depth.
    colour = colour * mix(0.72, 1.0, n.z);
    o.colour = vec4<f32>(colour, 1.0);
    return o;
}

@vertex
fn line_vs(@location(0) position: vec3<f32>) -> @builtin(position) vec4<f32> {
    return to_clip(frame.projection * frame.view * vec4<f32>(position, 1.0));
}

@fragment
fn line_fs() -> @location(0) vec4<f32> {
    return vec4<f32>(frame.line_colour.xyz, 1.0);
}

// Occlusion queries: a chunk's padded box (unit cube scaled to a..b).
@vertex
fn box_vs(@location(0) unit: vec3<f32>, @location(1) a: vec3<f32>, @location(2) b: vec3<f32>) -> @builtin(position) vec4<f32> {
    return to_clip(frame.projection * frame.view * vec4<f32>(a + unit * (b - a), 1.0));
}

@fragment
fn box_fs() -> @location(0) vec4<f32> {
    return vec4<f32>(0.0);
}
