// Ray-cast sphere impostors: one camera-facing quad per atom, the sphere
// surface and its depth computed per fragment (orthographic or perspective).

struct Camera {
    view: mat4x4<f32>,
    projection: mat4x4<f32>,
    light: vec4<f32>,       // direction to the light (view space), w unused
    params: vec4<f32>,      // x: 1 for perspective, 0 for orthographic
};
@group(0) @binding(0) var<uniform> camera: Camera;

struct Instance {
    @location(0) centre: vec3<f32>,
    @location(1) radius: f32,
    @location(2) colour: vec4<f32>,
};

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) centre_view: vec3<f32>,
    @location(2) radius: f32,
    @location(3) colour: vec4<f32>,
    @location(4) point_view: vec3<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vertex: u32, instance: Instance) -> VertexOut {
    var corners = array<vec2<f32>, 4>(vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, -1.0), vec2<f32>(-1.0, 1.0), vec2<f32>(1.0, 1.0));
    let corner = corners[vertex];
    let centre = (camera.view * vec4<f32>(instance.centre, 1.0)).xyz;
    // A quad slightly larger than the sphere's silhouette, in front of it.
    let grow = select(1.0, 1.5, camera.params.x > 0.5);
    let point = centre + vec3<f32>(corner * instance.radius * grow, instance.radius);
    var out: VertexOut;
    out.clip = camera.projection * vec4<f32>(point, 1.0);
    out.corner = corner;
    out.centre_view = centre;
    out.radius = instance.radius;
    out.colour = instance.colour;
    out.point_view = point;
    return out;
}

struct FragmentOut {
    @location(0) colour: vec4<f32>,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn fs_main(in: VertexOut) -> FragmentOut {
    // Ray from the eye (perspective) or along -z (orthographic) through the quad point.
    var origin = vec3<f32>(in.point_view.xy, 0.0);
    var direction = vec3<f32>(0.0, 0.0, -1.0);
    if (camera.params.x > 0.5) {
        origin = vec3<f32>(0.0);
        direction = normalize(in.point_view);
    }
    let oc = origin - in.centre_view;
    let b = dot(oc, direction);
    let c = dot(oc, oc) - in.radius * in.radius;
    let h = b * b - c;
    if (h < 0.0) {
        discard;
    }
    let t = -b - sqrt(h);
    let hit = origin + t * direction;
    let normal = normalize(hit - in.centre_view);
    let light = normalize(camera.light.xyz);
    let diffuse = max(dot(normal, light), 0.0);
    let half_vector = normalize(light + vec3<f32>(0.0, 0.0, 1.0));
    let specular = pow(max(dot(normal, half_vector), 0.0), 40.0);
    let rgb = in.colour.rgb * (0.25 + 0.75 * diffuse) + vec3<f32>(0.35 * specular);
    let clip = camera.projection * vec4<f32>(hit, 1.0);
    var out: FragmentOut;
    out.colour = vec4<f32>(rgb, 1.0);
    out.depth = clip.z / clip.w;
    return out;
}
