#version 450
// Atom-cloud point sprites (CloudRenderer.cpp's OpenGL shaders): sphere
// impostors from 8-byte packed atoms.
layout(set = 0, binding = 0) uniform Cloud
{
    mat4 view;
    mat4 projection;      // with the Vulkan clip correction
    vec4 light;           // view space, normalised; w = ambient
    vec4 params;          // pixels per unit, orthographic, radius scale, max point size
    vec4 colouring;       // x = colour axis (0 none, 1-3 x-z), yz = its range
} cloud;
layout(std430, set = 0, binding = 1) readonly buffer Species { vec4 species[]; };  // rgb colour, a radius
layout(push_constant) uniform Chunk { vec4 lower; vec4 scale; } chunk;  // lower.w = point enlargement

layout(location = 0) in uvec4 packedAtom;
layout(location = 0) out vec3 vColour;
layout(location = 1) out vec3 vCentre;
layout(location = 2) out float vRadius;

vec3 ramp(float t)
{
    t = clamp(t, 0.0, 1.0);
    return clamp(vec3(1.5 - abs(4.0 * t - 3.0), 1.5 - abs(4.0 * t - 2.0), 1.5 - abs(4.0 * t - 1.0)), 0.0, 1.0);
}

void main()
{
    vec3 p = chunk.lower.xyz + vec3(packedAtom.xyz) * chunk.scale.xyz;
    vec4 s = species[packedAtom.w];
    vColour = s.rgb;
    int axis = int(cloud.colouring.x + 0.5);
    if (axis > 0) vColour = ramp((p[axis - 1] - cloud.colouring.y) / max(cloud.colouring.z - cloud.colouring.y, 1e-6));
    vRadius = s.a * cloud.params.z;
    vec4 viewPosition = cloud.view * vec4(p, 1.0);
    vCentre = viewPosition.xyz;
    gl_Position = cloud.projection * viewPosition;
    // Rasterise at the sphere's front so that the depth written later is never
    // nearer than this (conservative depth keeps early depth rejection).
    vec4 front = cloud.projection * vec4(viewPosition.xy, viewPosition.z + vRadius, 1.0);
    gl_Position.z = front.z / front.w * gl_Position.w;
    float pixels = vRadius * cloud.params.x / (cloud.params.y > 0.5 ? 1.0 : max(-viewPosition.z, 1e-3));
    gl_PointSize = clamp(2.0 * pixels * chunk.lower.w, 1.0, cloud.params.w);
}
