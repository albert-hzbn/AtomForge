#version 450
// Isosurfaces of the electronic-structure viewport. The orthographic
// projection (corrected for Vulkan) is a per-axis scale plus a depth offset.
layout(push_constant) uniform Push
{
    mat4 view;
    vec4 clip;    // x, y, z scales, z offset
    vec4 params;  // low, high, opacity, specular
    float shininess;
    int palette;
} pc;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in float scalar;
layout(location = 0) out vec3 eyePosition;
layout(location = 1) out vec3 eyeNormal;
layout(location = 2) out float value;
void main()
{
    vec4 p = pc.view * vec4(position,1.0);
    eyePosition = p.xyz;
    eyeNormal = mat3(pc.view) * normal;
    value = scalar;
    gl_Position = vec4(p.x*pc.clip.x, p.y*pc.clip.y, p.z*pc.clip.z + pc.clip.w*p.w, p.w);
}
