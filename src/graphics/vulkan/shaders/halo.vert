#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) in vec2 corner;
layout(location = 1) in vec4 atom;   // centre xyz, radius w
layout(location = 0) out vec2 vLocal;
void main()
{
    // A camera-facing quad scaled to the sphere's silhouette (larger than the
    // radius under perspective), at the sphere's front surface.
    float extent = push.extra.x;
    mat4 view = frame.view;
    vec3 right = vec3(view[0][0], view[1][0], view[2][0]);
    vec3 up = vec3(view[0][1], view[1][1], view[2][1]);
    float silhouette = 1.0;
    if (frame.projection[3][3] < 0.5)
    {
        float distance = length((view * vec4(atom.xyz, 1.0)).xyz);
        silhouette = distance / sqrt(max(distance * distance - atom.w * atom.w, 1e-6 * distance * distance));
    }
    vec3 world = atom.xyz + (right * corner.x + up * corner.y) * atom.w * silhouette * extent;
    vLocal = corner * extent;
    vec4 eye = view * vec4(world, 1.0);
    eye.z += atom.w;
    gl_Position = frame.projection * eye;
}
