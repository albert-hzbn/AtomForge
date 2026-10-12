#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// shadow_atom.vert for GPU-culled atoms.
layout(location = 0) in vec3 position;
void main()
{
    vec4 atom = atoms[instance()];
    gl_Position = frame.lightClip * vec4(position * atom.w + atom.xyz, 1.0);
}
