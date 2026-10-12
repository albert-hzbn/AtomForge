#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// shadow_bond.vert for GPU-culled bonds.
layout(location = 0) in vec3 position;
void main()
{
    uint i = instance();
    vec4 s = bondStarts[i], e = bondEnds[i];
    vec3 axis = e.xyz - s.xyz;
    float lengthAxis = length(axis);
    vec3 dir = (lengthAxis > 1e-5) ? axis / lengthAxis : vec3(0.0, 0.0, 1.0);
    vec3 helper = (abs(dir.y) < 0.999) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(helper, dir));
    vec3 bitangent = cross(dir, tangent);
    vec3 local = vec3(position.xy * s.w, position.z * lengthAxis);
    vec3 center = 0.5 * (s.xyz + e.xyz);
    gl_Position = frame.lightClip * vec4(center + tangent * local.x + bitangent * local.y + dir * local.z, 1.0);
}
