#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 bondStart;
layout(location = 2) in vec3 bondEnd;
layout(location = 5) in float bondRadius;
void main()
{
    vec3 axis = bondEnd - bondStart;
    float lengthAxis = length(axis);
    vec3 dir = (lengthAxis > 1e-5) ? axis / lengthAxis : vec3(0.0, 0.0, 1.0);
    vec3 helper = (abs(dir.y) < 0.999) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(helper, dir));
    vec3 bitangent = cross(dir, tangent);
    vec3 local = vec3(position.xy * bondRadius, position.z * lengthAxis);
    vec3 center = 0.5 * (bondStart + bondEnd);
    gl_Position = frame.lightClip * vec4(center + tangent * local.x + bitangent * local.y + dir * local.z, 1.0);
}
