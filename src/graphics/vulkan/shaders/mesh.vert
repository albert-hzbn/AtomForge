#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
// Lit triangles of the dialog previews (custom shapes, Wulff and cutting planes).
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 0) out vec3 vNormal;
layout(location = 1) out vec3 vWorldPos;
void main()
{
    vWorldPos = position;
    vNormal = normal;
    gl_Position = frame.projection * frame.view * vec4(position, 1.0);
}
