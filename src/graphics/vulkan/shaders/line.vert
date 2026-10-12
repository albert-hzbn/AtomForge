#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) in vec3 position;
void main()
{
    gl_Position = frame.projection * frame.view * vec4(position, 1.0);
}
