#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 instancePos;
layout(location = 3) in float instanceScale;
void main()
{
    gl_Position = frame.lightClip * vec4(position * instanceScale + instancePos, 1.0);
}
