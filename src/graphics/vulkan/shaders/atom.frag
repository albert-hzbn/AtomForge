#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragWorldPos;
layout(location = 2) in vec3 fragNormal;
layout(location = 3) in float fragShininess;
layout(location = 4) in vec4 fragPosLight;
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(shade(normalize(fragNormal), fragWorldPos, fragColor, fragShininess, fragPosLight, true), 1.0);
}
