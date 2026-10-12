#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 fragColorA;
layout(location = 1) in vec3 fragColorB;
layout(location = 2) in float fragAxis;
layout(location = 3) in vec3 fragWorldPos;
layout(location = 4) in vec3 fragNormal;
layout(location = 5) in float fragShininessA;
layout(location = 6) in float fragShininessB;
layout(location = 0) out vec4 color;
void main()
{
    vec3 baseColor = (fragAxis < 0.0) ? fragColorA : fragColorB;
    float shininess = (fragAxis < 0.0) ? fragShininessA : fragShininessB;
    color = vec4(shade(normalize(fragNormal), fragWorldPos, baseColor, shininess, vec4(0.0), false), 1.0);
}
