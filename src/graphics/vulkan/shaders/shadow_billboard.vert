#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 instancePos;
layout(location = 3) in float instanceScale;
layout(location = 0) out vec2 quadUv;
void main()
{
    vec3 viewDir = normalize(instancePos - frame.viewPos.xyz);
    vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
    vec3 up = normalize(cross(viewDir, right));
    vec3 worldPos = instancePos + (right * position.x + up * position.y) * instanceScale;
    gl_Position = frame.lightClip * vec4(worldPos, 1.0);
    quadUv = position.xy;
}
