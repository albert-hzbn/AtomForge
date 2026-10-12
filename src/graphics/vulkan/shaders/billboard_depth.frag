#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
// The depth pre-pass of an impostor: billboard_gpu.frag's depth, no colour.
layout(depth_greater) out float gl_FragDepth;
layout(location = 1) in vec3 fragWorldCenter;
layout(location = 2) in float fragRadius;
layout(location = 4) in vec2 quadUv;

void main()
{
    float r2 = dot(quadUv, quadUv);
    if (r2 > 1.0) discard;
    float nz = sqrt(1.0 - r2);
    vec3 viewDir = normalize(frame.viewPos.xyz - fragWorldCenter);
    vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
    vec3 up = normalize(cross(viewDir, right));
    vec3 norm = normalize(right * quadUv.x + up * quadUv.y + viewDir * nz);
    vec3 worldPos = fragWorldCenter + norm * fragRadius;
    vec4 clipPos = frame.projection * frame.view * vec4(worldPos, 1.0);
    gl_FragDepth = clipPos.z / clipPos.w;
}
