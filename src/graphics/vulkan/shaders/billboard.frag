#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragWorldCenter;
layout(location = 2) in float fragRadius;
layout(location = 3) in float fragShininess;
layout(location = 4) in vec2 quadUv;
layout(location = 5) in vec4 fragPosLight;
layout(location = 0) out vec4 color;
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
    // The sphere's own depth (the clip correction already maps depth to [0, 1]).
    vec4 clipPos = frame.projection * frame.view * vec4(worldPos, 1.0);
    gl_FragDepth = clipPos.z / clipPos.w;
    color = vec4(shade(norm, worldPos, fragColor, fragShininess, fragPosLight, true), 1.0);
}
