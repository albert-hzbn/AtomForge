#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 instancePos;
layout(location = 2) in vec3 instanceColor;
layout(location = 3) in float instanceScale;
layout(location = 4) in float instanceShininess;
layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragWorldCenter;
layout(location = 2) out float fragRadius;
layout(location = 3) out float fragShininess;
layout(location = 4) out vec2 quadUv;
layout(location = 5) out vec4 fragPosLight;
void main()
{
    if (outsideFrustum(instancePos, instanceScale))
    {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        fragColor = vec3(0.0); fragWorldCenter = vec3(0.0); fragRadius = 0.0;
        fragShininess = instanceShininess; quadUv = vec2(0.0); fragPosLight = vec4(0.0);
        return;
    }
    vec3 viewDir = normalize(instancePos - frame.viewPos.xyz);
    vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
    vec3 up = normalize(cross(viewDir, right));
    vec3 worldPos = instancePos + (right * position.x + up * position.y) * instanceScale;
    gl_Position = frame.projection * frame.view * vec4(worldPos, 1.0);
    fragPosLight = frame.lightMVP * vec4(instancePos, 1.0);
    fragColor = instanceColor;
    fragWorldCenter = instancePos;
    fragRadius = instanceScale;
    fragShininess = instanceShininess;
    quadUv = position.xy;
}
