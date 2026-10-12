#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 instancePos;
layout(location = 2) in vec3 instanceColor;
layout(location = 3) in float instanceScale;
layout(location = 4) in float instanceShininess;
layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragWorldPos;
layout(location = 2) out vec3 fragNormal;
layout(location = 3) out float fragShininess;
layout(location = 4) out vec4 fragPosLight;
void main()
{
    if (outsideFrustum(instancePos, instanceScale))
    {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        fragColor = vec3(0.0); fragWorldPos = vec3(0.0); fragNormal = vec3(0.0, 0.0, 1.0);
        fragShininess = instanceShininess; fragPosLight = vec4(0.0);
        return;
    }
    vec3 worldPos = position * instanceScale + instancePos;
    gl_Position = frame.projection * frame.view * vec4(worldPos, 1.0);
    fragPosLight = frame.lightMVP * vec4(worldPos, 1.0);
    fragColor = instanceColor;
    fragWorldPos = worldPos;
    fragNormal = normalize(position);
    fragShininess = instanceShininess;
}
