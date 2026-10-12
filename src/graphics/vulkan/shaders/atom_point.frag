#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
// A point-sprite atom shaded as the sphere it stands for (the normal from the
// position in the point); its depth is the centre's, so early depth tests stay on.
layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragCentre;
layout(location = 2) in float fragRadius;
layout(location = 3) in float fragShininess;
layout(location = 4) in vec4 fragPosLight;
layout(location = 0) out vec4 color;

void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float r2 = dot(c, c);
    if (r2 > 1.0) discard;
    vec3 viewDir = normalize(frame.viewPos.xyz - fragCentre);
    vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
    vec3 up = normalize(cross(viewDir, right));
    vec3 norm = normalize(right * c.x + up * c.y + viewDir * sqrt(1.0 - r2));
    vec3 worldPos = fragCentre + norm * fragRadius;
    color = vec4(shadeTaps(norm, worldPos, fragColor, fragShininess, fragPosLight, 1), 1.0);
}
