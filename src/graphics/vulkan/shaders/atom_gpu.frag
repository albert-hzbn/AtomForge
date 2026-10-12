#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
// atom.frag for GPU-culled atoms: coarse levels (atoms a few pixels wide)
// take one shadow tap instead of the 5 x 5 kernel (push.extra.y = 1).
layout(push_constant) uniform Push { uvec4 list; vec4 extra; } push;
layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragWorldPos;
layout(location = 2) in vec3 fragNormal;
layout(location = 3) in float fragShininess;
layout(location = 4) in vec4 fragPosLight;
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(shadeTaps(normalize(fragNormal), fragWorldPos, fragColor, fragShininess, fragPosLight, push.extra.y > 0.5 ? 1 : 25), 1.0);
}
