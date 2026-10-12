#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// The farthest level of detail: an atom a few pixels wide is one point sprite.
layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragCentre;
layout(location = 2) out float fragRadius;
layout(location = 3) out float fragShininess;
layout(location = 4) out vec4 fragPosLight;

invariant gl_Position;  // the same depth in the depth pre-pass and the shading pass
void main()
{
    uint i = instance();
    vec4 atom = atoms[i];
    vec4 colour = atomColors[i];
    gl_Position = frame.projection * frame.view * vec4(atom.xyz, 1.0);
    gl_PointSize = max(1.0, 2.0 * atom.w * push.extra.x / max(gl_Position.w, 1e-4));
    fragPosLight = frame.lightMVP * vec4(atom.xyz, 1.0);
    fragColor = colour.rgb;
    fragCentre = atom.xyz;
    fragRadius = atom.w;
    fragShininess = colour.a;
}
