#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// An atom a few shadow-map texels wide casts its shadow as a disc (a point
// sprite) through its centre: the side facing the light lies in front of
// it, so the atom does not shadow itself. push.extra.x: texels per unit.
void main()
{
    vec4 atom = atoms[instance()];
    gl_Position = frame.lightClip * vec4(atom.xyz, 1.0);
    gl_PointSize = max(1.0, 2.0 * atom.w * push.extra.x / max(gl_Position.w, 1e-4));
}
