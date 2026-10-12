#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// atom.vert for GPU-culled atoms (the mesh of the atom's level of detail).
layout(location = 0) in vec3 position;
layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragWorldPos;
layout(location = 2) out vec3 fragNormal;
layout(location = 3) out float fragShininess;
layout(location = 4) out vec4 fragPosLight;

invariant gl_Position;  // the same depth in the depth pre-pass and the shading pass
void main()
{
    uint i = instance();
    vec4 atom = atoms[i];
    vec4 colour = atomColors[i];
    vec3 worldPos = position * atom.w + atom.xyz;
    gl_Position = frame.projection * frame.view * vec4(worldPos, 1.0);
    fragPosLight = frame.lightMVP * vec4(worldPos, 1.0);
    fragColor = colour.rgb;
    fragWorldPos = worldPos;
    fragNormal = normalize(position);
    fragShininess = colour.a;
}
