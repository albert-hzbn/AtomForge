#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// billboard.vert for GPU-culled atoms: ray-cast sphere impostors, the quad
// at the sphere's front (toward the camera; it still covers the outline).
layout(location = 0) in vec3 position;
layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragWorldCenter;
layout(location = 2) out float fragRadius;
layout(location = 3) out float fragShininess;
layout(location = 4) out vec2 quadUv;
layout(location = 5) out vec4 fragPosLight;

invariant gl_Position;  // the same depth in the depth pre-pass and the shading pass
void main()
{
    uint i = instance();
    vec4 atom = atoms[i];
    vec4 colour = atomColors[i];
    vec3 viewDir = normalize(atom.xyz - frame.viewPos.xyz);
    vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
    vec3 up = normalize(cross(viewDir, right));
    vec3 worldPos = atom.xyz + (right * position.x + up * position.y - viewDir) * atom.w;
    gl_Position = frame.projection * frame.view * vec4(worldPos, 1.0);
    fragPosLight = frame.lightMVP * vec4(atom.xyz, 1.0);
    fragColor = colour.rgb;
    fragWorldCenter = atom.xyz;
    fragRadius = atom.w;
    fragShininess = colour.a;
    quadUv = position.xy;
}
