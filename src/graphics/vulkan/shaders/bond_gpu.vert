#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
#include "gpu.glsl"
// bond.vert for GPU-culled bonds (the cylinder of the bond's level of detail).
layout(location = 0) in vec3 position;
layout(location = 0) out vec3 fragColorA;
layout(location = 1) out vec3 fragColorB;
layout(location = 2) out float fragAxis;
layout(location = 3) out vec3 fragWorldPos;
layout(location = 4) out vec3 fragNormal;
layout(location = 5) out float fragShininessA;
layout(location = 6) out float fragShininessB;

void main()
{
    uint i = instance();
    vec4 s = bondStarts[i], e = bondEnds[i], a = bondColorsA[i], b = bondColorsB[i];
    vec3 axis = e.xyz - s.xyz;
    float lengthAxis = length(axis);
    vec3 dir = (lengthAxis > 1e-5) ? axis / lengthAxis : vec3(0.0, 0.0, 1.0);
    vec3 helper = (abs(dir.y) < 0.999) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(helper, dir));
    vec3 bitangent = cross(dir, tangent);
    vec3 local = vec3(position.xy * s.w, position.z * lengthAxis);
    vec3 center = 0.5 * (s.xyz + e.xyz);
    vec3 worldPos = center + tangent * local.x + bitangent * local.y + dir * local.z;
    vec3 localNormal = normalize(vec3(position.xy, 0.0));
    gl_Position = frame.projection * frame.view * vec4(worldPos, 1.0);
    fragColorA = a.rgb;
    fragColorB = b.rgb;
    fragAxis = position.z;
    fragWorldPos = worldPos;
    fragNormal = normalize(tangent * localNormal.x + bitangent * localNormal.y + dir * localNormal.z);
    fragShininessA = e.w;
    fragShininessB = a.w;
}
