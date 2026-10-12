#version 450
#extension GL_GOOGLE_include_directive : require
#include "common.glsl"
// The shading of the OpenGL dialog shaders (CustomStructureShaders.cpp,
// CellSculptorShaders.cpp and the Wulff planes of the nanocrystal builder).
// push.extra: x = style (0 two-sided lit mesh, 1 Wulff plane, 2 cutting
// plane); for style 0, yzw = the light direction.
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) in vec3 vNormal;
layout(location = 1) in vec3 vWorldPos;
layout(location = 0) out vec4 outColor;

vec4 facet(float diffuseFloor, float shininess, float fresnelPower, float baseAmbient, float envMix,
           float specularWeight, float fresnelWeight)
{
    vec3 norm = normalize(vNormal);
    if (!gl_FrontFacing) norm = -norm;
    vec3 viewDir = normalize(frame.viewPos.xyz - vWorldPos);
    vec3 lightDir = normalize(frame.lightPos.xyz - vWorldPos);
    vec3 halfVec = normalize(lightDir + viewDir);
    float diffuse = max(max(dot(norm, lightDir), 0.0), diffuseFloor);
    float specular = pow(max(dot(norm, halfVec), 0.0), shininess);
    float fresnel = pow(1.0 - max(dot(norm, viewDir), 0.0), fresnelPower);
    vec3 reflected = reflect(-viewDir, norm);
    float skyMix = clamp(reflected.y * 0.5 + 0.5, 0.0, 1.0);
    vec3 envColor = mix(vec3(0.10, 0.11, 0.13), vec3(0.72, 0.78, 0.86), pow(skyMix, 1.35));
    vec3 baseLit = push.color.rgb * (baseAmbient + (1.0 - baseAmbient) * diffuse);
    vec3 lit = mix(baseLit, envColor, envMix * fresnel);
    return vec4(lit + vec3(specularWeight * specular + fresnelWeight * fresnel), push.color.a);
}

void main()
{
    int style = int(push.extra.x + 0.5);
    if (style == 1)
    {
        outColor = facet(0.34, 20.0, 2.6, 0.30, 0.14, 0.12, 0.08);
        return;
    }
    if (style == 2)
    {
        outColor = facet(0.30, 18.0, 2.5, 0.28, 0.12, 0.10, 0.07);
        return;
    }
    vec3 n = normalize(vNormal);
    vec3 lightDir = normalize(push.extra.yzw);
    float diff = max(dot(n, lightDir), 0.0);
    float diffBack = max(dot(-n, lightDir), 0.0);
    diff = max(diff, diffBack * 0.7);
    vec3 viewDir = normalize(frame.viewPos.xyz - vWorldPos);
    vec3 halfDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(n, halfDir), 0.0), 32.0) * 0.3;
    float specBack = pow(max(dot(-n, halfDir), 0.0), 32.0) * 0.15;
    spec = max(spec, specBack);
    outColor = vec4(push.color.rgb * (0.15 + diff * 0.75) + vec3(spec), push.color.a);
}
