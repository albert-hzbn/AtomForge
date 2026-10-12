#version 450
#extension GL_GOOGLE_include_directive : require
#include "electronic_palette.glsl"
layout(push_constant) uniform Push
{
    mat4 view;
    vec4 clip;
    vec4 params;  // low, high, opacity, specular
    float shininess;
    int palette;
} pc;

layout(location = 0) in vec3 eyePosition;
layout(location = 1) in vec3 eyeNormal;
layout(location = 2) in float value;
layout(location = 0) out vec4 outputColor;
void main()
{
    float low = pc.params.x, high = pc.params.y;
    float t = high > low ? clamp((value-low)/(high-low),0.0,1.0) : 0.5;
    vec3 n = normalize(eyeNormal);
    vec3 v = normalize(-eyePosition);
    if (dot(n,v) < 0.0) n = -n;
    vec3 key = normalize(vec3(-0.4,0.6,1.0));
    vec3 fill = normalize(vec3(0.7,0.1,0.6));
    float diffuse = 0.27 + 0.57*max(dot(n,key),0.0) + 0.16*max(dot(n,fill),0.0);
    float highlight = pc.params.w * pow(max(dot(n,normalize(key+v)),0.0),pc.shininess);
    outputColor = vec4(clamp(ramp(t,pc.palette)*diffuse + vec3(highlight),0.0,1.0),pc.params.z);
}
