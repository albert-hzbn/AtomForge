#version 450
#extension GL_GOOGLE_include_directive : require
#include "electronic_palette.glsl"
// Ray-marched density volume of the electronic-structure viewport.
layout(set = 0, binding = 0) uniform sampler3D density;
layout(push_constant) uniform Push
{
    mat4 transform;  // inverse cell * inverse rotation (upper 3x3)
    vec4 dimensions; // xyz, aspect
    vec4 view;       // pan, extent, opacity
    vec4 range;      // threshold, low, high
    ivec4 options;   // ray steps, palette, filter in the shader
} pc;

layout(location = 0) in vec2 screen;
layout(location = 0) out vec4 outputColor;

// GL_LINEAR with clamp to edge, where the device cannot filter 32-bit floats.
float sampleDensity(vec3 coord)
{
    if (pc.options.z == 0) return texture(density,coord).r;
    ivec3 size = ivec3(pc.dimensions.xyz);
    vec3 u = coord*pc.dimensions.xyz-0.5;
    vec3 f = fract(u);
    ivec3 a = clamp(ivec3(floor(u)),ivec3(0),size-1), b = clamp(ivec3(floor(u))+1,ivec3(0),size-1);
    float c000 = texelFetch(density,ivec3(a.x,a.y,a.z),0).r, c100 = texelFetch(density,ivec3(b.x,a.y,a.z),0).r;
    float c010 = texelFetch(density,ivec3(a.x,b.y,a.z),0).r, c110 = texelFetch(density,ivec3(b.x,b.y,a.z),0).r;
    float c001 = texelFetch(density,ivec3(a.x,a.y,b.z),0).r, c101 = texelFetch(density,ivec3(b.x,a.y,b.z),0).r;
    float c011 = texelFetch(density,ivec3(a.x,b.y,b.z),0).r, c111 = texelFetch(density,ivec3(b.x,b.y,b.z),0).r;
    float y0 = mix(mix(c000,c100,f.x),mix(c010,c110,f.x),f.y);
    float y1 = mix(mix(c001,c101,f.x),mix(c011,c111,f.x),f.y);
    return mix(y0,y1,f.z);
}

void main()
{
    vec3 background=vec3(.96,.97,.985);
    mat3 transform=mat3(pc.transform);
    vec3 dimensions=pc.dimensions.xyz;
    float aspect=pc.dimensions.w, extent=pc.view.z, opacity=pc.view.w;
    float threshold=pc.range.x, low=pc.range.y, high=pc.range.z;
    int raySteps=pc.options.x;
    vec3 origin=transform*vec3(screen.x*extent*aspect-pc.view.x,screen.y*extent-pc.view.y,2.0)+vec3(.5);
    vec3 direction=transform*vec3(0,0,-1);
    float entry=0.0, leave=1e20;
    for(int axis=0;axis<3;++axis) {
        if(abs(direction[axis])<1e-8) {
            if(origin[axis]<0.0 || origin[axis]>1.0) { outputColor=vec4(background,1); return; }
        } else {
            float a=-origin[axis]/direction[axis], b=(1.0-origin[axis])/direction[axis];
            entry=max(entry,min(a,b)); leave=min(leave,max(a,b));
        }
    }
    if(leave<=entry) { outputColor=vec4(background,1); return; }
    float stepLength=(leave-entry)/float(raySteps);
    vec4 accumulated=vec4(0);
    for(int i=0;i<1024;++i) {
        if(i>=raySteps) break;
        vec3 p=origin+direction*(entry+(float(i)+.5)*stepLength);
        float value=sampleDensity((p*(dimensions-1.0)+.5)/dimensions);
        if(value<threshold) continue;
        float t=high>low ? clamp((value-low)/(high-low),0.0,1.0) : .5;
        float opticalDepth=.65*opacity/max(1e-6,1.0-opacity);
        float alpha=opacity>=1.0 ? 1.0 : 1.0-exp(-opticalDepth*12.0*max(.05,t)*stepLength);
        accumulated.rgb+=(1.0-accumulated.a)*alpha*ramp(t,pc.options.y);
        accumulated.a+=(1.0-accumulated.a)*alpha;
        if(accumulated.a>.995) break;
    }
    outputColor=vec4(accumulated.rgb+(1.0-accumulated.a)*background,1);
}
