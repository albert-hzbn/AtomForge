// The transfer palettes of the electronic-structure viewport (as ElectronicViewport.cpp).
vec3 ramp(float t, int palette)
{
    if (palette == 1) return t < 0.5 ? mix(vec3(0.10,0.28,0.80),vec3(0.96),2.0*t)
                                     : mix(vec3(0.96),vec3(0.80,0.10,0.12),2.0*t-1.0);
    if (palette == 2) return mix(vec3(0.91,0.95,1.0),vec3(0.08,0.24,0.70),t);
    if (t < 0.25) return mix(vec3(0.10,0.20,0.85),vec3(0.0,0.75,0.95),4.0*t);
    if (t < 0.5) return mix(vec3(0.0,0.75,0.95),vec3(0.12,0.8,0.3),4.0*t-1.0);
    if (t < 0.75) return mix(vec3(0.12,0.8,0.3),vec3(1.0,0.85,0.08),4.0*t-2.0);
    return mix(vec3(1.0,0.85,0.08),vec3(0.85,0.10,0.08),4.0*t-3.0);
}
