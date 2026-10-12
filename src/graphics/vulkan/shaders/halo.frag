#version 450
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) in vec2 vLocal;
layout(location = 0) out vec4 fragColor;
void main()
{
    float extent = push.extra.x;
    float opacity = push.extra.y;
    float d = length(vLocal);
    if (d < 0.985 || d > extent) discard;
    float pixel = fwidth(d);
    float width = max(0.12, 2.0 * pixel);
    float inner = 1.0, outer = 1.0 + width;
    float ring = smoothstep(inner - pixel, inner + pixel, d) * (1.0 - smoothstep(outer - pixel, outer + pixel, d));
    float edge = smoothstep(outer - pixel, outer, d) * (1.0 - smoothstep(outer + pixel, outer + 2.0 * pixel, d));
    float glow = 0.38 * pow(1.0 - smoothstep(outer, extent, d), 2.0) * step(outer, d);
    vec3 color = mix(push.color.rgb, vec3(0.06), edge * 0.75);
    float alpha = max(max(ring, edge * 0.75), glow);
    if (alpha < 0.01) discard;
    fragColor = vec4(color, alpha * opacity);
}
