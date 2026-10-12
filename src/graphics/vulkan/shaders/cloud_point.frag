#version 450
layout(set = 0, binding = 0) uniform Cloud
{
    mat4 view;
    mat4 projection;
    vec4 light;
    vec4 params;
    vec4 colouring;
} cloud;

layout(depth_greater) out float gl_FragDepth;
layout(location = 0) in vec3 vColour;
layout(location = 1) in vec3 vCentre;
layout(location = 2) in float vRadius;
layout(location = 0) out vec4 fragColour;

void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float d2 = dot(c, c);
    if (d2 > 1.0) discard;
    vec3 n = vec3(c, sqrt(1.0 - d2));
    // Depth of the sphere surface, so that impostors intersect correctly.
    vec4 clip = cloud.projection * vec4(vCentre + n * vRadius, 1.0);
    gl_FragDepth = clip.z / clip.w;
    vec3 light = cloud.light.xyz;
    float ambient = cloud.light.w;
    float diffuse = max(dot(n, light), 0.0);
    float specular = pow(max(dot(reflect(-light, n), vec3(0.0, 0.0, 1.0)), 0.0), 40.0) * 0.35;
    vec3 colour = vColour * (ambient + (1.0 - ambient) * diffuse) + vec3(specular);
    colour *= mix(0.72, 1.0, n.z);  // darken the rim for depth
    fragColour = vec4(colour, 1.0);
}
