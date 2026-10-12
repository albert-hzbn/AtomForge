// Shared by the Vulkan scene shaders: the frame's uniforms, the shadow lookup
// and the six-light shading of the OpenGL renderer (Renderer.cpp), unchanged.

layout(set = 0, binding = 0) uniform Frame
{
    mat4 projection;      // OpenGL projection with the Vulkan clip correction
    mat4 view;
    mat4 lightMVP;        // OpenGL light matrix: the shadow lookup is the same as in OpenGL
    mat4 lightClip;       // lightMVP with depth mapped to [0, 1], for drawing the shadow map
    vec4 frustumPlanes[6];
    vec4 lightPos;
    vec4 viewPos;
    vec4 params0;         // ambient, saturation, contrast, specular intensity
    vec4 params1;         // shadow strength, shininess scale, shininess floor, shadow texel size
} frame;

layout(set = 0, binding = 1) uniform sampler2DShadow shadowMap;

float computeShadow(vec4 pos, float cosTheta)
{
    vec3 proj = pos.xyz / pos.w;
    proj = proj * 0.5 + 0.5;
    if (proj.x < 0.0 || proj.x > 1.0 || proj.y < 0.0 || proj.y > 1.0 || proj.z < 0.0 || proj.z > 1.0)
        return 0.0;
    // Slope-scale bias: larger at grazing angles to prevent self-shadowing acne.
    float bias = max(0.012 * (1.0 - cosTheta), 0.002);
    float shadow = 0.0;
    vec2 texelSize = vec2(frame.params1.w);
    for (int x = -2; x <= 2; ++x)
        for (int y = -2; y <= 2; ++y)
            shadow += 1.0 - texture(shadowMap, vec3(proj.xy + vec2(x, y) * texelSize, proj.z - bias));
    return shadow / 25.0;
}

// Six fill lights at symmetric cube-octant positions; with shadow from the key light.
// One tap (the compare sampler filters 2 x 2 texels): for atoms a few pixels
// wide, where the 5 x 5 kernel above is invisible and its cost is not.
float computeShadowFast(vec4 pos, float cosTheta)
{
    vec3 proj = pos.xyz / pos.w;
    proj = proj * 0.5 + 0.5;
    if (proj.x < 0.0 || proj.x > 1.0 || proj.y < 0.0 || proj.y > 1.0 || proj.z < 0.0 || proj.z > 1.0)
        return 0.0;
    float bias = max(0.012 * (1.0 - cosTheta), 0.002);
    return 1.0 - texture(shadowMap, vec3(proj.xy, proj.z - bias));
}

// shadowTaps: 0 no shadow, 1 one tap, else the 5 x 5 kernel.
vec3 shadeTaps(vec3 norm, vec3 worldPos, vec3 baseColor, float shininess, vec4 posLight, int shadowTaps)
{
    vec3 lightPos = frame.lightPos.xyz;
    vec3 lights[6] = vec3[6](lightPos,
                             vec3(-lightPos.x,  lightPos.y,  lightPos.z),
                             vec3( lightPos.x, -lightPos.y,  lightPos.z),
                             vec3(-lightPos.x, -lightPos.y,  lightPos.z),
                             vec3( lightPos.x,  lightPos.y, -lightPos.z),
                             vec3(-lightPos.x,  lightPos.y, -lightPos.z));
    float weights[6] = float[6](0.22, 0.17, 0.17, 0.15, 0.15, 0.14);
    vec3 viewDir = normalize(frame.viewPos.xyz - worldPos);
    float specPower = max(shininess * frame.params1.y, frame.params1.z);
    float diff = 0.0, spec = 0.0;
    for (int i = 0; i < 6; ++i)
    {
        vec3 lightDir = normalize(lights[i] - worldPos);
        vec3 halfDir = normalize(lightDir + viewDir);
        diff += weights[i] * max(dot(norm, lightDir), 0.0);
        spec += weights[i] * pow(max(dot(norm, halfDir), 0.0), specPower);
    }
    float shadow = 0.0;
    if (shadowTaps != 0)
    {
        float cosTheta = clamp(dot(norm, normalize(lightPos - worldPos)), 0.0, 1.0);
        shadow = shadowTaps == 1 ? computeShadowFast(posLight, cosTheta) : computeShadow(posLight, cosTheta);
    }
    float luma = dot(baseColor, vec3(0.2126, 0.7152, 0.0722));
    baseColor = mix(vec3(luma), baseColor, frame.params0.y);
    baseColor = clamp((baseColor - 0.5) * frame.params0.z + 0.5, 0.0, 1.0);
    float ambient = frame.params0.x;
    float shadowStrength = frame.params1.x;
    float litFactor = ambient + (1.0 - ambient) * diff * (1.0 - shadow * shadowStrength);
    return baseColor * litFactor + vec3(frame.params0.w * spec * (1.0 - shadow * shadowStrength));
}

vec3 shade(vec3 norm, vec3 worldPos, vec3 baseColor, float shininess, vec4 posLight, bool shadowed)
{
    return shadeTaps(norm, worldPos, baseColor, shininess, posLight, shadowed ? 25 : 0);
}

bool outsideFrustum(vec3 center, float radius)
{
    for (int i = 0; i < 6; ++i)
        if (dot(frame.frustumPlanes[i].xyz, center) + frame.frustumPlanes[i].w < -radius) return true;
    return false;
}
