// (Included by cull.comp, with subgroup aggregation, and cull_basic.comp.)
// GPU culling of atoms or bonds for one view: view-frustum culling, two-phase
// hierarchical-Z occlusion culling (Haar & Aaltonen, "GPU-Driven Rendering
// Pipelines", SIGGRAPH 2015), small-feature culling and a level of detail
// from each instance's radius on screen. Visible instances are appended to
// compact per-level lists and counted into indirect draw commands.
layout(local_size_x = 256) in;

layout(set = 0, binding = 0) uniform Cull
{
    mat4 viewProjection;  // Vulkan clip space (depth in [0, 1], y down)
    vec4 planes[6];       // world-space frustum planes, inward
    vec4 screen;          // target width, height (pixels), Hi-Z mip 0 width, height
    vec4 scale;           // pixels per unit at clip w = 1, Hi-Z levels, cull radius (px), -
    vec4 thresholds[2];   // radius (px) at and above which level i is used
    vec4 inflation[2];    // level i's mesh reaches out to inflation[i] radii
    vec4 eye;             // the camera's position; w = 1 for an orthographic projection
    vec4 forward;         // the camera's viewing direction
    uvec4 slots[8];       // per level: command, list start, forward (1) or backward (0)
    uvec4 counts;         // instances, levels, phase, Hi-Z test (1)
} cull;
layout(std430, set = 0, binding = 1) buffer Commands { uint commands[]; };  // VkDrawIndexedIndirectCommand x n
layout(set = 0, binding = 2) uniform sampler2D hiz;

layout(std430, set = 1, binding = 0) readonly buffer Atoms { vec4 atoms[]; };        // xyz, radius
layout(std430, set = 1, binding = 1) readonly buffer BondStarts { vec4 bondStarts[]; };  // xyz, radius
layout(std430, set = 1, binding = 2) readonly buffer BondEnds { vec4 bondEnds[]; };
layout(std430, set = 1, binding = 3) writeonly buffer Lists { uint lists[]; };
layout(std430, set = 1, binding = 4) buffer Visibility { uint visibility[]; };  // atoms drawn last frame

layout(push_constant) uniform Push { uint kind; } push;  // 0 atoms, 1 bonds

const uint kPhaseFirst = 0u;   // atoms visible last frame (no occlusion test)
const uint kPhaseSecond = 1u;  // the rest, against the Hi-Z of the first phase; updates the history
const uint kPhaseSingle = 2u;  // one pass (shadows, bonds)

float threshold(uint level)
{
    return cull.thresholds[level >> 2][level & 3u];
}

// True when the sphere is hidden behind the depth in the Hi-Z pyramid (the
// farthest depth per texel; level 0 is the window's resolution).
bool occluded(vec3 centre, float radius)
{
    // The box around the sphere against at most 2 x 2 texels of a coarse level.
    vec2 lo = vec2(1.0), hi = vec2(0.0);
    float nearest = 1.0;
    for (int k = 0; k < 8; ++k)
    {
        vec3 corner = centre + radius * vec3((k & 1) != 0 ? 1.0 : -1.0, (k & 2) != 0 ? 1.0 : -1.0, (k & 4) != 0 ? 1.0 : -1.0);
        vec4 clip = cull.viewProjection * vec4(corner, 1.0);
        if (clip.w <= 1e-4) return false;  // reaches behind the camera
        vec3 ndc = clip.xyz / clip.w;
        vec2 uv = ndc.xy * 0.5 + 0.5;
        lo = min(lo, uv);
        hi = max(hi, uv);
        nearest = min(nearest, ndc.z);
    }
    lo = clamp(lo, 0.0, 1.0);
    hi = clamp(hi, 0.0, 1.0);
    vec2 size = (hi - lo) * cull.screen.zw;  // in level-0 texels
    int levels = int(cull.scale.y);
    float level = clamp(ceil(log2(max(max(size.x, size.y), 1.0))), 0.0, float(levels - 1));
    ivec2 extent = textureSize(hiz, int(level));
    ivec2 a = clamp(ivec2(lo * vec2(extent)), ivec2(0), extent - 1);
    ivec2 b = clamp(ivec2(hi * vec2(extent)), ivec2(0), extent - 1);
    float farthest = 0.0;
    for (int y = a.y; y <= min(b.y, a.y + 1); ++y)
        for (int x = a.x; x <= min(b.x, a.x + 1); ++x)
            farthest = max(farthest, texelFetch(hiz, ivec2(x, y), int(level)).r);
    if (nearest > farthest) return true;

    // Else the texels of the sphere's disc on screen, at a level where it is at
    // most 8 texels in radius: the box's corners reach into the gaps between
    // atoms, but an atom behind another one along the view (columns of a
    // crystal seen along a lattice direction) has a smaller disc, inside it.
    vec4 clip = cull.viewProjection * vec4(centre, 1.0);
    if (clip.w <= 1e-4) return false;
    vec2 c = (clip.xy / clip.w * 0.5 + 0.5) * cull.screen.zw;
    // A radius around the projected centre containing the sphere's outline:
    // in perspective the outline is an ellipse, stretched away from the
    // view's centre; at angle theta from the view direction, a sphere seen
    // under the half-angle alpha reaches f (tan(theta + alpha) - tan theta).
    float f = cull.scale.x;
    float r;
    if (cull.eye.w > 0.5)
        r = radius * f;
    else
    {
        vec3 d = centre - cull.eye.xyz;
        float dist = length(d);
        if (dist <= radius * 1.001) return false;
        float theta = acos(clamp(dot(d / dist, cull.forward.xyz), -1.0, 1.0));
        float alpha = asin(radius / dist);
        if (theta + alpha >= 1.5) return false;
        r = f * max(tan(theta + alpha) - tan(theta), tan(alpha) / cos(theta));
    }
    r += 0.75;
    int disc = int(max(0.0, ceil(log2(r / 8.0))));
    if (disc >= levels) return false;
    float scale = exp2(float(disc));
    vec2 cl = c / scale;
    float rl = r / scale;
    extent = textureSize(hiz, disc);
    a = clamp(ivec2(floor(cl - rl)), ivec2(0), extent - 1);
    b = clamp(ivec2(floor(cl + rl)), ivec2(0), extent - 1);
    for (int y = a.y; y <= b.y; ++y)
        for (int x = a.x; x <= b.x; ++x)
        {
            vec2 q = clamp(cl, vec2(x, y), vec2(x + 1, y + 1));  // the texel's point nearest the centre
            if (distance(q, cl) > rl) continue;
            if (texelFetch(hiz, ivec2(x, y), disc).r >= nearest) return false;
        }
    return true;
}

// The level an instance is drawn at, or none (~0u); updates the history.
uint classify(uint i)
{
    vec3 centre;
    float radius;       // of the sphere (atoms) or around the bond
    float pixelRadius;  // whose size on screen sets the level
    if (push.kind == 0u)
    {
        vec4 atom = atoms[i];
        centre = atom.xyz;
        radius = atom.w;
        pixelRadius = atom.w;
    }
    else
    {
        vec4 s = bondStarts[i], e = bondEnds[i];
        centre = 0.5 * (s.xyz + e.xyz);
        radius = 0.5 * length(e.xyz - s.xyz) + s.w;
        pixelRadius = s.w;  // the level of detail follows the bond's thickness
    }
    // Level of detail from the radius on screen; coarse meshes enclose the
    // sphere, reaching out to cull.inflation[level] radii.
    vec4 clip = cull.viewProjection * vec4(centre, 1.0);
    float px = pixelRadius * cull.scale.x / max(clip.w, 1e-4);
    uint levels = cull.counts.y;
    uint level = 0u;
    while (level + 1u < levels && px < threshold(level)) ++level;
    float bound = radius * cull.inflation[level >> 2][level & 3u];

    uint phase = cull.counts.z;
    bool inside = true;
    for (int p = 0; p < 6; ++p)
        if (dot(cull.planes[p].xyz, centre) + cull.planes[p].w < -bound) inside = false;
    bool drawnBefore = phase != kPhaseSingle && visibility[i] != 0u;
    if (phase == kPhaseFirst && (!inside || !drawnBefore)) return ~0u;
    if (!inside)
    {
        if (phase == kPhaseSecond) visibility[i] = 0u;
        return ~0u;
    }
    if (phase != kPhaseFirst && cull.counts.w != 0u && occluded(centre, bound))
    {
        if (phase == kPhaseSecond) visibility[i] = 0u;
        return ~0u;
    }
    if (phase == kPhaseSecond)
    {
        visibility[i] = 1u;
        if (drawnBefore) return ~0u;  // drawn in the first phase
    }
    if (px < cull.scale.z) return ~0u;  // smaller than a feature worth drawing
    return level;
}

void append(uint i, uint level, uint n)
{
    uvec4 slot = cull.slots[level];
    lists[slot.z != 0u ? slot.y + n : slot.y - n] = i;
}
