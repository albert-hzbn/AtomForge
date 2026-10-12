// The instance data of GPU-culled draws: atoms and bonds in storage buffers,
// reached through the visible lists written by cull.comp.
layout(std430, set = 1, binding = 0) readonly buffer Atoms { vec4 atoms[]; };          // xyz, radius
layout(std430, set = 1, binding = 1) readonly buffer AtomColors { vec4 atomColors[]; };  // rgb, shininess
layout(std430, set = 1, binding = 2) readonly buffer Lists { uint lists[]; };
layout(std430, set = 1, binding = 3) readonly buffer BondStarts { vec4 bondStarts[]; };  // xyz, radius
layout(std430, set = 1, binding = 4) readonly buffer BondEnds { vec4 bondEnds[]; };      // xyz, shininess A
layout(std430, set = 1, binding = 5) readonly buffer BondColorsA { vec4 bondColorsA[]; };  // rgb, shininess B
layout(std430, set = 1, binding = 6) readonly buffer BondColorsB { vec4 bondColorsB[]; };

// list.x: first entry, list.y: 1 forward, 0 backward; extra.x: pixels per unit at clip w = 1
layout(push_constant) uniform Push { uvec4 list; vec4 extra; } push;

uint instance()
{
    uint k = uint(gl_InstanceIndex);
    return lists[push.list.y != 0u ? push.list.x + k : push.list.x - k];
}
