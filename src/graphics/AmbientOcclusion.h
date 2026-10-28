#pragma once

// Ambient occlusion per atom: the fraction of the sky an atom's neighbours
// hide from it, used to darken buried atoms (Display Settings > Lighting).

#include <glm/glm.hpp>

#include <cstddef>
#include <vector>

namespace atomforge::display
{
// Atoms are spheres (centres, radii). From each centre, rays along evenly
// spread directions count as blocked when they meet another sphere within
// `reach` of the atom's surface; the result is the blocked fraction (0 open,
// 1 enclosed) of each atom. Empty when there are more than maxAtoms atoms.
std::vector<float> ambientOcclusion(const std::vector<glm::vec3>& centres, const std::vector<float>& radii,
                                    float reach, std::size_t maxAtoms = 2000000);

// Darkens each colour by strength x its occlusion, with rays reaching four
// times the mean radius: past the first neighbours of ball-and-stick atoms
// (drawn at 0.6 of their element radii) to the next shells (nothing for
// strength 0 or mismatched sizes).
void darkenOccluded(const std::vector<glm::vec3>& centres, const std::vector<float>& radii,
                    std::vector<glm::vec3>& colours, float strength);

// The directions sampled (unit vectors spread over the sphere).
const std::vector<glm::vec3>& occlusionDirections();
}  // namespace atomforge::display
