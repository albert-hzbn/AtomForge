#pragma once
// A translucent triangle mesh drawn over the structure view, e.g. the surface
// from Analysis > Structure and Defects > Surface Mesh. Triangles are shaded
// by their orientation and drawn back to front.

#include <glm/glm.hpp>

#include <string>
#include <vector>

struct ImDrawList;

struct SurfaceMeshOverlay
{
    std::vector<glm::vec3> vertices;  // three per triangle (Angstrom)
    bool visible = false;
    float opacity = 0.45f;
    glm::vec3 colour{0.30f, 0.62f, 0.86f};
};

SurfaceMeshOverlay& surfaceMeshOverlay();

// Triangles of a Wavefront OBJ text ("v x y z" and "f a b c ..." lines;
// polygons are fanned). Throws on malformed input.
std::vector<glm::vec3> trianglesFromObj(const std::string& text);

void drawSurfaceMeshOverlay(ImDrawList* drawList, const glm::mat4& projection, const glm::mat4& view,
                            int framebufferWidth, int framebufferHeight);
