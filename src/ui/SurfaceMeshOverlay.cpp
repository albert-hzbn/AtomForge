#include "ui/SurfaceMeshOverlay.h"
#include "ui/OverlayProjection.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

SurfaceMeshOverlay& surfaceMeshOverlay()
{
    static SurfaceMeshOverlay overlay;
    return overlay;
}

std::vector<glm::vec3> trianglesFromObj(const std::string& text)
{
    std::vector<glm::vec3> points, triangles;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream words(line);
        std::string tag;
        words >> tag;
        if (tag == "v") {
            glm::vec3 p;
            if (!(words >> p.x >> p.y >> p.z)) throw std::runtime_error("Malformed OBJ vertex: " + line);
            points.push_back(p);
        } else if (tag == "f") {
            std::vector<long> face;
            std::string item;
            while (words >> item) {
                const long index = std::stol(item.substr(0, item.find('/')));
                face.push_back(index < 0 ? static_cast<long>(points.size()) + index : index - 1);
            }
            for (std::size_t k = 1; k + 1 < face.size(); ++k)
                for (long index : {face[0], face[k], face[k + 1]}) {
                    if (index < 0 || static_cast<std::size_t>(index) >= points.size()) throw std::runtime_error("OBJ face refers to a missing vertex");
                    triangles.push_back(points[static_cast<std::size_t>(index)]);
                }
        }
    }
    return triangles;
}

void drawSurfaceMeshOverlay(ImDrawList* drawList, const glm::mat4& projection, const glm::mat4& view,
                            int framebufferWidth, int framebufferHeight)
{
    const SurfaceMeshOverlay& overlay = surfaceMeshOverlay();
    if (!overlay.visible || overlay.vertices.size() < 3 || !drawList) return;
    struct Triangle { float depth; ImVec2 corners[3]; float shade; };
    std::vector<Triangle> triangles;
    triangles.reserve(overlay.vertices.size() / 3);
    const glm::mat4 viewProjection = projection * view;
    for (std::size_t t = 0; t + 2 < overlay.vertices.size(); t += 3) {
        Triangle tri{};
        bool visible = true;
        glm::vec3 eye[3];
        for (int c = 0; c < 3; ++c) {
            const glm::vec3& world = overlay.vertices[t + static_cast<std::size_t>(c)];
            if (!projectToView(viewProjection, world, framebufferWidth, framebufferHeight, tri.corners[c])) { visible = false; break; }
            eye[c] = glm::vec3(view * glm::vec4(world, 1.0f));
        }
        if (!visible) continue;
        const glm::vec3 normal = glm::cross(eye[1] - eye[0], eye[2] - eye[0]);
        const float length = glm::length(normal);
        if (length <= 0) continue;
        // Lambert-like shading from the viewer; both sides lit alike.
        tri.shade = 0.35f + 0.65f * std::abs(normal.z / length);
        tri.depth = (eye[0].z + eye[1].z + eye[2].z) / 3.0f;
        triangles.push_back(tri);
    }
    // Far (most negative eye z) first.
    std::sort(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) { return a.depth < b.depth; });
    const int alpha = static_cast<int>(std::clamp(overlay.opacity, 0.0f, 1.0f) * 255.0f);
    for (const Triangle& tri : triangles) {
        const ImU32 colour = IM_COL32(static_cast<int>(overlay.colour.r * tri.shade * 255.0f), static_cast<int>(overlay.colour.g * tri.shade * 255.0f),
                                      static_cast<int>(overlay.colour.b * tri.shade * 255.0f), alpha);
        drawList->AddTriangleFilled(tri.corners[0], tri.corners[1], tri.corners[2], colour);
    }
}
