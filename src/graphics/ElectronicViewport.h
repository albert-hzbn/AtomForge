#pragma once

#include "electronic/Volume.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <cstdint>
#include <memory>
#include <vector>

struct ElectronicVulkan;

// Offscreen surface renderer. GPU objects live only while a GL context exists
// (OpenGL) or in the Vulkan context (Vulkan).
class ElectronicViewport
{
public:
    // An ImTextureID: the GL texture name, or the Vulkan image's ImGui texture.
    using Texture = std::uint64_t;
    ElectronicViewport() = default;
    ~ElectronicViewport();
    ElectronicViewport(const ElectronicViewport&) = delete;
    ElectronicViewport& operator=(const ElectronicViewport&) = delete;
    void setMesh(const atomforge::electronic::Mesh& mesh);
    void setVolume(const atomforge::electronic::Grid& grid);
    void setQuality(int quality);
    // RGBA, top row first, of the last render.
    std::vector<unsigned char> pixels() const;
    Texture renderVolume(int width, int height, float yaw, float pitch, float zoom,
                         glm::vec2 pan, float opacity, float threshold, float low, float high, int palette);
    Texture render(int width, int height, float yaw, float pitch, float zoom,
                   glm::vec2 pan, float opacity, float low, float high, int palette,
                   float specular, float shininess);
    // The texture's rows are bottom first (OpenGL): show it with uv (0,1)-(1,0), else (0,0)-(1,1).
    static bool flippedRows();
    static glm::vec3 color(float t, int palette);
    // Normalized screen coordinates for world-space guides, using the same
    // framing and camera as the selected renderer. Y increases upwards.
    glm::vec2 project(glm::dvec3 point, bool volume, float aspect, float yaw,
                      float pitch, float zoom, glm::vec2 pan) const;

private:
    struct Vertex { glm::vec3 position; glm::vec3 normal; float scalar; };
    void release();
    // Back to front for the transparent surface; true when the order changed.
    bool sortTriangles(float yaw, float pitch);
    // Vulkan (ElectronicViewportVulkan.cpp).
    Texture renderVulkan(int width, int height, float yaw, float pitch, float zoom, glm::vec2 pan,
                         float opacity, float low, float high, int palette, float specular, float shininess);
    Texture renderVolumeVulkan(int width, int height, float yaw, float pitch, float zoom, glm::vec2 pan,
                               float opacity, float threshold, float low, float high, int palette);
    void setVolumeVulkan(const std::vector<float>& values);
    std::vector<unsigned char> pixelsVulkan() const;
    GLuint m_fbo = 0, m_texture = 0, m_depth = 0, m_program = 0, m_vao = 0, m_vbo = 0, m_ebo = 0;
    int m_width = 0, m_height = 0;
    bool m_dirty = false;
    float m_sortedYaw = 100, m_sortedPitch = 100;
    std::vector<Vertex> m_vertices;
    std::vector<unsigned int> m_indices;
    GLuint m_volumeTexture=0, m_volumeProgram=0;
    glm::mat3 m_inverseCell{1};
    glm::ivec3 m_volumeShape{0};
    glm::dvec3 m_meshCenter{0}, m_volumeCenter{0};
    double m_meshRadius=1, m_volumeRadius=1;
    int m_volumeResolution=96, m_raySteps=256;
    std::shared_ptr<ElectronicVulkan> m_vk;
};
