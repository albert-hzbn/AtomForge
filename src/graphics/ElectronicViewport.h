#pragma once

#include "electronic/Volume.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>

// Offscreen surface renderer. GPU objects live only while a GL context exists.
class ElectronicViewport
{
public:
    ElectronicViewport() = default;
    ~ElectronicViewport();
    ElectronicViewport(const ElectronicViewport&) = delete;
    ElectronicViewport& operator=(const ElectronicViewport&) = delete;
    void setMesh(const atomforge::electronic::Mesh& mesh);
    void setVolume(const atomforge::electronic::Grid& grid);
    GLuint renderVolume(int width, int height, float yaw, float pitch, float zoom,
                        glm::vec2 pan, float opacity, float threshold, float low, float high, int palette);
    GLuint render(int width, int height, float yaw, float pitch, float zoom,
                  glm::vec2 pan, float opacity, float low, float high, int palette,
                  float specular, float shininess);
    static glm::vec3 color(float t, int palette);

private:
    struct Vertex { glm::vec3 position; glm::vec3 normal; float scalar; };
    void release();
    GLuint m_fbo = 0, m_texture = 0, m_depth = 0, m_program = 0, m_vao = 0, m_vbo = 0, m_ebo = 0;
    int m_width = 0, m_height = 0;
    bool m_dirty = false;
    float m_sortedYaw = 100, m_sortedPitch = 100;
    std::vector<Vertex> m_vertices;
    std::vector<unsigned int> m_indices;
    GLuint m_volumeTexture=0, m_volumeProgram=0;
    glm::mat3 m_inverseCell{1};
    glm::ivec3 m_volumeShape{0};
};
