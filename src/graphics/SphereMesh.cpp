#include "SphereMesh.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cstring>

static void createSphereIndexed(std::vector<float>& vertices,
                                std::vector<unsigned int>& indices,
                                int stacks,
                                int slices)
{
    constexpr float kPi = 3.14159265358979323846f;
    
    std::unordered_map<uint64_t, unsigned int> vertexMap;
    unsigned int vertexIndex = 0;

    auto getOrCreateVertex = [&](float x, float y, float z) -> unsigned int {
        // Normalize position
        float len = std::sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;

        // Create hash for this vertex using memcpy (safe from aliasing issues)
        uint32_t ix, iy, iz;
        std::memcpy(&ix, &x, sizeof(uint32_t));
        std::memcpy(&iy, &y, sizeof(uint32_t));
        std::memcpy(&iz, &z, sizeof(uint32_t));
        uint64_t hash = ((uint64_t)ix << 32) | iy;
        hash ^= (uint64_t)iz;

        if (vertexMap.find(hash) != vertexMap.end())
        {
            return vertexMap[hash];
        }

        // Add new vertex
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);
        vertexMap[hash] = vertexIndex;
        return vertexIndex++;
    };

    for (int i = 0; i < stacks; ++i)
    {
        float lat0 = kPi * (-0.5f + i / (float)stacks);
        float lat1 = kPi * (-0.5f + (i + 1) / (float)stacks);

        for (int j = 0; j < slices; ++j)
        {
            float lng0 = 2 * kPi * j / (float)slices;
            float lng1 = 2 * kPi * (j + 1) / (float)slices;

            float x0 = std::cos(lng0);
            float y0 = std::sin(lng0);

            float x1 = std::cos(lng1);
            float y1 = std::sin(lng1);

            float z0 = std::sin(lat0);
            float zr0 = std::cos(lat0);

            float z1 = std::sin(lat1);
            float zr1 = std::cos(lat1);

            // Triangle 1
            unsigned int v0 = getOrCreateVertex(x0 * zr0, y0 * zr0, z0);
            unsigned int v1 = getOrCreateVertex(x0 * zr1, y0 * zr1, z1);
            unsigned int v2 = getOrCreateVertex(x1 * zr1, y1 * zr1, z1);

            indices.push_back(v0);
            indices.push_back(v1);
            indices.push_back(v2);

            // Triangle 2
            unsigned int v3 = getOrCreateVertex(x0 * zr0, y0 * zr0, z0);
            unsigned int v4 = getOrCreateVertex(x1 * zr1, y1 * zr1, z1);
            unsigned int v5 = getOrCreateVertex(x1 * zr0, y1 * zr0, z0);

            indices.push_back(v3);
            indices.push_back(v4);
            indices.push_back(v5);
        }
    }
}

SphereMesh::SphereMesh(int st, int sl)
    : ebo(0), indexCount(0)
{
    stacks = st;
    slices = sl;

    std::vector<float> vertices;
    std::vector<unsigned int> indices;

    createSphereIndexed(vertices, indices, stacks, slices);

    vertexCount = vertices.size() / 3;
    indexCount = indices.size();

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 vertices.size() * sizeof(float),
                 vertices.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 indices.size() * sizeof(unsigned int),
                 indices.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          3 * sizeof(float),
                          (void*)0);

    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

SphereMesh::~SphereMesh()
{
    glDeleteBuffers(1, &ebo);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

