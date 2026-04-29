#include "SphereMesh.h"
#include <cmath>
#include <vector>

static void createSphereIndexed(std::vector<float>& vertices,
                                std::vector<unsigned int>& indices,
                                int stacks,
                                int slices)
{
    constexpr float kPi = 3.14159265358979323846f;

    // Generate (stacks+1) x (slices+1) vertices analytically.
    // Vertex at grid position (i, j) has index i*(slices+1)+j.
    // This guarantees shared edges between adjacent quads reference exactly
    // the same index — no hashing, no floating-point comparison, no collisions.
    for (int i = 0; i <= stacks; ++i)
    {
        const float lat = kPi * (-0.5f + (float)i / (float)stacks);
        const float cosLat = std::cos(lat);
        const float sinLat = std::sin(lat);

        for (int j = 0; j <= slices; ++j)
        {
            const float lng = 2.0f * kPi * (float)j / (float)slices;
            vertices.push_back(std::cos(lng) * cosLat); // x
            vertices.push_back(std::sin(lng) * cosLat); // y
            vertices.push_back(sinLat);                  // z
        }
    }

    // Build two triangles per quad, winding CCW when viewed from outside.
    for (int i = 0; i < stacks; ++i)
    {
        for (int j = 0; j < slices; ++j)
        {
            const unsigned int row0 = (unsigned int)(i       * (slices + 1));
            const unsigned int row1 = (unsigned int)((i + 1) * (slices + 1));

            const unsigned int v00 = row0 + (unsigned int)j;
            const unsigned int v01 = row0 + (unsigned int)(j + 1);
            const unsigned int v10 = row1 + (unsigned int)j;
            const unsigned int v11 = row1 + (unsigned int)(j + 1);

            // Triangle 1
            indices.push_back(v00);
            indices.push_back(v10);
            indices.push_back(v01);

            // Triangle 2
            indices.push_back(v01);
            indices.push_back(v10);
            indices.push_back(v11);
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

