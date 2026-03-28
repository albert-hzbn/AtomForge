#include "LowPolyMesh.h"
#include <cmath>
#include <map>
#include <vector>

namespace
{
// Normalize a 3-component vertex to unit length.
void normalize(float& x, float& y, float& z)
{
    float len = std::sqrt(x * x + y * y + z * z);
    if (len > 1e-8f) { x /= len; y /= len; z /= len; }
}

// Return the index of the midpoint between two vertices, adding it if new.
typedef std::map<std::pair<unsigned int, unsigned int>, unsigned int> EdgeMap;

unsigned int midpoint(std::vector<float>& verts, EdgeMap& cache,
                      unsigned int a, unsigned int b)
{
    if (a > b) std::swap(a, b);
    std::pair<unsigned int, unsigned int> key(a, b);

    EdgeMap::iterator it = cache.find(key);
    if (it != cache.end())
        return it->second;

    float mx = (verts[3 * a] + verts[3 * b]) * 0.5f;
    float my = (verts[3 * a + 1] + verts[3 * b + 1]) * 0.5f;
    float mz = (verts[3 * a + 2] + verts[3 * b + 2]) * 0.5f;
    normalize(mx, my, mz);

    unsigned int idx = static_cast<unsigned int>(verts.size() / 3);
    verts.push_back(mx);
    verts.push_back(my);
    verts.push_back(mz);

    cache[key] = idx;
    return idx;
}
} // namespace

LowPolyMesh::LowPolyMesh()
{
    // Golden ratio
    const float PHI = 1.6180339887f;

    // Start with a base icosahedron (12 vertices, 20 faces).
    std::vector<float> verts = {
        -1,  PHI, 0,
         1,  PHI, 0,
        -1, -PHI, 0,
         1, -PHI, 0,

         0, -1,  PHI,
         0,  1,  PHI,
         0, -1, -PHI,
         0,  1, -PHI,

         PHI, 0, -1,
         PHI, 0,  1,
        -PHI, 0, -1,
        -PHI, 0,  1
    };

    // Normalize to unit sphere.
    for (size_t i = 0; i < verts.size(); i += 3)
        normalize(verts[i], verts[i + 1], verts[i + 2]);

    std::vector<unsigned int> tris = {
        0, 11, 5,   0, 5, 1,    0, 1, 7,    0, 7, 10,   0, 10, 11,
        1, 5, 9,    5, 11, 4,   11, 10, 2,  10, 7, 6,    7, 1, 8,
        3, 9, 4,    3, 4, 2,    3, 2, 6,    3, 6, 8,     3, 8, 9,
        4, 9, 5,    2, 4, 11,   6, 2, 10,   8, 6, 7,     9, 8, 1
    };

    // Subdivide twice: 20 → 80 → 320 faces.
    for (int sub = 0; sub < 2; ++sub)
    {
        EdgeMap cache;
        std::vector<unsigned int> newTris;
        newTris.reserve(tris.size() * 4);

        for (size_t i = 0; i < tris.size(); i += 3)
        {
            unsigned int a = tris[i];
            unsigned int b = tris[i + 1];
            unsigned int c = tris[i + 2];

            unsigned int ab = midpoint(verts, cache, a, b);
            unsigned int bc = midpoint(verts, cache, b, c);
            unsigned int ca = midpoint(verts, cache, c, a);

            newTris.push_back(a);  newTris.push_back(ab); newTris.push_back(ca);
            newTris.push_back(b);  newTris.push_back(bc); newTris.push_back(ab);
            newTris.push_back(c);  newTris.push_back(ca); newTris.push_back(bc);
            newTris.push_back(ab); newTris.push_back(bc); newTris.push_back(ca);
        }

        tris.swap(newTris);
    }

    vertexCount = static_cast<int>(verts.size() / 3);
    indexCount  = static_cast<int>(tris.size());

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(verts.size() * sizeof(float)),
                 verts.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(tris.size() * sizeof(unsigned int)),
                 tris.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

LowPolyMesh::~LowPolyMesh()
{
    glDeleteBuffers(1, &ebo);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}
