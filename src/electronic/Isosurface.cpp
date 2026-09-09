#include "electronic/Volume.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace atomforge::electronic
{
Mesh isosurface(const Grid& grid, double level, const Grid* color)
{
    grid.validate();
    if (!std::isfinite(level)) throw std::invalid_argument("Isovalue must be finite");
    if (color) color->validate();
    Mesh out;
    const glm::ivec3 offsets[8] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const int tetrahedra[6][4] = {{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
    const int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (int z = 0; z < grid.shape[2] - (grid.periodic ? 0 : 1); ++z)
        for (int y = 0; y < grid.shape[1] - (grid.periodic ? 0 : 1); ++y)
            for (int x = 0; x < grid.shape[0] - (grid.periodic ? 0 : 1); ++x)
            {
                glm::dvec3 p[8];
                double values[8];
                for (int i = 0; i < 8; ++i)
                {
                    const auto q = glm::ivec3(x,y,z) + offsets[i];
                    p[i] = grid.position(q.x,q.y,q.z);
                    values[i] = grid.values[grid.index(q.x % grid.shape[0],q.y % grid.shape[1],q.z % grid.shape[2])];
                }
                for (const auto& tetra : tetrahedra)
                {
                    std::vector<glm::dvec3> hits;
                    glm::dvec3 inside(0), outside(0);
                    int ni = 0, no = 0;
                    for (int i : tetra) if (values[i] < level) { inside += p[i]; ++ni; } else { outside += p[i]; ++no; }
                    if (!ni || !no) continue;
                    for (const auto& edge : edges)
                    {
                        const int a = tetra[edge[0]], b = tetra[edge[1]];
                        if ((values[a] < level) != (values[b] < level))
                        {
                            const auto hit = glm::mix(p[a], p[b], (level - values[a]) / (values[b] - values[a]));
                            bool duplicate = false;
                            for (auto v : hits) if (glm::length(hit - v) < 1e-12) duplicate = true;
                            if (!duplicate) hits.push_back(hit);
                        }
                    }
                    if (hits.size() < 3) continue;
                    const auto normal = glm::normalize(outside / static_cast<double>(no) - inside / static_cast<double>(ni));
                    glm::dvec3 center(0);
                    for (auto v : hits) center += v;
                    center /= hits.size();
                    const auto u = glm::normalize(hits[0] - center), v = glm::normalize(glm::cross(normal,u));
                    std::sort(hits.begin(), hits.end(), [&](glm::dvec3 a, glm::dvec3 b)
                    { return std::atan2(glm::dot(a-center,v),glm::dot(a-center,u)) < std::atan2(glm::dot(b-center,v),glm::dot(b-center,u)); });
                    for (std::size_t i = 1; i + 1 < hits.size(); ++i)
                    {
                        if (glm::length(glm::cross(hits[i] - hits[0], hits[i+1] - hits[0])) < 1e-16) continue;
                        for (auto vertex : {hits[0],hits[i],hits[i+1]})
                        {
                            out.vertices.push_back(vertex);
                            out.colors.push_back(color ? color->sample(vertex) : level);
                        }
                    }
                    if (out.vertices.size() > 12000000) throw std::runtime_error("Isosurface exceeds four million triangles; reduce resolution");
                }
            }
    return out;
}

void saveMesh(const Mesh& mesh, const std::string& path)
{
    if (mesh.vertices.size()%3 || mesh.colors.size() != mesh.vertices.size())
        throw std::invalid_argument("Mesh requires triangles and one scalar per vertex");
    std::ofstream out{std::filesystem::u8path(path)};
    if (!out) throw std::runtime_error("Cannot write OBJ mesh");
    if (std::filesystem::u8path(path).extension() == ".ply")
    {
        out << "ply\nformat ascii 1.0\nelement vertex " << mesh.vertices.size()
            << "\nproperty double x\nproperty double y\nproperty double z\nproperty double scalar\nelement face "
            << mesh.vertices.size()/3 << "\nproperty list uchar int vertex_indices\nend_header\n" << std::setprecision(17);
        for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
        {
            const auto p = mesh.vertices[i];
            out << p.x << ' ' << p.y << ' ' << p.z << ' ' << mesh.colors[i] << '\n';
        }
        for (std::size_t i = 0; i < mesh.vertices.size(); i += 3) out << "3 " << i << ' ' << i+1 << ' ' << i+2 << '\n';
        if (!out) throw std::runtime_error("Failed while writing PLY mesh");
        return;
    }
    out << std::setprecision(17) << "# AtomForge isosurface; scalar colors in vertex comments\n";
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        auto p = mesh.vertices[i];
        out << "v " << p.x << ' ' << p.y << ' ' << p.z << " # " << mesh.colors.at(i) << '\n';
    }
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) out << "f " << i+1 << ' ' << i+2 << ' ' << i+3 << '\n';
    if (!out) throw std::runtime_error("Failed while writing OBJ mesh");
}
}
