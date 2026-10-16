#pragma once
// Triangle meshes of the scientific tools (surface meshes, LAMMPS regions) as
// Wavefront OBJ text, for a result's "files" (shown in the view, or saved).

#include <glm/glm.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace atomforge::science
{
class ObjWriter
{
public:
    explicit ObjWriter(const std::string& comment)
    {
        m_text.precision(8);
        m_text << "# " << comment << '\n';
    }

    // Adds an object of triangles (three vertices each; shared vertices are not merged).
    void add(const std::string& name, const std::vector<glm::dvec3>& triangles)
    {
        if (!name.empty()) m_text << "o " << name << '\n';
        for (const auto& v : triangles) m_text << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';
        for (std::size_t t = 0; t + 2 < triangles.size(); t += 3)
            m_text << "f " << m_vertices + t + 1 << ' ' << m_vertices + t + 2 << ' ' << m_vertices + t + 3 << '\n';
        m_vertices += triangles.size();
    }

    std::string text() const { return m_text.str(); }

private:
    std::ostringstream m_text;
    std::size_t m_vertices = 0;
};
}
