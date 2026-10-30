#pragma once
// Shapes for carving crystals, modelled as in CAD tools: primitives (sphere,
// box, cylinder, cone, torus, capsule, ellipsoid, prisms, extruded polygons,
// half-spaces, slabs and crystallographic {hkl} facet polyhedra) combined by
// booleans (union, intersection, difference and their smooth versions),
// transforms (translate, rotate, scale, mirror), modifiers (shell, offset) and
// patterns (linear, polar and grid arrays). Regions assign an element to the
// shape's atoms inside them (core-shell particles, inclusions, Janus
// particles); they add no volume of their own.
//
// A shape is written as a small OpenSCAD-like script:
//
//     R = 20;                                  // variables and arithmetic
//     difference() {
//         sphere(r = R);
//         translate([0, 0, R]) cylinder(r = 5, h = 2 * R);
//     }
//     region("Au") sphere(r = R / 2);
//
// and evaluated as a signed distance function (negative inside, Angstrom).

#include "electronic/Volume.h"

#include <glm/glm.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace atomforge::shape
{
// A script error, with the 1-based line it was found on.
struct ShapeError : std::runtime_error
{
    int line;
    ShapeError(int at, const std::string& message)
        : std::runtime_error("line " + std::to_string(at) + ": " + message), line(at) {}
};

struct Bounds
{
    glm::dvec3 low{0.0}, high{0.0};
    bool finite = false;
    bool empty = false;  // no volume (regions)
};

struct Node;

class Model
{
public:
    // `cell` (columns a, b, c) gives Miller indices of facets() their plane
    // normals; without it they are Cartesian directions.
    static Model parse(const std::string& script, const glm::dmat3* cell = nullptr);

    // Signed distance (Angstrom): negative inside the shape.
    double distance(const glm::dvec3& p) const;
    bool contains(const glm::dvec3& p) const { return distance(p) <= 0.0; }
    // The element of the last region containing p (nested regions win), or
    // nullptr outside every region. "none" marks atoms to remove.
    const std::string* material(const glm::dvec3& p) const;
    Bounds bounds() const;
    // Region elements in script order, without repeats.
    std::vector<std::string> regionElements() const;
    int primitives() const { return m_primitives; }
    // Volume (Angstrom^3) counted on a grid of the given spacing.
    double volume(double spacing) const;
    // The surface as triangles (marching tetrahedra on a grid of the given spacing).
    electronic::Mesh mesh(double spacing) const;
    bool empty() const { return !m_root; }

private:
    std::shared_ptr<Node> m_root;
    int m_primitives = 0;
};

// Writes a mesh as Wavefront OBJ.
void saveObj(const electronic::Mesh& mesh, const std::string& path);

// Ready-made shapes: (title, script).
const std::vector<std::pair<std::string, std::string>>& templates();
// Snippets for each command, for inserting into a script: (command, snippet, description).
struct Command
{
    const char* group;
    const char* name;
    const char* snippet;
    const char* help;
};
const std::vector<Command>& commands();
}
