#include "science/SurfaceMesh.h"
#include "electronic/Volume.h"
#include "science/MeshText.h"
#include "util/TaskControl.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kInfinity = std::numeric_limits<double>::infinity();

struct Lattice
{
    std::array<int, 3> n{};
    bool periodic = false;
    glm::dvec3 origin{0.0};
    glm::dmat3 steps{1.0};    // columns: grid step vectors
    glm::dmat3 inverse{1.0};

    std::size_t size() const { return static_cast<std::size_t>(n[0]) * n[1] * n[2]; }
    std::size_t index(int x, int y, int z) const { return (static_cast<std::size_t>(z) * n[1] + y) * n[0] + x; }
    glm::dvec3 position(const glm::ivec3& i) const { return origin + steps * glm::dvec3(i); }
    // Wraps a periodic index or reports it outside the grid.
    bool wrap(glm::ivec3& i) const
    {
        for (int d = 0; d < 3; ++d) {
            if (periodic) i[d] = ((i[d] % n[d]) + n[d]) % n[d];
            else if (i[d] < 0 || i[d] >= n[d]) return false;
        }
        return true;
    }
};

// Integer offsets within `radius` (Angstrom) of a grid point.
std::vector<glm::ivec3> stencil(const Lattice& g, double radius)
{
    std::array<int, 3> reach{};
    for (int d = 0; d < 3; ++d) {
        // Row d of the inverse: how many steps along d a unit Cartesian length spans.
        const glm::dvec3 row(g.inverse[0][d], g.inverse[1][d], g.inverse[2][d]);
        reach[static_cast<std::size_t>(d)] = static_cast<int>(std::ceil(radius * glm::length(row))) + 1;
    }
    std::vector<glm::ivec3> out;
    for (int x = -reach[0]; x <= reach[0]; ++x)
        for (int y = -reach[1]; y <= reach[1]; ++y)
            for (int z = -reach[2]; z <= reach[2]; ++z)
                if (glm::length(g.steps * glm::dvec3(x, y, z)) <= radius) out.push_back({x, y, z});
    return out;
}

// Lowers each grid point's distance to the nearest of the source points
// (Cartesian), for points within `radius` of a source.
void stamp(const Lattice& g, const std::vector<glm::dvec3>& sources, double radius, std::vector<double>& distance)
{
    const double diagonal = 0.5 * (glm::length(g.steps[0]) + glm::length(g.steps[1]) + glm::length(g.steps[2]));
    const auto offsets = stencil(g, radius + diagonal);
    for (std::size_t s = 0; s < sources.size(); ++s) {
        if (s % 1024 == 0) taskCheckpoint();
        const glm::dvec3 u = g.inverse * (sources[s] - g.origin);
        const glm::ivec3 centre(static_cast<int>(std::lround(u.x)), static_cast<int>(std::lround(u.y)), static_cast<int>(std::lround(u.z)));
        for (const auto& o : offsets) {
            const glm::ivec3 unwrapped = centre + o;
            const double d = glm::length(g.position(unwrapped) - sources[s]);
            if (d > radius) continue;
            glm::ivec3 i = unwrapped;
            if (!g.wrap(i)) continue;
            double& slot = distance[g.index(i.x, i.y, i.z)];
            slot = std::min(slot, d);
        }
    }
}

// Connected components of the voxels where inside[k] is true (6-neighbour, periodic wrap).
std::vector<std::size_t> components(const Lattice& g, const std::vector<char>& inside)
{
    std::vector<int> label(g.size(), -1);
    std::vector<std::size_t> sizes;
    std::vector<glm::ivec3> stack;
    const glm::ivec3 steps[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (int z = 0; z < g.n[2]; ++z)
        for (int y = 0; y < g.n[1]; ++y)
            for (int x = 0; x < g.n[0]; ++x) {
                const std::size_t start = g.index(x, y, z);
                if (!inside[start] || label[start] >= 0) continue;
                const int id = static_cast<int>(sizes.size());
                std::size_t count = 0;
                label[start] = id;
                stack.push_back({x, y, z});
                while (!stack.empty()) {
                    const glm::ivec3 v = stack.back();
                    stack.pop_back();
                    ++count;
                    for (const auto& s : steps) {
                        glm::ivec3 w = v + s;
                        if (!g.wrap(w)) continue;
                        const std::size_t k = g.index(w.x, w.y, w.z);
                        if (inside[k] && label[k] < 0) { label[k] = id; stack.push_back(w); }
                    }
                }
                sizes.push_back(count);
            }
    std::sort(sizes.begin(), sizes.end(), std::greater<>());
    return sizes;
}
}

Json surfaceMesh(const Parameters& p)
{
    const auto all = points(p.array("positions"), "positions");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    if (periodic && !(pbc[0] && pbc[1] && pbc[2])) throw std::runtime_error("Use periodic boundaries along all three cell vectors (a slab: add vacuum) or none");
    if (periodic && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    const double probe = positive(p.number("probe_radius_A", 4.0), "probe_radius_A");
    const double spacing = positive(p.number("grid_spacing_A", 0.5), "grid_spacing_A");
    std::vector<double> mask(all.size(), 1.0);
    if (p.has("mask")) {
        mask = finiteArray(p.array("mask"), "mask", 1).values;
        if (mask.size() != all.size()) throw std::runtime_error("mask needs one value per atom");
    }
    std::vector<glm::dvec3> atoms;
    for (std::size_t i = 0; i < all.size(); ++i)
        if (mask[i] != 0) atoms.push_back({all[i][0], all[i][1], all[i][2]});
    if (atoms.empty()) throw std::runtime_error("No atoms define the solid");

    Lattice g;
    g.periodic = periodic;
    if (periodic) {
        const Mat3 cell = cellMatrix(p.array("cell"));
        glm::dmat3 lattice;
        for (int d = 0; d < 3; ++d) lattice[d] = glm::dvec3(cell[static_cast<std::size_t>(d)][0], cell[static_cast<std::size_t>(d)][1], cell[static_cast<std::size_t>(d)][2]);
        for (int d = 0; d < 3; ++d) {
            g.n[static_cast<std::size_t>(d)] = std::max(4, static_cast<int>(std::ceil(glm::length(lattice[d]) / spacing)));
            g.steps[d] = lattice[d] / static_cast<double>(g.n[static_cast<std::size_t>(d)]);
        }
    } else {
        glm::dvec3 low = atoms[0], high = atoms[0];
        for (const auto& x : atoms) { low = glm::min(low, x); high = glm::max(high, x); }
        const double pad = probe + 3 * spacing;
        g.origin = low - glm::dvec3(pad);
        const glm::dvec3 extent = high - low + glm::dvec3(2 * pad);
        g.steps = glm::dmat3(0.0);
        for (int d = 0; d < 3; ++d) {
            g.n[static_cast<std::size_t>(d)] = static_cast<int>(std::ceil(extent[d] / spacing)) + 1;
            g.steps[d][d] = extent[d] / static_cast<double>(g.n[static_cast<std::size_t>(d)] - 1);
        }
    }
    if (g.size() > 60'000'000) throw std::runtime_error("The grid is too fine for this structure; increase grid_spacing_A");
    g.inverse = glm::inverse(g.steps);
    const double voxel = std::abs(glm::determinant(g.steps));

    // 1. Distance to the nearest atom; probe centres may sit where it is >= the probe radius.
    std::vector<double> toAtom(g.size(), kInfinity);
    taskProgress(0.1);
    stamp(g, atoms, probe, toAtom);
    std::vector<char> allowed(g.size());
    for (std::size_t k = 0; k < g.size(); ++k) allowed[k] = toAtom[k] >= probe;
    // 2. Distance to the region of probe centres (from its boundary voxels).
    taskProgress(0.4);
    std::vector<glm::dvec3> boundary;
    const glm::ivec3 neighbours[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (int z = 0; z < g.n[2]; ++z)
        for (int y = 0; y < g.n[1]; ++y)
            for (int x = 0; x < g.n[0]; ++x) {
                if (!allowed[g.index(x, y, z)]) continue;
                for (const auto& s : neighbours) {
                    glm::ivec3 w = glm::ivec3(x, y, z) + s;
                    if (g.wrap(w) && !allowed[g.index(w.x, w.y, w.z)]) { boundary.push_back(g.position({x, y, z})); break; }
                }
            }
    std::vector<double> toProbe(g.size(), kInfinity);
    for (std::size_t k = 0; k < g.size(); ++k)
        if (allowed[k]) toProbe[k] = 0;
    stamp(g, boundary, probe + 2 * spacing, toProbe);
    // Solid where no probe sphere reaches: field = distance to probe centres - probe radius.
    atomforge::electronic::Grid field;
    field.shape = g.n;
    field.periodic = periodic;
    field.origin = g.origin;
    for (int d = 0; d < 3; ++d) field.cell[d] = g.steps[d] * static_cast<double>(periodic ? g.n[static_cast<std::size_t>(d)] : g.n[static_cast<std::size_t>(d)] - 1);
    field.values.resize(g.size());
    std::vector<char> solid(g.size()), empty(g.size());
    std::size_t solidCount = 0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        field.values[k] = std::min(toProbe[k], 2 * probe) - probe;
        solid[k] = field.values[k] > 0;
        empty[k] = !solid[k];
        solidCount += solid[k];
    }
    taskProgress(0.7);
    const auto mesh = atomforge::electronic::isosurface(field, 0.0);
    double area = 0;
    for (std::size_t t = 0; t + 2 < mesh.vertices.size(); t += 3)
        area += 0.5 * glm::length(glm::cross(mesh.vertices[t + 1] - mesh.vertices[t], mesh.vertices[t + 2] - mesh.vertices[t]));
    const auto solidRegions = components(g, solid);
    const auto emptyRegions = components(g, empty);
    // Surface atoms: a probe sphere touches them (their distance to the probe
    // centres is within a grid step of the probe radius).
    std::vector<bool> surfaceAtom(all.size(), false);
    {
        const auto offsets = stencil(g, probe + 3 * spacing);
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (mask[i] == 0) continue;
            const glm::dvec3 x(all[i][0], all[i][1], all[i][2]);
            const glm::dvec3 u = g.inverse * (x - g.origin);
            const glm::ivec3 centre(static_cast<int>(std::lround(u.x)), static_cast<int>(std::lround(u.y)), static_cast<int>(std::lround(u.z)));
            double nearest = kInfinity;
            for (const auto& o : offsets) {
                glm::ivec3 w = centre + o;
                const glm::dvec3 at = g.position(w);
                if (!g.wrap(w) || !allowed[g.index(w.x, w.y, w.z)]) continue;
                nearest = std::min(nearest, glm::length(at - x));
            }
            surfaceAtom[i] = nearest <= probe + 1.5 * spacing;
        }
    }
    std::size_t surfaceCount = 0;
    for (bool s : surfaceAtom) surfaceCount += s;
    const double total = voxel * static_cast<double>(g.size());
    Json result = Json::object();
    result["surface_area_A2"] = area;
    result["solid_volume_A3"] = voxel * static_cast<double>(solidCount);
    if (periodic) {
        result["cell_volume_A3"] = total;
        result["solid_fraction"] = static_cast<double>(solidCount) / static_cast<double>(g.size());
        result["specific_surface_area_per_A"] = area / total;
    }
    Json solidJson = Json::array(), emptyJson = Json::array();
    for (std::size_t c : solidRegions) solidJson.push(voxel * static_cast<double>(c));
    for (std::size_t c : emptyRegions) emptyJson.push(voxel * static_cast<double>(c));
    result["solid_region_volumes_A3"] = solidJson;
    result["empty_region_volumes_A3"] = emptyJson;
    result["solid_regions"] = static_cast<double>(solidRegions.size());
    result["empty_regions"] = static_cast<double>(emptyRegions.size());
    result["surface_atom"] = toJson(surfaceAtom);
    result["surface_atoms"] = static_cast<double>(surfaceCount);
    result["triangles"] = static_cast<double>(mesh.vertices.size() / 3);
    result["grid"] = Json::array({static_cast<double>(g.n[0]), static_cast<double>(g.n[1]), static_cast<double>(g.n[2])});
    result["probe_radius_A"] = probe;
    if (p.boolean("export_mesh", true) && mesh.vertices.size() <= 3'000'000) {
        std::ostringstream comment;
        comment << "AtomForge surface mesh: probe radius " << probe << " A, area " << area << " A^2";
        ObjWriter obj(comment.str());
        obj.add("", mesh.vertices);
        Json files = Json::object();
        files["surface.obj"] = obj.text();
        result["files"] = files;
    }
    return result;
}
}
