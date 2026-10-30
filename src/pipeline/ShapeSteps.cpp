// The Shape Modeller step: carves a crystal into a shape written as a CAD-like
// script (algorithms/ShapeModel), or keeps, deletes or selects the atoms of a
// structure inside such a shape. Regions give the atoms inside them another
// element (core-shell particles, inclusions, Janus particles).

#include "algorithms/ShapeModel.h"
#include "pipeline/ModifierSupport.h"
#include "pipeline/PipelineSteps.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace atomforge::pipeline
{
namespace
{
using namespace support;

glm::dvec3 toGlm(const Vec3& v) { return {v[0], v[1], v[2]}; }
Vec3 fromGlm(const glm::dvec3& v) { return {v.x, v.y, v.z}; }

void setElement(AtomSite& a, const std::string& symbol)
{
    if (a.symbol == symbol) return;
    a.symbol = symbol;
    a.atomicNumber = science::atomicNumber(symbol);
    getDefaultElementColor(a.atomicNumber, a.r, a.g, a.b);
}

// Region elements checked once ("none" removes atoms).
std::map<std::string, std::string> regionSymbols(const shape::Model& model)
{
    std::map<std::string, std::string> out;
    for (const std::string& name : model.regionElements())
        out[name] = name == "none" || name == "vacancy" ? std::string("none") : standardElement(name);
    return out;
}

void buildShape(const Json& p, PipelineData& data)
{
    const char* id = "build-shape";
    const std::string mode = text(p, id, "mode");
    if (mode != "carve" && mode != "keep" && mode != "delete" && mode != "select")
        throw std::runtime_error("build-shape: mode is carve, keep, delete or select");
    Structure& input = data.structure;
    glm::dmat3 lattice(1.0);
    const bool hasCell = input.hasUnitCell;
    if (hasCell) {
        const Mat3 cell = cellOf(input);
        lattice = glm::dmat3(toGlm(cell[0]), toGlm(cell[1]), toGlm(cell[2]));  // columns a, b, c
    }
    shape::Model model;
    try {
        model = shape::Model::parse(text(p, id, "script"), hasCell ? &lattice : nullptr);
    } catch (const shape::ShapeError& error) {
        throw std::runtime_error(std::string("Shape script, ") + error.what());
    }
    const shape::Bounds bounds = model.bounds();
    const auto symbols = regionSymbols(model);
    // Where the shape's origin sits.
    std::string anchor = text(p, id, "anchor");
    if (anchor == "auto") anchor = mode == "carve" ? "origin" : hasCell ? "cell-centre" : "centroid";
    Vec3 origin{0, 0, 0};
    if (anchor == "cell-centre") {
        requireCell(data, "The cell centre");
        const Mat3 cell = cellOf(input);
        origin = science::add(input.cellOffset, science::scale(science::add(science::add(cell[0], cell[1]), cell[2]), 0.5));
    } else if (anchor == "centroid") {
        for (const auto& a : input.atoms) origin = science::add(origin, positionOf(a));
        if (!input.atoms.empty()) origin = science::scale(origin, 1.0 / input.atoms.size());
    } else if (anchor != "origin")
        throw std::runtime_error("build-shape: anchor is auto, origin, cell-centre or centroid");
    origin = science::add(origin, vector(p, id, "position"));
    const glm::dvec3 shapeOrigin = toGlm(origin);
    std::map<std::string, std::size_t> regionCounts;
    std::size_t removedByRegion = 0;
    std::ostringstream note;
    note.precision(5);

    if (mode == "carve") {
        // The crystal (this structure's cell and atoms) repeated over the shape, keeping the atoms inside.
        requireCell(data, "Carving a shape (open a crystal first, e.g. Build > Bulk Crystal)");
        if (input.atoms.empty()) throw std::runtime_error("build-shape: the crystal has no atoms");
        if (!bounds.finite) throw std::runtime_error("build-shape: the shape is unbounded; intersect it with a box or other finite shape");
        const glm::dmat3 inverse = glm::inverse(lattice);
        const glm::dvec3 cellOrigin = toGlm(input.cellOffset);
        glm::dvec3 low(1e300), high(-1e300);
        for (int c = 0; c < 8; ++c) {
            const glm::dvec3 corner((c & 1) ? bounds.high.x : bounds.low.x, (c & 2) ? bounds.high.y : bounds.low.y, (c & 4) ? bounds.high.z : bounds.low.z);
            const glm::dvec3 f = inverse * (corner + shapeOrigin - cellOrigin);
            low = glm::min(low, f);
            high = glm::max(high, f);
        }
        long long range[3][2];
        double images = 1;
        for (int k = 0; k < 3; ++k) {
            range[k][0] = static_cast<long long>(std::floor(low[k])) - 1;
            range[k][1] = static_cast<long long>(std::ceil(high[k])) + 1;
            images *= double(range[k][1] - range[k][0] + 1);
        }
        const double maxAtoms = number(p, id, "max_atoms");
        const double cellVolume = std::abs(glm::determinant(lattice));
        const double expected = model.volume(std::max(0.5, std::cbrt((bounds.high.x - bounds.low.x) * (bounds.high.y - bounds.low.y) *
                                                                     (bounds.high.z - bounds.low.z)) / 80)) / cellVolume * input.atoms.size();
        if (expected > maxAtoms)
            throw std::runtime_error("build-shape: about " + std::to_string(static_cast<long long>(expected)) + " atoms, more than max_atoms (" +
                                     std::to_string(static_cast<long long>(maxAtoms)) + ")");
        if (images * input.atoms.size() > 4e9) throw std::runtime_error("build-shape: the shape is far too large for this cell");
        Structure out;
        out.atoms.reserve(static_cast<std::size_t>(std::max(0.0, expected * 1.1)));
        for (long long i = range[0][0]; i <= range[0][1]; ++i)
            for (long long j = range[1][0]; j <= range[1][1]; ++j)
                for (long long k = range[2][0]; k <= range[2][1]; ++k) {
                    const glm::dvec3 shift = lattice * glm::dvec3(double(i), double(j), double(k));
                    for (const AtomSite& a : input.atoms) {
                        const glm::dvec3 position = glm::dvec3(a.x, a.y, a.z) + shift;
                        const glm::dvec3 q = position - shapeOrigin;
                        if (!model.contains(q)) continue;
                        AtomSite atom = a;
                        place(atom, fromGlm(position));
                        if (const std::string* region = model.material(q)) {
                            const std::string& symbol = symbols.at(*region);
                            if (symbol == "none") { ++removedByRegion; continue; }
                            setElement(atom, symbol);
                            ++regionCounts[symbol];
                        }
                        out.atoms.push_back(atom);
                    }
                }
        if (out.atoms.empty()) throw std::runtime_error("build-shape: no atom of the crystal lies inside the shape");
        // An open box around the shape with vacuum on every side.
        const double vacuum = number(p, id, "vacuum");
        glm::dvec3 lo(1e300), hi(-1e300);
        for (const auto& a : out.atoms) {
            lo = glm::min(lo, glm::dvec3(a.x, a.y, a.z));
            hi = glm::max(hi, glm::dvec3(a.x, a.y, a.z));
        }
        const glm::dvec3 move = glm::dvec3(vacuum / 2) - lo;
        for (auto& a : out.atoms) { a.x += move.x; a.y += move.y; a.z += move.z; }
        const glm::dvec3 size = hi - lo + glm::dvec3(vacuum);
        out.hasUnitCell = true;
        out.cellVectors = {Vec3{size.x, 0, 0}, Vec3{0, size.y, 0}, Vec3{0, 0, size.z}};
        out.cellOffset = {0, 0, 0};
        out.pbc = {false, false, false};
        data.structure = std::move(out);
        data.clearSelection();
        note << data.structure.atoms.size() << " atoms carved from the crystal";
        const std::string meshFile = text(p, id, "mesh_file");
        if (!meshFile.empty()) {
            electronic::Mesh mesh = model.mesh(number(p, id, "mesh_spacing"));
            for (auto& v : mesh.vertices) v += shapeOrigin + move;
            shape::saveObj(mesh, meshFile);
            data.notes.push_back("Shape surface written to " + meshFile + " (" + std::to_string(mesh.vertices.size() / 3) + " triangles)");
        }
    } else {
        // The structure's own atoms against the shape.
        std::vector<char> inside(input.atoms.size(), 0);
        std::size_t count = 0;
        for (std::size_t i = 0; i < input.atoms.size(); ++i)
            if (model.contains(toGlm(positionOf(input.atoms[i])) - shapeOrigin)) { inside[i] = 1; ++count; }
        if (mode == "select") {
            select(data, inside, text(p, id, "combine"));
            note << count << " atoms inside the shape selected";
        } else {
            std::vector<std::size_t> keep;
            for (std::size_t i = 0; i < input.atoms.size(); ++i) {
                if ((mode == "keep") != static_cast<bool>(inside[i])) continue;
                if (mode == "keep") {
                    if (const std::string* region = model.material(toGlm(positionOf(input.atoms[i])) - shapeOrigin)) {
                        const std::string& symbol = symbols.at(*region);
                        if (symbol == "none") { ++removedByRegion; continue; }
                        setElement(input.atoms[i], symbol);
                        ++regionCounts[symbol];
                    }
                }
                keep.push_back(i);
            }
            rebuild(data, keep);
            note << (mode == "keep" ? std::to_string(keep.size()) + " atoms inside the shape kept" : std::to_string(count) + " atoms inside the shape deleted");
        }
        const std::string meshFile = text(p, id, "mesh_file");
        if (!meshFile.empty()) {
            electronic::Mesh mesh = model.mesh(number(p, id, "mesh_spacing"));
            for (auto& v : mesh.vertices) v += shapeOrigin;
            shape::saveObj(mesh, meshFile);
            data.notes.push_back("Shape surface written to " + meshFile + " (" + std::to_string(mesh.vertices.size() / 3) + " triangles)");
        }
    }
    note << " (" << model.primitives() << " primitive" << (model.primitives() == 1 ? "" : "s");
    if (bounds.finite)
        note << ", " << bounds.high.x - bounds.low.x << " x " << bounds.high.y - bounds.low.y << " x " << bounds.high.z - bounds.low.z << " Angstrom";
    note << ")";
    for (const auto& [symbol, n] : regionCounts) note << "; " << n << " " << symbol << " in regions";
    if (removedByRegion) note << "; " << removedByRegion << " removed by \"none\" regions";
    data.notes.insert(data.notes.begin(), note.str());
}
}  // namespace

std::vector<ModifierType> shapeStepTypes()
{
    return {
        {"build-shape", "Shape Modeller", "Build",
         "Carves a crystal into any shape modelled as in CAD: primitives (sphere, box, cylinder, cone, torus, capsule, ellipsoid, prism, "
         "extruded polygon, half-space, slab and crystallographic {hkl} facets), booleans (union, intersection, difference, smooth "
         "union), transforms (translate, rotate, scale, mirror), shells, offsets and linear, polar and grid arrays, written as an "
         "OpenSCAD-like script with variables. Regions give the atoms inside them another element (core-shell particles, inclusions, "
         "Janus particles; \"none\" removes them). Carve fills the shape with the crystal (this structure's cell, repeated) in an open box; "
         "keep, delete and select use the structure's own atoms (voids, notches, pores). The shape's surface can be written as OBJ.",
         {{"script", "Shape script", "string", "sphere(r = 15);", ""},
          {"mode", "Mode", "choice", "carve", "carve|keep|delete|select"},
          {"anchor", "Shape origin at", "choice", "auto", "auto|origin|cell-centre|centroid"},
          {"position", "Shape offset (Angstrom)", "vector", "0 0 0", ""},
          {"vacuum", "Vacuum around a carved shape (Angstrom)", "float", "10", ""},
          {"max_atoms", "At most this many atoms", "float", "20000000", ""},
          {"combine", "Combine with the selection (select)", "choice", "replace", kModes},
          {"mesh_file", "Write the shape surface to (OBJ; empty = no)", "string", "", ""},
          {"mesh_spacing", "Surface mesh spacing (Angstrom)", "float", "0.5", ""}},
         buildShape},
    };
}
}  // namespace atomforge::pipeline
