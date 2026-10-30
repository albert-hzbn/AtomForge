// Geometric and chemical edits: mirror, inversion, rotation and alignment,
// sorting, hydrogen saturation, setting and measuring bond geometry,
// duplicating the selection, cell reduction and scaling, swapping elements,
// making molecules whole and selecting fragments.

#include "pipeline/ModifierSupport.h"
#include "pipeline/PipelineSteps.h"
#include "pipeline/StructureEdits.h"
#include "science/ScienceData.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <filesystem>
#include <random>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <queue>
#include <sstream>

namespace atomforge::pipeline
{
namespace
{
using namespace support;
using science::Mat3;
using science::Vec3;
constexpr double kPi = 3.14159265358979323846;
constexpr double kDegree = kPi / 180.0;

Vec3 unit(const Vec3& v)
{
    const double n = science::norm(v);
    if (!(n > 1e-12)) throw std::runtime_error("A direction must not be zero");
    return science::scale(v, 1.0 / n);
}

// Wraps every atom into the cell along its periodic vectors.
void wrapAll(Structure& s)
{
    if (!s.hasUnitCell) return;
    const Mat3 cell = cellOf(s), inverse = science::inverse(cell);
    for (auto& a : s.atoms) {
        Vec3 f = science::rowTimes(science::sub(positionOf(a), s.cellOffset), inverse);
        for (std::size_t k = 0; k < 3; ++k)
            if (s.pbc[k]) f[k] -= std::floor(f[k]);
        place(a, science::add(science::rowTimes(f, cell), s.cellOffset));
    }
}

Vec3 centroid(const Structure& s, const std::vector<std::size_t>& atoms)
{
    Vec3 c{0, 0, 0};
    for (std::size_t i : atoms) c = science::add(c, positionOf(s.atoms[i]));
    return atoms.empty() ? c : science::scale(c, 1.0 / atoms.size());
}

Vec3 cellCentre(const Structure& s)
{
    const Mat3 cell = cellOf(s);
    return science::add(s.cellOffset, science::scale(science::add(science::add(cell[0], cell[1]), cell[2]), 0.5));
}

// The pivot of a mirror, inversion or rotation.
Vec3 pivotOf(const Json& p, const char* id, const PipelineData& data, const std::vector<std::size_t>& atoms)
{
    const std::string centre = text(p, id, "center");
    if (centre == "origin") return {0, 0, 0};
    if (centre == "point") return vector(p, id, "point");
    if (centre == "cell") {
        requireCell(data, "The cell centre");
        return cellCentre(data.structure);
    }
    return centroid(data.structure, atoms);
}

std::vector<std::size_t> parseIndices(const std::string& text, std::size_t atoms, const char* id)
{
    std::istringstream in(text);
    std::vector<std::size_t> out;
    for (std::string w; in >> w;) {
        char* end = nullptr;
        const long value = std::strtol(w.c_str(), &end, 10);
        if (*end || value < 0 || static_cast<std::size_t>(value) >= atoms)
            throw std::runtime_error(std::string(id) + ": '" + w + "' is not an atom index (0 to " + std::to_string(atoms ? atoms - 1 : 0) + ")");
        out.push_back(static_cast<std::size_t>(value));
    }
    return out;
}

// ---------------------------------------------------------------- symmetry edits
void mirror(const Json& p, PipelineData& data)
{
    const char* id = "mirror";
    const bool selectedOnly = text(p, id, "target") == "selected";
    const auto atoms = targets(data, selectedOnly, id);
    const Vec3 n = unit(vector(p, id, "normal"));
    const Vec3 point = pivotOf(p, id, data, atoms);
    const auto reflect = [&](const Vec3& x) { return science::sub(x, science::scale(n, 2 * science::dot(science::sub(x, point), n))); };
    for (std::size_t i : atoms) place(data.structure.atoms[i], reflect(positionOf(data.structure.atoms[i])));
    Structure& s = data.structure;
    if (!selectedOnly && s.hasUnitCell) {
        // The whole crystal is reflected with its lattice; the reflected vectors
        // a, b and -c span the same lattice and keep the cell right-handed.
        for (auto& row : s.cellVectors) row = science::sub(row, science::scale(n, 2 * science::dot(row, n)));
        s.cellVectors[2] = science::scale(s.cellVectors[2], -1);
        s.cellOffset = reflect(s.cellOffset);
        s.cellOffset = science::sub(s.cellOffset, s.cellVectors[2]);
        wrapAll(s);
    } else if (flag(p, id, "wrap"))
        wrapAll(s);
}

void invert(const Json& p, PipelineData& data)
{
    const char* id = "invert";
    const auto atoms = targets(data, text(p, id, "target") == "selected", id);
    const Vec3 point = pivotOf(p, id, data, atoms);
    for (std::size_t i : atoms)
        place(data.structure.atoms[i], science::sub(science::scale(point, 2), positionOf(data.structure.atoms[i])));
    if (flag(p, id, "wrap")) wrapAll(data.structure);
}

// A rotation from Euler angles (degrees): intrinsic z-y-z or z-x-z (alpha, beta, gamma),
// or extrinsic x-y-z (about the fixed x, then y, then z).
Mat3 eulerRotation(const Vec3& angles, const std::string& convention)
{
    const Vec3 x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1};
    const double a = angles[0] * kDegree, b = angles[1] * kDegree, c = angles[2] * kDegree;
    if (convention == "zyz") return science::multiply(rotationAbout(z, a), science::multiply(rotationAbout(y, b), rotationAbout(z, c)));
    if (convention == "zxz") return science::multiply(rotationAbout(z, a), science::multiply(rotationAbout(x, b), rotationAbout(z, c)));
    if (convention == "xyz") return science::multiply(rotationAbout(z, c), science::multiply(rotationAbout(y, b), rotationAbout(x, a)));
    throw std::runtime_error("Euler convention must be zyz, zxz or xyz");
}

void rotateStep(const Json& p, PipelineData& data)
{
    const char* id = "rotate";
    const bool selectedOnly = text(p, id, "target") == "selected";
    const auto atoms = targets(data, selectedOnly, id);
    const std::string mode = text(p, id, "mode");
    if (mode != "axis" && mode != "euler") throw std::runtime_error("rotate: mode is axis or euler");
    const Mat3 r = mode == "euler" ? eulerRotation(vector(p, id, "euler"), text(p, id, "convention"))
                                   : rotationAbout(vector(p, id, "axis"), number(p, id, "angle") * kDegree);
    Structure& crystal = data.structure;
    const bool withCell = !selectedOnly && flag(p, id, "rotate_cell") && crystal.hasUnitCell;
    // A crystal turned with its cell turns about the cell origin: fractional coordinates stay.
    const Vec3 point = withCell ? Vec3{crystal.cellOffset[0], crystal.cellOffset[1], crystal.cellOffset[2]} : pivotOf(p, id, data, atoms);
    for (std::size_t i : atoms)
        place(data.structure.atoms[i], science::add(point, rotate(r, science::sub(positionOf(data.structure.atoms[i]), point))));
    if (withCell)
        for (auto& row : crystal.cellVectors) row = rotate(r, row);
    if (mode == "euler") {
        // The equivalent single turn, for the record.
        const double trace = r[0][0] + r[1][1] + r[2][2];
        const double angle = std::acos(std::clamp((trace - 1) / 2, -1.0, 1.0)) / kDegree;
        std::ostringstream note;
        note.precision(6);
        Vec3 axis{r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]};
        if (science::norm(axis) > 1e-12) axis = science::scale(axis, 1 / science::norm(axis));
        note << "Equivalent to a turn of " << angle << " degrees about (" << axis[0] << ", " << axis[1] << ", " << axis[2] << ")";
        data.notes.push_back(note.str());
    }
}

void alignDirection(const Json& p, PipelineData& data)
{
    const char* id = "align-direction";
    Structure& s = data.structure;
    const bool lattice = flag(p, id, "lattice");
    const auto cartesian = [&](const Vec3& v) {
        if (!lattice) return v;
        requireCell(data, "A lattice direction");
        return science::rowTimes(v, cellOf(s));
    };
    const auto axisOf = [&](const std::string& name) -> Vec3 {
        return name == "x" ? Vec3{1, 0, 0} : name == "y" ? Vec3{0, 1, 0} : Vec3{0, 0, 1};
    };
    const Vec3 d = unit(cartesian(vector(p, id, "direction")));
    const Vec3 e = axisOf(text(p, id, "axis"));
    Mat3 r = rotationBetween(d, e);
    const std::string secondAxis = text(p, id, "second_axis");
    if (secondAxis != "none") {
        const Vec3 f = axisOf(secondAxis);
        if (science::norm(science::cross(e, f)) < 0.5) throw std::runtime_error("The second axis must differ from the first");
        Vec3 t = rotate(r, unit(cartesian(vector(p, id, "second_direction"))));
        t = science::sub(t, science::scale(e, science::dot(t, e)));
        if (science::norm(t) < 1e-6) throw std::runtime_error("The second direction is parallel to the first");
        // Turn about the first axis to bring the second direction onto the plane of e and f, towards f.
        const double angle = std::atan2(science::dot(science::cross(unit(t), f), e), science::dot(unit(t), f));
        r = science::multiply(rotationAbout(e, angle), r);
    }
    for (auto& a : s.atoms) place(a, rotate(r, positionOf(a)));
    if (s.hasUnitCell) {
        for (auto& row : s.cellVectors) row = rotate(r, row);
        s.cellOffset = rotate(r, s.cellOffset);
    }
}

// ---------------------------------------------------------------- ordering
void sortAtoms(const Json& p, PipelineData& data)
{
    const char* id = "sort-atoms";
    const Structure& s = data.structure;
    const Vec3 c = centroid(s, targets(data, false, id));
    // Species order (VASP POSCAR blocks): listed elements first, in that order.
    std::map<int, int> rank;
    {
        std::string order = text(p, id, "element_order");
        std::replace(order.begin(), order.end(), ',', ' ');
        std::istringstream in(order);
        int k = 0;
        for (std::string w; in >> w;) rank.emplace(science::atomicNumber(standardElement(w)), k++);
    }
    const auto values = [&](const std::string& key) {
        std::vector<double> value(s.atoms.size());
        for (std::size_t i = 0; i < s.atoms.size(); ++i) {
            const auto& a = s.atoms[i];
            if (key == "none") value[i] = 0;
            else if (key == "element") {
                const auto it = rank.find(a.atomicNumber);
                value[i] = it != rank.end() ? it->second - 1000.0 : a.atomicNumber;
            } else if (key == "x") value[i] = a.x;
            else if (key == "y") value[i] = a.y;
            else if (key == "z") value[i] = a.z;
            else if (key == "distance") value[i] = science::norm(science::sub(positionOf(a), c));
            else if (key == "magmom") value[i] = a.magmom;
            else if (key == "charge") value[i] = a.charge;
            else if (key == "mass") value[i] = elementAtomicMass(a.atomicNumber);
            else if (key == "selected") value[i] = i < data.selected.size() && data.selected[i] ? 0 : 1;
            else if (key == "property") {
                if (s.atomProperty.size() != s.atoms.size()) throw std::runtime_error("sort-atoms by property needs a per-atom property");
                value[i] = std::isnan(s.atomProperty[i]) ? HUGE_VAL : s.atomProperty[i];
            } else
                throw std::runtime_error("Unknown sort key " + key);
        }
        return value;
    };
    const auto first = values(text(p, id, "key")), second = values(text(p, id, "then"));
    const double tolerance = number(p, id, "tolerance");
    std::vector<std::size_t> order(s.atoms.size());
    std::iota(order.begin(), order.end(), 0);
    const bool descending = flag(p, id, "descending");
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
        // Values within the tolerance tie on the first key (atoms of one layer), then the second decides.
        if (std::abs(first[x] - first[y]) > tolerance) return descending ? first[x] > first[y] : first[x] < first[y];
        return second[x] < second[y];
    });
    rebuild(data, order);
}

// ---------------------------------------------------------------- chemistry
void saturateHydrogen(const Json& p, PipelineData& data)
{
    const char* id = "saturate-hydrogen";
    std::set<std::string> elements;
    std::istringstream in(text(p, id, "elements"));
    for (std::string w; in >> w;) elements.insert(standardElement(w));
    if (data.selected.size() != data.structure.atoms.size()) data.selected.assign(data.structure.atoms.size(), 0);
    const std::vector<char> mask = data.selected;
    if (text(p, id, "target") == "selected" && std::none_of(mask.begin(), mask.end(), [](char c) { return c; }))
        throw std::runtime_error("saturate-hydrogen: nothing is selected");
    SaturationOptions options;
    options.terminator = standardElement(text(p, id, "terminator"));
    options.bondLength = number(p, id, "bond_length");
    options.maxPerAtom = integerOf(p, id, "max_per_atom", 1);
    if (text(p, id, "target") == "selected") options.mask = &mask;
    if (flag(p, id, "select")) std::fill(data.selected.begin(), data.selected.end(), 0);
    const std::size_t added = saturate(data.structure, data.selected, elements, text(p, id, "hybridization") == "sp3", number(p, id, "bond_scale"),
                                       flag(p, id, "select"), options);
    data.notes.push_back(std::to_string(added) + " " + (options.terminator == "H" ? std::string("hydrogen") : options.terminator) + " atoms added");
}

void swapElements(const Json& p, PipelineData& data)
{
    const char* id = "swap-elements";
    const std::string a = standardElement(text(p, id, "first")), b = standardElement(text(p, id, "second"));
    std::size_t na = 0, nb = 0;
    for (auto& atom : data.structure.atoms) {
        if (atom.symbol == a) { atom.symbol = b; ++na; }
        else if (atom.symbol == b) { atom.symbol = a; ++nb; }
        else continue;
        atom.atomicNumber = science::atomicNumber(atom.symbol);
        getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
    }
    data.notes.push_back(std::to_string(na) + " " + a + " -> " + b + ", " + std::to_string(nb) + " " + b + " -> " + a);
}

// ---------------------------------------------------------------- bond geometry
// Vectors from `atom` to every atom bonded to it through paths that avoid the
// bond (blocked, barrier); the result maps atom -> vector from `atom`.
std::vector<std::pair<std::size_t, Vec3>> fragment(const Neighbours& nb, std::size_t start, std::size_t barrier)
{
    std::vector<std::pair<std::size_t, Vec3>> out{{start, {0, 0, 0}}};
    std::vector<char> seen(nb.atoms.size(), 0);
    seen[start] = 1;
    std::queue<std::pair<std::size_t, Vec3>> queue;
    queue.push({start, {0, 0, 0}});
    while (!queue.empty()) {
        const auto [atom, offset] = queue.front();
        queue.pop();
        for (std::size_t k = 0; k < nb.atoms[atom].size(); ++k) {
            const std::size_t next = nb.atoms[atom][k];
            if (atom == start && next == barrier) continue;
            if (next == barrier) throw std::runtime_error("The atoms are in a ring: the side to move is bonded to the fixed side (use move = atom)");
            if (seen[next]) continue;
            seen[next] = 1;
            const Vec3 to = science::add(offset, nb.vectors[atom][k]);
            out.push_back({next, to});
            queue.push({next, to});
        }
    }
    return out;
}

struct Geometry
{
    std::vector<Vec3> points;  // unwrapped positions of the listed atoms (minimum image from the first)
};

Geometry geometryOf(const Structure& s, const std::vector<std::size_t>& atoms)
{
    Geometry g;
    Vec3 previous = positionOf(s.atoms[atoms[0]]);
    g.points.push_back(previous);
    for (std::size_t k = 1; k < atoms.size(); ++k) {
        Vec3 d = science::sub(positionOf(s.atoms[atoms[k]]), positionOf(s.atoms[atoms[k - 1]]));
        if (s.hasUnitCell) d = minimumImage(s, d);
        previous = science::add(previous, d);
        g.points.push_back(previous);
    }
    return g;
}

double angleOf(const Vec3& a, const Vec3& b, const Vec3& c)
{
    const Vec3 u = science::sub(a, b), v = science::sub(c, b);
    return std::acos(std::clamp(science::dot(u, v) / (science::norm(u) * science::norm(v)), -1.0, 1.0)) / kDegree;
}

double dihedralOf(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
{
    const Vec3 b1 = science::sub(b, a), b2 = science::sub(c, b), b3 = science::sub(d, c);
    const Vec3 n1 = science::cross(b1, b2), n2 = science::cross(b2, b3);
    const Vec3 m = science::cross(n1, science::scale(b2, 1.0 / science::norm(b2)));
    return std::atan2(science::dot(m, n2), science::dot(n1, n2)) / kDegree;
}

double measure(const Geometry& g)
{
    const auto& q = g.points;
    if (q.size() == 2) return science::norm(science::sub(q[1], q[0]));
    if (q.size() == 3) return angleOf(q[0], q[1], q[2]);
    return dihedralOf(q[0], q[1], q[2], q[3]);
}

const char* quantity(std::size_t atoms) { return atoms == 2 ? "distance" : atoms == 3 ? "angle" : "dihedral"; }

void measureGeometry(const Json& p, PipelineData& data)
{
    const char* id = "measure-geometry";
    const auto atoms = parseIndices(text(p, id, "atoms"), data.structure.atoms.size(), id);
    if (atoms.size() < 2 || atoms.size() > 4) throw std::runtime_error("measure-geometry: give 2, 3 or 4 atom indices");
    const double value = measure(geometryOf(data.structure, atoms));
    std::ostringstream note;
    note.precision(6);
    note << quantity(atoms.size()) << " " << text(p, id, "atoms") << " = " << value << (atoms.size() == 2 ? " Angstrom" : " degrees");
    data.notes.push_back(note.str());
}

void setGeometry(const Json& p, PipelineData& data)
{
    const char* id = "set-geometry";
    Structure& s = data.structure;
    const auto atoms = parseIndices(text(p, id, "atoms"), s.atoms.size(), id);
    if (atoms.size() < 2 || atoms.size() > 4) throw std::runtime_error("set-geometry: give 2 (distance), 3 (angle) or 4 (dihedral) atom indices");
    for (std::size_t x = 0; x < atoms.size(); ++x)
        for (std::size_t y = x + 1; y < atoms.size(); ++y)
            if (atoms[x] == atoms[y]) throw std::runtime_error("set-geometry: the atoms must differ");
    const double target = number(p, id, "value");
    const Geometry g = geometryOf(s, atoms);
    const double current = measure(g);
    // The moving side: the last atom (distance, angle) or the third (dihedral), with its fragment.
    const std::size_t mover = atoms.size() == 4 ? atoms[2] : atoms.back();
    const std::size_t fixed = atoms.size() == 4 ? atoms[1] : atoms[atoms.size() - 2];
    const std::string move = text(p, id, "move");
    if (move != "fragment" && move != "atom" && move != "both") throw std::runtime_error("set-geometry: move is fragment, atom or both");
    // Side B (moves forward) and, for move = both, side A (moves back by the other half).
    std::vector<std::pair<std::size_t, Vec3>> moving, other;
    const Neighbours nb = move == "atom" ? Neighbours{} : bondNeighbours(s, number(p, id, "bond_scale"));
    if (move != "atom") moving = fragment(nb, mover, fixed);
    else moving = {{mover, {0, 0, 0}}};
    if (atoms.size() == 4 && move == "atom") moving = {{atoms[3], science::sub(g.points[3], g.points[2])}};
    const Vec3 moverPoint = atoms.size() == 4 ? g.points[2] : g.points.back();
    const std::size_t otherIndex = atoms.size() == 4 ? 1 : 0;
    if (move == "both") other = fragment(nb, atoms[otherIndex], atoms.size() == 4 ? atoms[2] : atoms[1]);
    const Vec3 otherPoint = g.points[otherIndex];
    for (const auto& [atom, offset] : moving) {
        for (std::size_t k = 0; k < atoms.size(); ++k)
            if (atoms[k] == atom && atom != mover && !(atoms.size() == 4 && k == 3))
                throw std::runtime_error("set-geometry: the moving side contains atom " + std::to_string(atom) + " of the fixed side");
        for (const auto& o : other)
            if (o.first == atom) throw std::runtime_error("set-geometry: the two sides are bonded to each other (a ring); use move = atom");
    }
    const double share = move == "both" ? 0.5 : 1.0;  // of the change made by side B
    std::vector<std::pair<std::size_t, Vec3>> result;  // atom -> new position
    if (atoms.size() == 2) {
        if (!(target > 0)) throw std::runtime_error("set-geometry: a distance must be positive");
        const Vec3 direction = unit(science::sub(g.points[1], g.points[0]));
        for (const auto& [atom, offset] : moving)
            result.push_back({atom, science::add(science::add(moverPoint, offset), science::scale(direction, share * (target - current)))});
        for (const auto& [atom, offset] : other)
            result.push_back({atom, science::add(science::add(otherPoint, offset), science::scale(direction, -(1 - share) * (target - current)))});
    } else {
        Vec3 pivot, axis;
        if (atoms.size() == 3) {
            pivot = g.points[1];
            axis = science::cross(science::sub(g.points[0], g.points[1]), science::sub(g.points[2], g.points[1]));
            if (science::norm(axis) < 1e-8) axis = science::cross(science::sub(g.points[0], g.points[1]), std::abs(g.points[0][0] - g.points[1][0]) < 0.9 * science::norm(science::sub(g.points[0], g.points[1])) ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
        } else {
            pivot = g.points[2];
            axis = science::sub(g.points[2], g.points[1]);
        }
        // Turn by the difference; the sign that reaches the target.
        double delta = target - current;
        if (atoms.size() == 4) delta = std::remainder(delta, 360.0);
        const auto turned = [&](double degrees) {
            const Mat3 forward = rotationAbout(axis, share * degrees * kDegree), back = rotationAbout(axis, -(1 - share) * degrees * kDegree);
            std::vector<std::pair<std::size_t, Vec3>> out;
            for (const auto& [atom, offset] : moving)
                out.push_back({atom, science::add(pivot, rotate(forward, science::sub(science::add(moverPoint, offset), pivot)))});
            for (const auto& [atom, offset] : other)
                out.push_back({atom, science::add(pivot, rotate(back, science::sub(science::add(otherPoint, offset), pivot)))});
            return out;
        };
        const auto check = [&](const std::vector<std::pair<std::size_t, Vec3>>& moved) {
            Geometry h = g;
            for (const auto& [atom, position] : moved)
                for (std::size_t k = 0; k < atoms.size(); ++k)
                    if (atoms[k] == atom) h.points[k] = position;
            const double value = measure(h);
            return atoms.size() == 4 ? std::abs(std::remainder(value - target, 360.0)) : std::abs(value - target);
        };
        auto forward = turned(delta), backward = turned(-delta);
        result = check(forward) <= check(backward) ? forward : backward;
        if (atoms.size() == 3 && (target < 0 || target > 180)) throw std::runtime_error("set-geometry: an angle lies between 0 and 180 degrees");
    }
    for (const auto& [atom, position] : result) place(s.atoms[atom], position);
    std::ostringstream note;
    note.precision(6);
    note << quantity(atoms.size()) << " " << current << " -> " << measure(geometryOf(s, atoms)) << "; moved " << result.size() << " atom(s)";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- duplication
void duplicateSelected(const Json& p, PipelineData& data)
{
    const char* id = "duplicate-selected";
    const auto atoms = targets(data, true, id);
    const Vec3 offset = vector(p, id, "offset");
    const int copies = static_cast<int>(number(p, id, "copies"));
    if (copies < 1) throw std::runtime_error("duplicate-selected: copies must be at least 1");
    const bool selectNew = flag(p, id, "select");
    if (selectNew) std::fill(data.selected.begin(), data.selected.end(), 0);
    Structure& s = data.structure;
    for (int c = 1; c <= copies; ++c)
        for (std::size_t i : atoms) {
            const AtomSite original = s.atoms[i];
            appendAtom(s, data.selected, original.symbol, science::add(positionOf(original), science::scale(offset, c)), selectNew);
            AtomSite& copy = s.atoms.back();
            copy.r = original.r; copy.g = original.g; copy.b = original.b;
            copy.fixed = original.fixed;
        }
    data.notes.push_back(std::to_string(atoms.size() * static_cast<std::size_t>(copies)) + " atoms added");
}

// ---------------------------------------------------------------- cell
// Niggli reduction (Krivy & Gruber 1976, with the tolerances of
// Grosse-Kunstleve, Sauter & Adams 2004) on the cell vectors.
Mat3 niggli(Mat3 v, double eps)
{
    const auto dotRows = [&](int i, int j) { return science::dot(v[i], v[j]); };
    const auto flipFor = [&](bool positive) {
        // Signs of a, b, c (an even number of flips) so that 2 b.c, 2 a.c and 2 a.b are all > 0 (or all <= 0).
        const double xi = 2 * dotRows(1, 2), eta = 2 * dotRows(0, 2), zeta = 2 * dotRows(0, 1);
        for (int pass = 0; pass < 2; ++pass)
            for (int mask = 0; mask < 8; ++mask) {
                const double sa = mask & 1 ? -1 : 1, sb = mask & 2 ? -1 : 1, sc = mask & 4 ? -1 : 1;
                if (pass == 0 && sa * sb * sc < 0) continue;
                const double x = sb * sc * xi, y = sa * sc * eta, z = sa * sb * zeta;
                const bool ok = positive ? (x > eps && y > eps && z > eps) : (x <= eps && y <= eps && z <= eps);
                if (!ok) continue;
                v[0] = science::scale(v[0], sa); v[1] = science::scale(v[1], sb); v[2] = science::scale(v[2], sc);
                if (sa * sb * sc < 0) for (auto& row : v) row = science::scale(row, -1);
                return;
            }
    };
    const auto sign = [](double x) { return x > 0 ? 1.0 : -1.0; };
    // Signs within the tolerance (0 for |x| < eps), as Grosse-Kunstleve et al.
    const auto signOf = [&](double x) { return x > eps ? 1 : x < -eps ? -1 : 0; };
    for (int iteration = 0; iteration < 1000; ++iteration) {
        double A = dotRows(0, 0), B = dotRows(1, 1), C = dotRows(2, 2);
        double xi = 2 * dotRows(1, 2), eta = 2 * dotRows(0, 2), zeta = 2 * dotRows(0, 1);
        if (A > B + eps || (std::abs(A - B) < eps && std::abs(xi) > std::abs(eta) + eps)) {
            std::swap(v[0], v[1]); v[2] = science::scale(v[2], -1); continue;
        }
        if (B > C + eps || (std::abs(B - C) < eps && std::abs(eta) > std::abs(zeta) + eps)) {
            std::swap(v[1], v[2]); v[0] = science::scale(v[0], -1); continue;
        }
        if (signOf(xi) * signOf(eta) * signOf(zeta) == 1) flipFor(true);
        else flipFor(false);
        xi = 2 * dotRows(1, 2); eta = 2 * dotRows(0, 2); zeta = 2 * dotRows(0, 1);
        if (std::abs(xi) > B + eps || (std::abs(xi - B) < eps && 2 * eta < zeta - eps) || (std::abs(xi + B) < eps && zeta < -eps)) {
            v[2] = science::sub(v[2], science::scale(v[1], sign(xi))); continue;
        }
        if (std::abs(eta) > A + eps || (std::abs(eta - A) < eps && 2 * xi < zeta - eps) || (std::abs(eta + A) < eps && zeta < -eps)) {
            v[2] = science::sub(v[2], science::scale(v[0], sign(eta))); continue;
        }
        if (std::abs(zeta) > A + eps || (std::abs(zeta - A) < eps && 2 * xi < eta - eps) || (std::abs(zeta + A) < eps && eta < -eps)) {
            v[1] = science::sub(v[1], science::scale(v[0], sign(zeta))); continue;
        }
        if (xi + eta + zeta + A + B < -eps || (std::abs(xi + eta + zeta + A + B) < eps && 2 * (A + eta) + zeta > eps)) {
            v[2] = science::add(v[2], science::add(v[0], v[1])); continue;
        }
        return v;
    }
    throw std::runtime_error("Niggli reduction did not converge");
}

void niggliReduce(const Json&, PipelineData& data)
{
    requireCell(data, "niggli-reduce");
    Structure& s = data.structure;
    if (!(s.pbc[0] && s.pbc[1] && s.pbc[2])) throw std::runtime_error("niggli-reduce needs a cell periodic along a, b and c");
    const Mat3 cell = cellOf(s);
    const double volume = std::abs(science::determinant(cell));
    const Mat3 reduced = niggli(cell, 1e-5 * std::pow(volume, 2.0 / 3.0));
    for (int r = 0; r < 3; ++r) s.cellVectors[r] = reduced[r];
    wrapAll(s);
    std::ostringstream note;
    note.precision(6);
    const auto len = [&](int k) { return science::norm(reduced[k]); };
    const auto ang = [&](int i, int j) { return std::acos(science::dot(reduced[i], reduced[j]) / (len(i) * len(j))) / kDegree; };
    note << "a b c = " << len(0) << " " << len(1) << " " << len(2) << " Angstrom; alpha beta gamma = " << ang(1, 2) << " " << ang(0, 2) << " " << ang(0, 1);
    data.notes.push_back(note.str());
}

void scaleCell(const Json& p, PipelineData& data)
{
    const char* id = "scale-cell";
    requireCell(data, id);
    Structure& s = data.structure;
    const std::string mode = text(p, id, "mode");
    const double value = number(p, id, "value");
    if (!(value > 0)) throw std::runtime_error("scale-cell: the value must be positive");
    const double volume = std::abs(science::determinant(cellOf(s)));
    double mass = 0;
    for (const auto& a : s.atoms) mass += elementAtomicMass(a.atomicNumber);
    if (mode == "axes") {
        // Lattice vectors scaled independently (uniaxial or biaxial strain); fractional coordinates kept.
        const Vec3 f = vector(p, id, "factors");
        if (!(f[0] > 0 && f[1] > 0 && f[2] > 0)) throw std::runtime_error("scale-cell: the factors must be positive");
        const Mat3 inverse = science::inverse(cellOf(s));
        if (flag(p, id, "scale_atoms"))
            for (auto& a : s.atoms) {
                Vec3 fr = science::rowTimes(science::sub(positionOf(a), s.cellOffset), inverse);
                for (int k = 0; k < 3; ++k) fr[k] *= f[k];
                place(a, science::add(s.cellOffset, science::rowTimes(fr, cellOf(s))));
            }
        for (int k = 0; k < 3; ++k) s.cellVectors[k] = science::scale(s.cellVectors[k], f[k]);
        const double v = std::abs(science::determinant(cellOf(s)));
        std::ostringstream note;
        note.precision(6);
        note << "a, b, c x " << f[0] << ", " << f[1] << ", " << f[2] << " (strains " << 100 * (f[0] - 1) << "%, " << 100 * (f[1] - 1) << "%, "
             << 100 * (f[2] - 1) << "%); volume " << v << " Angstrom^3";
        if (!s.atoms.empty()) note << ", density " << mass * 1.66053906660 / v << " g/cm3";
        data.notes.push_back(note.str());
        return;
    }
    double factor = value;
    if (mode == "volume") factor = std::cbrt(value / volume);
    else if (mode == "density") {
        if (s.atoms.empty()) throw std::runtime_error("scale-cell to a density needs atoms");
        factor = std::cbrt(mass * 1.66053906660 / value / volume);
    } else if (mode != "factor")
        throw std::runtime_error("scale-cell: mode is factor, volume, density or axes");
    if (flag(p, id, "scale_atoms"))
        for (auto& a : s.atoms) place(a, science::add(s.cellOffset, science::scale(science::sub(positionOf(a), s.cellOffset), factor)));
    for (auto& row : s.cellVectors) row = science::scale(row, factor);
    std::ostringstream note;
    note.precision(6);
    note << "Lengths x " << factor << "; volume " << volume * factor * factor * factor << " Angstrom^3";
    if (!s.atoms.empty()) note << ", density " << mass * 1.66053906660 / (volume * factor * factor * factor) << " g/cm3";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- molecules
std::vector<int> components(const Neighbours& nb)
{
    std::vector<int> label(nb.atoms.size(), -1);
    int next = 0;
    for (std::size_t start = 0; start < nb.atoms.size(); ++start) {
        if (label[start] >= 0) continue;
        std::queue<std::size_t> queue;
        queue.push(start);
        label[start] = next;
        while (!queue.empty()) {
            const std::size_t atom = queue.front();
            queue.pop();
            for (std::size_t other : nb.atoms[atom])
                if (label[other] < 0) { label[other] = next; queue.push(other); }
        }
        ++next;
    }
    return label;
}

void unwrapMolecules(const Json& p, PipelineData& data)
{
    const char* id = "unwrap-molecules";
    Structure& s = data.structure;
    const Neighbours nb = bondNeighbours(s, number(p, id, "bond_scale"));
    std::vector<char> done(s.atoms.size(), 0);
    int molecules = 0;
    std::size_t moved = 0;
    for (std::size_t start = 0; start < s.atoms.size(); ++start) {
        if (done[start]) continue;
        ++molecules;
        done[start] = 1;
        std::queue<std::size_t> queue;
        queue.push(start);
        while (!queue.empty()) {
            const std::size_t atom = queue.front();
            queue.pop();
            for (std::size_t k = 0; k < nb.atoms[atom].size(); ++k) {
                const std::size_t other = nb.atoms[atom][k];
                if (done[other]) continue;
                done[other] = 1;
                const Vec3 target = science::add(positionOf(s.atoms[atom]), nb.vectors[atom][k]);
                if (science::norm(science::sub(target, positionOf(s.atoms[other]))) > 1e-6) ++moved;
                place(s.atoms[other], target);
                queue.push(other);
            }
        }
    }
    data.notes.push_back(std::to_string(molecules) + " molecules or fragments; " + std::to_string(moved) + " atoms moved next to their bonded neighbours");
}

void selectFragment(const Json& p, PipelineData& data)
{
    const char* id = "select-fragment";
    const Neighbours nb = bondNeighbours(data.structure, number(p, id, "bond_scale"));
    const std::vector<int> label = components(nb);
    const int count = label.empty() ? 0 : *std::max_element(label.begin(), label.end()) + 1;
    std::vector<std::size_t> size(static_cast<std::size_t>(count), 0);
    for (int l : label) ++size[static_cast<std::size_t>(l)];
    std::vector<char> chosen(static_cast<std::size_t>(count), 0);
    if (text(p, id, "rule") == "touching") {
        for (std::size_t i = 0; i < label.size(); ++i)
            if (data.selected[i]) chosen[static_cast<std::size_t>(label[i])] = 1;
    } else {
        const double low = number(p, id, "min_atoms"), high = number(p, id, "max_atoms");
        for (int l = 0; l < count; ++l) {
            const double n = static_cast<double>(size[static_cast<std::size_t>(l)]);
            chosen[static_cast<std::size_t>(l)] = n >= low && (high <= 0 || n <= high);
        }
    }
    std::vector<char> mask(label.size(), 0);
    int fragments = 0;
    for (int l = 0; l < count; ++l) fragments += chosen[static_cast<std::size_t>(l)];
    for (std::size_t i = 0; i < label.size(); ++i) mask[i] = chosen[static_cast<std::size_t>(label[i])];
    select(data, mask, text(p, id, "mode"));
    data.notes.push_back(std::to_string(fragments) + " of " + std::to_string(count) + " fragments");
}
// ---------------------------------------------------------------- per-atom simulation data
// VASP's MAGMOM line: runs of equal moments as count*value.
std::string magmomLine(const Structure& s)
{
    std::ostringstream out;
    out.precision(6);
    for (std::size_t i = 0; i < s.atoms.size();) {
        std::size_t j = i;
        while (j < s.atoms.size() && s.atoms[j].magmom == s.atoms[i].magmom) ++j;
        out << (i ? " " : "") << (j - i) << "*" << s.atoms[i].magmom;
        i = j;
    }
    return out.str();
}

// "Fe 2.5 O 0": a value per element.
std::map<int, double> elementValues(std::string list, const char* what)
{
    std::replace(list.begin(), list.end(), ',', ' ');
    std::istringstream in(list);
    std::map<int, double> out;
    for (std::string symbol, value; in >> symbol;) {
        if (!(in >> value)) throw std::runtime_error(std::string(what) + ": give pairs of an element and a value, e.g. Fe 2.5 O 0");
        char* end = nullptr;
        const double v = std::strtod(value.c_str(), &end);
        if (*end) throw std::runtime_error(std::string(what) + ": '" + value + "' is not a number");
        out[science::atomicNumber(standardElement(symbol))] = v;
    }
    if (out.empty()) throw std::runtime_error(std::string(what) + ": no element values given");
    return out;
}

// Layer index of each atom along a Cartesian axis (layers: positions within the tolerance).
std::vector<int> layersAlong(const Structure& s, const std::vector<std::size_t>& atoms, int axis, double tolerance)
{
    std::vector<std::size_t> order(atoms);
    const auto coordinate = [&](std::size_t i) { return positionOf(s.atoms[i])[axis]; };
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return coordinate(a) < coordinate(b); });
    std::map<std::size_t, int> layer;
    int current = -1;
    double start = -HUGE_VAL;
    for (std::size_t i : order) {
        if (coordinate(i) > start + tolerance) { ++current; start = coordinate(i); }
        layer[i] = current;
    }
    std::vector<int> out;
    for (std::size_t i : atoms) out.push_back(layer[i]);
    return out;
}

int axisIndex(const std::string& name) { return name == "x" ? 0 : name == "y" ? 1 : 2; }

void setMagneticMoments(const Json& p, PipelineData& data)
{
    const char* id = "set-magnetic-moments";
    const auto atoms = targets(data, text(p, id, "target") == "selected", id);
    const std::string mode = text(p, id, "mode");
    const double value = number(p, id, "value");
    Structure& s = data.structure;
    std::ostringstream what;
    what.precision(6);
    if (mode == "value") {
        for (std::size_t i : atoms) s.atoms[i].magmom = value;
        what << "set to " << value << " Bohr magnetons";
    } else if (mode == "elements") {
        const auto values = elementValues(text(p, id, "moments"), id);
        for (std::size_t i : atoms)
            if (const auto it = values.find(s.atoms[i].atomicNumber); it != values.end()) s.atoms[i].magmom = it->second;
        what << "set per element";
    } else if (mode == "layered") {
        // Layered antiferromagnet (A-type, ...): the sign alternates every `period` layers along the axis.
        const auto layer = layersAlong(s, atoms, axisIndex(text(p, id, "axis")), number(p, id, "layer_tolerance"));
        const int period = integerOf(p, id, "period", 1);
        for (std::size_t k = 0; k < atoms.size(); ++k) s.atoms[atoms[k]].magmom = (layer[k] / period) % 2 ? -value : value;
        what << "+-" << value << " alternating every " << period << " layer(s) along " << text(p, id, "axis") << " ("
             << (atoms.empty() ? 0 : *std::max_element(layer.begin(), layer.end()) + 1) << " layers)";
    } else if (mode == "checkerboard") {
        // G-type: the sign of (-1)^(i+j+k) over layers along x, y and z.
        const double tolerance = number(p, id, "layer_tolerance");
        const auto lx = layersAlong(s, atoms, 0, tolerance), ly = layersAlong(s, atoms, 1, tolerance), lz = layersAlong(s, atoms, 2, tolerance);
        for (std::size_t k = 0; k < atoms.size(); ++k) s.atoms[atoms[k]].magmom = (lx[k] + ly[k] + lz[k]) % 2 ? -value : value;
        what << "+-" << value << " in a checkerboard (G-type)";
    } else if (mode == "random") {
        // Paramagnetic-like: random signs with zero net moment where possible (half up, half down).
        std::vector<std::size_t> order(atoms);
        std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
        std::shuffle(order.begin(), order.end(), rng);
        for (std::size_t k = 0; k < order.size(); ++k) s.atoms[order[k]].magmom = k < order.size() / 2 ? -value : value;
        what << "+-" << value << " at random (half up, half down)";
    } else
        throw std::runtime_error("set-magnetic-moments: mode is value, elements, layered, checkerboard or random");
    double total = 0, absolute = 0;
    for (const auto& a : s.atoms) { total += a.magmom; absolute += std::abs(a.magmom); }
    std::ostringstream note;
    note.precision(6);
    note << atoms.size() << " atoms " << what.str() << "; total " << total << ", absolute " << absolute << "; MAGMOM = " << magmomLine(s)
         << " (in this atom order)";
    data.notes.push_back(note.str());
}

void setCharges(const Json& p, PipelineData& data)
{
    const char* id = "set-charges";
    const auto atoms = targets(data, text(p, id, "target") == "selected", id);
    Structure& s = data.structure;
    const std::string mode = text(p, id, "mode");
    if (mode == "value") {
        const double value = number(p, id, "value");
        for (std::size_t i : atoms) s.atoms[i].charge = value;
    } else if (mode == "elements") {
        const auto values = elementValues(text(p, id, "charges"), id);
        for (std::size_t i : atoms)
            if (const auto it = values.find(s.atoms[i].atomicNumber); it != values.end()) s.atoms[i].charge = it->second;
    } else
        throw std::runtime_error("set-charges: mode is value or elements");
    double total = 0;
    for (const auto& a : s.atoms) total += a.charge;
    std::ostringstream note;
    note.precision(6);
    // Neutralising: the excess spread equally over the target atoms, or the target charges scaled.
    const std::string neutralize = text(p, id, "neutralize");
    if (neutralize == "shift" && !atoms.empty() && std::abs(total) > 1e-12) {
        const double shift = total / atoms.size();
        for (std::size_t i : atoms) s.atoms[i].charge -= shift;
        note << "shifted by " << -shift << " e to neutralise; ";
    } else if (neutralize == "scale" && std::abs(total) > 1e-12) {
        double positive = 0, negative = 0, other = total;
        for (std::size_t i : atoms) {
            (s.atoms[i].charge > 0 ? positive : negative) += s.atoms[i].charge;
            other -= s.atoms[i].charge;
        }
        // Scale the side in excess so the whole is neutral.
        if (total > 0 && positive > 0) {
            const double f = (-negative - other) / positive;
            for (std::size_t i : atoms) if (s.atoms[i].charge > 0) s.atoms[i].charge *= f;
            note << "positive charges x " << f << " to neutralise; ";
        } else if (total < 0 && negative < 0) {
            const double f = (-positive - other) / negative;
            for (std::size_t i : atoms) if (s.atoms[i].charge < 0) s.atoms[i].charge *= f;
            note << "negative charges x " << f << " to neutralise; ";
        }
    } else if (neutralize != "none" && neutralize != "shift" && neutralize != "scale")
        throw std::runtime_error("set-charges: neutralize is none, shift or scale");
    total = 0;
    for (const auto& a : s.atoms) total += a.charge;
    note << atoms.size() << " atoms set; total charge " << total << " e" << (std::abs(total) > 1e-9 ? " (not neutral)" : " (neutral)");
    data.notes.push_back(note.str());
}

// Maxwell-Boltzmann velocities at a temperature (Angstrom/fs), without net
// momentum, scaled so the kinetic temperature (3N - 3 degrees of freedom) is exact.
void setVelocities(const Json& p, PipelineData& data)
{
    const char* id = "set-velocities";
    const auto atoms = targets(data, text(p, id, "target") == "selected", id);
    const double temperature = number(p, id, "temperature");
    if (temperature < 0) throw std::runtime_error("set-velocities: the temperature must not be negative");
    constexpr double kB = 8.617333262e-5;          // eV/K
    constexpr double eVPerAmuA2fs2 = 103.6427;     // 1 u Angstrom^2/fs^2 in eV
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    std::normal_distribution<double> gauss(0.0, 1.0);
    Structure& s = data.structure;
    std::vector<double> mass;
    for (std::size_t i : atoms) mass.push_back(elementAtomicMass(s.atoms[i].atomicNumber));
    std::vector<Vec3> v(atoms.size());
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        const double sigma = std::sqrt(kB * temperature / (mass[k] * eVPerAmuA2fs2));
        v[k] = {gauss(rng) * sigma, gauss(rng) * sigma, gauss(rng) * sigma};
    }
    const bool zeroMomentum = flag(p, id, "zero_momentum") && atoms.size() > 1;
    if (zeroMomentum) {
        Vec3 momentum{0, 0, 0};
        double total = 0;
        for (std::size_t k = 0; k < atoms.size(); ++k) { momentum = science::add(momentum, science::scale(v[k], mass[k])); total += mass[k]; }
        const Vec3 drift = science::scale(momentum, 1.0 / total);
        for (auto& vk : v) vk = science::sub(vk, drift);
    }
    // No rotation (clusters, molecules in vacuum): the angular velocity w = I^-1 L about the centre of mass removed.
    const bool zeroRotation = flag(p, id, "zero_rotation") && atoms.size() > 2;
    if (zeroRotation) {
        Vec3 com{0, 0, 0};
        double total = 0;
        for (std::size_t k = 0; k < atoms.size(); ++k) { com = science::add(com, science::scale(positionOf(s.atoms[atoms[k]]), mass[k])); total += mass[k]; }
        com = science::scale(com, 1.0 / total);
        Vec3 angular{0, 0, 0};
        Mat3 inertia{};
        for (std::size_t k = 0; k < atoms.size(); ++k) {
            const Vec3 r = science::sub(positionOf(s.atoms[atoms[k]]), com);
            angular = science::add(angular, science::scale(science::cross(r, v[k]), mass[k]));
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) inertia[i][j] += mass[k] * ((i == j ? science::dot(r, r) : 0.0) - r[i] * r[j]);
        }
        if (std::abs(science::determinant(inertia)) > 1e-12) {
            const Vec3 omega = rotate(science::inverse(inertia), angular);
            for (std::size_t k = 0; k < atoms.size(); ++k)
                v[k] = science::sub(v[k], science::cross(omega, science::sub(positionOf(s.atoms[atoms[k]]), com)));
        }
    }
    double kinetic = 0;  // eV
    for (std::size_t k = 0; k < atoms.size(); ++k) kinetic += 0.5 * mass[k] * science::dot(v[k], v[k]) * eVPerAmuA2fs2;
    const double dof = 3.0 * atoms.size() - (zeroMomentum ? 3.0 : 0.0) - (zeroRotation ? 3.0 : 0.0);
    if (temperature > 0 && kinetic > 0 && dof > 0) {
        const double factor = std::sqrt(temperature / (2 * kinetic / (dof * kB)));
        for (auto& vk : v) vk = science::scale(vk, factor);
    }
    for (std::size_t k = 0; k < atoms.size(); ++k) s.atoms[atoms[k]].velocity = v[k];
    std::ostringstream note;
    note.precision(6);
    note << atoms.size() << " velocities at " << temperature << " K (Angstrom/fs" << (zeroMomentum ? ", no net momentum" : "")
         << (zeroRotation ? ", no rotation" : "") << ")";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- stacking
// Places another structure's atoms on top of this one along c (a van der
// Waals heterostructure, a multilayer): the gap is between the top of this
// structure and the bottom of the other, c grows to hold both, and the
// other's in-plane coordinates follow this cell (strained to it) or stay.
void stackLayers(const Json& p, PipelineData& data)
{
    const char* id = "stack-layers";
    requireCell(data, id);
    const std::string file = text(p, id, "file");
    if (file.empty()) throw std::runtime_error("Choose the structure file to stack");
    const auto frames = science::readFrames(std::filesystem::u8path(file));
    const Structure& top = frames.front().structure;
    if (top.atoms.empty()) throw std::runtime_error("The structure to stack has no atoms");
    Structure& s = data.structure;
    const double gap = number(p, id, "gap"), vacuum = number(p, id, "vacuum");
    const bool strain = flag(p, id, "strain") && top.hasUnitCell;
    const std::string match = text(p, id, "match");
    if (match != "layer" && match != "average" && match != "substrate") throw std::runtime_error("stack-layers: match is layer, average or substrate");
    const int copies = integerOf(p, id, "copies", 1);
    const Vec3 shift = vector(p, id, "shift");  // fractional in-plane offset per copy
    const double twist = number(p, id, "twist");
    Mat3 cell = cellOf(s);
    Vec3 normal = science::cross(cell[0], cell[1]);
    normal = science::scale(normal, 1.0 / science::norm(normal));
    if (science::dot(normal, cell[2]) < 0) normal = science::scale(normal, -1);
    std::vector<Vec3> placed;  // the other structure's atoms in this cell's frame, before lifting
    if (strain) {
        const Mat3 otherCell{Vec3{top.cellVectors[0][0], top.cellVectors[0][1], top.cellVectors[0][2]},
                             Vec3{top.cellVectors[1][0], top.cellVectors[1][1], top.cellVectors[1][2]},
                             Vec3{top.cellVectors[2][0], top.cellVectors[2][1], top.cellVectors[2][2]}};
        // The common in-plane lattice: this cell's (the layer strains), the layer's (this structure strains), or between.
        const double la = science::norm(otherCell[0]) / science::norm(cell[0]), lb = science::norm(otherCell[1]) / science::norm(cell[1]);
        const double fa = match == "layer" ? 1 : match == "substrate" ? la : 0.5 * (1 + la);
        const double fb = match == "layer" ? 1 : match == "substrate" ? lb : 0.5 * (1 + lb);
        if (fa != 1 || fb != 1) {
            const Mat3 inverse = science::inverse(cell);
            for (auto& a : s.atoms) {
                Vec3 f = science::rowTimes(science::sub(positionOf(a), s.cellOffset), inverse);
                place(a, science::add(s.cellOffset, science::add(science::add(science::scale(cell[0], f[0] * fa), science::scale(cell[1], f[1] * fb)),
                                                                 science::scale(cell[2], f[2]))));
            }
            s.cellVectors[0] = science::scale(cell[0], fa);
            s.cellVectors[1] = science::scale(cell[1], fb);
            cell = cellOf(s);
        }
        Vec3 otherNormal = science::cross(otherCell[0], otherCell[1]);
        otherNormal = science::scale(otherNormal, 1.0 / science::norm(otherNormal));
        const Mat3 inverse = science::inverse(otherCell);
        for (const auto& a : top.atoms) {
            const Vec3 f = science::rowTimes(positionOf(a), inverse);
            const Vec3 inPlane = science::add(science::scale(cell[0], f[0]), science::scale(cell[1], f[1]));
            placed.push_back(science::add(inPlane, science::scale(normal, science::dot(positionOf(a), otherNormal))));
        }
        std::ostringstream note;
        note.precision(4);
        note << "Strains: stacked layer " << 100 * (fa / la - 1) << "% along a and " << 100 * (fb / lb - 1) << "% along b; this structure "
             << 100 * (fa - 1) << "% and " << 100 * (fb - 1) << "%";
        data.notes.push_back(note.str());
    } else
        for (const auto& a : top.atoms) placed.push_back(positionOf(a));
    if (twist != 0) {
        // A turn about the normal through the layer's centroid (commensurate only for special angles).
        Vec3 c{0, 0, 0};
        for (const Vec3& q : placed) c = science::add(c, q);
        c = science::scale(c, 1.0 / placed.size());
        const Mat3 r = rotationAbout(normal, twist * kDegree);
        for (Vec3& q : placed) q = science::add(c, rotate(r, science::sub(q, c)));
        data.notes.push_back("The stacked layer is turned about the normal; unless the angle is commensurate with the cell, the stack is not periodic");
    }
    const std::size_t before = s.atoms.size();
    if (data.selected.size() != before) data.selected.assign(before, 0);
    std::fill(data.selected.begin(), data.selected.end(), 0);
    double low = HUGE_VAL;
    for (const Vec3& q : placed) low = std::min(low, science::dot(q, normal));
    for (int copy = 1; copy <= copies; ++copy) {
        double high = -HUGE_VAL;
        for (const auto& a : s.atoms) high = std::max(high, science::dot(positionOf(a), normal));
        if (s.atoms.empty()) high = 0;
        const Vec3 lateral = science::add(science::scale(cell[0], shift[0] * copy), science::scale(cell[1], shift[1] * copy));
        const Vec3 lift = science::add(science::scale(normal, high + gap - low), lateral);
        for (std::size_t k = 0; k < top.atoms.size(); ++k) {
            appendAtom(s, data.selected, top.atoms[k].symbol, science::add(placed[k], lift), flag(p, id, "select"));
            s.atoms.back().magmom = top.atoms[k].magmom;
            s.atoms.back().charge = top.atoms[k].charge;
            s.atoms.back().fixed = top.atoms[k].fixed;
        }
    }
    if (shift[0] != 0 || shift[1] != 0) {
        // The shifted copies back into the cell along a and b (not c, which is about to grow).
        const Mat3 inverse = science::inverse(cell);
        for (std::size_t i = before; i < s.atoms.size(); ++i) {
            Vec3 f = science::rowTimes(science::sub(positionOf(s.atoms[i]), s.cellOffset), inverse);
            f[0] -= std::floor(f[0]);
            f[1] -= std::floor(f[1]);
            place(s.atoms[i], science::add(science::rowTimes(f, cell), s.cellOffset));
        }
    }
    // c grows (along its own direction) so the stack plus vacuum fits.
    double newHigh = -HUGE_VAL, newLow = HUGE_VAL;
    for (const auto& a : s.atoms) {
        newHigh = std::max(newHigh, science::dot(positionOf(a), normal));
        newLow = std::min(newLow, science::dot(positionOf(a), normal));
    }
    const double needed = newHigh - newLow + vacuum;
    const double along = science::dot(cell[2], normal);
    if (needed > along) s.cellVectors[2] = science::scale(cell[2], needed / along);
    std::ostringstream note;
    note.precision(6);
    note << copies << " x " << top.atoms.size() << " atoms stacked " << gap << " Angstrom apart; c = " << science::norm(s.cellVectors[2]) << " Angstrom";
    data.notes.push_back(note.str());
}
}  // namespace

std::vector<ModifierType> editStepTypes()
{
    const ModifierParameter target = {"target", "Atoms", "choice", "all", "all|selected"};
    const ModifierParameter centre = {"center", "About", "choice", "centroid", "centroid|cell|origin|point"};
    const ModifierParameter point = {"point", "Point (Angstrom; for center = point)", "vector", "0 0 0", ""};
    return {
        {"mirror", "Mirror", "Modification",
         "Reflects atoms through the plane with the given normal through the atoms' centroid, the cell centre, the origin or a point. "
         "Reflecting all atoms reflects the cell too (made right-handed again), so a crystal stays periodic.",
         {{"normal", "Plane normal", "vector", "0 0 1", ""}, centre, point, target, {"wrap", "Wrap into the cell", "bool", "false", ""}}, mirror},
        {"invert", "Invert through a point", "Modification", "Inversion x -> 2p - x through the centroid, the cell centre, the origin or a point.",
         {centre, point, target, {"wrap", "Wrap into the cell", "bool", "false", ""}}, invert},
        {"rotate", "Rotate", "Modification",
         "Turns atoms by an angle (degrees, right-handed) about an axis, or by Euler angles (intrinsic z-y-z or z-x-z, or extrinsic x-y-z), "
         "about the centroid, the cell centre, the origin or a point. All atoms turn with the cell by default, about the cell origin, so a "
         "crystal keeps its periodicity and fractional coordinates.",
         {{"axis", "Axis", "vector", "0 0 1", ""}, {"angle", "Angle (degrees)", "float", "90", ""}, {"mode", "Rotation", "choice", "axis", "axis|euler"},
          {"euler", "Euler angles (degrees)", "vector", "0 0 0", ""}, {"convention", "Euler convention", "choice", "zyz", "zyz|zxz|xyz"}, centre, point, target,
          {"rotate_cell", "Rotate the cell too (all atoms)", "bool", "true", ""}}, rotateStep},
        {"align-direction", "Align a direction with an axis", "Cell",
         "Turns the whole structure and its cell so that a direction (lattice [u v w] or Cartesian) lies along x, y or z, and optionally a second "
         "direction in the plane of that axis and a second one (as Atomsk -orient).",
         {{"direction", "Direction", "vector", "1 1 1", ""}, {"lattice", "Lattice direction [u v w]", "bool", "true", ""},
          {"axis", "Along", "choice", "z", "x|y|z"}, {"second_direction", "Second direction", "vector", "1 -1 0", ""},
          {"second_axis", "Second direction towards", "choice", "none", "none|x|y|z"}}, alignDirection},
        {"sort-atoms", "Sort atoms", "Modification",
         "Reorders the atoms by element (atomic number, or a given species order such as O Ti for POSCAR blocks), a coordinate, distance from "
         "the centroid, magnetic moment, charge, mass, selection or the per-atom property; ties (within the tolerance, e.g. one layer) are "
         "broken by a second key.",
         {{"key", "Sort by", "choice", "element", "element|x|y|z|distance|magmom|charge|mass|selected|property"},
          {"then", "Then by", "choice", "none", "none|element|x|y|z|distance|magmom|charge|mass|selected|property"},
          {"element_order", "Species order (elements first; empty = atomic number)", "string", "", ""},
          {"tolerance", "Tie tolerance", "float", "0", ""}, {"descending", "Descending", "bool", "false", ""}}, sortAtoms},
        {"saturate-hydrogen", "Saturate with hydrogen", "Modification",
         "Adds hydrogen (or another terminator: F, Cl, OH-like O, ...) to under-coordinated atoms (bonds from covalent radii): tetrahedral "
         "(sp3) directions, or in-plane for planar carbon, boron and nitrogen with two bonds near 120 degrees; X-H lengths C 1.09, Si 1.48, "
         "N 1.01, O 0.97 Angstrom, other terminators at the sum of covalent radii unless a length is given. Only selected atoms, or at most "
         "a number per atom, on request.",
         {{"elements", "Elements to saturate", "elements", "C Si Ge N O B P S", ""}, {"hybridization", "Geometry", "choice", "auto", "auto|sp3"},
          {"terminator", "Terminating element", "element", "H", ""},
          {"bond_length", "Bond length (Angstrom; 0 = table or covalent radii)", "float", "0", ""},
          {"max_per_atom", "At most per atom", "int", "4", ""}, {"target", "Atoms to saturate", "choice", "all", "all|selected"},
          {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}, {"select", "Select the added atoms", "bool", "true", ""}}, saturateHydrogen},
        {"swap-elements", "Swap elements", "Modification", "Exchanges two elements everywhere (A becomes B and B becomes A).",
         {{"first", "Element A", "element", "Ni", ""}, {"second", "Element B", "element", "Al", ""}}, swapElements},
        {"set-geometry", "Set distance, angle or dihedral", "Modification",
         "Sets the distance (2 atoms), angle (3) or dihedral (4; atom indices from 0) by moving the last atom, or the third and everything beyond the "
         "middle bond of a dihedral, with its bonded fragment (as in Avogadro); or both sides by half each.",
         {{"atoms", "Atom indices", "list", "0 1", ""}, {"value", "Value (Angstrom or degrees)", "float", "1.5", ""},
          {"move", "Move", "choice", "fragment", "fragment|atom|both"}, {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}}, setGeometry},
        {"measure-geometry", "Measure distance, angle or dihedral", "Analysis",
         "Reports the distance (2 atoms), angle (3) or dihedral angle (4) between atoms given by index (from 0), with minimum images in periodic cells.",
         {{"atoms", "Atom indices", "list", "0 1 2 3", ""}}, measureGeometry},
        {"duplicate-selected", "Duplicate selection", "Modification", "Adds copies of the selected atoms shifted by k x offset (k = 1 ... copies).",
         {{"offset", "Offset (Angstrom)", "vector", "0 0 3", ""}, {"copies", "Copies", "int", "1", ""}, {"select", "Select the copies", "bool", "true", ""}}, duplicateSelected},
        {"niggli-reduce", "Niggli-reduced cell", "Cell",
         "Replaces the cell by its Niggli-reduced form (shortest, most orthogonal basis of the same lattice; Krivy-Gruber) and wraps the atoms into it.",
         {}, niggliReduce},
        {"scale-cell", "Scale cell", "Cell",
         "Scales the cell and atom positions uniformly by a factor, or to a volume (Angstrom^3) or a density (g/cm3); or a, b and c by separate "
         "factors (uniaxial and biaxial strain). Without scaling the atoms, their Cartesian positions stay.",
         {{"mode", "Scale to", "choice", "factor", "factor|volume|density|axes"}, {"value", "Value", "float", "1.01", ""},
          {"factors", "Factors for a, b, c (axes mode)", "vector", "1 1 1", ""}, {"scale_atoms", "Scale the atom positions", "bool", "true", ""}}, scaleCell},
        {"unwrap-molecules", "Make molecules whole", "Modification",
         "Moves atoms across periodic boundaries so that bonded atoms sit next to each other (molecules split by the cell become whole).",
         {{"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}}, unwrapMolecules},
        {"set-magnetic-moments", "Set magnetic moments", "Modification",
         "Sets the initial magnetic moment (Bohr magnetons) of the selected (or all) atoms, for spin-polarised calculations: written to "
         "extended XYZ (initial_magmoms), CASTEP (SPIN=) and Quantum ESPRESSO (species per moment, nspin = 2); the notes give VASP's MAGMOM line. "
         "Modes: one value; a value per element (Fe 2.5 O 0); layered antiferromagnets (sign alternating every `period` layers along an axis: "
         "A-type); checkerboard (G-type); random signs with zero net moment. Other orders: select a sublattice and set it.",
         {{"value", "Moment (Bohr magnetons)", "float", "2", ""}, {"mode", "Moments", "choice", "value", "value|elements|layered|checkerboard|random"},
          {"moments", "Per element (e.g. Fe 2.5 O 0)", "string", "", ""}, {"axis", "Layer axis", "choice", "z", "x|y|z"},
          {"period", "Layers per sign", "int", "1", ""}, {"layer_tolerance", "Layer tolerance (Angstrom)", "float", "0.5", ""},
          {"seed", "Seed (random)", "int", "1", ""}, {"target", "Atoms", "choice", "all", "all|selected"}}, setMagneticMoments},
        {"set-charges", "Set charges", "Modification",
         "Sets the charge (e) of the selected (or all) atoms, e.g. formal ionic charges for a LAMMPS charge or full model; written to "
         "LAMMPS data (atom style charge) and extended XYZ (initial_charges); the notes give the total charge. One value, or a value per element; "
         "neutralising shifts the target atoms' charges equally, or scales the side in excess.",
         {{"value", "Charge (e)", "float", "1", ""}, {"mode", "Charges", "choice", "value", "value|elements"},
          {"charges", "Per element (e.g. Ti 2.4 O -1.2)", "string", "", ""},
          {"neutralize", "Neutralise", "choice", "none", "none|shift|scale"}, {"target", "Atoms", "choice", "all", "all|selected"}}, setCharges},
        {"set-velocities", "Set velocities (temperature)", "Modification",
         "Gives the selected (or all) atoms Maxwell-Boltzmann velocities at a temperature (Angstrom/fs), without net momentum and scaled to "
         "that temperature exactly, reproducibly for a seed; written to LAMMPS data (Velocities, Angstrom/ps) and extended XYZ (velo). "
         "A temperature of 0 sets them to zero.",
         {{"temperature", "Temperature (K)", "float", "300", ""}, {"seed", "Seed", "int", "1", ""},
          {"zero_momentum", "Remove net momentum", "bool", "true", ""}, {"zero_rotation", "Remove rotation (clusters, molecules)", "bool", "false", ""},
          {"target", "Atoms", "choice", "all", "all|selected"}}, setVelocities},
        {"stack-layers", "Stack a layer", "Modification",
         "Places the atoms of another structure file on top of this one along c, `gap` above the top atom (van der Waals heterostructures, "
         "multilayers), copies times. With strain, the two share one in-plane lattice: this cell's (the layer strains), the layer's (this "
         "structure strains) or the average (both strain; the strains are reported). Each copy can be shifted laterally (fractional u v, "
         "for AB stacking) and turned about the normal; c grows to hold the stack and the vacuum.",
         {{"file", "Structure file to stack", "string", "", ""}, {"gap", "Gap (Angstrom)", "float", "3.35", ""},
          {"copies", "Copies", "int", "1", ""}, {"shift", "Lateral shift per copy (fractional u v)", "vector", "0 0 0", ""},
          {"twist", "Turn about the normal (degrees)", "float", "0", ""},
          {"strain", "Strain to a common lattice", "bool", "true", ""}, {"match", "Common lattice", "choice", "layer", "layer|average|substrate"},
          {"vacuum", "Vacuum (Angstrom)", "float", "15", ""},
          {"select", "Select the stacked atoms", "bool", "true", ""}}, stackLayers},
        {"select-fragment", "Select fragments", "Selection",
         "Selects whole bonded fragments (molecules, clusters): those touching the current selection, or those with a number of atoms in a range.",
         {{"rule", "Fragments", "choice", "touching", "touching|size"}, {"min_atoms", "Minimum atoms (size rule)", "int", "1", ""},
          {"max_atoms", "Maximum atoms (size rule; 0 = any)", "int", "0", ""}, {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""},
          kModeParameter}, selectFragment},
    };
}
}  // namespace atomforge::pipeline
