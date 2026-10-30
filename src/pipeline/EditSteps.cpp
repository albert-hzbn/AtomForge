// Geometric and chemical edits: mirror, inversion, rotation and alignment,
// sorting, hydrogen saturation, setting and measuring bond geometry,
// duplicating the selection, cell reduction and scaling, swapping elements,
// making molecules whole and selecting fragments.

#include "pipeline/ModifierSupport.h"
#include "pipeline/PipelineSteps.h"
#include "pipeline/StructureEdits.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <numeric>
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

void rotateStep(const Json& p, PipelineData& data)
{
    const char* id = "rotate";
    const bool selectedOnly = text(p, id, "target") == "selected";
    const auto atoms = targets(data, selectedOnly, id);
    const Mat3 r = rotationAbout(vector(p, id, "axis"), number(p, id, "angle") * kDegree);
    Structure& crystal = data.structure;
    const bool withCell = !selectedOnly && flag(p, id, "rotate_cell") && crystal.hasUnitCell;
    // A crystal turned with its cell turns about the cell origin: fractional coordinates stay.
    const Vec3 point = withCell ? Vec3{crystal.cellOffset[0], crystal.cellOffset[1], crystal.cellOffset[2]} : pivotOf(p, id, data, atoms);
    for (std::size_t i : atoms)
        place(data.structure.atoms[i], science::add(point, rotate(r, science::sub(positionOf(data.structure.atoms[i]), point))));
    if (withCell)
        for (auto& row : crystal.cellVectors) row = rotate(r, row);
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
    const std::string key = text(p, id, "key");
    const Structure& s = data.structure;
    std::vector<double> value(s.atoms.size());
    const Vec3 c = centroid(s, targets(data, false, id));
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const auto& a = s.atoms[i];
        if (key == "element") value[i] = a.atomicNumber;
        else if (key == "x") value[i] = a.x;
        else if (key == "y") value[i] = a.y;
        else if (key == "z") value[i] = a.z;
        else if (key == "distance") value[i] = science::norm(science::sub(positionOf(a), c));
        else if (key == "property") {
            if (s.atomProperty.size() != s.atoms.size()) throw std::runtime_error("sort-atoms by property needs a per-atom property");
            value[i] = std::isnan(s.atomProperty[i]) ? HUGE_VAL : s.atomProperty[i];
        } else
            throw std::runtime_error("Unknown sort key " + key);
    }
    std::vector<std::size_t> order(s.atoms.size());
    std::iota(order.begin(), order.end(), 0);
    const bool descending = flag(p, id, "descending");
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
        return descending ? value[x] > value[y] : value[x] < value[y];
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
    if (flag(p, id, "select")) std::fill(data.selected.begin(), data.selected.end(), 0);
    const std::size_t added = saturateWithHydrogen(data.structure, data.selected, elements, text(p, id, "hybridization") == "sp3",
                                                   number(p, id, "bond_scale"), flag(p, id, "select"));
    data.notes.push_back(std::to_string(added) + " hydrogen atoms added");
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
    std::vector<std::pair<std::size_t, Vec3>> moving;
    if (text(p, id, "move") == "fragment") moving = fragment(bondNeighbours(s, number(p, id, "bond_scale")), mover, fixed);
    else moving = {{mover, {0, 0, 0}}};
    if (atoms.size() == 4 && text(p, id, "move") == "atom") moving = {{atoms[3], science::sub(g.points[3], g.points[2])}};
    const Vec3 moverPoint = atoms.size() == 4 ? g.points[2] : g.points.back();
    for (const auto& [atom, offset] : moving)
        for (std::size_t k = 0; k < atoms.size(); ++k)
            if (atoms[k] == atom && atom != mover && !(atoms.size() == 4 && k == 3))
                throw std::runtime_error("set-geometry: the moving side contains atom " + std::to_string(atom) + " of the fixed side");
    std::vector<std::pair<std::size_t, Vec3>> result;  // atom -> new position
    if (atoms.size() == 2) {
        if (!(target > 0)) throw std::runtime_error("set-geometry: a distance must be positive");
        const Vec3 shiftBy = science::scale(unit(science::sub(g.points[1], g.points[0])), target - current);
        for (const auto& [atom, offset] : moving) result.push_back({atom, science::add(science::add(moverPoint, offset), shiftBy)});
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
            const Mat3 r = rotationAbout(axis, degrees * kDegree);
            std::vector<std::pair<std::size_t, Vec3>> out;
            for (const auto& [atom, offset] : moving)
                out.push_back({atom, science::add(pivot, rotate(r, science::sub(science::add(moverPoint, offset), pivot)))});
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
    double factor = value;
    if (mode == "volume") factor = std::cbrt(value / volume);
    else if (mode == "density") {
        if (s.atoms.empty()) throw std::runtime_error("scale-cell to a density needs atoms");
        factor = std::cbrt(mass * 1.66053906660 / value / volume);
    } else if (mode != "factor")
        throw std::runtime_error("scale-cell: mode is factor, volume or density");
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
         "Turns atoms by an angle (degrees, right-handed) about an axis through the centroid, the cell centre, the origin or a point. "
         "All atoms turn with the cell by default, about the cell origin, so a crystal keeps its periodicity and fractional coordinates.",
         {{"axis", "Axis", "vector", "0 0 1", ""}, {"angle", "Angle (degrees)", "float", "90", ""}, centre, point, target,
          {"rotate_cell", "Rotate the cell too (all atoms)", "bool", "true", ""}}, rotateStep},
        {"align-direction", "Align a direction with an axis", "Cell",
         "Turns the whole structure and its cell so that a direction (lattice [u v w] or Cartesian) lies along x, y or z, and optionally a second "
         "direction in the plane of that axis and a second one (as Atomsk -orient).",
         {{"direction", "Direction", "vector", "1 1 1", ""}, {"lattice", "Lattice direction [u v w]", "bool", "true", ""},
          {"axis", "Along", "choice", "z", "x|y|z"}, {"second_direction", "Second direction", "vector", "1 -1 0", ""},
          {"second_axis", "Second direction towards", "choice", "none", "none|x|y|z"}}, alignDirection},
        {"sort-atoms", "Sort atoms", "Modification", "Reorders the atoms by element (atomic number), a coordinate, distance from the centroid or the per-atom property.",
         {{"key", "Sort by", "choice", "element", "element|x|y|z|distance|property"}, {"descending", "Descending", "bool", "false", ""}}, sortAtoms},
        {"saturate-hydrogen", "Saturate with hydrogen", "Modification",
         "Adds hydrogen to under-coordinated atoms (bonds from covalent radii): tetrahedral (sp3) directions, or in-plane for planar carbon, boron "
         "and nitrogen with two bonds near 120 degrees; X-H lengths C 1.09, Si 1.48, N 1.01, O 0.97 Angstrom.",
         {{"elements", "Elements to saturate", "elements", "C Si Ge N O B P S", ""}, {"hybridization", "Geometry", "choice", "auto", "auto|sp3"},
          {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}, {"select", "Select the hydrogen", "bool", "true", ""}}, saturateHydrogen},
        {"swap-elements", "Swap elements", "Modification", "Exchanges two elements everywhere (A becomes B and B becomes A).",
         {{"first", "Element A", "element", "Ni", ""}, {"second", "Element B", "element", "Al", ""}}, swapElements},
        {"set-geometry", "Set distance, angle or dihedral", "Modification",
         "Sets the distance (2 atoms), angle (3) or dihedral (4; atom indices from 0) by moving the last atom, or the third and everything beyond the "
         "middle bond of a dihedral, with its bonded fragment (as in Avogadro).",
         {{"atoms", "Atom indices", "list", "0 1", ""}, {"value", "Value (Angstrom or degrees)", "float", "1.5", ""},
          {"move", "Move", "choice", "fragment", "fragment|atom"}, {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}}, setGeometry},
        {"measure-geometry", "Measure distance, angle or dihedral", "Analysis",
         "Reports the distance (2 atoms), angle (3) or dihedral angle (4) between atoms given by index (from 0), with minimum images in periodic cells.",
         {{"atoms", "Atom indices", "list", "0 1 2 3", ""}}, measureGeometry},
        {"duplicate-selected", "Duplicate selection", "Modification", "Adds copies of the selected atoms shifted by k x offset (k = 1 ... copies).",
         {{"offset", "Offset (Angstrom)", "vector", "0 0 3", ""}, {"copies", "Copies", "int", "1", ""}, {"select", "Select the copies", "bool", "true", ""}}, duplicateSelected},
        {"niggli-reduce", "Niggli-reduced cell", "Cell",
         "Replaces the cell by its Niggli-reduced form (shortest, most orthogonal basis of the same lattice; Krivy-Gruber) and wraps the atoms into it.",
         {}, niggliReduce},
        {"scale-cell", "Scale cell", "Cell", "Scales the cell and atom positions uniformly by a factor, or to a volume (Angstrom^3) or a density (g/cm3).",
         {{"mode", "Scale to", "choice", "factor", "factor|volume|density"}, {"value", "Value", "float", "1.01", ""}}, scaleCell},
        {"unwrap-molecules", "Make molecules whole", "Modification",
         "Moves atoms across periodic boundaries so that bonded atoms sit next to each other (molecules split by the cell become whole).",
         {{"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""}}, unwrapMolecules},
        {"select-fragment", "Select fragments", "Selection",
         "Selects whole bonded fragments (molecules, clusters): those touching the current selection, or those with a number of atoms in a range.",
         {{"rule", "Fragments", "choice", "touching", "touching|size"}, {"min_atoms", "Minimum atoms (size rule)", "int", "1", ""},
          {"max_atoms", "Maximum atoms (size rule; 0 = any)", "int", "0", ""}, {"bond_scale", "Bond cutoff (x sum of covalent radii)", "float", "1.2", ""},
          kModeParameter}, selectFragment},
    };
}
}  // namespace atomforge::pipeline
