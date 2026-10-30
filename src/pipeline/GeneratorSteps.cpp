// Builders that make a structure from parameters (nanotubes, 2D layers,
// ribbons, twisted bilayers, clusters, molecules, packings) or add to one
// (adsorbates, antisites, Frenkel pairs), and structure-wide constructions
// (near-cubic supercells, interpolation between two structures).

#include "pipeline/ModifierSupport.h"
#include "pipeline/MoleculeLibrary.h"
#include "pipeline/PipelineSteps.h"
#include "pipeline/StructureEdits.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>

namespace atomforge::pipeline
{
namespace
{
using namespace support;
using science::Mat3;
using science::Vec3;
constexpr double kPi = 3.14159265358979323846;

// ---------------------------------------------------------------- helpers
std::vector<std::string> words(const std::string& text)
{
    std::istringstream in(text);
    std::vector<std::string> out;
    for (std::string w; in >> w;) out.push_back(w);
    return out;
}

// Element symbols from a list ("C", "B N", "Mo S").
std::vector<std::string> elementList(const std::string& text, std::size_t count, const char* what)
{
    std::vector<std::string> out;
    for (const std::string& w : words(text)) out.push_back(standardElement(w));
    if (out.size() == 1 && count > 1) out.resize(count, out.front());
    if (out.size() != count)
        throw std::runtime_error(std::string(what) + ": give " + std::to_string(count) + " element symbol" + (count > 1 ? "s" : ""));
    return out;
}

int integerOf(const Json& p, const char* id, const char* name, int minimum)
{
    const double value = number(p, id, name);
    if (value != std::floor(value) || value < minimum)
        throw std::runtime_error(std::string(id) + ": " + name + " must be an integer >= " + std::to_string(minimum));
    return static_cast<int>(value);
}

// A finished new structure replaces the pipeline's.
void produce(PipelineData& data, Structure s)
{
    data.structure = std::move(s);
    data.clearSelection();
}

void addAtom(Structure& s, const std::string& symbol, const Vec3& position)
{
    std::vector<char> ignored(s.atoms.size(), 0);
    appendAtom(s, ignored, symbol, position, false);
}

void setCell(Structure& s, const Vec3& a, const Vec3& b, const Vec3& c, std::array<bool, 3> pbc)
{
    s.hasUnitCell = true;
    s.cellVectors = {a, b, c};
    s.cellOffset = {0, 0, 0};
    s.pbc = pbc;
}

// Moves every atom by d.
void shift(Structure& s, const Vec3& d)
{
    for (auto& a : s.atoms) place(a, science::add(positionOf(a), d));
}

Vec3 rotateZ(const Vec3& v, double angle)
{
    const double c = std::cos(angle), s = std::sin(angle);
    return {c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]};
}

Vec3 apply(const Mat3& r, const Vec3& v) { return rotate(r, v); }

// ---------------------------------------------------------------- 2D layers
// A hexagonal monolayer: cell a1 = (a, 0, 0), a2 = (-a/2, a sqrt3/2, 0) and
// the basis in fractional (in-plane) and Cartesian z coordinates.
struct Layer
{
    double a = 0;
    struct Site { std::string symbol; double u, v, z; int sublattice; };
    std::vector<Site> sites;
    double thickness = 0;
    double interlayer = 3.35;  // stacking distance in bilayers
};

Layer monolayer(const std::string& kind, const std::string& elementsText, double a, double thickness)
{
    Layer layer;
    if (kind == "graphene") {
        const auto e = elementList(elementsText.empty() ? "C" : elementsText, 1, "graphene");
        layer.a = a > 0 ? a : 2.4612;  // C-C 1.421 Angstrom
        layer.sites = {{e[0], 1.0 / 3, 2.0 / 3, 0, 0}, {e[0], 2.0 / 3, 1.0 / 3, 0, 1}};
        layer.interlayer = 3.35;
    } else if (kind == "hbn") {
        const auto e = elementList(elementsText.empty() ? "B N" : elementsText, 2, "h-BN");
        layer.a = a > 0 ? a : 2.504;
        layer.sites = {{e[0], 1.0 / 3, 2.0 / 3, 0, 0}, {e[1], 2.0 / 3, 1.0 / 3, 0, 1}};
        layer.interlayer = 3.33;
    } else if (kind == "mx2-2h" || kind == "mx2-1t") {
        const auto e = elementList(elementsText.empty() ? "Mo S" : elementsText, 2, "MX2");
        layer.a = a > 0 ? a : 3.18;
        layer.thickness = thickness > 0 ? thickness : 3.19;  // X-X distance across the layer (ASE mx2 defaults)
        const double h = layer.thickness / 2;
        layer.sites = {{e[0], 0, 0, 0, 0}, {e[1], 1.0 / 3, 2.0 / 3, h, 1}};
        if (kind == "mx2-2h") layer.sites.push_back({e[1], 1.0 / 3, 2.0 / 3, -h, 1});
        else layer.sites.push_back({e[1], 2.0 / 3, 1.0 / 3, -h, 1});
        layer.interlayer = 6.15;
    } else
        throw std::runtime_error("Unknown layer '" + kind + "': graphene, hbn, mx2-2h or mx2-1t");
    return layer;
}

Vec3 layerA1(const Layer& l) { return {l.a, 0, 0}; }
Vec3 layerA2(const Layer& l) { return {-l.a / 2, l.a * std::sqrt(3.0) / 2, 0}; }

void buildMonolayer(const Json& p, PipelineData& data)
{
    const char* id = "build-monolayer";
    const Layer layer = monolayer(text(p, id, "kind"), text(p, id, "elements"), number(p, id, "a"), number(p, id, "thickness"));
    const int na = integerOf(p, id, "na", 1), nb = integerOf(p, id, "nb", 1);
    const double vacuum = number(p, id, "vacuum");
    if (vacuum < 0) throw std::runtime_error("vacuum must be nonnegative");
    const Vec3 a1 = layerA1(layer), a2 = layerA2(layer);
    const double height = layer.thickness + vacuum;
    Structure s;
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (const auto& site : layer.sites) {
                const Vec3 xy = science::add(science::scale(a1, i + site.u), science::scale(a2, j + site.v));
                addAtom(s, site.symbol, {xy[0], xy[1], height / 2 + site.z});
            }
    setCell(s, science::scale(a1, na), science::scale(a2, nb), {0, 0, std::max(height, 1.0)}, {true, true, true});
    produce(data, std::move(s));
    data.notes.push_back(std::to_string(data.structure.atoms.size()) + " atoms, a = " + std::to_string(layer.a) + " Angstrom");
}

// ---------------------------------------------------------------- nanotubes
// Rolls the honeycomb lattice (bond d) along the chiral vector n a1 + m a2.
void addTubeWall(Structure& s, int n, int m, double bond, const std::vector<std::string>& species, int cells, double& radius, double& period)
{
    const double a = bond * std::sqrt(3.0);
    const Vec3 a1{a * std::sqrt(3.0) / 2, a / 2, 0}, a2{a * std::sqrt(3.0) / 2, -a / 2, 0};
    const Vec3 c = science::add(science::scale(a1, n), science::scale(a2, m));
    const int dR = std::gcd(2 * m + n, 2 * n + m);
    const int t1 = (2 * m + n) / dR, t2 = -(2 * n + m) / dR;
    const Vec3 t = science::add(science::scale(a1, t1), science::scale(a2, t2));
    const double lc = science::norm(c), lt = science::norm(t);
    radius = lc / (2 * kPi);
    period = lt;
    const Vec3 basis[2] = {{0, 0, 0}, science::scale(science::add(a1, a2), 1.0 / 3.0)};
    const int range = std::abs(n) + std::abs(m) + std::abs(t1) + std::abs(t2) + 2;
    const double eps = 1e-9;
    std::size_t found = 0;
    for (int i = -range; i <= range; ++i)
        for (int j = -range; j <= range; ++j)
            for (int b = 0; b < 2; ++b) {
                const Vec3 q = science::add(science::add(science::scale(a1, i), science::scale(a2, j)), basis[b]);
                const double u = science::dot(q, c) / (lc * lc), v = science::dot(q, t) / (lt * lt);
                if (u < -eps || u >= 1 - eps || v < -eps || v >= 1 - eps) continue;
                ++found;
                for (int k = 0; k < cells; ++k) {
                    const double theta = 2 * kPi * u;
                    addAtom(s, species[static_cast<std::size_t>(b)], {radius * std::cos(theta), radius * std::sin(theta), (v + k) * lt});
                }
            }
    const std::size_t expected = static_cast<std::size_t>(4 * (n * n + n * m + m * m) / dR);
    if (found != expected)
        throw std::runtime_error("Nanotube unit cell has " + std::to_string(found) + " atoms, expected " + std::to_string(expected));
}

void buildNanotube(const Json& p, PipelineData& data)
{
    const char* id = "build-nanotube";
    const int n = integerOf(p, id, "n", 0), m = integerOf(p, id, "m", 0);
    if (n == 0 && m == 0) throw std::runtime_error("The chiral indices (n, m) must not both be zero");
    if (m > n) throw std::runtime_error("Use n >= m (the tube (m, n) is the mirror image of (n, m))");
    const int cells = integerOf(p, id, "length", 1), walls = integerOf(p, id, "walls", 1);
    const double bond = number(p, id, "bond"), spacing = number(p, id, "wall_spacing"), vacuum = number(p, id, "vacuum");
    if (!(bond > 0) || !(spacing > 0) || vacuum < 0) throw std::runtime_error("bond and wall spacing must be positive, vacuum nonnegative");
    const auto species = elementList(text(p, id, "elements"), 2, "nanotube");
    if (walls > 1 && !(n == m || m == 0))
        throw std::runtime_error("Multi-wall tubes need armchair (n, n) or zigzag (n, 0) indices, which share the axial period");
    Structure s;
    double radius = 0, period = 0, outer = 0;
    std::ostringstream notes;
    for (int w = 0; w < walls; ++w) {
        // Each wall grows the radius by about the wall spacing.
        const double a = bond * std::sqrt(3.0);
        const int step = n == m ? static_cast<int>(std::lround(spacing * 2 * kPi / (a * std::sqrt(3.0))))
                                : static_cast<int>(std::lround(spacing * 2 * kPi / a));
        const int nw = n + w * step, mw = m == 0 ? 0 : m + w * step;
        addTubeWall(s, nw, mw, bond, species, cells, radius, period);
        outer = radius;
        notes << (w ? ", " : "") << "(" << nw << "," << mw << ") r = " << radius;
    }
    const double box = 2 * outer + vacuum + 2 * bond;
    shift(s, {box / 2, box / 2, 0});
    setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, period * cells}, {false, false, true});
    produce(data, std::move(s));
    data.notes.push_back(std::to_string(data.structure.atoms.size()) + " atoms; walls " + notes.str() + " Angstrom; period " + std::to_string(period));
}

// ---------------------------------------------------------------- nanoribbons
void buildNanoribbon(const Json& p, PipelineData& data)
{
    const char* id = "build-nanoribbon";
    const std::string edge = text(p, id, "edge");
    const int width = integerOf(p, id, "width", 1), length = integerOf(p, id, "length", 1);
    const double d = number(p, id, "bond"), vacuum = number(p, id, "vacuum");
    if (!(d > 0) || vacuum < 0) throw std::runtime_error("bond must be positive, vacuum nonnegative");
    const auto species = elementList(text(p, id, "elements"), 2, "nanoribbon");
    // Rectangular honeycomb cell, zigzag direction along x: period sqrt3 d by 3d.
    const double lx = std::sqrt(3.0) * d, ly = 3 * d;
    const Vec3 cellSites[4] = {{0, 0, 0}, {lx / 2, d / 2, 0}, {lx / 2, 1.5 * d, 0}, {0, 2 * d, 0}};
    const int sublattice[4] = {0, 1, 0, 1};
    const bool zigzag = edge == "zigzag";
    if (!zigzag && edge != "armchair") throw std::runtime_error("edge must be armchair or zigzag");
    // Levels across the ribbon: y for zigzag ribbons (periodic along x), x for armchair ones (periodic along y).
    struct Site { double along, across; int sublattice; };
    std::vector<Site> sites;
    const int repeats = width + 2;
    for (int r = 0; r < repeats; ++r)
        for (int k = 0; k < 4; ++k) {
            const Vec3 q = science::add(cellSites[k], zigzag ? Vec3{0, r * ly, 0} : Vec3{r * lx, 0, 0});
            sites.push_back(zigzag ? Site{q[0], q[1], sublattice[k]} : Site{q[1], q[0], sublattice[k]});
        }
    std::vector<double> levels;
    for (const Site& site : sites) {
        bool known = false;
        for (double l : levels) known = known || std::abs(l - site.across) < 1e-6;
        if (!known) levels.push_back(site.across);
    }
    std::sort(levels.begin(), levels.end());
    // A zigzag ribbon of N chains has 2N atom rows; an armchair ribbon N dimer lines.
    const std::size_t rows = zigzag ? 2 * static_cast<std::size_t>(width) : static_cast<std::size_t>(width);
    if (levels.size() < rows) throw std::runtime_error("ribbon too wide");
    const double top = levels[rows - 1];
    const double period = zigzag ? lx : ly;
    const double ribbonWidth = top - levels.front();
    Structure s;
    for (int k = 0; k < length; ++k)
        for (const Site& site : sites)
            if (site.across <= top + 1e-6)
                addAtom(s, species[static_cast<std::size_t>(site.sublattice)], {site.along + k * period, site.across, 0});
    const double across = ribbonWidth + vacuum + 2.2;  // room for edge hydrogen
    const double height = std::max(vacuum, 1.0);
    shift(s, {0, (across - ribbonWidth) / 2, height / 2});
    setCell(s, {period * length, 0, 0}, {0, across, 0}, {0, 0, height}, {true, false, false});
    std::vector<char> selected(s.atoms.size(), 0);
    std::size_t hydrogen = 0;
    if (flag(p, id, "passivate")) hydrogen = saturateWithHydrogen(s, selected, {species[0], species[1]}, false, 1.2, false);
    produce(data, std::move(s));
    data.notes.push_back(std::to_string(width) + (zigzag ? " chains, zigzag" : " dimer lines, armchair") + " edges; " +
                         std::to_string(data.structure.atoms.size()) + " atoms (" + std::to_string(hydrogen) + " H)");
}

// ---------------------------------------------------------------- twisted bilayers
// Fills the cell (L1, L2, c) with the in-plane lattice of a layer turned by `angle`.
void fillLayer(Structure& s, const Layer& layer, double angle, const Vec3& l1, const Vec3& l2, double z, std::size_t& count)
{
    const Vec3 a1 = rotateZ(layerA1(layer), angle), a2 = rotateZ(layerA2(layer), angle);
    Mat3 super{l1, l2, Vec3{0, 0, 1}};
    const Mat3 inverse = science::inverse(super);
    const double reach = (science::norm(l1) + science::norm(l2)) / layer.a + 3;
    const int range = static_cast<int>(std::ceil(reach)) + 2;
    const double eps = 1e-7;
    for (int i = -2 * range; i <= 2 * range; ++i)
        for (int j = -2 * range; j <= 2 * range; ++j)
            for (const auto& site : layer.sites) {
                const Vec3 q = science::add(science::scale(a1, i + site.u), science::scale(a2, j + site.v));
                const Vec3 f = science::rowTimes({q[0], q[1], 0}, inverse);
                if (f[0] < -eps || f[0] >= 1 - eps || f[1] < -eps || f[1] >= 1 - eps) continue;
                addAtom(s, site.symbol, {q[0], q[1], z + site.z});
                ++count;
            }
}

void buildTwistedBilayer(const Json& p, PipelineData& data)
{
    const char* id = "build-twisted-bilayer";
    const Layer layer = monolayer(text(p, id, "kind"), text(p, id, "elements"), number(p, id, "a"), number(p, id, "thickness"));
    // Commensurate angles: the layer vectors (i, i+1) and (i+1, i) of the 60-degree basis coincide.
    const auto angleOf = [](int i) {
        const double c = (3.0 * i * i + 3.0 * i + 0.5) / (3.0 * i * i + 3.0 * i + 1.0);
        return std::acos(std::min(1.0, c));
    };
    int index = integerOf(p, id, "index", 0);
    const double wanted = number(p, id, "angle");
    if (wanted > 0) {
        index = 1;
        for (int i = 1; i <= 300; ++i)
            if (std::abs(angleOf(i) * 180 / kPi - wanted) < std::abs(angleOf(index) * 180 / kPi - wanted)) index = i;
    }
    if (index < 1) throw std::runtime_error("Give a twist angle or a commensurate index >= 1");
    if (index > 300) throw std::runtime_error("index at most 300 (the cell grows as 3 i^2)");
    const double theta = angleOf(index);
    const Vec3 b1 = layerA1(layer), b2 = science::add(layerA1(layer), layerA2(layer));  // 60-degree basis
    const Vec3 l1 = science::add(science::scale(b1, index), science::scale(b2, index + 1));
    const Vec3 l2 = rotateZ(l1, kPi / 3);
    const double interlayer = number(p, id, "interlayer") > 0 ? number(p, id, "interlayer") : layer.interlayer;
    const double vacuum = number(p, id, "vacuum");
    const double height = interlayer + layer.thickness + vacuum;
    Structure s;
    std::size_t bottom = 0, top = 0;
    const double z0 = (height - interlayer) / 2;
    fillLayer(s, layer, 0, l1, l2, z0, bottom);
    // The second layer turned so that its vector (i+1, i) lies on L1.
    const Vec3 v2 = science::add(science::scale(b1, index + 1), science::scale(b2, index));
    const double turn = std::atan2(l1[1], l1[0]) - std::atan2(v2[1], v2[0]);
    fillLayer(s, layer, turn, l1, l2, z0 + interlayer, top);
    const std::size_t expected = layer.sites.size() * static_cast<std::size_t>(3 * index * index + 3 * index + 1);
    if (bottom != expected || top != expected)
        throw std::runtime_error("Twisted bilayer layers have " + std::to_string(bottom) + " and " + std::to_string(top) +
                                 " atoms, expected " + std::to_string(expected));
    setCell(s, l1, l2, {0, 0, height}, {true, true, true});
    produce(data, std::move(s));
    std::ostringstream note;
    note << "Twist " << theta * 180 / kPi << " degrees (index " << index << "), " << data.structure.atoms.size() << " atoms, |L| = "
         << science::norm(l1) << " Angstrom";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- clusters
void buildCluster(const Json& p, PipelineData& data)
{
    const char* id = "build-cluster";
    const std::string kind = text(p, id, "kind");
    const std::string symbol = standardElement(text(p, id, "element"));
    const int shells = integerOf(p, id, "shells", 1);
    const double d = number(p, id, "distance"), vacuum = number(p, id, "vacuum");
    if (!(d > 0) || vacuum < 0) throw std::runtime_error("distance must be positive, vacuum nonnegative");
    std::vector<Vec3> points{{0, 0, 0}};
    std::map<std::array<long long, 3>, bool> seen{{{0, 0, 0}, true}};
    const auto addUnique = [&](const Vec3& q) {
        const std::array<long long, 3> key{std::llround(q[0] / (1e-4 * d)), std::llround(q[1] / (1e-4 * d)), std::llround(q[2] / (1e-4 * d))};
        if (seen.emplace(key, true).second) points.push_back(q);
    };
    if (kind == "icosahedron") {
        // Mackay icosahedron (as ase.cluster.Icosahedron): shell k has its
        // vertices k d from the centre, edges and faces divided into k segments
        // (surface neighbours 5% further apart than radial ones).
        const double phi = (1 + std::sqrt(5.0)) / 2;
        std::vector<Vec3> v;
        for (int s1 = -1; s1 <= 1; s1 += 2)
            for (int s2 = -1; s2 <= 1; s2 += 2) {
                v.push_back({0, 1.0 * s1, phi * s2});
                v.push_back({1.0 * s1, phi * s2, 0});
                v.push_back({phi * s2, 0, 1.0 * s1});
            }
        std::vector<std::array<int, 3>> faces;
        for (int i = 0; i < 12; ++i)
            for (int j = i + 1; j < 12; ++j)
                for (int k = j + 1; k < 12; ++k) {
                    const auto edge = [&](int x, int y) { return std::abs(science::norm(science::sub(v[x], v[y])) - 2) < 1e-6; };
                    if (edge(i, j) && edge(j, k) && edge(i, k)) faces.push_back({i, j, k});
                }
        for (int shell = 1; shell < shells; ++shell)
            for (const auto& f : faces) {
                const double toVertex = shell * d / science::norm(v[0]);
                const Vec3 a = science::scale(v[f[0]], toVertex), b = science::scale(v[f[1]], toVertex), c = science::scale(v[f[2]], toVertex);
                for (int i = 0; i <= shell; ++i)
                    for (int j = 0; i + j <= shell; ++j)
                        addUnique(science::add(a, science::add(science::scale(science::sub(b, a), double(i) / shell),
                                                               science::scale(science::sub(c, a), double(j) / shell))));
            }
    } else if (kind == "decahedron") {
        // Regular (Ino-free) decahedron: a pentagonal bipyramid of five
        // tetrahedra sharing the five-fold axis, with `shells` atoms along each
        // edge (2: 7 atoms, 3: 23, 4: 54, ...); pentagon and slanted edges are
        // (shells - 1) d, the axis 5% longer.
        if (shells < 2) throw std::runtime_error("A decahedron needs at least 2 atoms along an edge (shells >= 2)");
        const int n = shells - 1;
        const double edge = n * d;
        const double r = edge / (2 * std::sin(kPi / 5));
        const double h = std::sqrt(std::max(0.0, edge * edge - r * r));
        const Vec3 top{0, 0, h}, bottom{0, 0, -h};
        points.clear();
        seen.clear();
        for (int k = 0; k < 5; ++k) {
            const Vec3 p1{r * std::cos(2 * kPi * k / 5), r * std::sin(2 * kPi * k / 5), 0};
            const Vec3 p2{r * std::cos(2 * kPi * (k + 1) / 5), r * std::sin(2 * kPi * (k + 1) / 5), 0};
            for (int i = 0; i <= n; ++i)
                for (int j = 0; i + j <= n; ++j)
                    for (int l = 0; i + j + l <= n; ++l) {
                        const int m = n - i - j - l;
                        addUnique(science::scale(science::add(science::add(science::scale(top, i), science::scale(bottom, j)),
                                                              science::add(science::scale(p1, l), science::scale(p2, m))), 1.0 / n));
                    }
        }
    } else if (kind == "cuboctahedron") {
        // fcc points (nearest-neighbour vectors (1,1,0) d/sqrt2) inside the
        // cuboctahedron with vertices k (1,1,0): |x|, |y|, |z| <= k and |x|+|y|+|z| <= 2k.
        const double unit = d / std::sqrt(2.0);
        const int k = shells - 1;
        points.clear();
        for (int x = -k; x <= k; ++x)
            for (int y = -k; y <= k; ++y)
                for (int z = -k; z <= k; ++z) {
                    if (((x + y + z) % 2 + 2) % 2 != 0) continue;
                    if (std::abs(x) + std::abs(y) + std::abs(z) > 2 * k) continue;
                    points.push_back({x * unit, y * unit, z * unit});
                }
    } else
        throw std::runtime_error("Unknown cluster '" + kind + "': icosahedron, decahedron or cuboctahedron");
    double extent = 0;
    for (const Vec3& q : points) extent = std::max(extent, science::norm(q));
    const double box = 2 * extent + vacuum;
    Structure s;
    for (const Vec3& q : points) addAtom(s, symbol, science::add(q, {box / 2, box / 2, box / 2}));
    setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, box}, {false, false, false});
    produce(data, std::move(s));
    data.notes.push_back(kind + " of " + std::to_string(shells) + " shells: " + std::to_string(points.size()) + " atoms");
}

// ---------------------------------------------------------------- molecules
void buildMolecule(const Json& p, PipelineData& data)
{
    const char* id = "build-molecule";
    const auto atoms = molecule(text(p, id, "name"));
    const double vacuum = number(p, id, "vacuum");
    Structure s;
    double extent = 0;
    for (const auto& a : atoms) extent = std::max(extent, science::norm(a.position));
    const double box = 2 * extent + vacuum;
    for (const auto& a : atoms) addAtom(s, a.symbol, vacuum > 0 ? science::add(a.position, {box / 2, box / 2, box / 2}) : a.position);
    if (vacuum > 0) setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, box}, {false, false, false});
    produce(data, std::move(s));
    data.notes.push_back(std::to_string(atoms.size()) + " atoms, " + std::to_string(moleculeMass(atoms)) + " u");
}

// Uniformly random rotation (from a random unit quaternion).
Mat3 randomRotation(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double u1 = u(rng), u2 = u(rng) * 2 * kPi, u3 = u(rng) * 2 * kPi;
    const double w = std::sqrt(1 - u1) * std::sin(u2), x = std::sqrt(1 - u1) * std::cos(u2);
    const double y = std::sqrt(u1) * std::sin(u3), z = std::sqrt(u1) * std::cos(u3);
    return {Vec3{1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)},
            Vec3{2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)},
            Vec3{2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)}};
}

void buildPacking(const Json& p, PipelineData& data)
{
    const char* id = "build-packing";
    // "H2O 100 CO2 10": molecules and counts.
    const auto list = words(text(p, id, "molecules"));
    if (list.empty() || list.size() % 2) throw std::runtime_error("molecules: give pairs of a molecule and a count, e.g. H2O 100");
    struct Kind { std::vector<MoleculeAtom> atoms; int count; std::string name; };
    std::vector<Kind> kinds;
    double mass = 0;
    for (std::size_t k = 0; k < list.size(); k += 2) {
        char* end = nullptr;
        const long count = std::strtol(list[k + 1].c_str(), &end, 10);
        if (*end || count < 0) throw std::runtime_error("molecules: '" + list[k + 1] + "' is not a count");
        kinds.push_back({molecule(list[k]), static_cast<int>(count), list[k]});
        mass += moleculeMass(kinds.back().atoms) * count;
    }
    Vec3 box = vector(p, id, "box");
    const double density = number(p, id, "density");
    if (box[0] <= 0 || box[1] <= 0 || box[2] <= 0) {
        if (!(density > 0)) throw std::runtime_error("Give a density or the box lengths");
        const double side = std::cbrt(mass * 1.66053906660 / density);  // u -> g, cm^3 -> A^3
        box = {side, side, side};
    }
    const double minimum = number(p, id, "min_distance");
    const int attempts = integerOf(p, id, "attempts", 1);
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    std::uniform_real_distribution<double> u(0.0, 1.0);
    // Placed atoms in a periodic grid of cells at least `minimum` wide.
    const int g[3] = {std::max(1, static_cast<int>(box[0] / std::max(minimum, 0.5))), std::max(1, static_cast<int>(box[1] / std::max(minimum, 0.5))),
                      std::max(1, static_cast<int>(box[2] / std::max(minimum, 0.5)))};
    std::unordered_map<long long, std::vector<Vec3>> grid;
    const auto key = [&](int x, int y, int z) {
        x = ((x % g[0]) + g[0]) % g[0]; y = ((y % g[1]) + g[1]) % g[1]; z = ((z % g[2]) + g[2]) % g[2];
        return (static_cast<long long>(x) * g[1] + y) * g[2] + z;
    };
    const auto wrapped = [&](Vec3 q) {
        for (int k = 0; k < 3; ++k) q[k] -= box[k] * std::floor(q[k] / box[k]);
        return q;
    };
    const auto cellOfPoint = [&](const Vec3& q, int k) { return std::min(g[k] - 1, static_cast<int>(q[k] / box[k] * g[k])); };
    const auto clear = [&](const Vec3& q) {
        const int cx = cellOfPoint(q, 0), cy = cellOfPoint(q, 1), cz = cellOfPoint(q, 2);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == grid.end()) continue;
                    for (const Vec3& r : it->second) {
                        Vec3 v = science::sub(q, r);
                        for (int k = 0; k < 3; ++k) v[k] -= box[k] * std::round(v[k] / box[k]);
                        if (science::norm(v) < minimum) return false;
                    }
                }
        return true;
    };
    // Larger molecules first: they are harder to place.
    std::stable_sort(kinds.begin(), kinds.end(), [](const Kind& a, const Kind& b) { return a.atoms.size() > b.atoms.size(); });
    Structure s;
    std::ostringstream summary;
    for (const Kind& kind : kinds) {
        for (int c = 0; c < kind.count; ++c) {
            bool placed = false;
            for (int attempt = 0; attempt < attempts && !placed; ++attempt) {
                const Mat3 r = randomRotation(rng);
                const Vec3 centre{u(rng) * box[0], u(rng) * box[1], u(rng) * box[2]};
                std::vector<Vec3> trial;
                bool ok = true;
                for (const auto& a : kind.atoms) {
                    const Vec3 q = wrapped(science::add(centre, apply(r, a.position)));
                    if (!clear(q)) { ok = false; break; }
                    trial.push_back(q);
                }
                if (!ok) continue;
                for (std::size_t k = 0; k < trial.size(); ++k) {
                    // Molecules stay whole: atoms are written unwrapped around their centre.
                    addAtom(s, kind.atoms[k].symbol, science::add(centre, apply(r, kind.atoms[k].position)));
                    grid[key(cellOfPoint(trial[k], 0), cellOfPoint(trial[k], 1), cellOfPoint(trial[k], 2))].push_back(trial[k]);
                }
                placed = true;
            }
            if (!placed)
                throw std::runtime_error("Placed " + std::to_string(c) + " of " + std::to_string(kind.count) + " " + kind.name +
                                         " molecules: lower the density or min_distance, or raise attempts");
        }
        summary << (summary.tellp() ? ", " : "") << kind.count << " " << kind.name;
    }
    setCell(s, {box[0], 0, 0}, {0, box[1], 0}, {0, 0, box[2]}, {true, true, true});
    produce(data, std::move(s));
    std::ostringstream note;
    note << summary.str() << " in " << box[0] << " x " << box[1] << " x " << box[2] << " Angstrom ("
         << mass * 1.66053906660 / (box[0] * box[1] * box[2]) << " g/cm3)";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- adsorbates
void addAdsorbate(const Json& p, PipelineData& data)
{
    const char* id = "add-adsorbate";
    requireCell(data, "add-adsorbate");
    Structure& s = data.structure;
    if (s.atoms.empty()) throw std::runtime_error("add-adsorbate needs a slab");
    const auto species = molecule(text(p, id, "species"));
    const std::string site = text(p, id, "site");
    const double height = number(p, id, "height");
    const Vec3 target = vector(p, id, "position");  // fractional (u, v) in the surface cell
    const Mat3 cell = cellOf(s);
    Vec3 normal = science::cross(cell[0], cell[1]);
    normal = science::scale(normal, 1 / science::norm(normal));
    if (science::dot(normal, cell[2]) < 0) normal = science::scale(normal, -1);
    // Layers along the normal.
    double top = -HUGE_VAL;
    for (const auto& a : s.atoms) top = std::max(top, science::dot(positionOf(a), normal));
    const double layerTolerance = 0.6;
    std::vector<std::size_t> surface, second;
    double secondHeight = -HUGE_VAL;
    for (const auto& a : s.atoms) {
        const double h = science::dot(positionOf(a), normal);
        if (h < top - layerTolerance) secondHeight = std::max(secondHeight, h);
    }
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const double h = science::dot(positionOf(s.atoms[i]), normal);
        if (h >= top - layerTolerance) surface.push_back(i);
        else if (h >= secondHeight - layerTolerance) second.push_back(i);
    }
    // In-plane displacements by the minimum image along a and b.
    const Mat3 inverse = science::inverse(cell);
    const auto inPlane = [&](Vec3 v) {
        v = science::sub(v, science::scale(normal, science::dot(v, normal)));
        Vec3 f = science::rowTimes(v, inverse);
        for (int k = 0; k < 2; ++k) f[k] -= std::round(f[k]);
        f[2] = 0;
        return science::rowTimes(f, cell);
    };
    // Nearest surface neighbour distance.
    double nearest = HUGE_VAL;
    for (std::size_t x = 0; x < surface.size(); ++x)
        for (std::size_t y = 0; y < surface.size(); ++y)
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j) {
                    Vec3 v = science::sub(positionOf(s.atoms[surface[y]]), positionOf(s.atoms[surface[x]]));
                    v = science::add(v, science::add(science::scale(cell[0], i), science::scale(cell[1], j)));
                    v = science::sub(v, science::scale(normal, science::dot(v, normal)));
                    if (science::norm(v) > 1e-6) nearest = std::min(nearest, science::norm(v));
                }
    if (site != "top" && !(nearest < HUGE_VAL)) throw std::runtime_error("A " + site + " site needs a periodic surface");
    std::vector<Vec3> candidates;  // in-plane site positions
    for (std::size_t x = 0; x < surface.size(); ++x) {
        const Vec3 px = positionOf(s.atoms[surface[x]]);
        if (site == "top") { candidates.push_back(px); continue; }
        std::vector<Vec3> near;  // vectors to neighbouring surface atoms (all images in a 3 x 3 block)
        for (std::size_t y = 0; y < surface.size(); ++y)
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j) {
                    Vec3 v = science::sub(positionOf(s.atoms[surface[y]]), px);
                    v = science::add(v, science::add(science::scale(cell[0], i), science::scale(cell[1], j)));
                    v = science::sub(v, science::scale(normal, science::dot(v, normal)));
                    const double len = science::norm(v);
                    if (len > 1e-6 && len < 1.2 * nearest) near.push_back(v);
                }
        if (site == "bridge") {
            for (const Vec3& v : near) candidates.push_back(science::add(px, science::scale(v, 0.5)));
            continue;
        }
        for (std::size_t m = 0; m < near.size(); ++m)
            for (std::size_t q = m + 1; q < near.size(); ++q) {
                if (science::norm(science::sub(near[m], near[q])) > 1.2 * nearest) continue;
                const Vec3 centroid = science::add(px, science::scale(science::add(near[m], near[q]), 1.0 / 3.0));
                // hcp hollows sit above a second-layer atom, fcc hollows do not.
                bool above = false;
                for (std::size_t k : second)
                    above = above || science::norm(inPlane(science::sub(positionOf(s.atoms[k]), centroid))) < 0.35 * nearest;
                if ((site == "hcp") == above) candidates.push_back(centroid);
            }
    }
    if (candidates.empty()) throw std::runtime_error("No " + site + " site found on this surface");
    const Vec3 wanted = science::rowTimes({target[0], target[1], 0}, cell);
    std::size_t best = 0;
    double bestDistance = HUGE_VAL;
    for (std::size_t k = 0; k < candidates.size(); ++k) {
        const double dist = science::norm(inPlane(science::sub(candidates[k], wanted)));
        if (dist < bestDistance - 1e-9) { bestDistance = dist; best = k; }
    }
    Vec3 anchor = candidates[best];
    anchor = science::add(science::sub(anchor, science::scale(normal, science::dot(anchor, normal))), science::scale(normal, top + height));
    // The molecule's z axis along the surface normal; its lowest atom at the height.
    const Mat3 r = rotationBetween({0, 0, 1}, normal);
    double lowest = HUGE_VAL;
    for (const auto& a : species) lowest = std::min(lowest, science::dot(apply(r, a.position), normal));
    if (data.selected.size() != s.atoms.size()) data.selected.assign(s.atoms.size(), 0);
    std::fill(data.selected.begin(), data.selected.end(), 0);
    for (const auto& a : species)
        appendAtom(s, data.selected, a.symbol, science::add(anchor, science::sub(apply(r, a.position), science::scale(normal, lowest))),
                   flag(p, id, "select"));
    std::ostringstream note;
    note << text(p, id, "species") << " on a " << site << " site, " << height << " Angstrom above the surface ("
         << candidates.size() << " such sites)";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- supercells
void nearCubicSupercell(const Json& p, PipelineData& data)
{
    const char* id = "near-cubic-supercell";
    requireCell(data, id);
    const Structure& s = data.structure;
    if (s.atoms.empty()) throw std::runtime_error("near-cubic-supercell needs atoms");
    const int target = integerOf(p, id, "target_atoms", 1);
    const int range = integerOf(p, id, "range", 1);
    if (range > 3) throw std::runtime_error("range at most 3 (the search grows as (2r+1)^9)");
    const long cells = std::lround(double(target) / double(s.atoms.size()));
    if (cells < 1) throw std::runtime_error("target_atoms is smaller than the cell");
    const Mat3 h = cellOf(s);
    const double volume = std::abs(science::determinant(h)) * cells;
    const double side2 = std::pow(volume, 2.0 / 3.0);
    // Deviation of the metric P H (P H)^T from a cube of the same volume.
    double bestScore = HUGE_VAL;
    int bestNegatives = 10, bestTrace = -100;
    Mat3 best{};
    std::array<int, 9> e{};
    const int span = 2 * range + 1;
    long total = 1;
    for (int k = 0; k < 9; ++k) total *= span;
    for (long code = 0; code < total; ++code) {
        long c = code;
        for (int k = 0; k < 9; ++k) { e[k] = static_cast<int>(c % span) - range; c /= span; }
        const long det = static_cast<long>(e[0]) * (e[4] * e[8] - e[5] * e[7]) - static_cast<long>(e[1]) * (e[3] * e[8] - e[5] * e[6]) +
                         static_cast<long>(e[2]) * (e[3] * e[7] - e[4] * e[6]);
        if (det != cells) continue;
        Mat3 q{};
        for (int r = 0; r < 3; ++r)
            for (int col = 0; col < 3; ++col) q[r][col] = e[3 * r] * h[0][col] + e[3 * r + 1] * h[1][col] + e[3 * r + 2] * h[2][col];
        double score = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const double g = science::dot(q[i], q[j]) - (i == j ? side2 : 0.0);
                score += g * g;
            }
        score = std::sqrt(score) / side2;
        int negatives = 0, trace = e[0] + e[4] + e[8];
        for (int k = 0; k < 9; ++k) negatives += e[k] < 0;
        const bool better = score < bestScore - 1e-9 ||
                            (score < bestScore + 1e-9 && (negatives < bestNegatives || (negatives == bestNegatives && trace > bestTrace)));
        if (better) {
            bestNegatives = negatives;
            bestTrace = trace;
            bestScore = score;
            for (int r = 0; r < 3; ++r)
                for (int col = 0; col < 3; ++col) best[r][col] = e[3 * r + col];
        }
    }
    if (!(bestScore < HUGE_VAL)) throw std::runtime_error("No supercell of " + std::to_string(cells) + " cells within the search range");
    std::ostringstream matrix;
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 3; ++col) matrix << (r || col ? " " : "") << best[r][col];
    const ModifierType* supercell = findModifierType("supercell");
    Json q = Json::object();
    q["matrix"] = matrix.str();
    supercell->apply(q, data);
    std::ostringstream note;
    note << "Supercell matrix " << matrix.str() << " (" << cells << " cells, " << data.structure.atoms.size()
         << " atoms); deviation from a cube " << bestScore;
    data.notes.insert(data.notes.begin(), note.str());
}

// ---------------------------------------------------------------- interpolation
void interpolate(const Json& p, PipelineData& data)
{
    const char* id = "interpolate";
    const std::string file = text(p, id, "file");
    if (file.empty()) throw std::runtime_error("Choose the final structure file");
    const double t = number(p, id, "fraction");
    const auto frames = science::readFrames(std::filesystem::u8path(file));
    const Structure& end = frames.front().structure;
    Structure& s = data.structure;
    if (end.atoms.size() != s.atoms.size())
        throw std::runtime_error("The final structure has " + std::to_string(end.atoms.size()) + " atoms, this one " + std::to_string(s.atoms.size()));
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (end.atoms[i].atomicNumber != s.atoms[i].atomicNumber)
            throw std::runtime_error("Atom " + std::to_string(i) + " is " + s.atoms[i].symbol + " here but " + end.atoms[i].symbol + " in the final structure");
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        Vec3 d = science::sub(positionOf(end.atoms[i]), positionOf(s.atoms[i]));
        if (s.hasUnitCell) d = minimumImage(s, d);  // atoms that crossed a periodic boundary move the short way
        place(s.atoms[i], science::add(positionOf(s.atoms[i]), science::scale(d, t)));
    }
    if (s.hasUnitCell && end.hasUnitCell)
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) s.cellVectors[r][c] += t * (end.cellVectors[r][c] - s.cellVectors[r][c]);
    data.notes.push_back("Interpolated " + std::to_string(t) + " of the way to " + std::filesystem::u8path(file).filename().u8string());
}

// ---------------------------------------------------------------- point defects
std::vector<std::size_t> atomsOf(const Structure& s, const std::string& symbol)
{
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (s.atoms[i].symbol == symbol) out.push_back(i);
    return out;
}

void antisite(const Json& p, PipelineData& data)
{
    const char* id = "create-antisite";
    const auto pair = elementList(text(p, id, "elements"), 2, "create-antisite");
    if (pair[0] == pair[1]) throw std::runtime_error("create-antisite: give two different elements");
    Structure& s = data.structure;
    const auto first = atomsOf(s, pair[0]), second = atomsOf(s, pair[1]);
    if (first.empty() || second.empty()) throw std::runtime_error("create-antisite: the structure needs both " + pair[0] + " and " + pair[1]);
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    const std::size_t a = first[std::uniform_int_distribution<std::size_t>(0, first.size() - 1)(rng)];
    std::size_t b = second.front();
    if (text(p, id, "partner") == "nearest") {
        double best = HUGE_VAL;
        for (std::size_t j : second) {
            Vec3 d = science::sub(positionOf(s.atoms[j]), positionOf(s.atoms[a]));
            if (s.hasUnitCell) d = minimumImage(s, d);
            if (science::norm(d) < best) { best = science::norm(d); b = j; }
        }
    } else
        b = second[std::uniform_int_distribution<std::size_t>(0, second.size() - 1)(rng)];
    std::swap(s.atoms[a].symbol, s.atoms[b].symbol);
    std::swap(s.atoms[a].atomicNumber, s.atoms[b].atomicNumber);
    getDefaultElementColor(s.atoms[a].atomicNumber, s.atoms[a].r, s.atoms[a].g, s.atoms[a].b);
    getDefaultElementColor(s.atoms[b].atomicNumber, s.atoms[b].r, s.atoms[b].g, s.atoms[b].b);
    data.selected.assign(s.atoms.size(), 0);
    data.selected[a] = data.selected[b] = 1;
    data.notes.push_back("Swapped atoms " + std::to_string(a) + " (" + pair[0] + ") and " + std::to_string(b) + " (" + pair[1] + ")");
}

void frenkelPair(const Json& p, PipelineData& data)
{
    const char* id = "create-frenkel-pair";
    const std::string symbol = standardElement(text(p, id, "element"));
    Structure& s = data.structure;
    if (atomsOf(s, symbol).empty()) throw std::runtime_error("create-frenkel-pair: no " + symbol + " atom");
    const std::size_t before = s.atoms.size();
    // The interstitial first (in the perfect structure, so the vacancy is not the void chosen).
    const ModifierType* insert = findModifierType("insert-interstitials");
    Json q = Json::object();
    q["element"] = symbol;
    q["count"] = 1.0;
    q["order"] = "random";
    q["seed"] = number(p, id, "seed");
    insert->apply(q, data);
    if (s.atoms.size() != before + 1) throw std::runtime_error("create-frenkel-pair: no interstitial void found");
    const Vec3 interstitial = positionOf(s.atoms.back());
    // The vacancy: an atom of the element at least min_distance away (the farthest when none is).
    const double minimum = number(p, id, "min_distance");
    std::vector<std::size_t> far;
    std::size_t farthest = 0;
    double farthestDistance = -1;
    for (std::size_t i = 0; i + 1 < s.atoms.size(); ++i) {
        if (s.atoms[i].symbol != symbol) continue;
        Vec3 d = science::sub(positionOf(s.atoms[i]), interstitial);
        if (s.hasUnitCell) d = minimumImage(s, d);
        const double dist = science::norm(d);
        if (dist >= minimum) far.push_back(i);
        if (dist > farthestDistance) { farthestDistance = dist; farthest = i; }
    }
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")) + 7);
    const std::size_t vacancy = far.empty() ? farthest : far[std::uniform_int_distribution<std::size_t>(0, far.size() - 1)(rng)];
    Vec3 d = science::sub(positionOf(s.atoms[vacancy]), interstitial);
    if (s.hasUnitCell) d = minimumImage(s, d);
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (i != vacancy) keep.push_back(i);
    rebuild(data, keep);
    data.selected.assign(data.structure.atoms.size(), 0);
    data.selected.back() = 1;
    data.notes.push_back("Frenkel pair: " + symbol + " atom " + std::to_string(vacancy) + " moved to an interstitial site " +
                         std::to_string(science::norm(d)) + " Angstrom away");
}
}  // namespace

std::vector<ModifierType> generatorStepTypes()
{
    const char* const layers = "graphene|hbn|mx2-2h|mx2-1t";
    return {
        {"build-nanotube", "Nanotube", "Build",
         "Rolls a honeycomb sheet into a tube with chiral indices (n, m): armchair (n, n), zigzag (n, 0) or chiral. Two elements make a BN-type tube "
         "(alternating sublattices). Multi-wall tubes (armchair or zigzag) add walls about wall_spacing apart. Periodic along the tube axis (z).",
         {{"n", "Chiral index n", "int", "6", ""}, {"m", "Chiral index m", "int", "6", ""}, {"length", "Unit cells along the axis", "int", "1", ""},
          {"elements", "Elements (one, or two for BN)", "elements", "C", ""}, {"bond", "Bond length (Angstrom)", "float", "1.421", ""},
          {"walls", "Walls", "int", "1", ""}, {"wall_spacing", "Wall spacing (Angstrom)", "float", "3.4", ""},
          {"vacuum", "Vacuum around the tube (Angstrom)", "float", "10", ""}}, buildNanotube},
        {"build-monolayer", "2D monolayer", "Build",
         "A hexagonal monolayer: graphene, hexagonal BN, or a 2H (trigonal prismatic) or 1T (octahedral) transition-metal dichalcogenide MX2, "
         "na x nb cells with vacuum along c. a = 0 takes the default lattice constant (graphene 2.461, h-BN 2.504, MoS2 3.18 Angstrom).",
         {{"kind", "Layer", "choice", "graphene", layers}, {"elements", "Elements (MX2: M X)", "elements", "", ""},
          {"a", "Lattice constant a (0 = default)", "float", "0", ""}, {"thickness", "MX2 X-X thickness (0 = 3.19)", "float", "0", ""},
          {"na", "Cells along a", "int", "1", ""}, {"nb", "Cells along b", "int", "1", ""}, {"vacuum", "Vacuum (Angstrom)", "float", "15", ""}}, buildMonolayer},
        {"build-nanoribbon", "Nanoribbon", "Build",
         "A graphene (or BN) ribbon with armchair edges (width in dimer lines, N-AGNR) or zigzag edges (width in zigzag chains, N-ZGNR), "
         "periodic along x, with hydrogen on the edge atoms.",
         {{"edge", "Edge", "choice", "armchair", "armchair|zigzag"}, {"width", "Width (dimer lines or chains)", "int", "7", ""},
          {"length", "Periods along the ribbon", "int", "1", ""}, {"elements", "Elements (one, or two for BN)", "elements", "C", ""},
          {"bond", "Bond length (Angstrom)", "float", "1.421", ""}, {"passivate", "Hydrogen on the edges", "bool", "true", ""},
          {"vacuum", "Vacuum (Angstrom)", "float", "10", ""}}, buildNanoribbon},
        {"build-twisted-bilayer", "Twisted bilayer", "Build",
         "Two hexagonal layers turned by a commensurate angle (cos t = (3i^2+3i+1/2)/(3i^2+3i+1): 21.79, 13.17, 9.43, ... degrees), in their "
         "moire supercell. Give the angle (the nearest commensurate one is used) or the index i.",
         {{"kind", "Layer", "choice", "graphene", layers}, {"elements", "Elements (MX2: M X)", "elements", "", ""},
          {"angle", "Twist angle (degrees; 0 = use the index)", "float", "0", ""}, {"index", "Commensurate index i", "int", "3", ""},
          {"a", "Lattice constant a (0 = default)", "float", "0", ""}, {"thickness", "MX2 X-X thickness (0 = 3.19)", "float", "0", ""},
          {"interlayer", "Interlayer distance (0 = default)", "float", "0", ""}, {"vacuum", "Vacuum (Angstrom)", "float", "15", ""}}, buildTwistedBilayer},
        {"build-cluster", "Magic-number cluster", "Build",
         "Closed-shell clusters of one element: Mackay icosahedra and fcc cuboctahedra (1, 13, 55, 147, ... atoms for 1, 2, 3, 4 shells) "
         "and regular decahedra (pentagonal bipyramids of five tetrahedra: 7, 23, 54, ... atoms for 2, 3, 4 atoms per edge), "
         "with the nearest-neighbour distance along the edges.",
         {{"kind", "Shape", "choice", "icosahedron", "icosahedron|decahedron|cuboctahedron"}, {"element", "Element", "element", "Cu", ""},
          {"shells", "Shells (1 = one atom)", "int", "3", ""}, {"distance", "Nearest-neighbour distance (Angstrom)", "float", "2.556", ""},
          {"vacuum", "Vacuum (Angstrom)", "float", "10", ""}}, buildCluster},
        {"build-molecule", "Molecule", "Build",
         "A molecule from the library (experimental geometries): H2, N2, O2, F2, Cl2, CO, NO, HF, HCl, H2O, H2S, CO2, SO2, NO2, O3, NH3, PH3, "
         "CH4, SiH4, CF4, C2H2, C2H4, C2H6, CH3OH, C6H6, C60 (or common names: water, benzene, ...). With vacuum, in an open box.",
         {{"name", "Molecule", "string", "H2O", ""}, {"vacuum", "Vacuum (Angstrom; 0 = no cell)", "float", "10", ""}}, buildMolecule},
        {"build-packing", "Molecular packing", "Build",
         "Packs molecules (from the molecule library, or single atoms) at random positions and orientations into a periodic box, "
         "keeping atoms of different molecules at least min_distance apart (as Packmol). The box follows from the density unless given.",
         {{"molecules", "Molecules and counts (e.g. H2O 100 CO2 10)", "list", "H2O 64", ""},
          {"density", "Density (g/cm3)", "float", "1.0", ""}, {"box", "Box lengths (Angstrom; 0 0 0 = from the density)", "vector", "0 0 0", ""},
          {"min_distance", "Minimum distance between molecules (Angstrom)", "float", "2.0", ""}, {"seed", "Seed", "int", "1", ""},
          {"attempts", "Attempts per molecule", "int", "20000", ""}}, buildPacking},
        {"add-adsorbate", "Add adsorbate", "Modification",
         "Places a molecule or atom on the top surface of a slab (normal a x b, towards c) at a top, bridge, fcc-hollow or hcp-hollow site, "
         "its lowest atom `height` above the top layer; the site nearest to the fractional surface position (u, v) is used.",
         {{"species", "Molecule or element", "string", "CO", ""}, {"site", "Site", "choice", "top", "top|bridge|fcc|hcp"},
          {"height", "Height above the surface (Angstrom)", "float", "1.8", ""}, {"position", "Near fractional position (u v, w ignored)", "vector", "0.5 0.5 0", ""},
          {"select", "Select the adsorbate", "bool", "true", ""}}, addAdsorbate},
        {"near-cubic-supercell", "Near-cubic supercell", "Cell",
         "The supercell (integer matrix with entries up to +-range) of about target_atoms atoms whose cell is closest to a cube, "
         "for point-defect calculations (as ASE's find_optimal_cell_shape).",
         {{"target_atoms", "Target number of atoms", "int", "64", ""}, {"range", "Matrix entries up to +-", "int", "2", ""}}, nearCubicSupercell},
        {"interpolate", "Interpolate to a structure", "Modification",
         "Moves every atom the fraction t of the way to its position in another structure with the same atoms in the same order "
         "(minimum image across periodic boundaries); the cell is interpolated too. Images for a nudged elastic band are t = k/(n+1).",
         {{"file", "Final structure file", "string", "", ""}, {"fraction", "Fraction t (0 = this, 1 = final)", "float", "0.5", ""}}, interpolate},
        {"create-antisite", "Antisite defect", "Modification",
         "Swaps an atom of the first element (chosen at random) with the nearest (or a random) atom of the second; both are selected.",
         {{"elements", "Elements A B", "elements", "Ni Al", ""}, {"partner", "Partner", "choice", "nearest", "nearest|random"},
          {"seed", "Seed", "int", "1", ""}}, antisite},
        {"create-frenkel-pair", "Frenkel pair", "Modification",
         "Moves an atom of the element into an interstitial void (as Insert interstitial atoms), leaving a vacancy at least min_distance away; "
         "the interstitial is selected.",
         {{"element", "Element", "element", "Cu", ""}, {"min_distance", "Minimum vacancy-interstitial distance (Angstrom)", "float", "5", ""},
          {"seed", "Seed", "int", "1", ""}}, frenkelPair},
    };
}
}  // namespace atomforge::pipeline
