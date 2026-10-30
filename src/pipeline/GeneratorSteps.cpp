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
// A monolayer: in-plane cell vectors a1, a2 and the basis in fractional
// (in-plane) and Cartesian z coordinates. Hexagonal layers have
// a1 = (a, 0, 0), a2 = (-a/2, a sqrt3/2, 0).
struct Layer
{
    double a = 0;
    Vec3 a1{}, a2{};
    bool hexagonal = true;
    struct Site { std::string symbol; double u, v, z; int sublattice; };
    std::vector<Site> sites;
    double thickness = 0;
    double interlayer = 3.35;  // distance between the mid-planes of stacked layers
    std::string kind;
};

Layer hexagonalLayer(double a)
{
    Layer l;
    l.a = a;
    l.a1 = {a, 0, 0};
    l.a2 = {-a / 2, a * std::sqrt(3.0) / 2, 0};
    return l;
}

// Buckled honeycombs (silicene, germanene, stanene): the sublattices `buckling` apart along z.
Layer xene(const std::string& element, double a, double buckling)
{
    Layer l = hexagonalLayer(a);
    l.sites = {{element, 1.0 / 3, 2.0 / 3, buckling / 2, 0}, {element, 2.0 / 3, 1.0 / 3, -buckling / 2, 1}};
    l.thickness = buckling;
    return l;
}

Layer monolayer(const std::string& kind, const std::string& elementsText, double a, double thickness)
{
    Layer layer;
    if (kind == "graphene") {
        const auto e = elementList(elementsText.empty() ? "C" : elementsText, 1, "graphene");
        layer = hexagonalLayer(a > 0 ? a : 2.4612);  // C-C 1.421 Angstrom
        layer.sites = {{e[0], 1.0 / 3, 2.0 / 3, 0, 0}, {e[0], 2.0 / 3, 1.0 / 3, 0, 1}};
        layer.interlayer = 3.35;
    } else if (kind == "hbn") {
        const auto e = elementList(elementsText.empty() ? "B N" : elementsText, 2, "h-BN");
        layer = hexagonalLayer(a > 0 ? a : 2.504);
        layer.sites = {{e[0], 1.0 / 3, 2.0 / 3, 0, 0}, {e[1], 2.0 / 3, 1.0 / 3, 0, 1}};
        layer.interlayer = 3.33;
    } else if (kind == "silicene" || kind == "germanene" || kind == "stanene") {
        // Low-buckled honeycombs (DFT-PBE lattice constants and bucklings, Cahangirov et al. 2009; Xu et al. 2013).
        const char* element = kind == "silicene" ? "Si" : kind == "germanene" ? "Ge" : "Sn";
        const double a0 = kind == "silicene" ? 3.866 : kind == "germanene" ? 4.02 : 4.67;
        const double b0 = kind == "silicene" ? 0.44 : kind == "germanene" ? 0.69 : 0.85;
        const auto e = elementList(elementsText.empty() ? element : elementsText, 1, kind.c_str());
        layer = xene(e[0], a > 0 ? a : a0, thickness > 0 ? thickness : b0);
        layer.interlayer = layer.thickness + 3.0;
    } else if (kind == "mx2-2h" || kind == "mx2-1t") {
        const auto e = elementList(elementsText.empty() ? "Mo S" : elementsText, 2, "MX2");
        layer = hexagonalLayer(a > 0 ? a : 3.18);
        layer.thickness = thickness > 0 ? thickness : 3.19;  // X-X distance across the layer (ASE mx2 defaults)
        const double h = layer.thickness / 2;
        layer.sites = {{e[0], 0, 0, 0, 0}, {e[1], 1.0 / 3, 2.0 / 3, h, 1}};
        if (kind == "mx2-2h") layer.sites.push_back({e[1], 1.0 / 3, 2.0 / 3, -h, 1});
        else layer.sites.push_back({e[1], 2.0 / 3, 1.0 / 3, -h, 1});
        layer.interlayer = 6.15;
    } else if (kind == "phosphorene") {
        // Black phosphorus layer (bulk Cmce, Brown & Rundqvist 1965): a along the
        // zigzag direction, b along the armchair direction, puckered 2.13 A thick.
        const auto e = elementList(elementsText.empty() ? "P" : elementsText, 1, "phosphorene");
        const double za = a > 0 ? a : 3.3136, zb = 4.3763 * (a > 0 ? a / 3.3136 : 1.0);
        layer.a = za;
        layer.a1 = {za, 0, 0};
        layer.a2 = {0, zb, 0};
        layer.hexagonal = false;
        const double h = thickness > 0 ? thickness / 2 : 0.10168 * 10.478;
        layer.sites = {{e[0], 0, 0.08056, h, 0}, {e[0], 0, 0.91944, -h, 1}, {e[0], 0.5, 0.58056, -h, 1}, {e[0], 0.5, 0.41944, h, 0}};
        layer.thickness = 2 * h;
        layer.interlayer = 5.239;  // bulk b / 2
    } else
        throw std::runtime_error("Unknown layer '" + kind + "': graphene, hbn, silicene, germanene, stanene, mx2-2h, mx2-1t or phosphorene");
    layer.kind = kind;
    return layer;
}

Vec3 layerA1(const Layer& l) { return l.a1; }
Vec3 layerA2(const Layer& l) { return l.a2; }

// Biaxial strain (fraction) of a layer's in-plane lattice.
void strainLayer(Layer& l, double strain)
{
    l.a *= 1 + strain;
    l.a1 = science::scale(l.a1, 1 + strain);
    l.a2 = science::scale(l.a2, 1 + strain);
}

// The in-plane registry of layer k of a stack: AA (none), AB (Bernal:
// alternate layers shifted by (2/3, 1/3)), ABC (rhombohedral: k (2/3, 1/3)),
// 2H (MX2: alternate layers turned by 180 degrees with M over X).
struct Registry { double du = 0, dv = 0; bool inverted = false; };
Registry registryOf(const std::string& stacking, int k, const Layer& layer)
{
    if (stacking == "AA" || k == 0) return {};
    if (!layer.hexagonal) throw std::runtime_error("Only AA stacking is available for " + layer.kind);
    if (stacking == "AB") return k % 2 ? Registry{2.0 / 3, 1.0 / 3, false} : Registry{};
    if (stacking == "ABC") return {std::fmod(k * 2.0 / 3, 1.0), std::fmod(k * 1.0 / 3, 1.0), false};
    if (stacking == "2H") return k % 2 ? Registry{1.0 / 3, 2.0 / 3, true} : Registry{};
    throw std::runtime_error("stacking must be AA, AB, ABC or 2H");
}

// Fills the in-plane cell (L1, L2) with a layer turned by `angle` (and shifted
// by its registry) at height z; the number of atoms placed is added to count.
void fillLayer(Structure& s, const Layer& layer, double angle, const Vec3& l1, const Vec3& l2, double z, std::size_t& count,
               const Registry& registry = {})
{
    const Vec3 a1 = rotateZ(layerA1(layer), angle), a2 = rotateZ(layerA2(layer), angle);
    Mat3 super{l1, l2, Vec3{0, 0, 1}};
    const Mat3 inverse = science::inverse(super);
    const double shortest = std::min(science::norm(layer.a1), science::norm(layer.a2));
    const double reach = (science::norm(l1) + science::norm(l2)) / shortest + 3;
    const int range = static_cast<int>(std::ceil(reach)) + 2;
    const double eps = 1e-7;
    for (int i = -2 * range; i <= 2 * range; ++i)
        for (int j = -2 * range; j <= 2 * range; ++j)
            for (const auto& site : layer.sites) {
                const double u = (registry.inverted ? -site.u : site.u) + registry.du, v = (registry.inverted ? -site.v : site.v) + registry.dv;
                const Vec3 q = science::add(science::scale(a1, i + u), science::scale(a2, j + v));
                const Vec3 f = science::rowTimes({q[0], q[1], 0}, inverse);
                if (f[0] < -eps || f[0] >= 1 - eps || f[1] < -eps || f[1] >= 1 - eps) continue;
                addAtom(s, site.symbol, {q[0], q[1], z + (registry.inverted ? -site.z : site.z)});
                ++count;
            }
}

void buildMonolayer(const Json& p, PipelineData& data)
{
    const char* id = "build-monolayer";
    Layer layer = monolayer(text(p, id, "kind"), text(p, id, "elements"), number(p, id, "a"), number(p, id, "thickness"));
    strainLayer(layer, number(p, id, "strain") / 100.0);
    const int na = integerOf(p, id, "na", 1), nb = integerOf(p, id, "nb", 1), layers = integerOf(p, id, "layers", 1);
    const double vacuum = number(p, id, "vacuum");
    if (vacuum < 0) throw std::runtime_error("vacuum must be nonnegative");
    const double interlayer = number(p, id, "interlayer") > 0 ? number(p, id, "interlayer") : layer.interlayer;
    const std::string stacking = text(p, id, "stacking");
    // The in-plane cell: the primitive one, or for hexagonal layers the rectangular a x a sqrt3.
    Vec3 l1 = layer.a1, l2 = layer.a2;
    const bool rectangular = flag(p, id, "rectangular") && layer.hexagonal;
    if (rectangular) l2 = science::add(layer.a1, science::scale(layer.a2, 2));
    l1 = science::scale(l1, na);
    l2 = science::scale(l2, nb);
    const double stack = (layers - 1) * interlayer;
    const double height = stack + layer.thickness + vacuum;
    Structure s;
    std::size_t count = 0;
    for (int k = 0; k < layers; ++k)
        fillLayer(s, layer, 0, l1, l2, (height - stack) / 2 + k * interlayer, count, registryOf(stacking, k, layer));
    setCell(s, l1, l2, {0, 0, std::max(height, 1.0)}, {true, true, flag(p, id, "periodic_c")});
    produce(data, std::move(s));
    std::ostringstream note;
    note.precision(6);
    note << data.structure.atoms.size() << " atoms, " << layers << " layer" << (layers > 1 ? "s " + stacking : std::string()) << ", a = " << layer.a
         << " Angstrom" << (rectangular ? ", rectangular cell" : "");
    data.notes.push_back(note.str());
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

// Removes atoms with fewer than two bonds, repeatedly (the dangling atoms a cut leaves).
std::size_t removeDangling(Structure& s, double bondScale = 1.2)
{
    std::size_t removed = 0;
    for (int pass = 0; pass < 20; ++pass) {
        const Neighbours nb = bondNeighbours(s, bondScale);
        std::vector<AtomSite> kept;
        for (std::size_t i = 0; i < s.atoms.size(); ++i)
            if (nb.atoms[i].size() >= 2) kept.push_back(s.atoms[i]);
        if (kept.size() == s.atoms.size()) break;
        removed += s.atoms.size() - kept.size();
        s.atoms = std::move(kept);
    }
    return removed;
}

// Turns a structure built along z so that its axis lies along x or y (the cell follows).
void axisTo(Structure& s, const std::string& axis)
{
    if (axis == "z") return;
    const auto permute = [&](std::array<double, 3> v) {
        return axis == "x" ? std::array<double, 3>{v[2], v[0], v[1]} : std::array<double, 3>{v[1], v[2], v[0]};
    };
    for (auto& a : s.atoms) place(a, permute({a.x, a.y, a.z}));
    if (s.hasUnitCell) {
        const auto old = s.cellVectors;
        const auto oldPbc = s.pbc;
        // New cell vector k is the old vector whose direction became axis k.
        for (int k = 0; k < 3; ++k) {
            const int from = axis == "x" ? (k + 2) % 3 : (k + 1) % 3;
            s.cellVectors[static_cast<std::size_t>(k)] = permute(old[static_cast<std::size_t>(from)]);
            s.pbc[static_cast<std::size_t>(k)] = oldPbc[static_cast<std::size_t>(from)];
        }
    }
}

void buildNanotube(const Json& p, PipelineData& data)
{
    const char* id = "build-nanotube";
    const int n = integerOf(p, id, "n", 0), m = integerOf(p, id, "m", 0);
    if (n == 0 && m == 0) throw std::runtime_error("The chiral indices (n, m) must not both be zero");
    if (m > n) throw std::runtime_error("Use n >= m (the tube (m, n) is the mirror image of (n, m))");
    int cells = integerOf(p, id, "length", 1);
    const int walls = integerOf(p, id, "walls", 1);
    const auto species = elementList(text(p, id, "elements"), 2, "nanotube");
    // Bond length: given, or C-C 1.421, B-N 1.446, else the sum of covalent radii.
    double bond = number(p, id, "bond");
    if (!(bond > 0)) {
        if (species[0] == "C" && species[1] == "C") bond = 1.421;
        else if ((species[0] == "B" && species[1] == "N") || (species[0] == "N" && species[1] == "B")) bond = 1.446;
        else bond = covalentRadius(science::atomicNumber(species[0])) + covalentRadius(science::atomicNumber(species[1]));
    }
    const double spacing = number(p, id, "wall_spacing"), vacuum = number(p, id, "vacuum");
    if (!(spacing > 0) || vacuum < 0) throw std::runtime_error("wall spacing must be positive, vacuum nonnegative");
    if (walls > 1 && !(n == m || m == 0))
        throw std::runtime_error("Multi-wall tubes need armchair (n, n) or zigzag (n, 0) indices, which share the axial period");
    // A length in Angstrom sets the number of unit cells.
    {
        double radius = 0, period = 0;
        Structure probe;
        addTubeWall(probe, n, m, bond, species, 1, radius, period);
        const double wanted = number(p, id, "length_A");
        if (wanted > 0) cells = std::max(1, static_cast<int>(std::lround(wanted / period)));
    }
    Structure s;
    double radius = 0, period = 0, outer = 0;
    std::ostringstream notes;
    notes.precision(5);
    for (int w = 0; w < walls; ++w) {
        // Each wall grows the radius by about the wall spacing.
        const double a = bond * std::sqrt(3.0);
        const int step = n == m ? static_cast<int>(std::lround(spacing * 2 * kPi / (a * std::sqrt(3.0))))
                                : static_cast<int>(std::lround(spacing * 2 * kPi / a));
        const int nw = n + w * step, mw = m == 0 ? 0 : m + w * step;
        addTubeWall(s, nw, mw, bond, species, cells, radius, period);
        outer = radius;
        notes << (w ? ", " : "") << "(" << nw << "," << mw << ") d = " << 2 * radius;
    }
    const bool finite = flag(p, id, "finite");
    const double box = 2 * outer + vacuum + 2 * bond;
    shift(s, {box / 2, box / 2, finite ? vacuum / 2 : 0});
    std::size_t capped = 0, dropped = 0;
    if (finite) {
        // A finite tube: open ends without dangling atoms, capped with the terminating atom.
        setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, period * cells + vacuum}, {false, false, false});
        dropped = removeDangling(s);
        const std::string cap = text(p, id, "end_cap");
        if (cap != "none") {
            std::vector<char> selected(s.atoms.size(), 0);
            SaturationOptions options;
            options.terminator = standardElement(cap);
            capped = saturate(s, selected, {species[0], species[1]}, false, 1.2, false, options);
        }
    } else
        setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, period * cells}, {false, false, true});
    axisTo(s, text(p, id, "axis"));
    produce(data, std::move(s));
    // Chiral angle, and for carbon the band character ((n - m) divisible by 3: metallic).
    const double chiral = std::atan2(std::sqrt(3.0) * m, 2.0 * n + m) * 180 / kPi;
    std::ostringstream note;
    note.precision(5);
    note << data.structure.atoms.size() << " atoms; walls " << notes.str() << " Angstrom; chiral angle " << chiral << " degrees; period " << period
         << " Angstrom x " << cells;
    if (species[0] == "C" && species[1] == "C")
        note << "; " << ((n - m) % 3 == 0 ? "metallic" : "semiconducting") << " (zone folding)";
    if (finite) note << "; finite, " << capped << " end atoms " << text(p, id, "end_cap") << (dropped ? ", " + std::to_string(dropped) + " dangling atoms removed" : "");
    data.notes.push_back(note.str());
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
    const double across = ribbonWidth + vacuum + 2.2;  // room for the edge atoms
    const double height = std::max(vacuum, 1.0);
    const bool finite = flag(p, id, "finite");
    shift(s, {finite ? vacuum / 2 : 0, (across - ribbonWidth) / 2, height / 2});
    if (finite) setCell(s, {period * length + vacuum, 0, 0}, {0, across, 0}, {0, 0, height}, {false, false, false});
    else setCell(s, {period * length, 0, 0}, {0, across, 0}, {0, 0, height}, {true, false, false});
    const std::size_t dropped = finite ? removeDangling(s) : 0;
    std::vector<char> selected(s.atoms.size(), 0);
    std::size_t terminated = 0;
    const std::string termination = text(p, id, "termination");
    if (flag(p, id, "passivate") && termination != "none") {
        SaturationOptions options;
        options.terminator = standardElement(termination);
        options.bondLength = number(p, id, "termination_bond");
        terminated = saturate(s, selected, {species[0], species[1]}, false, 1.2, false, options);
    }
    produce(data, std::move(s));
    std::ostringstream note;
    note.precision(5);
    note << width << (zigzag ? " chains, zigzag" : " dimer lines, armchair") << " edges, width " << ribbonWidth << " Angstrom; "
         << data.structure.atoms.size() << " atoms (" << terminated << " " << termination << ")";
    if (finite) note << "; finite flake" << (dropped ? ", " + std::to_string(dropped) + " dangling atoms removed" : "");
    if (!zigzag) note << "; family " << width % 3 << " (3p, 3p+1, 3p+2 differ in band gap)";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- twisted bilayers
void buildTwistedBilayer(const Json& p, PipelineData& data)
{
    const char* id = "build-twisted-bilayer";
    Layer layer = monolayer(text(p, id, "kind"), text(p, id, "elements"), number(p, id, "a"), number(p, id, "thickness"));
    if (!layer.hexagonal) throw std::runtime_error("Twisted bilayers need hexagonal layers");
    // The top layer: the same, or another hexagonal layer strained to the bottom lattice (a heterobilayer).
    const std::string topKind = text(p, id, "top_kind");
    Layer top = topKind == "same" ? layer : monolayer(topKind, text(p, id, "top_elements"), 0, 0);
    if (!top.hexagonal) throw std::runtime_error("Twisted bilayers need hexagonal layers");
    const double naturalTop = top.a;
    const double mismatch = layer.a / naturalTop - 1;
    strainLayer(top, mismatch);
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
    const double interlayer = number(p, id, "interlayer") > 0 ? number(p, id, "interlayer") : (layer.interlayer + top.interlayer) / 2;
    const double vacuum = number(p, id, "vacuum");
    const double height = interlayer + (layer.thickness + top.thickness) / 2 + vacuum;
    Structure s;
    std::size_t bottomCount = 0, topCount = 0;
    const double z0 = (height - interlayer) / 2;
    fillLayer(s, layer, 0, l1, l2, z0, bottomCount);
    // The second layer turned so that its vector (i+1, i) lies on L1; its registry
    // before the twist: AA (atoms on atoms at the rotation centre) or AB (Bernal).
    const Vec3 v2 = science::add(science::scale(b1, index + 1), science::scale(b2, index));
    const double turn = std::atan2(l1[1], l1[0]) - std::atan2(v2[1], v2[0]);
    const std::string stacking = text(p, id, "stacking");
    const Registry registry = stacking == "AB" ? Registry{2.0 / 3, 1.0 / 3, false} : Registry{};
    if (stacking != "AA" && stacking != "AB") throw std::runtime_error("stacking must be AA or AB");
    fillLayer(s, top, turn, l1, l2, z0 + interlayer, topCount, registry);
    const std::size_t expectedBottom = layer.sites.size() * static_cast<std::size_t>(3 * index * index + 3 * index + 1);
    const std::size_t expectedTop = top.sites.size() * static_cast<std::size_t>(3 * index * index + 3 * index + 1);
    if (bottomCount != expectedBottom || topCount != expectedTop)
        throw std::runtime_error("Twisted bilayer layers have " + std::to_string(bottomCount) + " and " + std::to_string(topCount) +
                                 " atoms, expected " + std::to_string(expectedBottom) + " and " + std::to_string(expectedTop));
    setCell(s, l1, l2, {0, 0, height}, {true, true, true});
    produce(data, std::move(s));
    std::ostringstream note;
    note << "Twist " << theta * 180 / kPi << " degrees (index " << index << ", " << stacking << " before the twist), " << data.structure.atoms.size()
         << " atoms, |L| = " << science::norm(l1) << " Angstrom, moire period " << layer.a / (2 * std::sin(theta / 2)) << " Angstrom";
    if (topKind != "same") note << "; top " << topKind << " strained by " << 100 * mismatch << "% from a = " << naturalTop;
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
    } else if (kind == "octahedron") {
        // `shells` atoms along each edge: fcc points (x + y + z of the parity of shells - 1) inside
        // |x| + |y| + |z| <= shells - 1, in units of d / sqrt 2.
        const double unit = d / std::sqrt(2.0);
        const int k = shells - 1;
        points.clear();
        for (int x = -k; x <= k; ++x)
            for (int y = -k; y <= k; ++y)
                for (int z = -k; z <= k; ++z)
                    if (((x + y + z - k) % 2 + 2) % 2 == 0 && std::abs(x) + std::abs(y) + std::abs(z) <= k)
                        points.push_back({x * unit, y * unit, z * unit});
    } else
        throw std::runtime_error("Unknown cluster '" + kind + "': icosahedron, decahedron, cuboctahedron or octahedron");
    double extent = 0;
    for (const Vec3& q : points) extent = std::max(extent, science::norm(q));
    const double box = 2 * extent + vacuum;
    // Nanoalloys: a fraction of the atoms becomes the second element, at random,
    // in the core, in the outer shell(s), or on one side (Janus).
    std::vector<std::string> symbols(points.size(), symbol);
    const std::string second = text(p, id, "second_element");
    std::size_t substituted = 0;
    if (!second.empty()) {
        const std::string other = standardElement(second);
        const double fraction = number(p, id, "fraction");
        if (fraction < 0 || fraction > 1) throw std::runtime_error("fraction must lie between 0 and 1");
        const std::string distribution = text(p, id, "distribution");
        std::vector<std::size_t> order(points.size());
        std::iota(order.begin(), order.end(), 0);
        if (distribution == "random") {
            std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
            std::shuffle(order.begin(), order.end(), rng);
        } else {
            std::vector<double> key(points.size());
            for (std::size_t i = 0; i < points.size(); ++i) {
                const double r = science::norm(points[i]);
                key[i] = distribution == "core" ? r : distribution == "shell" ? -r : distribution == "janus" ? -points[i][2] : 0;
            }
            if (distribution != "core" && distribution != "shell" && distribution != "janus")
                throw std::runtime_error("distribution must be random, core, shell or janus");
            std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return key[a] < key[b] - 1e-9; });
        }
        substituted = static_cast<std::size_t>(std::lround(fraction * points.size()));
        for (std::size_t k = 0; k < substituted; ++k) symbols[order[k]] = other;
    }
    Structure s;
    for (std::size_t i = 0; i < points.size(); ++i) addAtom(s, symbols[i], science::add(points[i], {box / 2, box / 2, box / 2}));
    setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, box}, {false, false, false});
    produce(data, std::move(s));
    std::ostringstream note;
    note.precision(5);
    note << kind << " of " << shells << " shells: " << points.size() << " atoms, diameter " << 2 * extent << " Angstrom";
    if (substituted) note << "; " << substituted << " " << standardElement(second) << " (" << text(p, id, "distribution") << ")";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- molecules
// Rotates molecule atoms by R (positions about the centre of mass).
std::vector<MoleculeAtom> rotated(std::vector<MoleculeAtom> atoms, const Mat3& r)
{
    for (auto& a : atoms) a.position = rotate(r, a.position);
    return atoms;
}

// The rotation taking the principal axes (smallest moment first) onto z, x, y:
// long molecules lie along z, planar ones in the xz plane.
Mat3 principalFrame(const std::vector<MoleculeAtom>& atoms)
{
    const auto axes = principalAxes(atoms);
    // Rows: the new x, y, z expressed in the old frame.
    return Mat3{axes[1], axes[2], axes[0]};
}

Mat3 randomRotation(std::mt19937_64& rng);

void buildMolecule(const Json& p, PipelineData& data)
{
    const char* id = "build-molecule";
    auto atoms = molecule(text(p, id, "name"));
    const std::string orientation = text(p, id, "orientation");
    if (orientation == "principal" && atoms.size() > 1) atoms = rotated(atoms, principalFrame(atoms));
    else if (orientation == "random") {
        std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
        atoms = rotated(atoms, randomRotation(rng));
    } else if (orientation != "library" && orientation != "principal")
        throw std::runtime_error("orientation must be library, principal or random");
    const double vacuum = number(p, id, "vacuum");
    Structure s;
    double extent = 0;
    for (const auto& a : atoms) extent = std::max(extent, science::norm(a.position));
    const double box = 2 * extent + vacuum;
    for (const auto& a : atoms) addAtom(s, a.symbol, vacuum > 0 ? science::add(a.position, {box / 2, box / 2, box / 2}) : a.position);
    if (vacuum > 0) setCell(s, {box, 0, 0}, {0, box, 0}, {0, 0, box}, {flag(p, id, "periodic"), flag(p, id, "periodic"), flag(p, id, "periodic")});
    produce(data, std::move(s));
    std::ostringstream note;
    note.precision(6);
    note << atoms.size() << " atoms, " << moleculeMass(atoms) << " u, extent " << 2 * extent << " Angstrom";
    data.notes.push_back(note.str());
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

// A Packmol-style constraint: every atom of a placed molecule lies inside it.
struct Region
{
    std::string kind = "all";
    Vec3 low{}, high{}, centre{};
    double radius = 0;
    int axis = 2;
    bool contains(const Vec3& q) const
    {
        if (kind == "all") return true;
        if (kind == "box") return q[0] >= low[0] && q[0] <= high[0] && q[1] >= low[1] && q[1] <= high[1] && q[2] >= low[2] && q[2] <= high[2];
        if (kind == "sphere") return science::norm(science::sub(q, centre)) <= radius;
        if (kind == "outside-sphere") return science::norm(science::sub(q, centre)) >= radius;
        if (kind == "slab") return q[axis] >= low[axis] && q[axis] <= high[axis];
        if (kind == "cylinder") {
            Vec3 d = science::sub(q, centre);
            d[axis] = 0;
            return science::norm(d) <= radius;
        }
        return true;
    }
    // Volume inside the box lengths (for automatic counts).
    double volume(const Vec3& box) const
    {
        if (kind == "box") return std::max(0.0, high[0] - low[0]) * std::max(0.0, high[1] - low[1]) * std::max(0.0, high[2] - low[2]);
        if (kind == "sphere") return 4.0 / 3.0 * kPi * radius * radius * radius;
        if (kind == "outside-sphere") return box[0] * box[1] * box[2] - 4.0 / 3.0 * kPi * radius * radius * radius;
        if (kind == "slab") return (high[axis] - low[axis]) * box[0] * box[1] * box[2] / box[axis];
        if (kind == "cylinder") return kPi * radius * radius * box[axis];
        return box[0] * box[1] * box[2];
    }
};

void buildPacking(const Json& p, PipelineData& data)
{
    const char* id = "build-packing";
    // "H2O 100 CO2 10": molecules (library names or files) and counts (0: from the density).
    const auto list = words(text(p, id, "molecules"));
    if (list.empty() || list.size() % 2) throw std::runtime_error("molecules: give pairs of a molecule and a count, e.g. H2O 100");
    struct Kind { std::vector<MoleculeAtom> atoms; int count; std::string name; bool automatic = false; };
    std::vector<Kind> kinds;
    for (std::size_t k = 0; k < list.size(); k += 2) {
        char* end = nullptr;
        const long count = std::strtol(list[k + 1].c_str(), &end, 10);
        if (*end || count < 0) throw std::runtime_error("molecules: '" + list[k + 1] + "' is not a count");
        kinds.push_back({molecule(list[k]), static_cast<int>(count), list[k], count == 0});
    }
    const bool solvate = text(p, id, "mode") == "solvate";
    const double density = number(p, id, "density");
    Vec3 box = vector(p, id, "box");
    Structure s;
    double soluteMass = 0;
    if (solvate) {
        // Around the input structure (the solute) in its orthogonal cell.
        requireCell(data, "Solvation");
        s = data.structure;
        const Mat3 cell = cellOf(s);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                if (r != c && std::abs(cell[r][c]) > 1e-6) throw std::runtime_error("Solvation needs an orthogonal cell (use shrink-wrap or set-cell)");
        box = {cell[0][0], cell[1][1], cell[2][2]};
        for (const auto& a : s.atoms) soluteMass += elementAtomicMass(a.atomicNumber);
    }
    Region region;
    region.kind = text(p, id, "region");
    region.low = vector(p, id, "region_low");
    region.high = vector(p, id, "region_high");
    region.centre = vector(p, id, "region_center");
    region.radius = number(p, id, "region_radius");
    const std::string axisName = text(p, id, "region_axis");
    region.axis = axisName == "x" ? 0 : axisName == "y" ? 1 : 2;
    // Box lengths from the density (all counts given) when not set.
    double fixedMass = 0;
    bool automatic = false;
    for (const Kind& k : kinds) {
        fixedMass += moleculeMass(k.atoms) * k.count;
        automatic = automatic || k.count == 0;
    }
    if (box[0] <= 0 || box[1] <= 0 || box[2] <= 0) {
        if (!(density > 0) || automatic) throw std::runtime_error("Give the box lengths, or a density with every count");
        const double side = std::cbrt(fixedMass * 1.66053906660 / density);
        box = {side, side, side};
    }
    // Counts of 0: fill the region (less the solute's share) to the density, shared equally by those molecules.
    if (automatic) {
        if (!(density > 0)) throw std::runtime_error("A count of 0 needs a density to fill to");
        // The region's free volume: voxels (about 0.4 Angstrom) not within the covalent radius plus half
        // solute_distance of a solute atom.
        int n[3];
        for (int k = 0; k < 3; ++k) n[k] = std::clamp(static_cast<int>(std::ceil(box[k] / 0.4)), 1, 400);
        std::vector<char> taken(static_cast<std::size_t>(n[0]) * n[1] * n[2], 0);
        for (const auto& a : s.atoms) {
            const double r = covalentRadius(a.atomicNumber) + number(p, id, "solute_distance") / 2;
            const Vec3 q{a.x, a.y, a.z};
            int lo[3], hi[3];
            for (int k = 0; k < 3; ++k) {
                lo[k] = static_cast<int>(std::floor((q[k] - r) / box[k] * n[k]));
                hi[k] = static_cast<int>(std::ceil((q[k] + r) / box[k] * n[k]));
            }
            for (int x = lo[0]; x <= hi[0]; ++x)
                for (int y = lo[1]; y <= hi[1]; ++y)
                    for (int z = lo[2]; z <= hi[2]; ++z) {
                        const Vec3 c{(x + 0.5) * box[0] / n[0], (y + 0.5) * box[1] / n[1], (z + 0.5) * box[2] / n[2]};
                        if (science::norm(science::sub(c, q)) > r) continue;
                        const int wx = ((x % n[0]) + n[0]) % n[0], wy = ((y % n[1]) + n[1]) % n[1], wz = ((z % n[2]) + n[2]) % n[2];
                        taken[(static_cast<std::size_t>(wx) * n[1] + wy) * n[2] + wz] = 1;
                    }
        }
        std::size_t open = 0, inside = 0;
        for (int x = 0; x < n[0]; ++x)
            for (int y = 0; y < n[1]; ++y)
                for (int z = 0; z < n[2]; ++z)
                    if (region.contains({(x + 0.5) * box[0] / n[0], (y + 0.5) * box[1] / n[1], (z + 0.5) * box[2] / n[2]})) {
                        ++inside;
                        if (!taken[(static_cast<std::size_t>(x) * n[1] + y) * n[2] + z]) ++open;
                    }
        const double free = inside ? region.volume(box) * static_cast<double>(open) / inside : 0.0;
        const double massToPlace = std::max(0.0, density * free / 1.66053906660 - fixedMass);
        double perSet = 0;
        for (const Kind& k : kinds)
            if (k.count == 0) perSet += moleculeMass(k.atoms);
        const int sets = static_cast<int>(std::floor(massToPlace / perSet));
        for (Kind& k : kinds)
            if (k.count == 0) k.count = sets;
    }
    const double minimum = number(p, id, "min_distance"), soluteMinimum = number(p, id, "solute_distance");
    const int attempts = integerOf(p, id, "attempts", 1);
    const bool fixedOrientation = flag(p, id, "fixed_orientation");
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double reach = std::max({minimum, soluteMinimum, 0.5});
    const int g[3] = {std::max(1, static_cast<int>(box[0] / reach)), std::max(1, static_cast<int>(box[1] / reach)), std::max(1, static_cast<int>(box[2] / reach))};
    std::unordered_map<long long, std::vector<std::pair<Vec3, bool>>> grid;  // position, solute
    const auto key = [&](int x, int y, int z) {
        x = ((x % g[0]) + g[0]) % g[0]; y = ((y % g[1]) + g[1]) % g[1]; z = ((z % g[2]) + g[2]) % g[2];
        return (static_cast<long long>(x) * g[1] + y) * g[2] + z;
    };
    const auto wrapped = [&](Vec3 q) {
        for (int k = 0; k < 3; ++k) q[k] -= box[k] * std::floor(q[k] / box[k]);
        return q;
    };
    const auto cellOfPoint = [&](const Vec3& q, int k) { return std::min(g[k] - 1, static_cast<int>(q[k] / box[k] * g[k])); };
    const auto insert = [&](const Vec3& q, bool solute) {
        const Vec3 w = wrapped(q);
        grid[key(cellOfPoint(w, 0), cellOfPoint(w, 1), cellOfPoint(w, 2))].push_back({w, solute});
    };
    const auto clear = [&](const Vec3& q) {
        const int cx = cellOfPoint(q, 0), cy = cellOfPoint(q, 1), cz = cellOfPoint(q, 2);
        const int span = static_cast<int>(std::ceil(std::max(minimum, soluteMinimum) / (box[0] / g[0]))) + 1;
        for (int dx = -span; dx <= span; ++dx)
            for (int dy = -span; dy <= span; ++dy)
                for (int dz = -span; dz <= span; ++dz) {
                    const auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == grid.end()) continue;
                    for (const auto& [r, solute] : it->second) {
                        Vec3 v = science::sub(q, r);
                        for (int k = 0; k < 3; ++k) v[k] -= box[k] * std::round(v[k] / box[k]);
                        if (science::norm(v) < (solute ? soluteMinimum : minimum)) return false;
                    }
                }
        return true;
    };
    for (const auto& a : s.atoms) insert({a.x, a.y, a.z}, true);
    // Larger molecules first: they are harder to place.
    std::stable_sort(kinds.begin(), kinds.end(), [](const Kind& a, const Kind& b) { return a.atoms.size() > b.atoms.size(); });
    const std::size_t soluteAtoms = s.atoms.size();
    std::ostringstream summary;
    double placedMass = 0;
    for (Kind& kind : kinds) {
        for (int c = 0; c < kind.count; ++c) {
            bool placed = false;
            for (int attempt = 0; attempt < attempts && !placed; ++attempt) {
                const Mat3 r = fixedOrientation ? science::identity() : randomRotation(rng);
                const Vec3 centre{u(rng) * box[0], u(rng) * box[1], u(rng) * box[2]};
                std::vector<Vec3> trial;
                bool ok = true;
                for (const auto& a : kind.atoms) {
                    const Vec3 raw = science::add(centre, apply(r, a.position));
                    if (!region.contains(raw) || !clear(wrapped(raw))) { ok = false; break; }
                    trial.push_back(raw);
                }
                if (!ok) continue;
                for (std::size_t k = 0; k < trial.size(); ++k) {
                    // Molecules stay whole: atoms are written unwrapped around their centre.
                    addAtom(s, kind.atoms[k].symbol, trial[k]);
                    insert(trial[k], false);
                }
                placed = true;
            }
            if (!placed && kind.automatic) {
                // Filling to the density is best effort: random insertion jams short of liquid densities.
                data.notes.push_back("Placed " + std::to_string(c) + " of the " + std::to_string(kind.count) + " " + kind.name +
                                     " molecules the density asks for (random insertion jammed; lower min_distance or relax the result)");
                kind.count = c;
                break;
            }
            if (!placed)
                throw std::runtime_error("Placed " + std::to_string(c) + " of " + std::to_string(kind.count) + " " + kind.name +
                                         " molecules: lower the density or min_distance, enlarge the region, or raise attempts");
        }
        placedMass += moleculeMass(kind.atoms) * kind.count;
        summary << (summary.tellp() ? ", " : "") << kind.count << " " << kind.name;
    }
    if (!solvate) setCell(s, {box[0], 0, 0}, {0, box[1], 0}, {0, 0, box[2]}, {true, true, true});
    produce(data, std::move(s));
    if (solvate)
        for (std::size_t i = soluteAtoms; i < data.selected.size(); ++i) data.selected[i] = 1;
    std::ostringstream note;
    note.precision(5);
    note << summary.str() << (solvate ? " around the " + std::to_string(soluteAtoms) + "-atom solute" : std::string()) << " in " << box[0] << " x "
         << box[1] << " x " << box[2] << " Angstrom";
    if (region.kind != "all") note << ", inside the " << region.kind << " region";
    note << " (" << placedMass * 1.66053906660 / region.volume(box) << " g/cm3 of molecules in the region)";
    data.notes.push_back(note.str());
}

// ---------------------------------------------------------------- adsorbates
void addAdsorbate(const Json& p, PipelineData& data)
{
    const char* id = "add-adsorbate";
    requireCell(data, "add-adsorbate");
    Structure& s = data.structure;
    if (s.atoms.empty()) throw std::runtime_error("add-adsorbate needs a slab");
    auto species = molecule(text(p, id, "species"));
    const std::string site = text(p, id, "site");
    if (site != "top" && site != "bridge" && site != "fcc" && site != "hcp") throw std::runtime_error("site must be top, bridge, fcc or hcp");
    const double height = number(p, id, "height");
    const Vec3 target = vector(p, id, "position");  // fractional (u, v) in the surface cell
    const std::string coverage = text(p, id, "coverage"), side = text(p, id, "surface");
    if (coverage != "one" && coverage != "all" && coverage != "fraction") throw std::runtime_error("coverage must be one, all or fraction");
    if (side != "top" && side != "bottom" && side != "both") throw std::runtime_error("surface must be top, bottom or both");
    const int anchorIndex = integerOf(p, id, "anchor", -1);
    if (anchorIndex >= static_cast<int>(species.size())) throw std::runtime_error("anchor: the molecule has " + std::to_string(species.size()) + " atoms");
    // Orientation in the molecule's frame (z: away from the surface): optional principal axes, tilt about y, then azimuth about z.
    if (flag(p, id, "principal") && species.size() > 1) species = rotated(species, principalFrame(species));
    const double tilt = number(p, id, "tilt") * kPi / 180, azimuth = number(p, id, "azimuth") * kPi / 180;
    species = rotated(species, rotationAbout({0, 1, 0}, tilt));
    species = rotated(species, rotationAbout({0, 0, 1}, azimuth));
    const Mat3 cell = cellOf(s);
    Vec3 up = science::cross(cell[0], cell[1]);
    up = science::scale(up, 1 / science::norm(up));
    if (science::dot(up, cell[2]) < 0) up = science::scale(up, -1);
    const Mat3 inverse = science::inverse(cell);
    // In-plane displacements by the minimum image along a and b.
    const auto inPlane = [&](Vec3 v) {
        v = science::sub(v, science::scale(up, science::dot(v, up)));
        Vec3 f = science::rowTimes(v, inverse);
        for (int k = 0; k < 2; ++k) f[k] -= std::round(f[k]);
        f[2] = 0;
        return science::rowTimes(f, cell);
    };
    const double layerTolerance = 0.6;
    struct Placement { Vec3 site; Vec3 normal; double surface; };
    std::vector<Placement> placements;
    std::size_t siteCount = 0;
    for (const double sign : {1.0, -1.0}) {
        if ((sign > 0 && side == "bottom") || (sign < 0 && side == "top")) continue;
        const Vec3 normal = science::scale(up, sign);
        // Layers along the outward normal.
        double top = -HUGE_VAL;
        for (const auto& a : s.atoms) top = std::max(top, science::dot(positionOf(a), normal));
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
        // Nearest surface neighbour distance.
        double nearest = HUGE_VAL;
        for (std::size_t x = 0; x < surface.size(); ++x)
            for (std::size_t y = 0; y < surface.size(); ++y)
                for (int i = -1; i <= 1; ++i)
                    for (int j = -1; j <= 1; ++j) {
                        Vec3 v = science::sub(positionOf(s.atoms[surface[y]]), positionOf(s.atoms[surface[x]]));
                        v = science::add(v, science::add(science::scale(cell[0], i), science::scale(cell[1], j)));
                        v = science::sub(v, science::scale(up, science::dot(v, up)));
                        if (science::norm(v) > 1e-6) nearest = std::min(nearest, science::norm(v));
                    }
        if (site != "top" && !(nearest < HUGE_VAL)) throw std::runtime_error("A " + site + " site needs a periodic surface");
        std::vector<Vec3> candidates;  // in-plane site positions
        const auto addCandidate = [&](const Vec3& q) {
            for (const Vec3& c : candidates)
                if (science::norm(inPlane(science::sub(c, q))) < 0.1) return;
            candidates.push_back(q);
        };
        for (std::size_t x = 0; x < surface.size(); ++x) {
            const Vec3 px = positionOf(s.atoms[surface[x]]);
            if (site == "top") { addCandidate(px); continue; }
            std::vector<Vec3> near;  // vectors to neighbouring surface atoms (all images in a 3 x 3 block)
            for (std::size_t y = 0; y < surface.size(); ++y)
                for (int i = -1; i <= 1; ++i)
                    for (int j = -1; j <= 1; ++j) {
                        Vec3 v = science::sub(positionOf(s.atoms[surface[y]]), px);
                        v = science::add(v, science::add(science::scale(cell[0], i), science::scale(cell[1], j)));
                        v = science::sub(v, science::scale(up, science::dot(v, up)));
                        const double len = science::norm(v);
                        if (len > 1e-6 && len < 1.2 * nearest) near.push_back(v);
                    }
            if (site == "bridge") {
                for (const Vec3& v : near) addCandidate(science::add(px, science::scale(v, 0.5)));
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
                    if ((site == "hcp") == above) addCandidate(centroid);
                }
        }
        if (candidates.empty()) throw std::runtime_error("No " + site + " site found on the " + (sign > 0 ? "top" : "bottom") + " surface");
        siteCount += candidates.size();
        if (coverage == "one") {
            const Vec3 wanted = science::rowTimes({target[0], target[1], 0}, cell);
            std::size_t best = 0;
            double bestDistance = HUGE_VAL;
            for (std::size_t k = 0; k < candidates.size(); ++k) {
                const double dist = science::norm(inPlane(science::sub(candidates[k], wanted)));
                if (dist < bestDistance - 1e-9) { bestDistance = dist; best = k; }
            }
            candidates = {candidates[best]};
        } else if (coverage == "fraction") {
            // A random subset, keeping adsorbates at least min_separation apart (in the plane, periodically).
            const double fraction = number(p, id, "fraction"), separation = number(p, id, "min_separation");
            if (fraction < 0 || fraction > 1) throw std::runtime_error("fraction must lie between 0 and 1");
            std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")) + (sign > 0 ? 0 : 7919));
            std::shuffle(candidates.begin(), candidates.end(), rng);
            const std::size_t wanted = static_cast<std::size_t>(std::lround(fraction * candidates.size()));
            std::vector<Vec3> chosen;
            for (const Vec3& c : candidates) {
                if (chosen.size() == wanted) break;
                bool ok = true;
                for (const Vec3& o : chosen) ok = ok && science::norm(inPlane(science::sub(c, o))) >= separation;
                if (ok) chosen.push_back(c);
            }
            if (chosen.size() < wanted)
                data.notes.push_back("Only " + std::to_string(chosen.size()) + " of " + std::to_string(wanted) +
                                     " sites could be covered at the minimum separation");
            candidates = chosen;
        }
        for (const Vec3& c : candidates) placements.push_back({c, normal, top});
    }
    if (data.selected.size() != s.atoms.size()) data.selected.assign(s.atoms.size(), 0);
    std::fill(data.selected.begin(), data.selected.end(), 0);
    for (const Placement& place : placements) {
        Vec3 anchor = science::add(science::sub(place.site, science::scale(place.normal, science::dot(place.site, place.normal))),
                                   science::scale(place.normal, place.surface + height));
        // The molecule's z axis along the outward normal; the anchor atom (or the lowest atom) at the height.
        const Mat3 r = rotationBetween({0, 0, 1}, place.normal);
        Vec3 reference{0, 0, 0};
        if (anchorIndex >= 0) {
            reference = rotate(r, species[static_cast<std::size_t>(anchorIndex)].position);
        } else {
            double lowest = HUGE_VAL;
            for (const auto& a : species) lowest = std::min(lowest, science::dot(apply(r, a.position), place.normal));
            reference = science::scale(place.normal, lowest);
        }
        for (const auto& a : species)
            appendAtom(s, data.selected, a.symbol, science::add(anchor, science::sub(apply(r, a.position), reference)), flag(p, id, "select"));
    }
    // Adsorbates closer than the slab's own spacing to one another or to the slab are worth a warning.
    std::ostringstream note;
    note.precision(5);
    note << placements.size() << " " << text(p, id, "species") << " on " << site << " sites (" << siteCount << " available), " << height
         << " Angstrom above the " << (side == "both" ? "two surfaces" : side + " surface");
    if (coverage != "one") {
        const double area = science::norm(science::cross(cell[0], cell[1]));
        note << "; coverage " << placements.size() / (side == "both" ? 2.0 : 1.0) / area * 100 << " per nm2";
    }
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
    const double tolerance = number(p, id, "tolerance") / 100;
    if (tolerance < 0 || tolerance >= 1) throw std::runtime_error("tolerance is a percentage from 0 to 100");
    const std::string shape = text(p, id, "shape");
    // Cell counts allowed: the nearest to the target, or all within the tolerance.
    const long nearest = std::max(1L, std::lround(double(target) / double(s.atoms.size())));
    long lowest = nearest, highest = nearest;
    if (tolerance > 0) {
        lowest = std::max(1L, static_cast<long>(std::ceil(target * (1 - tolerance) / s.atoms.size())));
        highest = std::max(lowest, static_cast<long>(std::floor(target * (1 + tolerance) / s.atoms.size())));
    }
    const Mat3 h = cellOf(s);
    const double cellVolume = std::abs(science::determinant(h));
    // The ideal metric of each shape (sc: a cube; fcc and bcc: the primitive cells of those lattices).
    double ideal[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (shape == "fcc") {
        const double m[3][3] = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
        std::copy(&m[0][0], &m[0][0] + 9, &ideal[0][0]);
    } else if (shape == "bcc") {
        const double m[3][3] = {{3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}};
        std::copy(&m[0][0], &m[0][0] + 9, &ideal[0][0]);
    } else if (shape != "sc" && shape != "image")
        throw std::runtime_error("shape is sc, fcc, bcc or image");
    const double idealDet = ideal[0][0] * (ideal[1][1] * ideal[2][2] - ideal[1][2] * ideal[2][1]) -
                            ideal[0][1] * (ideal[1][0] * ideal[2][2] - ideal[1][2] * ideal[2][0]) +
                            ideal[0][2] * (ideal[1][0] * ideal[2][1] - ideal[1][1] * ideal[2][0]);
    // Shortest lattice vector (combinations up to +-2 of the cell vectors).
    const auto shortest = [](const Mat3& q) {
        double best = HUGE_VAL;
        for (int i = -2; i <= 2; ++i)
            for (int j = -2; j <= 2; ++j)
                for (int k = -2; k <= 2; ++k) {
                    if (!i && !j && !k) continue;
                    const Vec3 v = science::add(science::add(science::scale(q[0], i), science::scale(q[1], j)), science::scale(q[2], k));
                    best = std::min(best, science::norm(v));
                }
        return best;
    };
    double bestScore = HUGE_VAL;
    int bestNegatives = 10, bestTrace = -100;
    long bestCells = 0;
    Mat3 best{}, bestCell{};
    std::array<int, 9> e{};
    const int span = 2 * range + 1;
    long total = 1;
    for (int k = 0; k < 9; ++k) total *= span;
    for (long code = 0; code < total; ++code) {
        long c = code;
        for (int k = 0; k < 9; ++k) { e[k] = static_cast<int>(c % span) - range; c /= span; }
        const long det = static_cast<long>(e[0]) * (e[4] * e[8] - e[5] * e[7]) - static_cast<long>(e[1]) * (e[3] * e[8] - e[5] * e[6]) +
                         static_cast<long>(e[2]) * (e[3] * e[7] - e[4] * e[6]);
        if (det < lowest || det > highest) continue;
        Mat3 q{};
        for (int r = 0; r < 3; ++r)
            for (int col = 0; col < 3; ++col) q[r][col] = e[3 * r] * h[0][col] + e[3 * r + 1] * h[1][col] + e[3 * r + 2] * h[2][col];
        const double volume = cellVolume * det;
        double score = 0;
        if (shape == "image") {
            // The largest shortest image distance, relative to that of an fcc cell of the volume (the densest packing).
            score = 1 - shortest(q) / std::cbrt(std::sqrt(2.0) * volume);
        } else {
            // Deviation of the metric P H (P H)^T from the ideal one of the same volume.
            const double scale = std::cbrt(volume * volume / idealDet);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    const double g = science::dot(q[i], q[j]) - ideal[i][j] * scale;
                    score += g * g;
                }
            score = std::sqrt(score) / std::cbrt(volume * volume);
        }
        // Ties: fewer cells off the target, then fewer negative entries, then a larger trace.
        score += 1e-6 * std::abs(det - double(target) / s.atoms.size());
        int negatives = 0, trace = e[0] + e[4] + e[8];
        for (int k = 0; k < 9; ++k) negatives += e[k] < 0;
        const bool better = score < bestScore - 1e-9 ||
                            (score < bestScore + 1e-9 && (negatives < bestNegatives || (negatives == bestNegatives && trace > bestTrace)));
        if (better) {
            bestNegatives = negatives;
            bestTrace = trace;
            bestScore = score;
            bestCells = det;
            bestCell = q;
            for (int r = 0; r < 3; ++r)
                for (int col = 0; col < 3; ++col) best[r][col] = e[3 * r + col];
        }
    }
    if (!(bestScore < HUGE_VAL))
        throw std::runtime_error("No supercell of " + std::to_string(lowest) + (highest > lowest ? " to " + std::to_string(highest) : std::string()) +
                                 " cells within the search range");
    std::ostringstream matrix;
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 3; ++col) matrix << (r || col ? " " : "") << best[r][col];
    const ModifierType* supercell = findModifierType("supercell");
    Json q = Json::object();
    q["matrix"] = matrix.str();
    supercell->apply(q, data);
    // The inscribed sphere: the smallest distance between opposite faces.
    const double volume = std::abs(science::determinant(bestCell));
    double height = HUGE_VAL;
    for (int k = 0; k < 3; ++k) height = std::min(height, volume / science::norm(science::cross(bestCell[(k + 1) % 3], bestCell[(k + 2) % 3])));
    std::ostringstream note;
    note.precision(5);
    note << "Supercell matrix " << matrix.str() << " (" << bestCells << " cells, " << data.structure.atoms.size() << " atoms); "
         << (shape == "image" ? "shortest image distance " : shape == "sc" ? std::string("deviation from a cube ") : "deviation from an " + shape + " cell ")
         << (shape == "image" ? shortest(bestCell) : bestScore) << "; shortest image " << shortest(bestCell) << " Angstrom, inscribed sphere diameter "
         << height << " Angstrom";
    data.notes.insert(data.notes.begin(), note.str());
}

// ---------------------------------------------------------------- interpolation
void interpolate(const Json& p, PipelineData& data)
{
    const char* id = "interpolate";
    const std::string file = text(p, id, "file");
    if (file.empty()) throw std::runtime_error("Choose the final structure file");
    const double t = number(p, id, "fraction");
    const std::string method = text(p, id, "method");
    if (method != "linear" && method != "idpp") throw std::runtime_error("method is linear or idpp");
    const auto frames = science::readFrames(std::filesystem::u8path(file));
    const Structure& end = frames.front().structure;
    Structure& s = data.structure;
    if (end.atoms.size() != s.atoms.size())
        throw std::runtime_error("The final structure has " + std::to_string(end.atoms.size()) + " atoms, this one " + std::to_string(s.atoms.size()));
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (end.atoms[i].atomicNumber != s.atoms[i].atomicNumber)
            throw std::runtime_error("Atom " + std::to_string(i) + " is " + s.atoms[i].symbol + " here but " + end.atoms[i].symbol + " in the final structure");
    const Structure start = s;
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        Vec3 d = science::sub(positionOf(end.atoms[i]), positionOf(s.atoms[i]));
        if (s.hasUnitCell) d = minimumImage(s, d);  // atoms that crossed a periodic boundary move the short way
        place(s.atoms[i], science::add(positionOf(s.atoms[i]), science::scale(d, t)));
    }
    if (s.hasUnitCell && end.hasUnitCell)
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) s.cellVectors[r][c] += t * (end.cellVectors[r][c] - s.cellVectors[r][c]);
    std::ostringstream note;
    note.precision(5);
    note << (method == "idpp" ? "IDPP" : "Linear") << " interpolation " << t << " of the way to " << std::filesystem::u8path(file).filename().u8string();
    if (method == "idpp") {
        // Image-dependent pair potential (Smidstrup et al. 2014): from the linear guess, minimise
        // sum w(d) (d_ij - d_ij(t))^2, w = 1/d^4, with d_ij(t) interpolating the end points' pair distances.
        const std::size_t n = s.atoms.size();
        if (n > 3000) throw std::runtime_error("IDPP: at most 3000 atoms");
        const double cutoff = number(p, id, "idpp_cutoff");
        const auto image = [&](const Structure& at, Vec3 d) { return at.hasUnitCell ? minimumImage(at, d) : d; };
        struct Pair { std::size_t i, j; double target; };
        std::vector<Pair> pairs;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) {
                const double d0 = science::norm(image(start, science::sub(positionOf(start.atoms[j]), positionOf(start.atoms[i]))));
                const double d1 = science::norm(image(end, science::sub(positionOf(end.atoms[j]), positionOf(end.atoms[i]))));
                if (cutoff > 0 && std::min(d0, d1) > cutoff) continue;
                pairs.push_back({i, j, (1 - t) * d0 + t * d1});
            }
        std::vector<Vec3> x(n);
        for (std::size_t i = 0; i < n; ++i) x[i] = positionOf(s.atoms[i]);
        const auto evaluate = [&](const std::vector<Vec3>& q, std::vector<Vec3>* force) {
            double objective = 0;
            if (force) force->assign(n, Vec3{0, 0, 0});
            for (const Pair& pr : pairs) {
                const Vec3 v = image(s, science::sub(q[pr.j], q[pr.i]));
                const double d = science::norm(v);
                if (d < 1e-8) continue;
                const double w = 1 / (d * d * d * d), diff = d - pr.target;
                objective += w * diff * diff;
                if (force) {
                    // -d/dd [w (d - target)^2] = -(2 w diff - 4 w diff^2 / d), along v/d for atom j.
                    const double g = 2 * w * diff - 4 * w * diff * diff / d;
                    const Vec3 f = science::scale(v, -g / d);
                    (*force)[pr.j] = science::add((*force)[pr.j], f);
                    (*force)[pr.i] = science::sub((*force)[pr.i], f);
                }
            }
            return objective;
        };
        // FIRE minimisation.
        std::vector<Vec3> velocity(n, Vec3{0, 0, 0}), force;
        double dt = 0.1, alpha = 0.1, before = evaluate(x, &force), now = before;
        int positive = 0, steps = 0;
        const int maxSteps = integerOf(p, id, "idpp_steps", 1);
        for (; steps < maxSteps; ++steps) {
            double power = 0, fnorm = 0, vnorm = 0, fmax = 0;
            for (std::size_t i = 0; i < n; ++i) {
                power += science::dot(force[i], velocity[i]);
                fnorm += science::dot(force[i], force[i]);
                vnorm += science::dot(velocity[i], velocity[i]);
                fmax = std::max(fmax, science::norm(force[i]));
            }
            if (fmax < 1e-5) break;
            fnorm = std::sqrt(fnorm);
            vnorm = std::sqrt(vnorm);
            for (std::size_t i = 0; i < n; ++i)
                velocity[i] = science::add(science::scale(velocity[i], 1 - alpha), science::scale(force[i], fnorm > 0 ? alpha * vnorm / fnorm : 0));
            if (power > 0) {
                if (++positive > 5) { dt = std::min(dt * 1.1, 1.0); alpha *= 0.99; }
            } else {
                positive = 0;
                dt *= 0.5;
                alpha = 0.1;
                std::fill(velocity.begin(), velocity.end(), Vec3{0, 0, 0});
            }
            for (std::size_t i = 0; i < n; ++i) {
                velocity[i] = science::add(velocity[i], science::scale(force[i], dt));
                Vec3 step = science::scale(velocity[i], dt);
                const double len = science::norm(step);
                if (len > 0.1) step = science::scale(step, 0.1 / len);  // at most 0.1 Angstrom per step
                x[i] = science::add(x[i], step);
            }
            now = evaluate(x, &force);
        }
        for (std::size_t i = 0; i < n; ++i) place(s.atoms[i], x[i]);
        note << "; objective " << before << " -> " << now << " in " << steps << " steps (" << pairs.size() << " pairs)";
    }
    data.notes.push_back(note.str());
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
    auto first = atomsOf(s, pair[0]), second = atomsOf(s, pair[1]);
    if (first.empty() || second.empty()) throw std::runtime_error("create-antisite: the structure needs both " + pair[0] + " and " + pair[1]);
    const int count = integerOf(p, id, "count", 1);
    const double separation = number(p, id, "min_separation");
    const bool nearestPartner = text(p, id, "partner") == "nearest";
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    const auto distance = [&](std::size_t i, std::size_t j) {
        Vec3 d = science::sub(positionOf(s.atoms[j]), positionOf(s.atoms[i]));
        if (s.hasUnitCell) d = minimumImage(s, d);
        return science::norm(d);
    };
    std::shuffle(first.begin(), first.end(), rng);
    std::vector<std::size_t> used;  // atoms already swapped
    const auto farFromUsed = [&](std::size_t i) {
        for (std::size_t u : used)
            if (i == u || distance(i, u) < separation) return false;
        return true;
    };
    std::ostringstream swaps;
    int made = 0;
    for (std::size_t a : first) {
        if (made == count) break;
        if (!farFromUsed(a)) continue;
        std::vector<std::size_t> partners;
        for (std::size_t j : second)
            if (farFromUsed(j)) partners.push_back(j);
        if (partners.empty()) continue;
        std::size_t b = partners[std::uniform_int_distribution<std::size_t>(0, partners.size() - 1)(rng)];
        if (nearestPartner) {
            double best = HUGE_VAL;
            for (std::size_t j : partners)
                if (distance(a, j) < best) { best = distance(a, j); b = j; }
        }
        std::swap(s.atoms[a].symbol, s.atoms[b].symbol);
        std::swap(s.atoms[a].atomicNumber, s.atoms[b].atomicNumber);
        getDefaultElementColor(s.atoms[a].atomicNumber, s.atoms[a].r, s.atoms[a].g, s.atoms[a].b);
        getDefaultElementColor(s.atoms[b].atomicNumber, s.atoms[b].r, s.atoms[b].g, s.atoms[b].b);
        used.push_back(a);
        used.push_back(b);
        swaps << (made ? ", " : "") << a << "<->" << b;
        ++made;
    }
    if (made == 0) throw std::runtime_error("create-antisite: no pair satisfies the minimum separation");
    data.selected.assign(s.atoms.size(), 0);
    for (std::size_t u : used) data.selected[u] = 1;
    std::string note = std::to_string(made) + " antisite pair(s) " + pair[0] + "<->" + pair[1] + ": atoms " + swaps.str();
    if (made < count) note += " (only " + std::to_string(made) + " of " + std::to_string(count) + " fit at the minimum separation)";
    data.notes.push_back(note);
}

void frenkelPair(const Json& p, PipelineData& data)
{
    const char* id = "create-frenkel-pair";
    const std::string symbol = standardElement(text(p, id, "element"));
    Structure& s = data.structure;
    if (atomsOf(s, symbol).empty()) throw std::runtime_error("create-frenkel-pair: no " + symbol + " atom");
    const std::string kind = text(p, id, "interstitial");
    if (kind != "void" && kind != "dumbbell") throw std::runtime_error("interstitial is void or dumbbell");
    const double minimum = number(p, id, "min_distance"), maximum = number(p, id, "max_distance");
    const std::size_t before = s.atoms.size();
    std::mt19937_64 rng(static_cast<std::uint64_t>(number(p, id, "seed")));
    const auto separation = [&](std::size_t i, const Vec3& q) {
        Vec3 d = science::sub(positionOf(s.atoms[i]), q);
        if (s.hasUnitCell) d = minimumImage(s, d);
        return science::norm(d);
    };
    std::vector<std::size_t> selected;
    Vec3 interstitial{};
    std::size_t host = before;  // the dumbbell's lattice site
    if (kind == "void") {
        // The interstitial first (in the perfect structure, so the vacancy is not the void chosen).
        const ModifierType* insert = findModifierType("insert-interstitials");
        Json q = Json::object();
        q["element"] = symbol;
        q["count"] = 1.0;
        q["order"] = "random";
        q["seed"] = number(p, id, "seed");
        insert->apply(q, data);
        if (s.atoms.size() != before + 1) throw std::runtime_error("create-frenkel-pair: no interstitial void found");
        interstitial = positionOf(s.atoms.back());
    } else {
        // A split (dumbbell) interstitial: a lattice atom replaced by two atoms along a direction, centred on its site.
        const auto sites = atomsOf(s, symbol);
        host = sites[std::uniform_int_distribution<std::size_t>(0, sites.size() - 1)(rng)];
        interstitial = positionOf(s.atoms[host]);
    }
    // The vacancy: an atom of the element between min_distance and max_distance away (the farthest when none is).
    std::vector<std::size_t> candidates;
    std::size_t farthest = before;
    double farthestDistance = -1;
    for (std::size_t i = 0; i < before; ++i) {
        if (s.atoms[i].symbol != symbol || i == host) continue;
        const double dist = separation(i, interstitial);
        if (dist >= minimum && (maximum <= 0 || dist <= maximum)) candidates.push_back(i);
        if (dist > farthestDistance && (maximum <= 0 || dist <= maximum)) { farthestDistance = dist; farthest = i; }
    }
    if (candidates.empty() && farthest == before) throw std::runtime_error("create-frenkel-pair: no " + symbol + " atom within max_distance");
    std::mt19937_64 pick(static_cast<std::uint64_t>(number(p, id, "seed")) + 7);
    const std::size_t vacancy = candidates.empty() ? farthest : candidates[std::uniform_int_distribution<std::size_t>(0, candidates.size() - 1)(pick)];
    const double distance = separation(vacancy, interstitial);
    std::ostringstream note;
    note.precision(5);
    if (kind == "dumbbell") {
        Vec3 direction = vector(p, id, "direction");
        if (flag(p, id, "lattice_direction")) {
            requireCell(data, "A lattice direction");
            direction = science::rowTimes(direction, cellOf(s));
        }
        const double length = science::norm(direction);
        if (!(length > 1e-12)) throw std::runtime_error("The dumbbell direction must not be zero");
        direction = science::scale(direction, 1 / length);
        double bond = number(p, id, "dumbbell_length");
        if (!(bond > 0)) {
            // Default: 0.9 x the nearest-neighbour distance of the host site.
            double nearest = HUGE_VAL;
            for (std::size_t i = 0; i < before; ++i)
                if (i != host) nearest = std::min(nearest, separation(i, interstitial));
            bond = 0.9 * nearest;
        }
        place(s.atoms[host], science::sub(interstitial, science::scale(direction, bond / 2)));
        if (data.selected.size() != s.atoms.size()) data.selected.assign(s.atoms.size(), 0);
        appendAtom(s, data.selected, symbol, science::add(interstitial, science::scale(direction, bond / 2)), false);
        note << "Frenkel pair: " << symbol << " dumbbell (length " << bond << " Angstrom) at atom " << host << ", vacancy at atom " << vacancy << " "
             << distance << " Angstrom away";
    } else
        note << "Frenkel pair: " << symbol << " atom " << vacancy << " moved to an interstitial void " << distance << " Angstrom away";
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if (i != vacancy) keep.push_back(i);
    rebuild(data, keep);
    data.selected.assign(data.structure.atoms.size(), 0);
    data.selected.back() = 1;
    if (kind == "dumbbell") data.selected[host > vacancy ? host - 1 : host] = 1;
    data.notes.push_back(note.str());
}
}  // namespace

std::vector<ModifierType> generatorStepTypes()
{
    static const std::string layerChoices = "graphene|hbn|silicene|germanene|stanene|mx2-2h|mx2-1t|phosphorene";
    static const std::string hexagonalChoices = "graphene|hbn|silicene|germanene|stanene|mx2-2h|mx2-1t";
    static const std::string topChoices = "same|" + hexagonalChoices;
    const char* const layers = layerChoices.c_str();
    return {
        {"build-nanotube", "Nanotube", "Build",
         "Rolls a honeycomb sheet into a tube with chiral indices (n, m): armchair (n, n), zigzag (n, 0) or chiral. Two elements make a BN-type tube "
         "(alternating sublattices). Multi-wall tubes (armchair or zigzag) add walls about wall_spacing apart. Periodic along the tube axis, or "
         "finite with the ends capped (H, F, ...). The notes give diameters, the chiral angle and, for carbon, metallic or semiconducting.",
         {{"n", "Chiral index n", "int", "6", ""}, {"m", "Chiral index m", "int", "6", ""}, {"length", "Unit cells along the axis", "int", "1", ""},
          {"length_A", "Length in Angstrom (0 = use the unit cells)", "float", "0", ""},
          {"elements", "Elements (one, or two for BN)", "elements", "C", ""}, {"bond", "Bond length (0 = C-C 1.421, B-N 1.446)", "float", "0", ""},
          {"walls", "Walls", "int", "1", ""}, {"wall_spacing", "Wall spacing (Angstrom)", "float", "3.4", ""},
          {"finite", "Finite tube (open ends)", "bool", "false", ""}, {"end_cap", "End atoms of a finite tube", "choice", "H", "H|F|Cl|none"},
          {"axis", "Tube axis", "choice", "z", "x|y|z"}, {"vacuum", "Vacuum around the tube (Angstrom)", "float", "10", ""}}, buildNanotube},
        {"build-monolayer", "2D layer", "Build",
         "A 2D layer or few-layer stack: graphene, hexagonal BN, buckled silicene, germanene and stanene, 2H or 1T MX2 dichalcogenides, or "
         "phosphorene; na x nb primitive cells (or the rectangular a x a sqrt3 cell of hexagonal layers), biaxial strain in %, and stacking of "
         "several layers (AA, AB Bernal, ABC rhombohedral, or 2H for MX2). a = 0 and thickness = 0 take the defaults.",
         {{"kind", "Layer", "choice", "graphene", layers}, {"elements", "Elements (MX2: M X)", "elements", "", ""},
          {"a", "Lattice constant a (0 = default)", "float", "0", ""},
          {"thickness", "MX2 X-X thickness or X-ene buckling (0 = default)", "float", "0", ""},
          {"strain", "Biaxial strain (%)", "float", "0", ""}, {"na", "Cells along a", "int", "1", ""}, {"nb", "Cells along b", "int", "1", ""},
          {"rectangular", "Rectangular cell (hexagonal layers)", "bool", "false", ""}, {"layers", "Layers", "int", "1", ""},
          {"stacking", "Stacking", "choice", "AB", "AA|AB|ABC|2H"}, {"interlayer", "Interlayer distance (0 = default)", "float", "0", ""},
          {"periodic_c", "Periodic along c (bulk-like stack)", "bool", "true", ""}, {"vacuum", "Vacuum (Angstrom)", "float", "15", ""}}, buildMonolayer},
        {"build-nanoribbon", "Nanoribbon", "Build",
         "A graphene (or BN) ribbon with armchair edges (width in dimer lines, N-AGNR) or zigzag edges (width in zigzag chains, N-ZGNR), "
         "periodic along x or a finite flake, with H, F or Cl on the edge atoms.",
         {{"edge", "Edge", "choice", "armchair", "armchair|zigzag"}, {"width", "Width (dimer lines or chains)", "int", "7", ""},
          {"length", "Periods along the ribbon", "int", "1", ""}, {"elements", "Elements (one, or two for BN)", "elements", "C", ""},
          {"bond", "Bond length (Angstrom)", "float", "1.421", ""}, {"passivate", "Terminate the edges", "bool", "true", ""},
          {"termination", "Edge atoms", "choice", "H", "H|F|Cl|none"}, {"termination_bond", "Edge bond length (0 = default)", "float", "0", ""},
          {"finite", "Finite flake (open ends)", "bool", "false", ""}, {"vacuum", "Vacuum (Angstrom)", "float", "10", ""}}, buildNanoribbon},
        {"build-twisted-bilayer", "Twisted bilayer", "Build",
         "Two hexagonal layers turned by a commensurate angle (cos t = (3i^2+3i+1/2)/(3i^2+3i+1): 21.79, 13.17, 9.43, ... degrees), in their "
         "moire supercell. Give the angle (the nearest commensurate one is used) or the index i. The top layer may be another layer, strained "
         "to the bottom lattice (a heterobilayer, the strain reported), and start in AA or AB registry.",
         {{"kind", "Layer", "choice", "graphene", layers}, {"elements", "Elements (MX2: M X)", "elements", "", ""},
          {"angle", "Twist angle (degrees; 0 = use the index)", "float", "0", ""}, {"index", "Commensurate index i", "int", "3", ""},
          {"a", "Lattice constant a (0 = default)", "float", "0", ""}, {"thickness", "MX2 X-X thickness (0 = 3.19)", "float", "0", ""},
          {"top_kind", "Top layer", "choice", "same", topChoices.c_str()},
          {"top_elements", "Top layer elements", "elements", "", ""}, {"stacking", "Registry before the twist", "choice", "AA", "AA|AB"},
          {"interlayer", "Interlayer distance (0 = default)", "float", "0", ""}, {"vacuum", "Vacuum (Angstrom)", "float", "15", ""}}, buildTwistedBilayer},
        {"build-cluster", "Magic-number cluster", "Build",
         "Closed-shell clusters: Mackay icosahedra and fcc cuboctahedra (1, 13, 55, 147, ... atoms for 1, 2, 3, 4 shells), "
         "regular decahedra (pentagonal bipyramids of five tetrahedra: 7, 23, 54, ... atoms for 2, 3, 4 atoms per edge) and fcc octahedra "
        "(shells = atoms per edge: 6, 19, 44, 85 atoms), "
         "with the nearest-neighbour distance along the edges. Nanoalloys: a fraction of the atoms becomes a second element, at random, "
         "in the core, in the shell (outermost atoms first) or on one side (Janus, along z).",
         {{"kind", "Shape", "choice", "icosahedron", "icosahedron|decahedron|cuboctahedron|octahedron"}, {"element", "Element", "element", "Cu", ""},
          {"shells", "Shells (1 = one atom)", "int", "3", ""}, {"distance", "Nearest-neighbour distance (Angstrom)", "float", "2.556", ""},
          {"second_element", "Alloying element (empty = none)", "string", "", ""}, {"fraction", "Alloying fraction", "float", "0.5", ""},
          {"distribution", "Alloying element in", "choice", "random", "random|core|shell|janus"}, {"seed", "Seed (random)", "int", "1", ""},
          {"vacuum", "Vacuum (Angstrom)", "float", "10", ""}}, buildCluster},
        {"build-molecule", "Molecule", "Build",
         "A molecule from the library (experimental geometries): H2, N2, O2, F2, Cl2, CO, NO, HF, HCl, H2O, H2S, CO2, SO2, NO2, O3, NH3, PH3, "
         "CH4, SiH4, CF4, C2H2, C2H4, C2H6, CH3OH, HCN, H2CO, C6H6, C60 (or common names: water, benzene, ...), or any structure file "
         "(XYZ, POSCAR, CIF, ...). Orientation: as in the library or file, along the principal axes of inertia (long axis on z, a planar "
         "molecule in the xz plane), or random. With vacuum, in a box (periodic or open).",
         {{"name", "Molecule (name or file)", "string", "H2O", ""}, {"orientation", "Orientation", "choice", "library", "library|principal|random"},
          {"seed", "Seed (random orientation)", "int", "1", ""}, {"vacuum", "Vacuum (Angstrom; 0 = no cell)", "float", "10", ""},
          {"periodic", "Periodic box", "bool", "false", ""}}, buildMolecule},
        {"build-packing", "Molecular packing", "Build",
         "Packs molecules (library names, single atoms or structure files) at random positions and orientations into a periodic box, "
         "keeping atoms of different molecules at least min_distance apart (as Packmol). The box follows from the density unless given. "
         "Solvate mode fills the input structure's orthogonal cell around it, keeping solute_distance from its atoms. A count of 0 fills "
         "the region to the density. Regions (all atoms of a molecule inside): box (low, high corners), sphere or outside-sphere (centre, "
         "radius), slab (low, high along the axis) or cylinder (centre, radius, along the axis).",
         {{"molecules", "Molecules and counts (e.g. H2O 100 CO2 10; 0 = to the density)", "list", "H2O 64", ""},
          {"mode", "Mode", "choice", "new", "new|solvate"},
          {"density", "Density (g/cm3)", "float", "1.0", ""}, {"box", "Box lengths (Angstrom; 0 0 0 = from the density)", "vector", "0 0 0", ""},
          {"min_distance", "Minimum distance between molecules (Angstrom)", "float", "2.0", ""},
          {"solute_distance", "Minimum distance to the solute (Angstrom)", "float", "2.5", ""},
          {"region", "Region", "choice", "all", "all|box|sphere|outside-sphere|slab|cylinder"},
          {"region_low", "Region low corner (Angstrom)", "vector", "0 0 0", ""}, {"region_high", "Region high corner (Angstrom)", "vector", "0 0 0", ""},
          {"region_center", "Region centre (Angstrom)", "vector", "0 0 0", ""}, {"region_radius", "Region radius (Angstrom)", "float", "0", ""},
          {"region_axis", "Region axis (slab, cylinder)", "choice", "z", "x|y|z"},
          {"fixed_orientation", "Keep the library orientation", "bool", "false", ""}, {"seed", "Seed", "int", "1", ""},
          {"attempts", "Attempts per molecule", "int", "20000", ""}}, buildPacking},
        {"add-adsorbate", "Add adsorbate", "Modification",
         "Places a molecule (library name, element or structure file) on a slab surface (normal a x b; top towards c, bottom, or both) at "
         "top, bridge, fcc-hollow or hcp-hollow sites, its anchor atom (or lowest atom) `height` above the outer layer. Coverage: one "
         "adsorbate at the site nearest the fractional position (u, v), all symmetry-distinct sites of the kind, or a random fraction of them "
         "kept min_separation apart. The molecule's z axis points away from the surface, after an optional principal-axes alignment, a tilt "
         "about y and an azimuthal turn about z.",
         {{"species", "Molecule, element or file", "string", "CO", ""}, {"site", "Site", "choice", "top", "top|bridge|fcc|hcp"},
          {"height", "Height above the surface (Angstrom)", "float", "1.8", ""}, {"position", "Near fractional position (u v, w ignored)", "vector", "0.5 0.5 0", ""},
          {"surface", "Surface", "choice", "top", "top|bottom|both"}, {"coverage", "Coverage", "choice", "one", "one|all|fraction"},
          {"fraction", "Fraction of the sites", "float", "0.5", ""}, {"min_separation", "Minimum adsorbate separation (Angstrom)", "float", "0", ""},
          {"seed", "Seed", "int", "1", ""}, {"anchor", "Anchor atom of the molecule (0-based; -1 = lowest)", "int", "-1", ""},
          {"principal", "Align the principal axes first", "bool", "false", ""}, {"tilt", "Tilt from the normal (degrees)", "float", "0", ""},
          {"azimuth", "Azimuth (degrees)", "float", "0", ""},
          {"select", "Select the adsorbate", "bool", "true", ""}}, addAdsorbate},
        {"near-cubic-supercell", "Near-cubic supercell", "Cell",
         "The supercell (integer matrix with entries up to +-range) of about target_atoms atoms (or any count within the tolerance) whose cell "
         "is closest to a cube, an fcc or bcc primitive cell (as ASE's find_optimal_cell_shape), or that has the largest shortest image "
         "distance, for point-defect calculations. The notes give the shortest image distance and the inscribed sphere.",
         {{"target_atoms", "Target number of atoms", "int", "64", ""}, {"shape", "Target shape", "choice", "sc", "sc|fcc|bcc|image"},
          {"tolerance", "Atom-count tolerance (%)", "float", "0", ""}, {"range", "Matrix entries up to +-", "int", "2", ""}}, nearCubicSupercell},
        {"interpolate", "Interpolate to a structure", "Modification",
         "Moves every atom the fraction t of the way to its position in another structure with the same atoms in the same order "
         "(minimum image across periodic boundaries); the cell is interpolated too. Images for a nudged elastic band are t = k/(n+1). "
         "IDPP (image-dependent pair potential, Smidstrup et al. 2014) then relaxes the linear guess towards interpolated interatomic "
         "distances, avoiding atoms passing too close (better NEB starting paths).",
         {{"file", "Final structure file", "string", "", ""}, {"fraction", "Fraction t (0 = this, 1 = final)", "float", "0.5", ""},
          {"method", "Method", "choice", "linear", "linear|idpp"},
          {"idpp_cutoff", "IDPP pair cutoff (Angstrom; 0 = all pairs)", "float", "0", ""}, {"idpp_steps", "IDPP steps", "int", "2000", ""}}, interpolate},
        {"create-antisite", "Antisite defect", "Modification",
         "Swaps atoms of the first element (chosen at random) with the nearest (or a random) atom of the second, count times, every swapped "
         "atom at least min_separation from the others (isolated defects); the swapped atoms are selected.",
         {{"elements", "Elements A B", "elements", "Ni Al", ""}, {"partner", "Partner", "choice", "nearest", "nearest|random"},
          {"count", "Antisite pairs", "int", "1", ""}, {"min_separation", "Minimum separation between defects (Angstrom)", "float", "0", ""},
          {"seed", "Seed", "int", "1", ""}}, antisite},
        {"create-frenkel-pair", "Frenkel pair", "Modification",
         "Moves an atom of the element into an interstitial void (as Insert interstitial atoms) or into a split (dumbbell) interstitial "
         "on a lattice site, e.g. <100> in fcc or <110> in bcc metals, leaving a vacancy between min_distance and max_distance away; the "
         "interstitial atoms are selected.",
         {{"element", "Element", "element", "Cu", ""}, {"interstitial", "Interstitial", "choice", "void", "void|dumbbell"},
          {"direction", "Dumbbell direction", "vector", "1 0 0", ""}, {"lattice_direction", "Lattice direction [u v w]", "bool", "true", ""},
          {"dumbbell_length", "Dumbbell length (Angstrom; 0 = 0.9 x nearest neighbour)", "float", "0", ""},
          {"min_distance", "Minimum vacancy-interstitial distance (Angstrom)", "float", "5", ""},
          {"max_distance", "Maximum vacancy-interstitial distance (Angstrom; 0 = any)", "float", "0", ""},
          {"seed", "Seed", "int", "1", ""}}, frenkelPair},
    };
}
}  // namespace atomforge::pipeline
