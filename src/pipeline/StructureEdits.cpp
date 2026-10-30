#include "pipeline/StructureEdits.h"

#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace atomforge::pipeline
{
using science::Vec3;

namespace
{
constexpr double kTetrahedral = 109.4712206;  // degrees
constexpr double kDegree = 3.14159265358979323846 / 180.0;

Vec3 unit(const Vec3& v)
{
    const double n = science::norm(v);
    return n > 1e-12 ? science::scale(v, 1.0 / n) : Vec3{0, 0, 0};
}

// Some unit vector perpendicular to v.
Vec3 perpendicular(const Vec3& v)
{
    const Vec3 trial = std::abs(v[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    return unit(science::cross(v, trial));
}

science::Mat3 cellMatrix(const Structure& s)
{
    science::Mat3 m{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) m[r][c] = s.cellVectors[r][c];
    return m;
}
}  // namespace

science::Mat3 rotationBetween(const Vec3& from, const Vec3& to)
{
    const Vec3 a = unit(from), b = unit(to);
    const Vec3 v = science::cross(a, b);
    const double c = science::dot(a, b);
    if (c < -1 + 1e-12) return rotationAbout(perpendicular(a), 3.14159265358979323846);
    if (science::norm(v) < 1e-14) return science::identity();
    return rotationAbout(v, std::atan2(science::norm(v), c));
}

science::Mat3 rotationAbout(const Vec3& axis, double angle)
{
    const Vec3 k = unit(axis);
    const double c = std::cos(angle), s = std::sin(angle);
    science::Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r[i][j] = (i == j ? c : 0.0) + (1 - c) * k[i] * k[j];
    r[0][1] -= s * k[2]; r[0][2] += s * k[1];
    r[1][0] += s * k[2]; r[1][2] -= s * k[0];
    r[2][0] -= s * k[1]; r[2][1] += s * k[0];
    return r;
}

Vec3 rotate(const science::Mat3& r, const Vec3& v)
{
    return {science::dot(r[0], v), science::dot(r[1], v), science::dot(r[2], v)};
}

double covalentRadius(int atomicNumber)
{
    static const std::vector<float> radii = makeLiteratureCovalentRadii();
    if (atomicNumber > 0 && static_cast<std::size_t>(atomicNumber) < radii.size() && radii[atomicNumber] > 0)
        return radii[atomicNumber];
    return 1.5;
}

std::vector<Bond> covalentBonds(const Structure& s, double scale)
{
    std::vector<Bond> bonds;
    const std::size_t n = s.atoms.size();
    if (n < 2 && !(n == 1 && s.hasUnitCell)) return bonds;
    std::vector<Vec3> positions(n);
    double largest = 0;
    for (std::size_t i = 0; i < n; ++i) {
        positions[i] = {s.atoms[i].x, s.atoms[i].y, s.atoms[i].z};
        largest = std::max(largest, covalentRadius(s.atoms[i].atomicNumber));
    }
    science::Pbc pbc{false, false, false};
    science::Mat3 cell = science::identity();
    if (s.hasUnitCell) {
        cell = cellMatrix(s);
        for (int k = 0; k < 3; ++k) pbc[static_cast<std::size_t>(k)] = s.pbc[static_cast<std::size_t>(k)];
    }
    for (const auto& link : science::neighborList(positions, cell, pbc, scale * 2.0 * largest)) {
        const auto i = static_cast<std::size_t>(link.i), j = static_cast<std::size_t>(link.j);
        if (j < i) continue;
        if (i == j) {
            // An own image: count each +-shift pair once.
            const auto& t = link.shift;
            if (t[0] < 0 || (t[0] == 0 && (t[1] < 0 || (t[1] == 0 && t[2] < 0)))) continue;
        }
        const double limit = scale * (covalentRadius(s.atoms[i].atomicNumber) + covalentRadius(s.atoms[j].atomicNumber));
        if (science::norm(link.vector) <= limit && science::norm(link.vector) > 1e-6) bonds.push_back({i, j, link.vector});
    }
    return bonds;
}

Neighbours bondNeighbours(const Structure& s, double scale)
{
    Neighbours out;
    out.atoms.assign(s.atoms.size(), {});
    out.vectors.assign(s.atoms.size(), {});
    for (const Bond& b : covalentBonds(s, scale)) {
        out.atoms[b.i].push_back(b.j);
        out.vectors[b.i].push_back(b.vector);
        out.atoms[b.j].push_back(b.i);
        out.vectors[b.j].push_back(science::scale(b.vector, -1.0));
    }
    return out;
}

void appendAtom(Structure& s, std::vector<char>& selected, const std::string& symbol, const Vec3& position, bool select)
{
    AtomSite atom;
    atom.symbol = symbol;
    atom.atomicNumber = science::atomicNumber(symbol);
    atom.x = position[0];
    atom.y = position[1];
    atom.z = position[2];
    getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
    const std::size_t n = s.atoms.size();
    if (s.grainColors.size() == n && n) s.grainColors.push_back({0.6f, 0.6f, 0.6f});
    if (s.grainRegionIds.size() == n && n) s.grainRegionIds.push_back(-1);
    if (s.atomProperty.size() == n && n) s.atomProperty.push_back(std::numeric_limits<double>::quiet_NaN());
    if (s.frozenProperty.size() == n && n) s.frozenProperty.push_back(std::numeric_limits<double>::quiet_NaN());
    s.atoms.push_back(atom);
    selected.resize(n, 0);
    selected.push_back(select ? 1 : 0);
}

std::size_t saturateWithHydrogen(Structure& s, std::vector<char>& selected, const std::set<std::string>& elements,
                                 bool forceSp3, double bondScale, bool selectNew)
{
    return saturate(s, selected, elements, forceSp3, bondScale, selectNew, SaturationOptions());
}

std::size_t saturate(Structure& s, std::vector<char>& selected, const std::set<std::string>& elements, bool forceSp3,
                     double bondScale, bool selectNew, const SaturationOptions& options)
{
    static const std::map<std::string, int> valence = {{"C", 4}, {"Si", 4}, {"Ge", 4}, {"Sn", 4}, {"N", 3}, {"P", 3},
                                                       {"B", 3}, {"Al", 3}, {"O", 2}, {"S", 2}, {"Se", 2}};
    static const std::map<std::string, double> hydrogenBond = {{"C", 1.09}, {"Si", 1.48}, {"Ge", 1.53}, {"Sn", 1.71}, {"N", 1.01},
                                                               {"P", 1.42}, {"B", 1.19}, {"Al", 1.58}, {"O", 0.97}, {"S", 1.34}, {"Se", 1.46}};
    const Neighbours nb = bondNeighbours(s, bondScale);
    std::vector<std::pair<std::size_t, Vec3>> additions;  // atom, hydrogen position
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const std::string& symbol = s.atoms[i].symbol;
        if (!elements.count(symbol)) continue;
        if (options.mask && (i >= options.mask->size() || !(*options.mask)[i])) continue;
        const auto v = valence.find(symbol);
        if (v == valence.end()) continue;
        std::vector<Vec3> u;
        for (const Vec3& d : nb.vectors[i]) u.push_back(unit(d));
        const int k = static_cast<int>(u.size());
        int missing = v->second - k;
        const bool planarCapable = !forceSp3 && (symbol == "C" || symbol == "B" || symbol == "N");
        // Planar (sp2) carbon: three bonds in a plane need nothing, two at ~120 degrees one hydrogen.
        if (planarCapable && k == 3 && science::norm(science::add(science::add(u[0], u[1]), u[2])) < 0.5) missing = 0;
        const bool sp2Pair = planarCapable && k == 2 && science::dot(u[0], u[1]) < std::cos(115.0 * kDegree);
        if (sp2Pair) missing = std::min(missing, 1);
        missing = std::min(missing, options.maxPerAtom);
        if (missing <= 0) continue;
        const auto length = hydrogenBond.find(symbol);
        const int terminatorZ = science::atomicNumber(options.terminator);
        const double bond = options.bondLength > 0 ? options.bondLength
                          : options.terminator == "H" && length != hydrogenBond.end() ? length->second
                          : covalentRadius(s.atoms[i].atomicNumber) + covalentRadius(terminatorZ);
        const Vec3 origin{s.atoms[i].x, s.atoms[i].y, s.atoms[i].z};
        std::vector<Vec3> directions;
        if (k == 0) {
            // A lone atom: tetrahedral directions, or bent water-like for divalent ones.
            const double c = 1.0 / std::sqrt(3.0);
            const std::vector<Vec3> tetra = {{c, c, c}, {c, -c, -c}, {-c, c, -c}, {-c, -c, c}};
            if (v->second == 2) {
                const double half = 0.5 * 104.5 * kDegree;
                directions = {{std::sin(half), 0, std::cos(half)}, {-std::sin(half), 0, std::cos(half)}};
            } else
                directions.assign(tetra.begin(), tetra.begin() + std::min(4, missing));
        } else if (k == 1) {
            // Around the one bond: at the tetrahedral angle (104.5 degrees for O, S).
            const double alpha = (v->second == 2 ? 104.5 : kTetrahedral) * kDegree;
            const Vec3 e1 = perpendicular(u[0]), e2 = unit(science::cross(u[0], e1));
            for (int m = 0; m < missing; ++m) {
                const double phi = 2.0 * 3.14159265358979323846 * m / 3.0;
                directions.push_back(science::add(science::scale(u[0], std::cos(alpha)),
                                                  science::scale(science::add(science::scale(e1, std::cos(phi)), science::scale(e2, std::sin(phi))),
                                                                 std::sin(alpha))));
            }
        } else if (k == 2) {
            const Vec3 b = unit(science::scale(science::add(u[0], u[1]), -1.0));
            if (sp2Pair || missing == 1) {
                if (sp2Pair) directions.push_back(b);
                else {
                    const Vec3 p = unit(science::cross(u[0], u[1]));
                    const double half = 0.5 * kTetrahedral * kDegree;
                    directions.push_back(unit(science::add(science::scale(b, std::cos(half)), science::scale(p, std::sin(half)))));
                }
            } else {
                const Vec3 p = unit(science::cross(u[0], u[1]));
                const double half = 0.5 * kTetrahedral * kDegree;
                directions.push_back(unit(science::add(science::scale(b, std::cos(half)), science::scale(p, std::sin(half)))));
                directions.push_back(unit(science::sub(science::scale(b, std::cos(half)), science::scale(p, std::sin(half)))));
            }
        } else {
            Vec3 d = science::scale(science::add(science::add(u[0], u[1]), u[2]), -1.0);
            if (science::norm(d) < 1e-3) d = science::cross(science::sub(u[1], u[0]), science::sub(u[2], u[0]));
            directions.push_back(unit(d));
        }
        for (const Vec3& d : directions) additions.push_back({i, science::add(origin, science::scale(unit(d), bond))});
    }
    for (const auto& [atom, position] : additions) appendAtom(s, selected, options.terminator, position, selectNew);
    return additions.size();
}
}  // namespace atomforge::pipeline
