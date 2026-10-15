#include "science/TemplateMatching.h"
#include "util/TaskControl.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using Quat = Eigen::Quaterniond;

const double kPi = std::acos(-1.0);
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kIdealCa = std::sqrt(8.0 / 3.0);

enum StructureKind { Other = 0, Fcc, Hcp, Bcc, Ico, Sc, CubicDiamond, HexDiamond, KindCount };
const char* const kNames[KindCount] = {"other", "fcc", "hcp", "bcc", "icosahedral", "simple cubic", "cubic diamond", "hexagonal diamond"};
const char* const kKeys[KindCount] = {"other", "fcc", "hcp", "bcc", "ico", "sc", "cubic-diamond", "hex-diamond"};

// ---------------------------------------------------------------- templates
using Cell = std::array<V3, 3>;

V3 cartesian(const Cell& cell, const V3& f) { return f[0] * cell[0] + f[1] * cell[1] + f[2] * cell[2]; }

// The `count` vectors from basis atom 0 to its nearest lattice neighbours,
// ordered by length in the `ideal` lattice and evaluated in `actual` (the same
// lattice with another axial ratio), so the order does not depend on c/a.
std::vector<V3> shell(const Cell& ideal, const Cell& actual, const std::vector<V3>& basis, std::size_t count)
{
    struct Candidate { V3 ideal, actual; long long key; };
    std::vector<Candidate> candidates;
    for (int i = -3; i <= 3; ++i)
        for (int j = -3; j <= 3; ++j)
            for (int k = -3; k <= 3; ++k)
                for (const V3& b : basis) {
                    const V3 f = b + V3(i, j, k) - basis[0];
                    const V3 v = cartesian(ideal, f);
                    if (v.norm() < 1e-9) continue;
                    candidates.push_back({v, cartesian(actual, f), std::llround(v.norm() * 1e6)});
                }
    std::stable_sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) { return a.key < b.key; });
    if (candidates.size() <= count || candidates[count].key == candidates[count - 1].key)
        throw std::logic_error("template shell is not closed");
    std::vector<V3> out;
    for (std::size_t n = 0; n < count; ++n) out.push_back(candidates[n].actual);
    return out;
}

Cell hexagonalCell(double ca) { return {V3(1, 0, 0), V3(-0.5, std::sqrt(3.0) / 2, 0), V3(0, 0, ca)}; }
const Cell kCubic = {V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)};

// Neighbour vectors of each structure in units of its lattice constant
// (cubic a; hexagonal a; the icosahedral and simple-cubic bond length).
std::vector<V3> templateVectors(int kind, double ca = kIdealCa)
{
    const std::vector<V3> fccBasis = {V3(0, 0, 0), V3(0, 0.5, 0.5), V3(0.5, 0, 0.5), V3(0.5, 0.5, 0)};
    switch (kind) {
    case Fcc: return shell(kCubic, kCubic, fccBasis, 12);
    case Hcp: return shell(hexagonalCell(kIdealCa), hexagonalCell(ca), {V3(0, 0, 0), V3(1.0 / 3, 2.0 / 3, 0.5)}, 12);
    case Bcc: return shell(kCubic, kCubic, {V3(0, 0, 0), V3(0.5, 0.5, 0.5)}, 14);
    case Sc: return shell(kCubic, kCubic, {V3(0, 0, 0)}, 6);
    case CubicDiamond: {
        std::vector<V3> basis = fccBasis;
        for (const V3& b : fccBasis) basis.push_back(b + V3(0.25, 0.25, 0.25));
        return shell(kCubic, kCubic, basis, 16);
    }
    case HexDiamond: {
        const double z = 1.0 / 16;
        const std::vector<V3> basis = {V3(1.0 / 3, 2.0 / 3, z), V3(2.0 / 3, 1.0 / 3, 0.5 + z), V3(2.0 / 3, 1.0 / 3, 1 - z), V3(1.0 / 3, 2.0 / 3, 0.5 - z)};
        return shell(hexagonalCell(kIdealCa), hexagonalCell(ca), basis, 16);
    }
    case Ico: {
        const double phi = (1 + std::sqrt(5.0)) / 2;
        std::vector<V3> out;
        for (double s : {-1.0, 1.0})
            for (double t : {-phi, phi}) {
                out.push_back(V3(0, s, t));
                out.push_back(V3(s, t, 0));
                out.push_back(V3(t, 0, s));
            }
        for (V3& v : out) v /= std::sqrt(1 + phi * phi);
        return out;
    }
    default: throw std::logic_error("unknown template");
    }
}

// Orthonormal frame with e1 along u and e2 in the (u, w) plane.
M3 frame(const V3& u, const V3& w)
{
    const V3 e1 = u.normalized();
    const V3 e2 = (w - w.dot(e1) * e1).normalized();
    M3 m;
    m.col(0) = e1;
    m.col(1) = e2;
    m.col(2) = e1.cross(e2);
    return m;
}

bool parallel(const V3& a, const V3& b) { return a.cross(b).norm() < 0.05 * a.norm() * b.norm(); }

struct Template
{
    std::vector<V3> unit;                         // scaled to mean length 1
    double meanLength = 0;                        // mean length in lattice-constant units
    std::vector<M3> symmetry;                     // proper rotations permuting the vectors
    std::vector<Quat> symmetryQuats;
    std::vector<std::pair<int, int>> anchors;     // ordered pairs, one per symmetry orbit
    std::vector<double> anchorRatio, anchorCos;   // |t2|/|t1| and cos(t1, t2)
};

Template makeTemplate(int kind)
{
    Template t;
    const auto vectors = templateVectors(kind);
    for (const V3& v : vectors) t.meanLength += v.norm();
    t.meanLength /= static_cast<double>(vectors.size());
    for (const V3& v : vectors) t.unit.push_back(v / t.meanLength);
    const auto& u = t.unit;
    const int n = static_cast<int>(u.size());
    // Symmetry: rotations taking (u0, ub) to an equivalent pair that map the set onto itself.
    int b = 1;
    while (parallel(u[0], u[static_cast<std::size_t>(b)])) ++b;
    const M3 reference = frame(u[0], u[static_cast<std::size_t>(b)]);
    std::vector<std::vector<int>> permutations;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const V3 &a = u[static_cast<std::size_t>(i)], &c = u[static_cast<std::size_t>(j)];
            if (i == j || std::abs(a.norm() - u[0].norm()) > 1e-6 || std::abs(c.norm() - u[static_cast<std::size_t>(b)].norm()) > 1e-6 ||
                std::abs(a.dot(c) - u[0].dot(u[static_cast<std::size_t>(b)])) > 1e-6)
                continue;
            const M3 rotation = frame(a, c) * reference.transpose();
            std::vector<int> permutation(static_cast<std::size_t>(n), -1);
            bool ok = true;
            for (int k = 0; k < n && ok; ++k) {
                const V3 image = rotation * u[static_cast<std::size_t>(k)];
                for (int m = 0; m < n; ++m)
                    if ((image - u[static_cast<std::size_t>(m)]).norm() < 1e-6) permutation[static_cast<std::size_t>(k)] = m;
                ok = permutation[static_cast<std::size_t>(k)] >= 0;
            }
            if (!ok) continue;
            t.symmetry.push_back(rotation);
            t.symmetryQuats.push_back(Quat(rotation).normalized());
            permutations.push_back(permutation);
        }
    // Anchor pairs: the lexicographically smallest member of each orbit.
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (i == j || parallel(u[static_cast<std::size_t>(i)], u[static_cast<std::size_t>(j)])) continue;
            bool smallest = true;
            for (const auto& perm : permutations) {
                const std::pair<int, int> image(perm[static_cast<std::size_t>(i)], perm[static_cast<std::size_t>(j)]);
                if (image < std::make_pair(i, j)) { smallest = false; break; }
            }
            if (!smallest) continue;
            const V3 &a = u[static_cast<std::size_t>(i)], &c = u[static_cast<std::size_t>(j)];
            t.anchors.push_back({i, j});
            t.anchorRatio.push_back(c.norm() / a.norm());
            t.anchorCos.push_back(a.dot(c) / (a.norm() * c.norm()));
        }
    return t;
}

const Template& templateOf(int kind)
{
    static const std::vector<Template> templates = [] {
        std::vector<Template> all(KindCount);
        for (int k = 1; k < KindCount; ++k) all[static_cast<std::size_t>(k)] = makeTemplate(k);
        return all;
    }();
    return templates[static_cast<std::size_t>(kind)];
}

std::size_t shellSize(int kind) { return templateOf(kind).unit.size(); }

// ---------------------------------------------------------------- fitting
// Minimum-cost perfect matching (Hungarian algorithm) of rows to columns.
std::vector<int> hungarian(const std::vector<std::vector<double>>& cost)
{
    const std::size_t n = cost.size();
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> u(n + 1), v(n + 1), minv(n + 1);
    std::vector<std::size_t> p(n + 1), way(n + 1);
    std::vector<char> used(n + 1);
    for (std::size_t i = 1; i <= n; ++i) {
        p[0] = i;
        std::size_t j0 = 0;
        std::fill(minv.begin(), minv.end(), inf);
        std::fill(used.begin(), used.end(), 0);
        do {
            used[j0] = 1;
            const std::size_t i0 = p[j0];
            double delta = inf;
            std::size_t j1 = 0;
            for (std::size_t j = 1; j <= n; ++j) {
                if (used[j]) continue;
                const double current = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if (current < minv[j]) { minv[j] = current; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (std::size_t j = 0; j <= n; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const std::size_t j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    std::vector<int> rowToColumn(n);
    for (std::size_t j = 1; j <= n; ++j) rowToColumn[p[j] - 1] = static_cast<int>(j - 1);
    return rowToColumn;
}

// Template index for each data vector after rotating the template by R.
std::vector<int> assign(const std::vector<V3>& t, const std::vector<V3>& b, const M3& rotation)
{
    const std::size_t n = b.size();
    std::vector<V3> rotated(n);
    for (std::size_t j = 0; j < n; ++j) rotated[j] = rotation * t[j];
    std::vector<int> map(n);
    std::vector<char> used(n, 0);
    bool bijective = true;
    for (std::size_t k = 0; k < n; ++k) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < n; ++j) {
            const double d = (rotated[j] - b[k]).squaredNorm();
            if (d < best) { best = d; map[k] = static_cast<int>(j); }
        }
        if (used[static_cast<std::size_t>(map[k])]++) bijective = false;
    }
    if (bijective) return map;
    std::vector<std::vector<double>> cost(n, std::vector<double>(n));
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j) cost[k][j] = (rotated[j] - b[k]).squaredNorm();
    return hungarian(cost);
}

// Rotation minimising sum |R t_map(k) - b_k|^2 (Kabsch).
M3 kabsch(const std::vector<V3>& t, const std::vector<V3>& b, const std::vector<int>& map)
{
    M3 a = M3::Zero();
    for (std::size_t k = 0; k < b.size(); ++k) a += b[k] * t[static_cast<std::size_t>(map[k])].transpose();
    const Eigen::JacobiSVD<M3> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
    M3 d = M3::Identity();
    if ((svd.matrixU() * svd.matrixV().transpose()).determinant() < 0) d(2, 2) = -1;
    return svd.matrixU() * d * svd.matrixV().transpose();
}

struct Fit
{
    int kind = Other;
    double rmsd = std::numeric_limits<double>::infinity();
    double scale = 0;          // data mean length / template mean length (Angstrom per lattice constant)
    M3 rotation = M3::Identity();
    std::vector<int> map;      // data neighbour k -> template vector
};

// Best rotation of template `kind` onto the data vectors (Angstrom).
Fit fitTemplate(int kind, const std::vector<V3>& data)
{
    const Template& tpl = templateOf(kind);
    const std::size_t n = data.size();
    double mean = 0;
    for (const V3& v : data) mean += v.norm();
    mean /= static_cast<double>(n);
    std::vector<V3> b(n);
    for (std::size_t k = 0; k < n; ++k) b[k] = data[k] / mean;
    // Anchor: the nearest neighbour and the neighbours at the smallest angles to it.
    std::size_t first = 0;
    for (std::size_t k = 1; k < n; ++k)
        if (b[k].norm() < b[first].norm()) first = k;
    std::vector<std::pair<double, std::size_t>> byAngle;
    for (std::size_t k = 0; k < n; ++k)
        if (k != first && !parallel(b[first], b[k])) byAngle.push_back({-b[first].dot(b[k]) / (b[first].norm() * b[k].norm()), k});
    std::sort(byAngle.begin(), byAngle.end());
    Fit best;
    best.kind = kind;
    for (std::size_t c = 0; c < std::min<std::size_t>(2, byAngle.size()); ++c) {
        const std::size_t second = byAngle[c].second;
        const double ratio = b[second].norm() / b[first].norm();
        const double cosine = -byAngle[c].first;
        const M3 dataFrame = frame(b[first], b[second]);
        for (std::size_t a = 0; a < tpl.anchors.size(); ++a) {
            if (std::abs(tpl.anchorRatio[a] - ratio) > 0.12 * tpl.anchorRatio[a]) continue;
            if (std::abs(std::acos(std::clamp(tpl.anchorCos[a], -1.0, 1.0)) - std::acos(std::clamp(cosine, -1.0, 1.0))) > 15.0 * kPi / 180) continue;
            const auto [i, j] = tpl.anchors[a];
            M3 rotation = dataFrame * frame(tpl.unit[static_cast<std::size_t>(i)], tpl.unit[static_cast<std::size_t>(j)]).transpose();
            std::vector<int> map = assign(tpl.unit, b, rotation);
            rotation = kabsch(tpl.unit, b, map);
            map = assign(tpl.unit, b, rotation);
            rotation = kabsch(tpl.unit, b, map);
            double cross = 0, templateNorm = 0;
            for (std::size_t k = 0; k < n; ++k) {
                const V3 rt = rotation * tpl.unit[static_cast<std::size_t>(map[k])];
                cross += b[k].dot(rt);
                templateNorm += rt.squaredNorm();
            }
            const double s = cross / templateNorm;
            double residual = 0;
            for (std::size_t k = 0; k < n; ++k)
                residual += (b[k] - s * (rotation * tpl.unit[static_cast<std::size_t>(map[k])])).squaredNorm();
            const double rmsd = std::sqrt(residual / static_cast<double>(n));
            if (rmsd < best.rmsd) {
                best.rmsd = rmsd;
                best.rotation = rotation;
                best.map = std::move(map);
                best.scale = mean * s / tpl.meanLength;
            }
        }
        if (best.rmsd < 0.02) break;
    }
    return best;
}

// Proper rotations of the lattice (not of one site): the two sublattices of
// hcp and diamond have shells rotated by a lattice symmetry that is no site
// symmetry, so orientations are reduced with the lattice point group (cubic
// 432, hexagonal 622 with c along z, icosahedral 532).
const std::vector<Quat>& latticeSymmetry(int kind)
{
    static const std::vector<Quat> hexagonal = [] {
        std::vector<Quat> out;
        for (int k = 0; k < 6; ++k) {
            out.push_back(Quat(Eigen::AngleAxisd(k * kPi / 3, V3::UnitZ())));
            out.push_back(Quat(Eigen::AngleAxisd(kPi, V3(std::cos(k * kPi / 6), std::sin(k * kPi / 6), 0))));
        }
        return out;
    }();
    if (kind == Hcp || kind == HexDiamond) return hexagonal;
    if (kind == Ico) return templateOf(Ico).symmetryQuats;
    return templateOf(Fcc).symmetryQuats;
}

// The symmetry-equivalent orientation closest to the identity (w >= 0).
Quat reducedOrientation(int kind, const M3& rotation)
{
    const Quat q(rotation);
    Quat best = q;
    double bestW = -1;
    for (const Quat& s : latticeSymmetry(kind)) {
        const Quat candidate = (q * s).normalized();
        if (std::abs(candidate.w()) > bestW) { bestW = std::abs(candidate.w()); best = candidate; }
    }
    if (best.w() < 0) best.coeffs() *= -1;
    return best;
}

// Smallest rotation angle (radians) between two orientations of the same structure.
double misorientation(int kind, const Quat& a, const Quat& b)
{
    const Quat delta = a.conjugate() * b;
    double bestW = 0;
    for (const Quat& s : latticeSymmetry(kind)) bestW = std::max(bestW, std::abs((delta * s).w()));
    return 2 * std::acos(std::min(1.0, bestW));
}

// ---------------------------------------------------------------- neighbours
struct Neighbors
{
    std::vector<std::vector<int>> index;   // nearest neighbours of each atom, by distance
    std::vector<std::vector<V3>> vectors;
};

int kindFromKey(const std::string& key)
{
    for (int k = 1; k < KindCount; ++k)
        if (key == kKeys[k]) return k;
    std::string known;
    for (int k = 1; k < KindCount; ++k) known += std::string(k > 1 ? ", " : "") + kKeys[k];
    throw std::runtime_error("Unknown structure '" + key + "'; choose from " + known);
}

std::vector<int> enabledKinds(const Parameters& p)
{
    if (!p.has("structures")) return {Fcc, Hcp, Bcc, Ico};
    std::vector<int> kinds;
    const Json& list = p.json("structures");
    if (list.isString()) kinds.push_back(kindFromKey(list.string()));
    else
        for (const Json& item : list.items()) kinds.push_back(kindFromKey(item.string()));
    if (kinds.empty()) throw std::runtime_error("structures must name at least one structure");
    return kinds;
}

struct Geometry
{
    std::vector<Vec3> positions;
    Mat3 cell{};
    Pbc pbc{};
    double cutoff = 0;
    bool automatic = false;  // cutoff chosen here, so it may be widened
};

Geometry geometry(const Parameters& p)
{
    Geometry g;
    g.positions = points(p.array("positions"), "positions");
    g.pbc = p.pbc("pbc", {false, false, false});
    if ((g.pbc[0] || g.pbc[1] || g.pbc[2]) && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    g.cell = p.has("cell") ? cellMatrix(p.array("cell")) : identity();
    g.cutoff = p.number("cutoff_A", 0.0);
    if (g.cutoff < 0) throw std::runtime_error("cutoff_A must be nonnegative");
    if (g.cutoff == 0 && !g.positions.empty()) {
        // Automatic: 1.8 times the mean atomic spacing reaches the second bcc and diamond shells.
        double volume = 0;
        if (p.has("cell")) volume = cellVolume(g.cell);
        else {
            Vec3 low = g.positions[0], high = g.positions[0];
            for (const auto& x : g.positions)
                for (int d = 0; d < 3; ++d) { low[static_cast<std::size_t>(d)] = std::min(low[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]); high[static_cast<std::size_t>(d)] = std::max(high[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]); }
            volume = 1;
            for (int d = 0; d < 3; ++d) volume *= high[static_cast<std::size_t>(d)] - low[static_cast<std::size_t>(d)] + 2.0;
        }
        g.cutoff = 1.8 * std::cbrt(volume / static_cast<double>(g.positions.size()));
        g.automatic = true;
    }
    return g;
}

Neighbors nearestNeighbors(Geometry& g, std::size_t needed)
{
    const std::size_t n = g.positions.size();
    for (int attempt = 0;; ++attempt) {
        Neighbors result;
        result.index.assign(n, {});
        result.vectors.assign(n, {});
        std::vector<std::vector<std::pair<double, std::pair<int, V3>>>> lists(n);
        for (const auto& link : neighborList(g.positions, g.cell, g.pbc, g.cutoff)) {
            const V3 v(link.vector[0], link.vector[1], link.vector[2]);
            if (v.norm() < 1e-10) throw std::runtime_error("Coincident atoms make local environment analysis undefined");
            lists[static_cast<std::size_t>(link.i)].push_back({v.norm(), {link.j, v}});
        }
        std::size_t enough = 0;
        for (std::size_t i = 0; i < n; ++i) {
            auto& list = lists[i];
            const std::size_t keep = std::min<std::size_t>(list.size(), 16);
            std::partial_sort(list.begin(), list.begin() + static_cast<std::ptrdiff_t>(keep), list.end(),
                              [](const auto& a, const auto& b) { return a.first < b.first; });
            for (std::size_t k = 0; k < keep; ++k) {
                result.index[i].push_back(list[k].second.first);
                result.vectors[i].push_back(list[k].second.second);
            }
            enough += list.size() >= needed;
        }
        // Too short a cutoff for most atoms: widen it (a few times at most).
        if (!g.automatic || attempt >= 3 || n == 0 || 2 * enough >= n) return result;
        g.cutoff *= 1.3;
    }
}

// Neighbour shell of atom i for a structure: the nearest atoms, or for diamond
// the four nearest and the three further neighbours of each of them.
bool environment(int kind, std::size_t i, const Neighbors& nb, std::vector<V3>& out, std::vector<int>& atoms)
{
    out.clear();
    atoms.clear();
    const auto& vectors = nb.vectors[i];
    if (kind != CubicDiamond && kind != HexDiamond) {
        const std::size_t count = shellSize(kind);
        if (vectors.size() < count) return false;
        out.assign(vectors.begin(), vectors.begin() + static_cast<std::ptrdiff_t>(count));
        atoms.assign(nb.index[i].begin(), nb.index[i].begin() + static_cast<std::ptrdiff_t>(count));
        return true;
    }
    if (vectors.size() < 4) return false;
    for (std::size_t k = 0; k < 4; ++k) { out.push_back(vectors[k]); atoms.push_back(nb.index[i][k]); }
    for (std::size_t k = 0; k < 4; ++k) {
        const std::size_t j = static_cast<std::size_t>(nb.index[i][k]);
        const auto& second = nb.vectors[j];
        if (second.size() < 4) return false;
        // The back bond is the one closest to -v; keep the other three.
        std::size_t back = 0;
        for (std::size_t m = 1; m < 4; ++m)
            if ((second[m] + vectors[k]).norm() < (second[back] + vectors[k]).norm()) back = m;
        for (std::size_t m = 0; m < 4; ++m)
            if (m != back) { out.push_back(vectors[k] + second[m]); atoms.push_back(nb.index[j][m]); }
    }
    return true;
}

struct Matching
{
    Geometry geometry;
    Neighbors neighbors;
    std::vector<Fit> fits;   // per atom; kind Other when unidentified
    std::vector<int> kinds;  // enabled
};

Matching matchAll(const Parameters& p, const std::vector<int>& kinds, double rmsdCutoff)
{
    std::size_t needed = 4;
    for (int k : kinds)
        if (k != CubicDiamond && k != HexDiamond) needed = std::max(needed, shellSize(k));
    Matching m;
    m.kinds = kinds;
    m.geometry = geometry(p);
    m.neighbors = nearestNeighbors(m.geometry, needed);
    const std::size_t n = m.geometry.positions.size();
    m.fits.assign(n, Fit());
    std::vector<V3> env;
    std::vector<int> atoms;
    for (std::size_t i = 0; i < n; ++i) {
        if (i % 256 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
        Fit best;
        for (int kind : kinds) {
            if (!environment(kind, i, m.neighbors, env, atoms)) continue;
            Fit fit = fitTemplate(kind, env);
            if (fit.rmsd < best.rmsd) best = std::move(fit);
        }
        if (best.kind == Other || !(best.rmsd <= rmsdCutoff || rmsdCutoff <= 0)) {
            const double rmsd = best.rmsd;
            best = Fit();
            best.rmsd = rmsd;
        }
        m.fits[i] = std::move(best);
    }
    return m;
}

double rmsdCutoffOf(const Parameters& p)
{
    const double cutoff = p.number("rmsd_cutoff", 0.1);
    if (cutoff < 0) throw std::runtime_error("rmsd_cutoff must be nonnegative (0 = no cutoff)");
    return cutoff;
}

Json quaternionJson(const Quat& q) { return Json::array({q.w(), q.x(), q.y(), q.z()}); }

Json countsJson(const std::vector<int>& types, const std::vector<std::string>& names, Json& fractions)
{
    std::vector<int> counts(names.size(), 0);
    for (int t : types) ++counts[static_cast<std::size_t>(t)];
    Json result = Json::object();
    fractions = Json::object();
    for (std::size_t k = 0; k < names.size(); ++k) {
        result[names[k]] = counts[k];
        fractions[names[k]] = types.empty() ? 0.0 : static_cast<double>(counts[k]) / static_cast<double>(types.size());
    }
    return result;
}

std::vector<std::string> kindNames()
{
    return std::vector<std::string>(kNames, kNames + KindCount);
}

// Union-find with path halving.
std::size_t root(std::vector<std::size_t>& parent, std::size_t i)
{
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}
}

Json polyhedralTemplateMatching(const Parameters& p)
{
    const Matching m = matchAll(p, enabledKinds(p), rmsdCutoffOf(p));
    const std::size_t n = m.fits.size();
    std::vector<int> types(n);
    std::vector<double> rmsd(n, kNaN), lattice(n, kNaN), angle(n, kNaN);
    Json orientations = Json::array();
    for (std::size_t i = 0; i < n; ++i) {
        const Fit& fit = m.fits[i];
        types[i] = fit.kind;
        if (std::isfinite(fit.rmsd)) rmsd[i] = fit.rmsd;
        if (fit.kind == Other) {
            orientations.push(Json::array({kNaN, kNaN, kNaN, kNaN}));
            continue;
        }
        lattice[i] = fit.scale;
        const Quat q = reducedOrientation(fit.kind, fit.rotation);
        angle[i] = 2 * std::acos(std::min(1.0, q.w())) * 180 / kPi;
        orientations.push(quaternionJson(q));
    }
    Json result = Json::object();
    result["structure_type"] = toJson(types);
    Json names = Json::array();
    for (const auto& name : kindNames()) names.push(name);
    result["type_names"] = names;
    Json fractions;
    result["counts"] = countsJson(types, kindNames(), fractions);
    result["fractions"] = fractions;
    result["rmsd"] = toJson(rmsd);
    result["orientation"] = orientations;
    result["orientation_angle_deg"] = toJson(angle);
    result["lattice_constant_A"] = toJson(lattice);
    result["cutoff_A"] = m.geometry.cutoff;
    return result;
}

Json grainSegmentation(const Parameters& p)
{
    const double threshold = positive(p.number("merge_threshold_deg", 4.0), "merge_threshold_deg") * kPi / 180;
    const long long minimum = integer(p.number("min_grain_size", 100), "min_grain_size");
    const bool adopt = p.boolean("adopt_orphans", true);
    const Matching m = matchAll(p, enabledKinds(p), rmsdCutoffOf(p));
    const std::size_t n = m.fits.size();
    std::vector<Quat> orientation(n, Quat::Identity());
    for (std::size_t i = 0; i < n; ++i)
        if (m.fits[i].kind != Other) orientation[i] = reducedOrientation(m.fits[i].kind, m.fits[i].rotation);
    // Bonds between crystalline neighbours of the same structure, by misorientation.
    struct Edge { double angle; std::size_t a, b; };
    std::vector<Edge> edges;
    std::vector<V3> env;
    std::vector<int> atoms;
    for (std::size_t i = 0; i < n; ++i) {
        const int kind = m.fits[i].kind;
        if (kind == Other || !environment(kind, i, m.neighbors, env, atoms)) continue;
        const std::size_t first = (kind == CubicDiamond || kind == HexDiamond) ? 4 : atoms.size();
        for (std::size_t k = 0; k < first; ++k) {
            const std::size_t j = static_cast<std::size_t>(atoms[k]);
            if (j <= i || m.fits[j].kind != kind) continue;
            const double angle = misorientation(kind, orientation[i], orientation[j]);
            if (angle < threshold) edges.push_back({angle, i, j});
        }
    }
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) { return a.angle < b.angle; });
    // Merge clusters while their mean orientations agree (keeps slowly
    // rotating chains of atoms from bridging distinct grains).
    std::vector<std::size_t> parent(n), size(n, 1);
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<Eigen::Vector4d> sum(n);
    for (std::size_t i = 0; i < n; ++i) sum[i] = orientation[i].coeffs();
    const auto meanOf = [&](std::size_t r) { Quat q; q.coeffs() = sum[r].normalized(); return q; };
    for (std::size_t e = 0; e < edges.size(); ++e) {
        if (e % 4096 == 0) taskCheckpoint();
        std::size_t a = root(parent, edges[e].a), b = root(parent, edges[e].b);
        if (a == b) continue;
        const int kind = m.fits[edges[e].a].kind;
        const Quat qa = meanOf(a), qb = meanOf(b);
        if (misorientation(kind, qa, qb) >= threshold) continue;
        if (size[a] < size[b]) std::swap(a, b);
        // Bring b's orientation sum into a's symmetry branch before adding it.
        const Quat delta = qa.conjugate() * meanOf(b);
        Quat bestS = Quat::Identity();
        double bestW = -1;
        for (const Quat& s : latticeSymmetry(kind)) {
            const double w = std::abs((delta * s).w());
            if (w > bestW) { bestW = w; bestS = s; }
        }
        Quat aligned;
        aligned.coeffs() = sum[b];
        aligned = aligned * bestS;
        Eigen::Vector4d add = aligned.coeffs();
        if (add.dot(sum[a]) < 0) add = -add;
        sum[a] += add;
        size[a] += size[b];
        parent[b] = a;
    }
    std::map<std::size_t, std::vector<std::size_t>> clusters;
    for (std::size_t i = 0; i < n; ++i)
        if (m.fits[i].kind != Other) clusters[root(parent, i)].push_back(i);
    std::vector<std::pair<std::size_t, std::vector<std::size_t>>> grains;
    for (auto& [r, members] : clusters)
        if (static_cast<long long>(members.size()) >= minimum) grains.push_back({r, std::move(members)});
    std::sort(grains.begin(), grains.end(), [](const auto& a, const auto& b) { return a.second.size() > b.second.size(); });
    std::vector<int> grain(n, 0);
    for (std::size_t g = 0; g < grains.size(); ++g)
        for (std::size_t i : grains[g].second) grain[i] = static_cast<int>(g + 1);
    // Orphans (non-crystalline atoms and small clusters) join the grain of
    // their nearest assigned neighbour, breadth first.
    std::size_t adopted = 0;
    if (adopt && !grains.empty())
        for (;;) {
            std::vector<std::pair<std::size_t, int>> claims;
            for (std::size_t i = 0; i < n; ++i) {
                if (grain[i]) continue;
                for (int j : m.neighbors.index[i])
                    if (grain[static_cast<std::size_t>(j)]) { claims.push_back({i, grain[static_cast<std::size_t>(j)]}); break; }
            }
            if (claims.empty()) break;
            for (const auto& [i, g] : claims) grain[i] = g;
            adopted += claims.size();
        }
    // Grain table and boundaries between grains.
    std::vector<long long> atomsPerGrain(grains.size() + 1, 0);
    for (int g : grain) ++atomsPerGrain[static_cast<std::size_t>(g)];
    Json table = Json::array();
    for (std::size_t g = 0; g < grains.size(); ++g) {
        const std::size_t r = grains[g].first;
        const int kind = m.fits[grains[g].second.front()].kind;
        Json item = Json::object();
        item["id"] = static_cast<double>(g + 1);
        item["atoms"] = static_cast<double>(atomsPerGrain[g + 1]);
        item["crystalline_atoms"] = static_cast<double>(grains[g].second.size());
        item["structure"] = kNames[kind];
        item["orientation"] = quaternionJson(reducedOrientation(kind, meanOf(r).toRotationMatrix()));
        table.push(item);
    }
    std::map<std::pair<int, int>, long long> contacts;
    for (std::size_t i = 0; i < n; ++i)
        for (int j : m.neighbors.index[i])
            if (grain[i] && grain[static_cast<std::size_t>(j)] && grain[i] < grain[static_cast<std::size_t>(j)])
                ++contacts[{grain[i], grain[static_cast<std::size_t>(j)]}];
    Json boundaries = Json::array();
    for (const auto& [pair, count] : contacts) {
        const auto& ga = grains[static_cast<std::size_t>(pair.first - 1)];
        const auto& gb = grains[static_cast<std::size_t>(pair.second - 1)];
        const int ka = m.fits[ga.second.front()].kind, kb = m.fits[gb.second.front()].kind;
        Json item = Json::object();
        item["grains"] = Json::array({static_cast<double>(pair.first), static_cast<double>(pair.second)});
        item["neighbor_pairs"] = static_cast<double>(count);
        item["misorientation_deg"] = ka == kb ? misorientation(ka, meanOf(ga.first), meanOf(gb.first)) * 180 / kPi : kNaN;
        boundaries.push(item);
    }
    std::vector<int> types(n);
    for (std::size_t i = 0; i < n; ++i) types[i] = m.fits[i].kind;
    Json result = Json::object();
    result["grain_id"] = toJson(grain);
    result["structure_type"] = toJson(types);
    result["grain_count"] = static_cast<double>(grains.size());
    result["grains"] = table;
    result["boundaries"] = boundaries;
    result["adopted_atoms"] = static_cast<double>(adopted);
    result["unassigned_atoms"] = static_cast<double>(atomsPerGrain[0]);
    result["cutoff_A"] = m.geometry.cutoff;
    return result;
}

Json elasticStrain(const Parameters& p)
{
    const double lattice = positive(p.number("lattice_constant_A"), "lattice_constant_A");
    const double ca = positive(p.number("c_over_a", kIdealCa), "c_over_a");
    const int kind = kindFromKey(p.has("lattice") ? p.json("lattice").string() : std::string("fcc"));
    if (kind == Ico) throw std::runtime_error("Elastic strain needs a crystal lattice, not icosahedral order");
    const Matching m = matchAll(p, {kind}, rmsdCutoffOf(p));
    const auto reference = templateVectors(kind, ca);
    const std::size_t n = m.fits.size();
    std::vector<double> volumetric(n, kNaN), shear(n, kNaN);
    Json strain = Json::array(), gradient = Json::array();
    std::vector<V3> env;
    std::vector<int> atoms;
    std::size_t identified = 0;
    double volumetricSum = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Fit& fit = m.fits[i];
        if (fit.kind != kind || !environment(kind, i, m.neighbors, env, atoms)) {
            strain.push(Json::array({kNaN, kNaN, kNaN, kNaN, kNaN, kNaN}));
            gradient.push(Json::array({kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN}));
            continue;
        }
        M3 bt = M3::Zero(), tt = M3::Zero();
        for (std::size_t k = 0; k < env.size(); ++k) {
            const V3 t = lattice * reference[static_cast<std::size_t>(fit.map[k])];
            bt += env[k] * t.transpose();
            tt += t * t.transpose();
        }
        const M3 f = bt * tt.inverse();
        // Green-Lagrange strain pushed to the spatial frame: (F F^T - I) / 2.
        const M3 e = 0.5 * (f * f.transpose() - M3::Identity());
        strain.push(Json::array({e(0, 0), e(1, 1), e(2, 2), e(1, 2), e(0, 2), e(0, 1)}));
        // F relative to the crystal's own rotation: R^T F, the lattice-frame stretch.
        const M3 local = fit.rotation.transpose() * f;
        gradient.push(Json::array({local(0, 0), local(0, 1), local(0, 2), local(1, 0), local(1, 1), local(1, 2), local(2, 0), local(2, 1), local(2, 2)}));
        volumetric[i] = e.trace() / 3;
        shear[i] = std::sqrt(e(0, 1) * e(0, 1) + e(1, 2) * e(1, 2) + e(0, 2) * e(0, 2) +
                             ((e(0, 0) - e(1, 1)) * (e(0, 0) - e(1, 1)) + (e(1, 1) - e(2, 2)) * (e(1, 1) - e(2, 2)) +
                              (e(0, 0) - e(2, 2)) * (e(0, 0) - e(2, 2))) / 6);
        ++identified;
        volumetricSum += volumetric[i];
    }
    Json result = Json::object();
    result["strain_voigt"] = strain;
    result["strain_components"] = Json::array({"xx", "yy", "zz", "yz", "xz", "xy"});
    result["volumetric_strain"] = toJson(volumetric);
    result["shear_strain"] = toJson(shear);
    result["deformation_gradient_lattice"] = gradient;
    result["structure"] = kNames[kind];
    result["identified_atoms"] = static_cast<double>(identified);
    result["mean_volumetric_strain"] = identified ? volumetricSum / static_cast<double>(identified) : kNaN;
    result["cutoff_A"] = m.geometry.cutoff;
    return result;
}

Json fccPlanarFaults(const Parameters& p)
{
    const Matching m = matchAll(p, {Fcc, Hcp}, rmsdCutoffOf(p));
    const std::size_t n = m.fits.size();
    static const std::vector<std::string> names = {"non-hcp", "intrinsic stacking fault", "coherent twin boundary", "multilayer stacking fault", "other hcp"};
    std::vector<int> fault(n, 0);
    std::vector<V3> env;
    std::vector<int> atoms;
    double spacing = 0;
    std::size_t spacingCount = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (m.fits[i].kind == Fcc || m.fits[i].kind == Hcp) {
            // Nearest-neighbour distance: a/sqrt 2 in fcc, a in hcp.
            spacing += m.fits[i].kind == Fcc ? m.fits[i].scale / std::sqrt(2.0) : m.fits[i].scale;
            ++spacingCount;
        }
        if (m.fits[i].kind != Hcp || !environment(Hcp, i, m.neighbors, env, atoms)) continue;
        const V3 normal = m.fits[i].rotation * V3(0, 0, 1);
        int up = 0, down = 0;
        for (std::size_t k = 0; k < env.size(); ++k) {
            const double c = env[k].normalized().dot(normal);
            if (m.fits[static_cast<std::size_t>(atoms[k])].kind != Hcp) continue;
            if (c > 0.5) ++up;
            else if (c < -0.5) ++down;
        }
        if (up == 0 && down == 0) fault[i] = 2;
        else if ((up == 3 && down == 0) || (up == 0 && down == 3)) fault[i] = 1;
        else if (up == 3 && down == 3) fault[i] = 3;
        else fault[i] = 4;
    }
    Json fractions;
    Json counts = countsJson(fault, names, fractions);
    std::vector<long long> count(names.size(), 0);
    for (int f : fault) ++count[static_cast<std::size_t>(f)];
    const double d = spacingCount ? spacing / static_cast<double>(spacingCount) : kNaN;
    const double areaPerAtom = std::sqrt(3.0) / 2 * d * d;  // one atom of a close-packed plane
    Json result = Json::object();
    result["fault_type"] = toJson(fault);
    Json typeNames = Json::array();
    for (const auto& name : names) typeNames.push(name);
    result["type_names"] = typeNames;
    result["counts"] = counts;
    result["nearest_neighbor_distance_A"] = d;
    // An intrinsic stacking fault is two hcp planes, a twin boundary one.
    result["stacking_fault_area_A2"] = static_cast<double>(count[1]) / 2 * areaPerAtom;
    result["twin_boundary_area_A2"] = static_cast<double>(count[2]) * areaPerAtom;
    result["cutoff_A"] = m.geometry.cutoff;
    return result;
}
}
