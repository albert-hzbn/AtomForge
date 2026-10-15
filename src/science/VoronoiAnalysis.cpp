#include "science/VoronoiAnalysis.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

namespace atomforge::science
{
namespace
{
const double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Atoms
{
    std::vector<Vec3> positions;
    Mat3 cell{};
    Pbc pbc{};
    bool hasCell = false;
};

Atoms atomsOf(const Parameters& p)
{
    Atoms a;
    a.positions = points(p.array("positions"), "positions");
    a.pbc = p.pbc("pbc", {false, false, false});
    a.hasCell = p.has("cell");
    if ((a.pbc[0] || a.pbc[1] || a.pbc[2]) && !a.hasCell) throw std::runtime_error("Periodic analysis requires a cell");
    a.cell = a.hasCell ? cellMatrix(p.array("cell")) : identity();
    return a;
}

// Mean spacing (V/N)^(1/3) from the cell or the padded bounding box.
double meanSpacing(const Atoms& a)
{
    double volume = 0;
    if (a.hasCell) volume = cellVolume(a.cell);
    else {
        Vec3 low = a.positions[0], high = a.positions[0];
        for (const auto& x : a.positions)
            for (std::size_t d = 0; d < 3; ++d) { low[d] = std::min(low[d], x[d]); high[d] = std::max(high[d], x[d]); }
        volume = 1;
        for (std::size_t d = 0; d < 3; ++d) volume *= high[d] - low[d] + 2.0;
    }
    return std::cbrt(volume / static_cast<double>(a.positions.size()));
}

// Neighbour vectors of every atom within the cutoff, nearest first.
std::vector<std::vector<Vec3>> sortedNeighbors(const Atoms& a, double cutoff)
{
    std::vector<std::vector<Vec3>> lists(a.positions.size());
    for (const auto& link : neighborList(a.positions, a.cell, a.pbc, cutoff)) {
        if (norm(link.vector) < 1e-10) throw std::runtime_error("Coincident atoms make the analysis undefined");
        lists[static_cast<std::size_t>(link.i)].push_back(link.vector);
    }
    for (auto& list : lists) std::sort(list.begin(), list.end(), [](const Vec3& x, const Vec3& y) { return dot(x, x) < dot(y, y); });
    return lists;
}

// ---------------------------------------------------------------- Voronoi cell
struct Face
{
    int neighbor = -1;  // -1: a face of the initial bounding cube
    std::vector<Vec3> points;
};

double polygonArea(const std::vector<Vec3>& pts)
{
    Vec3 sum{};
    for (std::size_t k = 1; k + 1 < pts.size(); ++k) sum = add(sum, cross(sub(pts[k], pts[0]), sub(pts[k + 1], pts[0])));
    return 0.5 * norm(sum);
}

// Cuts the convex cell with the half-space n.x <= d. Returns false when the plane misses it.
bool clip(std::vector<Face>& cell, const Vec3& n, double d, int neighbor, double eps)
{
    bool cuts = false;
    for (const auto& face : cell) {
        for (const auto& point : face.points)
            if (dot(n, point) - d > eps) { cuts = true; break; }
        if (cuts) break;
    }
    if (!cuts) return false;
    std::vector<Face> out;
    std::vector<Vec3> cap;
    for (const auto& face : cell) {
        Face clipped;
        clipped.neighbor = face.neighbor;
        const std::size_t m = face.points.size();
        for (std::size_t k = 0; k < m; ++k) {
            const Vec3& current = face.points[k];
            const Vec3& next = face.points[(k + 1) % m];
            const double dc = dot(n, current) - d, dn = dot(n, next) - d;
            const bool inside = dc <= eps, nextInside = dn <= eps;
            if (inside) {
                clipped.points.push_back(current);
                if (std::abs(dc) <= eps) cap.push_back(current);
            }
            if (inside != nextInside && std::abs(dc - dn) > 0) {
                const Vec3 crossing = add(current, scale(sub(next, current), dc / (dc - dn)));
                clipped.points.push_back(crossing);
                cap.push_back(crossing);
            }
        }
        // Drop repeated points.
        std::vector<Vec3> unique;
        for (const auto& point : clipped.points)
            if (unique.empty() || norm(sub(point, unique.back())) > eps) unique.push_back(point);
        while (unique.size() > 1 && norm(sub(unique.front(), unique.back())) <= eps) unique.pop_back();
        if (unique.size() >= 3) { clipped.points = std::move(unique); out.push_back(std::move(clipped)); }
    }
    std::vector<Vec3> capUnique;
    for (const auto& point : cap) {
        bool seen = false;
        for (const auto& q : capUnique) if (norm(sub(point, q)) <= eps) { seen = true; break; }
        if (!seen) capUnique.push_back(point);
    }
    if (capUnique.size() >= 3) {
        Vec3 centre{};
        for (const auto& point : capUnique) centre = add(centre, scale(point, 1.0 / static_cast<double>(capUnique.size())));
        const Vec3 helper = std::abs(n[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        const Vec3 u = scale(cross(n, helper), 1 / norm(cross(n, helper)));
        const Vec3 w = cross(n, u);
        std::sort(capUnique.begin(), capUnique.end(), [&](const Vec3& x, const Vec3& y) {
            const Vec3 dx = sub(x, centre), dy = sub(y, centre);
            return std::atan2(dot(dx, w), dot(dx, u)) < std::atan2(dot(dy, w), dot(dy, u));
        });
        out.push_back({neighbor, std::move(capUnique)});
    }
    cell = std::move(out);
    return true;
}

std::vector<Face> cube(double half)
{
    const double h = half;
    const auto v = [h](double x, double y, double z) { return Vec3{x * h, y * h, z * h}; };
    return {
        {-1, {v(-1, -1, -1), v(-1, 1, -1), v(1, 1, -1), v(1, -1, -1)}},
        {-1, {v(-1, -1, 1), v(1, -1, 1), v(1, 1, 1), v(-1, 1, 1)}},
        {-1, {v(-1, -1, -1), v(1, -1, -1), v(1, -1, 1), v(-1, -1, 1)}},
        {-1, {v(-1, 1, -1), v(-1, 1, 1), v(1, 1, 1), v(1, 1, -1)}},
        {-1, {v(-1, -1, -1), v(-1, -1, 1), v(-1, 1, 1), v(-1, 1, -1)}},
        {-1, {v(1, -1, -1), v(1, 1, -1), v(1, 1, 1), v(1, -1, 1)}},
    };
}

struct CellResult
{
    std::vector<Face> faces;
    double radius = 0;   // largest vertex distance
    bool open = false;   // still touches the bounding cube
};

CellResult voronoiCell(const std::vector<Vec3>& neighbors, double half, double eps)
{
    CellResult result;
    result.faces = cube(half);
    double radius = std::sqrt(3.0) * half;
    for (std::size_t k = 0; k < neighbors.size(); ++k) {
        const double r = norm(neighbors[k]);
        // Atoms farther than twice the cell radius cannot cut it.
        if (r > 2 * radius) break;
        clip(result.faces, scale(neighbors[k], 1 / r), r / 2, static_cast<int>(k), eps);
        if (k % 8 == 7 || k + 1 == neighbors.size()) {
            radius = 0;
            for (const auto& face : result.faces)
                for (const auto& point : face.points) radius = std::max(radius, norm(point));
        }
    }
    result.radius = 0;
    for (const auto& face : result.faces) {
        if (face.neighbor < 0) result.open = true;
        for (const auto& point : face.points) result.radius = std::max(result.radius, norm(point));
    }
    return result;
}

// ---------------------------------------------------------------- CNA
struct Signature { int common, bonds, chain; };

// Common-neighbour signature of bond j in a set of neighbour vectors.
Signature signature(const std::vector<Vec3>& v, std::size_t j, double cutoff)
{
    std::vector<std::size_t> common;
    for (std::size_t k = 0; k < v.size(); ++k)
        if (k != j && norm(sub(v[k], v[j])) < cutoff) common.push_back(k);
    std::vector<std::pair<std::size_t, std::size_t>> bonds;
    for (std::size_t a = 0; a < common.size(); ++a)
        for (std::size_t b = a + 1; b < common.size(); ++b)
            if (norm(sub(v[common[a]], v[common[b]])) < cutoff) bonds.push_back({a, b});
    // Bonds in the largest connected cluster of bonds.
    std::vector<std::size_t> parent(common.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](std::size_t x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
    for (const auto& [a, b] : bonds) parent[find(a)] = find(b);
    std::map<std::size_t, int> perCluster;
    for (const auto& bond : bonds) ++perCluster[find(bond.first)];
    int chain = 0;
    for (const auto& entry : perCluster) chain = std::max(chain, entry.second);
    return {static_cast<int>(common.size()), static_cast<int>(bonds.size()), chain};
}

// 1 fcc, 2 hcp, 3 bcc, 4 icosahedral, 0 other from a neighbour set and its cutoff.
int cnaType(const std::vector<Vec3>& v, double cutoff)
{
    int n421 = 0, n422 = 0, n444 = 0, n555 = 0, n666 = 0;
    for (std::size_t j = 0; j < v.size(); ++j) {
        const Signature s = signature(v, j, cutoff);
        if (s.common == 4 && s.bonds == 2 && s.chain == 1) ++n421;
        else if (s.common == 4 && s.bonds == 2 && s.chain == 2) ++n422;
        else if (s.common == 4 && s.bonds == 4 && s.chain == 4) ++n444;
        else if (s.common == 5 && s.bonds == 5 && s.chain == 5) ++n555;
        else if (s.common == 6 && s.bonds == 6 && s.chain == 6) ++n666;
    }
    if (v.size() == 12 && n421 == 12) return 1;
    if (v.size() == 12 && n421 == 6 && n422 == 6) return 2;
    if (v.size() == 12 && n555 == 12) return 4;
    if (v.size() == 14 && n444 == 6 && n666 == 8) return 3;
    return 0;
}
}

Json voronoiAnalysis(const Parameters& p)
{
    const Atoms a = atomsOf(p);
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    const double areaThreshold = p.number("face_area_threshold_A2", 0.0);
    const double relativeThreshold = p.number("relative_face_area_threshold", 0.0);
    const double edgeThreshold = p.number("edge_length_threshold_A", 0.0);
    if (areaThreshold < 0 || relativeThreshold < 0 || edgeThreshold < 0) throw std::runtime_error("Thresholds must be nonnegative");
    const bool periodic = a.pbc[0] && a.pbc[1] && a.pbc[2];
    double cutoff = p.number("cutoff_A", 0.0);
    if (cutoff <= 0) cutoff = 2.5 * meanSpacing(a);
    const double eps = 1e-9 * meanSpacing(a);
    std::vector<CellResult> cells(n);
    std::vector<char> done(n, 0);
    // Cells whose vertices reach beyond half the search radius may be cut by
    // farther atoms: widen the search for them.
    for (int attempt = 0; attempt < 6; ++attempt) {
        const auto neighbors = sortedNeighbors(a, cutoff);
        std::size_t pending = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (done[i]) continue;
            if (i % 256 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
            cells[i] = voronoiCell(neighbors[i], cutoff, eps);
            const bool complete = !cells[i].open && 2 * cells[i].radius <= cutoff;
            if (complete || (cells[i].open && !periodic && attempt > 0)) done[i] = 1;
            else ++pending;
        }
        if (!pending) break;
        cutoff *= 1.6;
    }
    std::vector<double> volume(n, kNaN), maxOrder(n, kNaN);
    std::vector<int> coordination(n, 0);
    std::vector<bool> open(n, false);
    Json indices = Json::array();
    std::map<std::string, long long> indexCounts;
    double total = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& cell = cells[i];
        open[i] = cell.open;
        std::vector<int> index(10, 0);
        double surface = 0;
        for (const auto& face : cell.faces) surface += polygonArea(face.points);
        double v = 0;
        int largest = 0;
        for (const auto& face : cell.faces) {
            const double area = polygonArea(face.points);
            double signed6 = 0;
            for (std::size_t k = 1; k + 1 < face.points.size(); ++k) signed6 += dot(face.points[0], cross(face.points[k], face.points[k + 1]));
            v += std::abs(signed6) / 6;
            if (face.neighbor < 0 || area <= areaThreshold || area <= relativeThreshold * surface) continue;
            ++coordination[i];
            int edges = 0;
            for (std::size_t k = 0; k < face.points.size(); ++k)
                edges += norm(sub(face.points[(k + 1) % face.points.size()], face.points[k])) >= edgeThreshold;
            if (edges >= 3 && edges <= 12) ++index[static_cast<std::size_t>(edges - 3)];
            largest = std::max(largest, edges);
        }
        if (!cell.open) { volume[i] = v; total += v; }
        maxOrder[i] = largest;
        Json row = Json::array();
        for (int c : index) row.push(static_cast<double>(c));
        indices.push(row);
        const std::string key = "<" + std::to_string(index[0]) + "," + std::to_string(index[1]) + "," + std::to_string(index[2]) + "," + std::to_string(index[3]) + ">";
        ++indexCounts[key];
    }
    std::vector<std::pair<long long, std::string>> ranked;
    for (const auto& [key, count] : indexCounts) ranked.push_back({count, key});
    std::sort(ranked.begin(), ranked.end(), [](const auto& x, const auto& y) { return x.first != y.first ? x.first > y.first : x.second < y.second; });
    Json common = Json::array();
    for (std::size_t k = 0; k < std::min<std::size_t>(20, ranked.size()); ++k) {
        Json item = Json::object();
        item["index"] = ranked[k].second;
        item["atoms"] = static_cast<double>(ranked[k].first);
        common.push(item);
    }
    Json result = Json::object();
    result["atomic_volume_A3"] = toJson(volume);
    result["coordination"] = toJson(coordination);
    result["voronoi_index"] = indices;
    result["index_columns"] = "faces with 3, 4, ..., 12 edges";
    result["max_face_order"] = toJson(maxOrder);
    result["open_cell"] = toJson(open);
    result["most_common_indices"] = common;
    result["total_volume_A3"] = total;
    if (a.hasCell) result["cell_volume_A3"] = cellVolume(a.cell);
    result["search_radius_A"] = cutoff;
    return result;
}

Json adaptiveCna(const Parameters& p)
{
    const Atoms a = atomsOf(p);
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    const double fixed = p.number("cutoff_A", 0.0);
    if (fixed < 0) throw std::runtime_error("cutoff_A must be nonnegative");
    const double factor = (1 + std::sqrt(2.0)) / 2;
    const auto neighbors = sortedNeighbors(a, fixed > 0 ? fixed : 1.8 * meanSpacing(a));
    static const std::vector<std::string> names = {"other", "fcc", "hcp", "bcc", "icosahedral"};
    std::vector<int> type(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (i % 512 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
        const auto& v = neighbors[i];
        if (fixed > 0) { type[i] = cnaType(v, fixed); continue; }
        if (v.size() >= 12) {
            double mean = 0;
            for (std::size_t k = 0; k < 12; ++k) mean += norm(v[k]) / 12;
            const std::vector<Vec3> shell(v.begin(), v.begin() + 12);
            type[i] = cnaType(shell, factor * mean);
        }
        if (type[i] == 0 && v.size() >= 14) {
            double first = 0, second = 0;
            for (std::size_t k = 0; k < 8; ++k) first += norm(v[k]) / 8;
            for (std::size_t k = 8; k < 14; ++k) second += norm(v[k]) / 6;
            const double cutoff = factor * (2 / std::sqrt(3.0) * first + second) / 2;
            const std::vector<Vec3> shell(v.begin(), v.begin() + 14);
            if (cnaType(shell, cutoff) == 3) type[i] = 3;
        }
    }
    std::vector<long long> counts(names.size(), 0);
    for (int t : type) ++counts[static_cast<std::size_t>(t)];
    Json result = Json::object(), countJson = Json::object(), typeNames = Json::array(), fractions = Json::object();
    for (std::size_t t = 0; t < names.size(); ++t) {
        countJson[names[t]] = static_cast<double>(counts[t]);
        fractions[names[t]] = static_cast<double>(counts[t]) / static_cast<double>(n);
        typeNames.push(names[t]);
    }
    result["structure_type"] = toJson(type);
    result["type_names"] = typeNames;
    result["counts"] = countJson;
    result["fractions"] = fractions;
    result["mode"] = fixed > 0 ? "fixed cutoff" : "adaptive";
    return result;
}
}
