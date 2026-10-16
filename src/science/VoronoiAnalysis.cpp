#include "science/VoronoiAnalysis.h"
#include "science/LocalGeometry.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cctype>
#include <filesystem>
#include <sstream>
#include <random>
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

// Every atom's cell; cells whose vertices reach beyond half the search
// radius may be cut by farther atoms, so the search widens for them.
// cutoff: the initial search radius (0 = automatic); returns the last one used.
std::vector<CellResult> allCells(const AtomSet& a, double& cutoff)
{
    const std::size_t n = a.positions.size();
    const bool periodic = a.fullyPeriodic();
    if (cutoff <= 0) cutoff = 2.5 * meanAtomicSpacing(a);
    const double eps = 1e-9 * meanAtomicSpacing(a);
    std::vector<CellResult> cells(n);
    std::vector<char> done(n, 0);
    for (int attempt = 0; attempt < 6; ++attempt) {
        const auto neighbors = sortedNeighborVectors(a, cutoff);
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
    return cells;
}

// ---------------------------------------------------------------- Voronoi topology
// Faces of a cell as cycles of vertex ids, counter-clockwise seen from outside.
std::vector<std::vector<int>> topologyOf(const CellResult& cell, double eps, int& vertexCount)
{
    std::vector<Vec3> vertices;
    std::vector<std::vector<int>> faces;
    for (const auto& face : cell.faces) {
        std::vector<int> ids;
        for (const auto& point : face.points) {
            int id = -1;
            for (std::size_t k = 0; k < vertices.size(); ++k)
                if (norm(sub(vertices[k], point)) <= eps) { id = static_cast<int>(k); break; }
            if (id < 0) { id = static_cast<int>(vertices.size()); vertices.push_back(point); }
            if (ids.empty() || ids.back() != id) ids.push_back(id);
        }
        while (ids.size() > 1 && ids.front() == ids.back()) ids.pop_back();
        if (ids.size() >= 3) faces.push_back(std::move(ids));
    }
    vertexCount = static_cast<int>(vertices.size());
    return faces;
}

// Canonical Weinberg code of a polyhedral graph (Weinberg, IEEE Trans.
// Circuit Theory 13, 142 (1966)), as in VoroTop (Lazar, Han and Srolovitz,
// PNAS 112, 13212 (2015)): the lexicographically smallest vertex-label
// sequence of the Weinberg traversal over all starting directed edges and
// both orientations, so mirror images share a code. Empty when the faces do
// not form a closed polyhedral surface.
std::vector<int> weinbergCode(const std::vector<std::vector<int>>& faces, int n)
{
    if (n < 4 || faces.size() < 4) return {};
    std::vector<int> best;
    for (int orientation = 0; orientation < 2; ++orientation) {
        // turn[u * n + v] = w: the face left of u->v continues v->w.
        std::vector<int> turn(static_cast<std::size_t>(n) * n, -1);
        std::vector<std::pair<int, int>> edges;
        for (const auto& original : faces) {
            std::vector<int> face = original;
            if (orientation) std::reverse(face.begin(), face.end());
            const std::size_t m = face.size();
            for (std::size_t k = 0; k < m; ++k) {
                const int u = face[k], v = face[(k + 1) % m], w = face[(k + 2) % m];
                int& slot = turn[static_cast<std::size_t>(u) * n + v];
                if (slot >= 0) return {};  // a directed edge in two faces: not a closed surface
                slot = w;
                edges.push_back({u, v});
            }
        }
        for (const auto& [u, v] : edges)
            if (turn[static_cast<std::size_t>(v) * n + u] < 0) return {};
        const std::size_t total = edges.size();
        std::vector<char> used(static_cast<std::size_t>(n) * n);
        std::vector<int> label(static_cast<std::size_t>(n)), code;
        for (const auto& [a, b] : edges) {
            std::fill(used.begin(), used.end(), 0);
            std::fill(label.begin(), label.end(), 0);
            code.clear();
            int next = 1;
            label[static_cast<std::size_t>(a)] = next;
            code.push_back(next);
            int u = a, v = b;
            std::size_t count = 0;
            bool worse = false;
            while (true) {
                used[static_cast<std::size_t>(u) * n + v] = 1;
                ++count;
                int w;
                if (!label[static_cast<std::size_t>(v)]) {
                    label[static_cast<std::size_t>(v)] = ++next;
                    code.push_back(next);
                    w = turn[static_cast<std::size_t>(u) * n + v];
                } else {
                    code.push_back(label[static_cast<std::size_t>(v)]);
                    if (count == total) break;
                    if (!used[static_cast<std::size_t>(v) * n + u]) w = u;
                    else {
                        w = turn[static_cast<std::size_t>(u) * n + v];
                        int guard = 0;
                        while (used[static_cast<std::size_t>(v) * n + w] && guard++ < n) w = turn[static_cast<std::size_t>(w) * n + v];
                        if (used[static_cast<std::size_t>(v) * n + w]) break;
                    }
                }
                // Abandon a traversal that is already larger than the best code.
                if (!best.empty()) {
                    const std::size_t k = code.size() - 1;
                    if (k < best.size() && code[k] != best[k]) {
                        if (code[k] > best[k] && std::equal(code.begin(), code.begin() + static_cast<std::ptrdiff_t>(k), best.begin())) { worse = true; break; }
                    }
                }
                u = v;
                v = w;
            }
            if (!worse && count == total && (best.empty() || code < best)) best = code;
        }
    }
    return best;
}

std::string codeKey(const std::vector<int>& code)
{
    std::string key;
    for (std::size_t k = 0; k < code.size(); ++k) key += (k ? "," : "") + std::to_string(code[k]);
    return key;
}

// Canonical codes of every atom's cell (empty for open cells).
std::vector<std::vector<int>> cellCodes(const AtomSet& a)
{
    double cutoff = 0;
    const auto cells = allCells(a, cutoff);
    const double eps = 1e-6 * meanAtomicSpacing(a);
    std::vector<std::vector<int>> codes(cells.size());
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (cells[i].open) continue;
        int n = 0;
        const auto faces = topologyOf(cells[i], eps, n);
        codes[i] = weinbergCode(faces, n);
    }
    return codes;
}

struct Filter
{
    std::vector<std::string> names = {"other"};
    std::map<std::string, int> types;  // code key -> type
};

// Topologies of thermally perturbed fcc, bcc and hcp crystals (Gaussian
// displacements of `noise` times the nearest-neighbour distance), the way
// VoroTop filters are sampled. Noise levels spread from half to one and a
// half times the given one; a code shared by lattices goes to the one in
// which it is most frequent.
const Filter& builtInFilter(double noise, int samples)
{
    static std::map<std::pair<long long, int>, Filter> cache;
    const std::pair<long long, int> key(std::llround(noise * 1e6), samples);
    if (const auto found = cache.find(key); found != cache.end()) return found->second;
    Filter filter;
    filter.names = {"other", "fcc", "bcc", "hcp"};
    std::map<std::string, std::array<long long, 4>> seen;  // code -> occurrences per lattice
    std::array<long long, 4> sampled{};                     // cells sampled per lattice
    std::mt19937 random(12345);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (int type = 1; type <= 3; ++type) {
        AtomSet crystal;
        std::vector<Vec3> basis;
        Mat3 unit{};
        int repeat[3] = {4, 4, 4};
        double nearest = 1;
        if (type == 1) { unit = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; basis = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}}; nearest = std::sqrt(0.5); }
        else if (type == 2) { unit = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; basis = {{0, 0, 0}, {.5, .5, .5}}; nearest = std::sqrt(0.75); repeat[0] = repeat[1] = repeat[2] = 5; }
        else { unit = {{{1, 0, 0}, {-0.5, std::sqrt(3.0) / 2, 0}, {0, 0, std::sqrt(8.0 / 3)}}}; basis = {{0, 0, 0}, {1.0 / 3, 2.0 / 3, .5}}; repeat[2] = 3; repeat[0] = repeat[1] = 5; }
        for (int i = 0; i < repeat[0]; ++i)
            for (int j = 0; j < repeat[1]; ++j)
                for (int k = 0; k < repeat[2]; ++k)
                    for (const auto& b : basis) crystal.positions.push_back(rowTimes({b[0] + i, b[1] + j, b[2] + k}, unit));
        for (std::size_t r = 0; r < 3; ++r) crystal.cell[r] = scale(unit[r], repeat[r]);
        crystal.pbc = {true, true, true};
        crystal.hasCell = true;
        const auto ideal = crystal.positions;
        for (int sample = 0; sample < samples; ++sample) {
            const double level = noise * (0.5 + static_cast<double>(sample) / (samples - 1));
            taskCheckpoint();
            for (std::size_t i = 0; i < ideal.size(); ++i)
                for (std::size_t d = 0; d < 3; ++d) crystal.positions[i][d] = ideal[i][d] + level * nearest * gauss(random);
            for (const auto& code : cellCodes(crystal)) {
                if (code.empty()) continue;
                ++seen[codeKey(code)][static_cast<std::size_t>(type)];
                ++sampled[static_cast<std::size_t>(type)];
            }
        }
    }
    // Each code goes to the lattice in which it is most frequent (per sampled cell).
    for (const auto& [k, counts] : seen) {
        int best = 0;
        double bestFrequency = 0;
        for (int type = 1; type <= 3; ++type) {
            const double frequency = static_cast<double>(counts[static_cast<std::size_t>(type)]) / static_cast<double>(std::max(1LL, sampled[static_cast<std::size_t>(type)]));
            if (frequency > bestFrequency) { bestFrequency = frequency; best = type; }
        }
        if (best) filter.types[k] = best;
    }
    return cache[key] = filter;
}

// A VoroTop filter file: "* id name" lines name the types, "id (1,2,3,...)" lines list codes.
Filter readFilter(const std::filesystem::path& path)
{
    const std::string text = readTextFile(path);
    Filter filter;
    std::map<int, std::string> names;
    std::vector<std::pair<int, std::string>> codes;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        const auto start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos || line[start] == '#') continue;
        if (line[start] == '*') {
            std::istringstream words(line.substr(start + 1));
            int id = 0;
            std::string name;
            if (words >> id >> name) names[id] = name;
            continue;
        }
        const auto open = line.find('('), close = line.find(')');
        if (open == std::string::npos || close == std::string::npos || close < open) continue;
        const int id = std::atoi(line.substr(start, open - start).c_str());
        std::string inner = line.substr(open + 1, close - open - 1);
        inner.erase(std::remove_if(inner.begin(), inner.end(), [](char c) { return std::isspace(static_cast<unsigned char>(c)); }), inner.end());
        codes.push_back({id, inner});
    }
    if (codes.empty()) throw std::runtime_error("No codes in the VoroTop filter " + path.u8string());
    std::map<int, int> index;
    for (const auto& [id, name] : names) { index[id] = static_cast<int>(filter.names.size()); filter.names.push_back(name); }
    for (const auto& [id, code] : codes) {
        if (!index.count(id)) { index[id] = static_cast<int>(filter.names.size()); filter.names.push_back("type " + std::to_string(id)); }
        filter.types[code] = index[id];
    }
    return filter;
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
    const AtomSet a = readAtomSet(p);
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    const double areaThreshold = p.number("face_area_threshold_A2", 0.0);
    const double relativeThreshold = p.number("relative_face_area_threshold", 0.0);
    const double edgeThreshold = p.number("edge_length_threshold_A", 0.0);
    if (areaThreshold < 0 || relativeThreshold < 0 || edgeThreshold < 0) throw std::runtime_error("Thresholds must be nonnegative");
    double cutoff = p.number("cutoff_A", 0.0);
    const std::vector<CellResult> cells = allCells(a, cutoff);
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

Json voroTop(const Parameters& p)
{
    const AtomSet a = readAtomSet(p);
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    const double noise = p.number("noise", 0.06);
    if (!(noise > 0 && noise < 0.3)) throw std::runtime_error("noise is a fraction of the nearest-neighbour distance between 0 and 0.3");
    const Filter filter = p.has("filter_file") ? readFilter(std::filesystem::u8path(p.json("filter_file").string())) : builtInFilter(noise, static_cast<int>(integer(p.number("samples", 60), "samples", 2)));
    const auto codes = cellCodes(a);
    std::vector<int> type(n, 0), topology(n, -1);
    std::map<std::string, int> distinct;
    std::vector<std::pair<std::string, long long>> topologies;
    for (std::size_t i = 0; i < n; ++i) {
        if (codes[i].empty()) continue;
        const std::string k = codeKey(codes[i]);
        if (const auto found = filter.types.find(k); found != filter.types.end()) type[i] = found->second;
        auto [entry, added] = distinct.insert({k, static_cast<int>(topologies.size())});
        if (added) topologies.push_back({k, 0});
        topology[i] = entry->second;
        ++topologies[static_cast<std::size_t>(entry->second)].second;
    }
    Json result = Json::object();
    addClassification(result, "structure_type", type, filter.names);
    std::vector<std::size_t> order(topologies.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return topologies[x].second > topologies[y].second; });
    Json common = Json::array();
    for (std::size_t k = 0; k < std::min<std::size_t>(10, order.size()); ++k) {
        const auto& [code, count] = topologies[order[k]];
        Json item = Json::object();
        item["atoms"] = static_cast<double>(count);
        item["weinberg_code"] = "(" + code + ")";
        const auto found = filter.types.find(code);
        item["type"] = filter.names[static_cast<std::size_t>(found == filter.types.end() ? 0 : found->second)];
        common.push(item);
    }
    result["topology_id"] = toJson(topology);
    result["distinct_topologies"] = static_cast<double>(topologies.size());
    std::size_t withoutCode = 0;
    for (const auto& code : codes) withoutCode += code.empty();
    result["cells_without_topology"] = static_cast<double>(withoutCode);
    result["most_common_topologies"] = common;
    result["filter"] = p.has("filter_file") ? p.json("filter_file").string() : std::string("built-in (perturbed fcc, bcc, hcp)");
    result["filter_codes"] = static_cast<double>(filter.types.size());
    if (p.boolean("export_filter", false)) {
        // The filter in VoroTop's format, to reuse with filter_file.
        std::string text = "# AtomForge VoroTop filter\n";
        for (std::size_t t = 1; t < filter.names.size(); ++t) text += "*\t" + std::to_string(t) + "\t" + filter.names[t] + "\n";
        for (const auto& [code, t] : filter.types) text += std::to_string(t) + "\t(" + code + ")\n";
        Json files = Json::object();
        files["vorotop.filter"] = text;
        result["files"] = files;
    }
    return result;
}

Json adaptiveCna(const Parameters& p)
{
    const AtomSet a = readAtomSet(p);
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    const double fixed = p.number("cutoff_A", 0.0);
    if (fixed < 0) throw std::runtime_error("cutoff_A must be nonnegative");
    const double factor = (1 + std::sqrt(2.0)) / 2;
    const auto neighbors = sortedNeighborVectors(a, fixed > 0 ? fixed : 1.8 * meanAtomicSpacing(a));
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
    Json result = Json::object();
    addClassification(result, "structure_type", type, names);
    result["mode"] = fixed > 0 ? "fixed cutoff" : "adaptive";
    return result;
}
}
