#include "science/Clusters.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);

std::size_t root(std::vector<std::size_t>& parent, std::size_t i)
{
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}

double medianNearest(const std::vector<Vec3>& positions, const Mat3& cell, const Pbc& pbc)
{
    std::vector<double> nearest(positions.size(), HUGE_VAL);
    for (const auto& n : neighborList(positions, cell, pbc, 6.0))
        nearest[static_cast<std::size_t>(n.i)] = std::min(nearest[static_cast<std::size_t>(n.i)], norm(n.vector));
    std::vector<double> finite;
    for (double d : nearest) if (std::isfinite(d)) finite.push_back(d);
    if (finite.empty()) throw std::runtime_error("No neighbours within 6 Angstrom; set atomic_radius_A");
    std::nth_element(finite.begin(), finite.begin() + static_cast<std::ptrdiff_t>(finite.size() / 2), finite.end());
    return finite[finite.size() / 2];
}
}

Json clusterAnalysis(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    if (periodic && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    const Mat3 cell = p.has("cell") ? cellMatrix(p.array("cell")) : identity();
    const double cutoff = positive(p.number("cutoff_A", 3.0), "cutoff_A");
    std::vector<std::size_t> selected;
    if (p.has("mask")) {
        const auto mask = finiteArray(p.array("mask"), "mask", 1).values;
        if (mask.size() != positions.size()) throw std::runtime_error("mask needs one value per atom");
        for (std::size_t i = 0; i < mask.size(); ++i) if (mask[i] != 0) selected.push_back(i);
    } else {
        selected.resize(positions.size());
        std::iota(selected.begin(), selected.end(), 0);
    }
    const long long minimum = integer(p.number("min_size", 1), "min_size");
    std::vector<Vec3> chosen;
    for (std::size_t i : selected) chosen.push_back(positions[i]);
    const auto links = neighborList(chosen, cell, pbc, cutoff);
    std::vector<std::size_t> parent(chosen.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<std::vector<std::pair<std::size_t, Vec3>>> adjacency(chosen.size());
    for (const auto& link : links) {
        parent[root(parent, static_cast<std::size_t>(link.i))] = root(parent, static_cast<std::size_t>(link.j));
        adjacency[static_cast<std::size_t>(link.i)].push_back({static_cast<std::size_t>(link.j), link.vector});
    }
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t k = 0; k < chosen.size(); ++k) groups[root(parent, k)].push_back(k);
    std::vector<std::vector<std::size_t>> clusters;
    for (auto& [r, members] : groups)
        if (static_cast<long long>(members.size()) >= minimum) clusters.push_back(members);
    std::sort(clusters.begin(), clusters.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
    std::vector<int> clusterId(positions.size(), -1);
    Json list = Json::array();
    std::map<std::size_t, int> sizeCounts;
    for (std::size_t c = 0; c < clusters.size(); ++c) {
        const auto& members = clusters[c];
        // Unwrap across periodic boundaries along the bond network.
        std::map<std::size_t, Vec3> unwrapped;
        unwrapped[members.front()] = chosen[members.front()];
        std::queue<std::size_t> queue;
        queue.push(members.front());
        bool percolates = false;
        while (!queue.empty()) {
            const std::size_t k = queue.front();
            queue.pop();
            for (const auto& [j, vector] : adjacency[k]) {
                const Vec3 reached = add(unwrapped[k], vector);
                if (!unwrapped.count(j)) { unwrapped[j] = reached; queue.push(j); }
                else if (norm(sub(reached, unwrapped[j])) > 1e-6) percolates = true;
            }
        }
        Vec3 centroid{0, 0, 0};
        for (std::size_t k : members) { centroid = add(centroid, unwrapped[k]); clusterId[selected[k]] = static_cast<int>(c); }
        centroid = scale(centroid, 1.0 / static_cast<double>(members.size()));
        double gyration = 0;
        for (std::size_t k : members) { const Vec3 d = sub(unwrapped[k], centroid); gyration += dot(d, d); }
        gyration = std::sqrt(gyration / static_cast<double>(members.size()));
        Json entry = Json::object();
        entry["size"] = members.size();
        entry["centroid_A"] = toJson(centroid);
        entry["radius_of_gyration_A"] = gyration;
        entry["percolating"] = percolates;
        list.push(entry);
        ++sizeCounts[members.size()];
    }
    Json histogram = Json::array();
    for (const auto& [size, count] : sizeCounts) histogram.push(Json::array({Json(size), Json(count)}));
    Json result = Json::object();
    result["cluster_count"] = clusters.size();
    result["selected_atoms"] = selected.size();
    result["clusters"] = list;
    result["size_histogram"] = histogram;
    result["cluster_id"] = toJson(clusterId);
    return result;
}

Json voidAnalysis(const Parameters& p)
{
    std::vector<Vec3> positions = points(p.array("positions"), "positions");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool anyPeriodic = pbc[0] || pbc[1] || pbc[2];
    if (anyPeriodic && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    Mat3 cell = p.has("cell") ? cellMatrix(p.array("cell")) : identity();
    double radius = p.number("atomic_radius_A", 0.0);
    if (radius <= 0) radius = 0.5 * medianNearest(positions, cell, pbc);
    const double probe = p.number("probe_radius_A", 1.0);
    if (!(probe >= 0)) throw std::runtime_error("probe_radius_A must be nonnegative");
    const double spacing = positive(p.number("grid_spacing_A", 0.3), "grid_spacing_A");
    const double reach = radius + probe;
    // Domain: the cell (periodic axes) or a padded bounding box.
    Mat3 domain = cell;
    Vec3 origin{0, 0, 0};
    if (!p.has("cell")) {
        Vec3 low = positions[0], high = positions[0];
        for (const auto& x : positions) for (int k = 0; k < 3; ++k) { low[k] = std::min(low[k], x[k]); high[k] = std::max(high[k], x[k]); }
        const double margin = reach + 2 * spacing;
        domain = {{{high[0] - low[0] + 2 * margin, 0, 0}, {0, high[1] - low[1] + 2 * margin, 0}, {0, 0, high[2] - low[2] + 2 * margin}}};
        origin = {low[0] - margin, low[1] - margin, low[2] - margin};
    }
    const Mat3 inverseDomain = inverse(domain);
    std::array<int, 3> size{};
    for (int k = 0; k < 3; ++k) size[static_cast<std::size_t>(k)] = std::max(4, static_cast<int>(std::ceil(norm(domain[k]) / spacing)));
    const long long total = static_cast<long long>(size[0]) * size[1] * size[2];
    if (total > 60000000) throw std::runtime_error("Void grid too large; increase grid_spacing_A");
    // Atom images near the domain, binned by the exclusion reach.
    struct Image { Vec3 x; std::size_t atom; };
    std::vector<Image> images;
    int range[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        if (pbc[k]) range[k] = static_cast<int>(std::ceil(reach * norm({inverseDomain[0][k], inverseDomain[1][k], inverseDomain[2][k]}))) + 1;
    for (std::size_t a = 0; a < positions.size(); ++a) {
        Vec3 f = rowTimes(sub(positions[a], origin), inverseDomain);
        for (int k = 0; k < 3; ++k) if (pbc[k]) f[k] -= std::floor(f[k]);
        for (int i = -range[0]; i <= range[0]; ++i)
            for (int j = -range[1]; j <= range[1]; ++j)
                for (int k = -range[2]; k <= range[2]; ++k)
                    images.push_back({add(origin, rowTimes({f[0] + i, f[1] + j, f[2] + k}, domain)), a});
    }
    auto key = [&](const Vec3& x) {
        return std::array<long long, 3>{static_cast<long long>(std::floor(x[0] / reach)), static_cast<long long>(std::floor(x[1] / reach)), static_cast<long long>(std::floor(x[2] / reach))};
    };
    auto hash = [](const std::array<long long, 3>& k) { return (k[0] * 73856093LL) ^ (k[1] * 19349663LL) ^ (k[2] * 83492791LL); };
    std::unordered_map<long long, std::vector<std::size_t>> bins;
    for (std::size_t i = 0; i < images.size(); ++i) bins[hash(key(images[i].x))].push_back(i);
    auto pointAt = [&](int i, int j, int k) {
        return add(origin, rowTimes({(i + 0.5) / size[0], (j + 0.5) / size[1], (k + 0.5) / size[2]}, domain));
    };
    auto nearAtoms = [&](const Vec3& x, double distance, const std::function<void(std::size_t)>& visit) {
        const auto c = key(x);
        const int span = static_cast<int>(std::ceil(distance / reach));
        for (long long dx = -span; dx <= span; ++dx)
            for (long long dy = -span; dy <= span; ++dy)
                for (long long dz = -span; dz <= span; ++dz) {
                    const auto found = bins.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                    if (found == bins.end()) continue;
                    for (std::size_t index : found->second)
                        if (norm(sub(images[index].x, x)) < distance) visit(index);
                }
    };
    std::vector<char> accessible(static_cast<std::size_t>(total), 0);
    auto index = [&](int i, int j, int k) { return (static_cast<std::size_t>(i) * size[1] + j) * size[2] + k; };
    for (int i = 0; i < size[0]; ++i) {
        taskProgress(0.7 * i / size[0]);
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k) {
                bool blocked = false;
                const Vec3 x = pointAt(i, j, k);
                const auto c = key(x);
                for (long long dx = -1; dx <= 1 && !blocked; ++dx)
                    for (long long dy = -1; dy <= 1 && !blocked; ++dy)
                        for (long long dz = -1; dz <= 1 && !blocked; ++dz) {
                            const auto found = bins.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                            if (found == bins.end()) continue;
                            for (std::size_t index2 : found->second)
                                if (dot(sub(images[index2].x, x), sub(images[index2].x, x)) < reach * reach) { blocked = true; break; }
                        }
                accessible[index(i, j, k)] = !blocked;
            }
    }
    // Connected components with periodic wrapping; open boundaries are exterior.
    std::vector<int> label(static_cast<std::size_t>(total), -1);
    const double voxel = cellVolume(domain) / static_cast<double>(total);
    struct Void { std::size_t points = 0; Vec3 sum{0, 0, 0}; bool exterior = false; };
    std::vector<Void> voids;
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k) {
                if (!accessible[index(i, j, k)] || label[index(i, j, k)] >= 0) continue;
                Void v;
                const int id = static_cast<int>(voids.size());
                std::queue<std::array<int, 6>> queue;  // grid index and unwrapped image offsets
                queue.push({i, j, k, 0, 0, 0});
                label[index(i, j, k)] = id;
                while (!queue.empty()) {
                    const auto q = queue.front();
                    queue.pop();
                    ++v.points;
                    v.sum = add(v.sum, {(q[0] + 0.5) / size[0] + q[3], (q[1] + 0.5) / size[1] + q[4], (q[2] + 0.5) / size[2] + q[5]});
                    static const int steps[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
                    for (const auto& s : steps) {
                        int g[3] = {q[0] + s[0], q[1] + s[1], q[2] + s[2]};
                        int offset[3] = {q[3], q[4], q[5]};
                        bool outside = false;
                        for (int axis = 0; axis < 3; ++axis) {
                            if (g[axis] >= 0 && g[axis] < size[static_cast<std::size_t>(axis)]) continue;
                            if (!pbc[axis]) { outside = true; continue; }
                            const int wrap = g[axis] < 0 ? -1 : 1;
                            g[axis] -= wrap * size[static_cast<std::size_t>(axis)];
                            offset[axis] += wrap;
                        }
                        if (outside) { v.exterior = true; continue; }
                        const std::size_t n = index(g[0], g[1], g[2]);
                        if (!accessible[n] || label[n] >= 0) continue;
                        label[n] = id;
                        queue.push({g[0], g[1], g[2], offset[0], offset[1], offset[2]});
                    }
                }
                voids.push_back(v);
            }
    taskProgress(0.85);
    // Keep interior voids, largest first, and record the atoms lining them.
    std::vector<int> order;
    for (std::size_t v = 0; v < voids.size(); ++v) if (!voids[v].exterior) order.push_back(static_cast<int>(v));
    std::sort(order.begin(), order.end(), [&](int a, int b) { return voids[static_cast<std::size_t>(a)].points > voids[static_cast<std::size_t>(b)].points; });
    std::map<int, int> rank;
    for (std::size_t r = 0; r < order.size(); ++r) rank[order[r]] = static_cast<int>(r);
    std::vector<int> lining(positions.size(), -1);
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k) {
                const int id = label[index(i, j, k)];
                if (id < 0 || !rank.count(id)) continue;
                nearAtoms(pointAt(i, j, k), reach + 1.5 * spacing, [&](std::size_t image) {
                    int& owner = lining[images[image].atom];
                    if (owner < 0) owner = rank[id];
                });
            }
    Json list = Json::array();
    double totalVolume = 0;
    for (int id : order) {
        const Void& v = voids[static_cast<std::size_t>(id)];
        const double volume = static_cast<double>(v.points) * voxel;
        totalVolume += volume;
        const Vec3 centroid = add(origin, rowTimes(scale(v.sum, 1.0 / static_cast<double>(v.points)), domain));
        const double equivalent = std::cbrt(3 * volume / (4 * kPi));
        Json entry = Json::object();
        entry["accessible_volume_A3"] = volume;
        entry["centroid_A"] = toJson(centroid);
        entry["equivalent_radius_A"] = equivalent;
        entry["pore_radius_estimate_A"] = equivalent + reach;
        entry["grid_points"] = v.points;
        list.push(entry);
    }
    Json result = Json::object();
    result["void_count"] = order.size();
    result["voids"] = list;
    result["accessible_volume_A3"] = totalVolume;
    result["void_fraction"] = p.has("cell") ? totalVolume / cellVolume(cell) : Json();
    result["atomic_radius_A"] = radius;
    result["probe_radius_A"] = probe;
    result["grid"] = Json::array({size[0], size[1], size[2]});
    result["lining_void_id"] = toJson(lining);
    return result;
}
}
