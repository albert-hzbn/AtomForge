#include "science/LocalEnvironment.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);
const double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Atoms
{
    std::vector<Vec3> positions;
    Mat3 cell{};
    Pbc pbc{};
};

Atoms atomsOf(const Parameters& p)
{
    Atoms a;
    a.positions = points(p.array("positions"), "positions");
    a.pbc = p.pbc("pbc", {false, false, false});
    if ((a.pbc[0] || a.pbc[1] || a.pbc[2]) && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    a.cell = p.has("cell") ? cellMatrix(p.array("cell")) : identity();
    return a;
}

std::vector<std::vector<std::size_t>> byCenter(const std::vector<Neighbor>& neighbors, std::size_t count)
{
    std::vector<std::vector<std::size_t>> result(count);
    for (std::size_t n = 0; n < neighbors.size(); ++n) result[static_cast<std::size_t>(neighbors[n].i)].push_back(n);
    return result;
}

// Orthonormal associated Legendre factor of Y_lm (Condon-Shortley phase).
double legendre(int l, int m, double x)
{
    double pmm = 1.0;
    if (m > 0) {
        const double omx2 = (1.0 - x) * (1.0 + x);
        double fact = 1.0;
        for (int i = 1; i <= m; ++i) { pmm *= omx2 * fact / (fact + 1.0); fact += 2.0; }
    }
    pmm = std::sqrt((2 * m + 1) * pmm / (4 * kPi));
    if (m & 1) pmm = -pmm;
    if (l == m) return pmm;
    double pmmp1 = x * std::sqrt(2.0 * m + 3.0) * pmm;
    if (l == m + 1) return pmmp1;
    double oldFactor = std::sqrt(2.0 * m + 3.0), pll = 0;
    for (int ll = m + 2; ll <= l; ++ll) {
        const double factor = std::sqrt((4.0 * ll * ll - 1.0) / (static_cast<double>(ll) * ll - static_cast<double>(m) * m));
        pll = (x * pmmp1 - pmm / oldFactor) * factor;
        oldFactor = factor;
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

}

Json chillPlus(const Parameters& p)
{
    const Atoms a = atomsOf(p);
    const double cutoff = positive(p.number("cutoff_A", 3.5), "cutoff_A");
    const auto neighbors = neighborList(a.positions, a.cell, a.pbc, cutoff);
    const auto groups = byCenter(neighbors, a.positions.size());
    const std::size_t n = a.positions.size();
    // q_3m for m = 0..3 (negative m follow by symmetry).
    std::vector<std::array<std::complex<double>, 4>> q(n);
    for (std::size_t i = 0; i < n; ++i) {
        q[i].fill(0.0);
        for (std::size_t g : groups[i]) {
            const Vec3& v = neighbors[g].vector;
            const double r = norm(v);
            if (r < 1e-10) throw std::runtime_error("Coincident atoms make Chill+ undefined");
            const double cosine = std::clamp(v[2] / r, -1.0, 1.0), azimuth = std::atan2(v[1], v[0]);
            for (int m = 0; m <= 3; ++m) q[i][static_cast<std::size_t>(m)] += legendre(3, m, cosine) * std::polar(1.0, m * azimuth);
        }
    }
    const auto correlation = [&](std::size_t i, std::size_t j) {
        double dotProduct = std::real(q[i][0] * std::conj(q[j][0]));
        double ni = std::norm(q[i][0]), nj = std::norm(q[j][0]);
        for (std::size_t m = 1; m <= 3; ++m) {
            dotProduct += 2 * std::real(q[i][m] * std::conj(q[j][m]));
            ni += 2 * std::norm(q[i][m]);
            nj += 2 * std::norm(q[j][m]);
        }
        return ni > 0 && nj > 0 ? dotProduct / std::sqrt(ni * nj) : kNaN;
    };
    static const std::vector<std::string> names = {"other", "hexagonal ice", "cubic ice", "interfacial ice", "hydrate", "interfacial hydrate"};
    std::vector<int> type(n, 0), staggeredCount(n, 0), eclipsedCount(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (i % 1024 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
        if (groups[i].size() != 4) continue;
        int staggered = 0, eclipsed = 0;
        for (std::size_t g : groups[i]) {
            const double c = correlation(i, static_cast<std::size_t>(neighbors[g].j));
            if (c <= -0.8) ++staggered;
            else if (c >= -0.35 && c <= 0.25) ++eclipsed;
        }
        staggeredCount[i] = staggered;
        eclipsedCount[i] = eclipsed;
        if (eclipsed == 4) type[i] = 4;
        else if (eclipsed == 3) type[i] = 5;
        else if (staggered == 4) type[i] = 2;
        else if (staggered == 3 && eclipsed == 1) type[i] = 1;
        else if (staggered == 3 || staggered == 2) type[i] = 3;
    }
    std::vector<long long> counts(names.size(), 0);
    for (int t : type) ++counts[static_cast<std::size_t>(t)];
    Json result = Json::object(), countJson = Json::object(), typeNames = Json::array();
    for (std::size_t t = 0; t < names.size(); ++t) {
        countJson[names[t]] = static_cast<double>(counts[t]);
        typeNames.push(names[t]);
    }
    result["chill_type"] = toJson(type);
    result["type_names"] = typeNames;
    result["counts"] = countJson;
    result["staggered_bonds"] = toJson(staggeredCount);
    result["eclipsed_bonds"] = toJson(eclipsedCount);
    result["cutoff_A"] = cutoff;
    return result;
}

Json localEntropy(const Parameters& p)
{
    const Atoms a = atomsOf(p);
    const double cutoff = positive(p.number("cutoff_A", 5.0), "cutoff_A");
    const double sigma = positive(p.number("sigma_A", 0.2), "sigma_A");
    const bool localDensity = p.boolean("local_density", false);
    const bool average = p.boolean("average", true);
    const double averageCutoff = positive(p.number("average_cutoff_A", cutoff), "average_cutoff_A");
    const std::size_t n = a.positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    double globalDensity = kNaN;
    if (p.has("cell") && a.pbc[0] && a.pbc[1] && a.pbc[2]) globalDensity = static_cast<double>(n) / cellVolume(a.cell);
    if (!localDensity && !std::isfinite(globalDensity))
        throw std::runtime_error("The global density needs a fully periodic cell; set local_density for other systems");
    const auto neighbors = neighborList(a.positions, a.cell, a.pbc, std::max(cutoff + 3 * sigma, averageCutoff));
    const auto groups = byCenter(neighbors, n);
    // Radial grid for the integral, from dr/2 to the cutoff.
    const double dr = sigma / 10;
    const std::size_t bins = static_cast<std::size_t>(std::ceil(cutoff / dr));
    std::vector<double> radius(bins);
    for (std::size_t b = 0; b < bins; ++b) radius[b] = (static_cast<double>(b) + 0.5) * cutoff / static_cast<double>(bins);
    const double width = cutoff / static_cast<double>(bins);
    const double gaussianNorm = 1 / std::sqrt(2 * kPi * sigma * sigma);
    std::vector<double> entropy(n, 0.0);
    std::vector<double> g(bins);
    for (std::size_t i = 0; i < n; ++i) {
        if (i % 256 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
        double density = globalDensity;
        if (localDensity) {
            std::size_t inside = 0;
            for (std::size_t k : groups[i]) inside += norm(neighbors[k].vector) <= cutoff;
            density = static_cast<double>(inside) / (4.0 / 3.0 * kPi * cutoff * cutoff * cutoff);
        }
        if (!(density > 0)) { entropy[i] = 0; continue; }
        std::fill(g.begin(), g.end(), 0.0);
        for (std::size_t k : groups[i]) {
            const double rij = norm(neighbors[k].vector);
            // Only bins within 5 sigma of r_ij contribute measurably.
            const std::size_t low = static_cast<std::size_t>(std::max(0.0, (rij - 5 * sigma) / width));
            const std::size_t high = std::min(bins, static_cast<std::size_t>(std::max(0.0, (rij + 5 * sigma) / width)) + 1);
            for (std::size_t b = low; b < high; ++b) {
                const double d = radius[b] - rij;
                g[b] += gaussianNorm * std::exp(-d * d / (2 * sigma * sigma));
            }
        }
        double integral = 0;
        for (std::size_t b = 0; b < bins; ++b) {
            const double r = radius[b];
            const double gr = g[b] / (4 * kPi * density * r * r);
            const double integrand = (gr > 0 ? gr * std::log(gr) : 0.0) - gr + 1;
            integral += integrand * r * r * width;
        }
        entropy[i] = -2 * kPi * density * integral;
    }
    std::vector<double> averaged = entropy;
    if (average) {
        for (std::size_t i = 0; i < n; ++i) {
            double sum = entropy[i];
            double count = 1;
            for (std::size_t k : groups[i])
                if (norm(neighbors[k].vector) <= averageCutoff) { sum += entropy[static_cast<std::size_t>(neighbors[k].j)]; count += 1; }
            averaged[i] = sum / count;
        }
    }
    double mean = 0;
    for (double s : averaged) mean += s;
    Json result = Json::object();
    result["entropy_kB"] = toJson(averaged);
    if (average) result["unaveraged_entropy_kB"] = toJson(entropy);
    result["mean_entropy_kB"] = mean / static_cast<double>(n);
    result["cutoff_A"] = cutoff;
    result["sigma_A"] = sigma;
    return result;
}

Json findRings(const Parameters& p)
{
    const Atoms a = atomsOf(p);
    const double cutoff = positive(p.number("cutoff_A", 2.0), "cutoff_A");
    const long long maxSize = integer(p.number("max_size", 10), "max_size", 3);
    if (maxSize > 24) throw std::runtime_error("max_size must be at most 24");
    const std::size_t n = a.positions.size();
    struct Bond { int j; std::array<int, 3> shift; };
    std::vector<std::vector<Bond>> bonds(n);
    for (const auto& link : neighborList(a.positions, a.cell, a.pbc, cutoff))
        bonds[static_cast<std::size_t>(link.i)].push_back({link.j, link.shift});
    // A node of the unwrapped network: atom and periodic image.
    using Node = std::array<int, 4>;
    struct NodeHash
    {
        std::size_t operator()(const Node& k) const
        {
            std::size_t h = static_cast<std::size_t>(k[0]) * 1000003u;
            for (int d = 1; d < 4; ++d) h = (h ^ static_cast<std::size_t>(k[static_cast<std::size_t>(d)] + 512)) * 1000003u;
            return h;
        }
    };
    const auto step = [](const Node& from, const Bond& b) {
        return Node{b.j, from[1] + b.shift[0], from[2] + b.shift[1], from[3] + b.shift[2]};
    };
    // Canonical key: rotation and direction that give the smallest sequence,
    // with images taken relative to the first node.
    const auto canonical = [](const std::vector<Node>& ring) {
        const std::size_t m = ring.size();
        std::vector<Node> best;
        for (std::size_t start = 0; start < m; ++start)
            for (int direction : {1, -1}) {
                std::vector<Node> sequence(m);
                const Node& origin = ring[start];
                for (std::size_t k = 0; k < m; ++k) {
                    const std::size_t index = direction > 0 ? (start + k) % m : (start + m - k % m) % m;
                    const Node& node = ring[index];
                    sequence[k] = {node[0], node[1] - origin[1], node[2] - origin[2], node[3] - origin[3]};
                }
                if (best.empty() || sequence < best) best = std::move(sequence);
            }
        return best;
    };
    std::set<std::vector<Node>> rings;
    const std::size_t pathLimit = 64;
    for (std::size_t i = 0; i < n; ++i) {
        if (i % 128 == 0) taskProgress(static_cast<double>(i) / static_cast<double>(n));
        for (const Bond& bond : bonds[i]) {
            const std::size_t j = static_cast<std::size_t>(bond.j);
            // Each undirected bond once.
            if (j < i || (j == i && bond.shift < std::array<int, 3>{0, 0, 0})) continue;
            const Node source{static_cast<int>(i), 0, 0, 0};
            const Node target{bond.j, bond.shift[0], bond.shift[1], bond.shift[2]};
            // Breadth-first search from i to j without the bond itself.
            std::unordered_map<Node, int, NodeHash> distance;
            std::unordered_map<Node, std::vector<Node>, NodeHash> parents;
            distance[source] = 0;
            std::vector<Node> level = {source};
            bool found = false;
            for (int depth = 0; depth + 1 < maxSize && !level.empty() && !found; ++depth) {
                std::vector<Node> next;
                for (const Node& node : level)
                    for (const Bond& b : bonds[static_cast<std::size_t>(node[0])]) {
                        const Node to = step(node, b);
                        if ((node == source && to == target) || (node == target && to == source)) continue;
                        const auto known = distance.find(to);
                        if (known == distance.end()) {
                            distance[to] = depth + 1;
                            parents[to].push_back(node);
                            next.push_back(to);
                        } else if (known->second == depth + 1) parents[to].push_back(node);
                        if (to == target) found = true;
                    }
                level = std::move(next);
            }
            if (!found) continue;
            // Every shortest path back from the target (bounded in number).
            std::vector<std::vector<Node>> paths = {{target}};
            std::vector<std::vector<Node>> complete;
            while (!paths.empty() && complete.size() < pathLimit) {
                std::vector<Node> path = std::move(paths.back());
                paths.pop_back();
                const Node& last = path.back();
                if (last == source) { complete.push_back(std::move(path)); continue; }
                for (const Node& parent : parents[last]) {
                    if (paths.size() + complete.size() > 4 * pathLimit) break;
                    std::vector<Node> extended = path;
                    extended.push_back(parent);
                    paths.push_back(std::move(extended));
                }
            }
            for (const auto& path : complete) rings.insert(canonical(path));
        }
    }
    std::map<std::size_t, long long> sizes;
    std::vector<int> perAtom(n, 0);
    Json list = Json::array();
    const std::size_t listLimit = 10000;
    for (const auto& ring : rings) {
        ++sizes[ring.size()];
        Json members = Json::array();
        for (const Node& node : ring) {
            ++perAtom[static_cast<std::size_t>(node[0])];
            members.push(static_cast<double>(node[0]));
        }
        if (list.size() < listLimit) list.push(members);
    }
    Json sizeJson = Json::array(), countJson = Json::array();
    for (const auto& [size, count] : sizes) {
        sizeJson.push(static_cast<double>(size));
        countJson.push(static_cast<double>(count));
    }
    Json result = Json::object();
    result["ring_sizes"] = sizeJson;
    result["ring_counts"] = countJson;
    result["total_rings"] = static_cast<double>(rings.size());
    result["rings_per_atom"] = toJson(perAtom);
    result["rings"] = list;
    result["rings_listed"] = static_cast<double>(list.size());
    result["cutoff_A"] = cutoff;
    result["max_size"] = static_cast<double>(maxSize);
    return result;
}

Json displacementVectors(const Parameters& p)
{
    const auto reference = points(p.array("reference"), "reference");
    const auto current = points(p.array("current"), "current");
    if (reference.size() != current.size()) throw std::runtime_error("reference and current need the same number of atoms");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    if (periodic && !p.has("reference_cell")) throw std::runtime_error("Periodic displacements need reference_cell");
    const Mat3 referenceCell = p.has("reference_cell") ? cellMatrix(p.array("reference_cell")) : identity();
    const Mat3 currentCell = p.has("current_cell") ? cellMatrix(p.array("current_cell")) : referenceCell;
    const std::string mapping = p.has("affine_mapping") ? p.json("affine_mapping").string() : std::string("off");
    if (mapping != "off" && mapping != "to_reference" && mapping != "to_current")
        throw std::runtime_error("affine_mapping must be off, to_reference or to_current");
    if (mapping != "off" && !(p.has("reference_cell") && p.has("current_cell")))
        throw std::runtime_error("Affine mapping needs reference_cell and current_cell");
    const Mat3& imageCell = mapping == "to_current" ? currentCell : referenceCell;
    const MinimumImage mic(imageCell, pbc);
    const std::size_t n = reference.size();
    Json vectors = Json::array();
    std::vector<double> magnitude(n);
    double sum = 0, sumSquares = 0, largest = 0;
    for (std::size_t i = 0; i < n; ++i) {
        Vec3 from = reference[i], to = current[i];
        if (mapping == "to_reference") to = rowTimes(fractional(to, currentCell), referenceCell);
        else if (mapping == "to_current") from = rowTimes(fractional(from, referenceCell), currentCell);
        Vec3 d = sub(to, from);
        if (periodic) d = mic(d);
        vectors.push(toJson(d));
        magnitude[i] = norm(d);
        sum += magnitude[i];
        sumSquares += magnitude[i] * magnitude[i];
        largest = std::max(largest, magnitude[i]);
    }
    Json result = Json::object();
    result["displacement_A"] = vectors;
    result["magnitude_A"] = toJson(magnitude);
    result["mean_magnitude_A"] = n ? sum / static_cast<double>(n) : kNaN;
    result["rms_displacement_A"] = n ? std::sqrt(sumSquares / static_cast<double>(n)) : kNaN;
    result["max_magnitude_A"] = largest;
    return result;
}
}
