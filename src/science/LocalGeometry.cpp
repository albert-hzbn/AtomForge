#include "science/LocalGeometry.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atomforge::science
{
AtomSet readAtomSet(const Parameters& p, const char* positions)
{
    AtomSet a;
    a.positions = points(p.array(positions), positions);
    a.pbc = p.pbc("pbc", {false, false, false});
    a.hasCell = p.has("cell");
    if (a.periodic() && !a.hasCell) throw std::runtime_error("Periodic analysis requires a cell");
    a.cell = a.hasCell ? cellMatrix(p.array("cell")) : identity();
    return a;
}

double meanAtomicSpacing(const AtomSet& atoms)
{
    if (atoms.positions.empty()) throw std::runtime_error("No atoms");
    double volume = 1;
    if (atoms.hasCell) volume = cellVolume(atoms.cell);
    else {
        Vec3 low = atoms.positions[0], high = atoms.positions[0];
        for (const auto& x : atoms.positions)
            for (std::size_t d = 0; d < 3; ++d) {
                low[d] = std::min(low[d], x[d]);
                high[d] = std::max(high[d], x[d]);
            }
        for (std::size_t d = 0; d < 3; ++d) volume *= high[d] - low[d] + 2.0;
    }
    return std::cbrt(volume / static_cast<double>(atoms.positions.size()));
}

std::vector<std::vector<std::size_t>> neighborsByCenter(const std::vector<Neighbor>& neighbors, std::size_t atoms)
{
    std::vector<std::vector<std::size_t>> result(atoms);
    for (std::size_t n = 0; n < neighbors.size(); ++n) result[static_cast<std::size_t>(neighbors[n].i)].push_back(n);
    return result;
}

std::vector<std::vector<Vec3>> sortedNeighborVectors(const AtomSet& atoms, double cutoff)
{
    std::vector<std::vector<Vec3>> lists(atoms.positions.size());
    for (const auto& link : neighborList(atoms.positions, atoms.cell, atoms.pbc, cutoff)) {
        if (norm(link.vector) < 1e-10) throw std::runtime_error("Coincident atoms make local environment analysis undefined");
        lists[static_cast<std::size_t>(link.i)].push_back(link.vector);
    }
    for (auto& list : lists) std::sort(list.begin(), list.end(), [](const Vec3& x, const Vec3& y) { return dot(x, x) < dot(y, y); });
    return lists;
}

double normalizedLegendre(int l, int m, double x)
{
    const double pi = std::acos(-1.0);
    double pmm = 1.0;
    if (m > 0) {
        const double omx2 = (1.0 - x) * (1.0 + x);
        double fact = 1.0;
        for (int i = 1; i <= m; ++i) { pmm *= omx2 * fact / (fact + 1.0); fact += 2.0; }
    }
    pmm = std::sqrt((2 * m + 1) * pmm / (4 * pi));
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

std::size_t unionRoot(std::vector<std::size_t>& parent, std::size_t i)
{
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}

void addClassification(Json& result, const std::string& key, const std::vector<int>& types,
                       const std::vector<std::string>& names, const std::string& prefix)
{
    std::vector<long long> counts(names.size(), 0);
    for (int t : types) ++counts.at(static_cast<std::size_t>(t));
    Json nameList = Json::array(), countJson = Json::object(), fractions = Json::object();
    for (std::size_t t = 0; t < names.size(); ++t) {
        nameList.push(names[t]);
        countJson[names[t]] = static_cast<double>(counts[t]);
        fractions[names[t]] = types.empty() ? 0.0 : static_cast<double>(counts[t]) / static_cast<double>(types.size());
    }
    result[key] = toJson(types);
    result[prefix.empty() ? "type_names" : prefix + "names"] = nameList;
    result[prefix + "counts"] = countJson;
    result[prefix + "fractions"] = fractions;
}
}
