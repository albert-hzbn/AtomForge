#include "science/BondsAndRegions.h"
#include "electronic/Volume.h"
#include "util/ElementData.h"
#include "util/TaskControl.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kHuge = 1e30;

std::vector<double> speciesOf(const Parameters& p, std::size_t n)
{
    const auto values = finiteArray(p.array("types"), "types", 1).values;
    if (values.size() != n) throw std::runtime_error("types needs one atomic number per atom");
    return values;
}

// ---------------------------------------------------------------- regions
struct Region
{
    std::string id, style;
    std::vector<double> args;
    int axis = 0;          // cylinder and cone: 0 x, 1 y, 2 z
    bool outside = false;  // side out
    Vec3 origin{}, a{}, b{}, c{};  // prism

    // Signed distance, positive inside (approximate for cones and prisms).
    double distance(const Vec3& x) const
    {
        double d = 0;
        if (style == "block") {
            d = kHuge;
            for (std::size_t k = 0; k < 3; ++k) d = std::min({d, x[k] - args[2 * k], args[2 * k + 1] - x[k]});
        } else if (style == "sphere") {
            d = args[3] - norm(sub(x, {args[0], args[1], args[2]}));
        } else if (style == "cylinder" || style == "cone") {
            const std::size_t ax = static_cast<std::size_t>(axis), u = (ax + 1) % 3, v = (ax + 2) % 3;
            // LAMMPS lists the two other coordinates in x, y, z order.
            const std::size_t first = std::min(u, v), second = std::max(u, v);
            const double radial = std::hypot(x[first] - args[0], x[second] - args[1]);
            double radius = args[2], low = args[3], high = args[4];
            if (style == "cone") {
                low = args[4];
                high = args[5];
                const double t = high > low ? std::clamp((x[ax] - low) / (high - low), 0.0, 1.0) : 0.0;
                radius = args[2] + t * (args[3] - args[2]);
            }
            d = std::min({radius - radial, x[ax] - low, high - x[ax]});
        } else if (style == "plane") {
            const Vec3 normal = {args[3], args[4], args[5]};
            d = dot(sub(x, {args[0], args[1], args[2]}), scale(normal, 1 / norm(normal)));
        } else if (style == "prism") {
            const Mat3 cell = {a, b, c};
            const Vec3 f = fractional(sub(x, origin), cell);
            d = kHuge;
            for (std::size_t k = 0; k < 3; ++k) {
                // Fractional distance to the faces, scaled by the face spacing.
                const double height = cellVolume(cell) / norm(cross(cell[(k + 1) % 3], cell[(k + 2) % 3]));
                d = std::min({d, f[k] * height, (1 - f[k]) * height});
            }
        }
        return outside ? -d : d;
    }
};

// A bound of a region: a number, INF or -INF (infinite on the bound's own
// side, as in LAMMPS), or EDGE (the domain boundary). `fixed` marks bounds
// that lattice units must not scale.
double parseBound(const std::string& word, bool lower, double edge, bool& fixed)
{
    fixed = true;
    if (word == "INF" || word == "-INF") return lower ? -kHuge : kHuge;
    if (word == "EDGE") return edge;
    fixed = false;
    std::size_t used = 0;
    const double value = std::stod(word, &used);
    if (used != word.size()) throw std::runtime_error("Not a number in a region command: " + word);
    return value;
}

std::vector<Region> parseRegions(const std::string& text, double latticeScale, const Vec3& low, const Vec3& high)
{
    std::vector<Region> regions;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        if (const auto hash = line.find('#'); hash != std::string::npos) line.erase(hash);
        std::istringstream in(line);
        std::vector<std::string> words;
        for (std::string w; in >> w;) words.push_back(w);
        if (words.empty()) continue;
        if (words[0] != "region") throw std::runtime_error("Each line must be a LAMMPS region command: " + line);
        if (words.size() < 4) throw std::runtime_error("Incomplete region command: " + line);
        Region r;
        r.id = words[1];
        r.style = words[2];
        static const std::map<std::string, std::size_t> counts = {{"block", 6}, {"sphere", 4}, {"cylinder", 6}, {"cone", 7}, {"plane", 6}, {"prism", 9}};
        const auto found = counts.find(r.style);
        if (found == counts.end()) throw std::runtime_error("Unsupported region style " + r.style + " (block, sphere, cylinder, cone, plane, prism)");
        std::size_t k = 3;
        if (r.style == "cylinder" || r.style == "cone") {
            if (words[3] != "x" && words[3] != "y" && words[3] != "z") throw std::runtime_error("The axis of a " + r.style + " is x, y or z");
            r.axis = words[3][0] - 'x';
            ++k;
        }
        const std::size_t needed = found->second - ((r.style == "cylinder" || r.style == "cone") ? 1 : 0);
        if (words.size() < k + needed) throw std::runtime_error("Too few values for region " + r.id);
        // Each argument's coordinate and whether it is a lower bound.
        std::vector<char> fixed(needed, 0);
        for (std::size_t m = 0; m < needed; ++m) {
            std::size_t coordinate = 0;
            bool lower = true;
            if ((r.style == "block" || r.style == "prism") && m < 6) { coordinate = m / 2; lower = m % 2 == 0; }
            else if (r.style == "cylinder" && m >= 3) { coordinate = static_cast<std::size_t>(r.axis); lower = m == 3; }
            else if (r.style == "cone" && m >= 4) { coordinate = static_cast<std::size_t>(r.axis); lower = m == 4; }
            bool isFixed = false;
            r.args.push_back(parseBound(words[k + m], lower, lower ? low[coordinate] : high[coordinate], isFixed));
            fixed[m] = isFixed;
        }
        bool lattice = true;  // LAMMPS's default
        for (std::size_t m = k + needed; m < words.size(); ++m) {
            if (words[m] == "side" && m + 1 < words.size()) { r.outside = words[++m] == "out"; }
            else if (words[m] == "units" && m + 1 < words.size()) { lattice = words[++m] == "lattice"; }
            else if (words[m] == "open" && m + 1 < words.size()) { ++m; }
            else throw std::runtime_error("Unsupported region keyword " + words[m] + " (side, units and open are understood)");
        }
        if (lattice) {
            for (std::size_t m = 0; m < r.args.size(); ++m)
                if (!fixed[m] && !(r.style == "plane" && m >= 3)) r.args[m] *= latticeScale;
        }
        if (r.style == "prism") {
            r.origin = {r.args[0], r.args[2], r.args[4]};
            r.a = {r.args[1] - r.args[0], 0, 0};
            r.b = {r.args[6], r.args[3] - r.args[2], 0};
            r.c = {r.args[7], r.args[8], r.args[5] - r.args[4]};
        }
        regions.push_back(r);
    }
    if (regions.empty()) throw std::runtime_error("No region commands");
    return regions;
}
}

Json bondOrders(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    const std::size_t n = positions.size();
    const auto species = speciesOf(p, n);
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    if ((pbc[0] || pbc[1] || pbc[2]) && !p.has("cell")) throw std::runtime_error("Periodic analysis requires a cell");
    const Mat3 cell = p.has("cell") ? cellMatrix(p.array("cell")) : identity();
    const double tolerance = positive(p.number("tolerance", 1.25), "tolerance");
    static const std::vector<float> covalent = makeLiteratureCovalentRadii();
    const auto radius = [](double z) {
        const std::size_t k = static_cast<std::size_t>(std::llround(z));
        return k < covalent.size() && covalent[k] > 0 ? static_cast<double>(covalent[k]) : 1.5;
    };
    double largest = 0;
    for (double z : species) largest = std::max(largest, radius(z));
    std::vector<double> valence(n, 0.0);
    std::vector<int> coordination(n, 0);
    Json bonds = Json::array();
    std::map<double, long long> rounded;
    for (const auto& link : neighborList(positions, cell, pbc, 2 * largest * tolerance)) {
        if (link.j < link.i || (link.j == link.i && link.shift < std::array<int, 3>{0, 0, 0})) continue;
        const double r = norm(link.vector);
        const double single = radius(species[static_cast<std::size_t>(link.i)]) + radius(species[static_cast<std::size_t>(link.j)]);
        if (r > tolerance * single || r <= 0) continue;
        const double order = std::pow(10.0, (single - r) / 0.6);
        // The nearest conventional order.
        double nearest = 1;
        for (double candidate : {1.0, 1.5, 2.0, 3.0})
            if (std::abs(order - candidate) < std::abs(order - nearest)) nearest = candidate;
        ++rounded[nearest];
        for (int atom : {link.i, link.j}) { valence[static_cast<std::size_t>(atom)] += order; ++coordination[static_cast<std::size_t>(atom)]; }
        Json bond = Json::object();
        bond["atoms"] = Json::array({static_cast<double>(link.i), static_cast<double>(link.j)});
        bond["length_A"] = r;
        bond["order"] = order;
        bond["rounded_order"] = nearest;
        bonds.push(bond);
    }
    Json counts = Json::object();
    for (const auto& [order, count] : rounded) {
        char key[16];
        std::snprintf(key, sizeof key, "%g", order);
        counts[key] = static_cast<double>(count);
    }
    Json result = Json::object();
    result["bonds"] = bonds;
    result["bond_count"] = static_cast<double>(bonds.size());
    result["rounded_order_counts"] = counts;
    result["valence_sum"] = toJson(valence);
    result["coordination"] = toJson(coordination);
    return result;
}

Json lammpsRegions(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    if (!p.has("regions")) throw std::runtime_error("Give the region commands in regions");
    std::string text = p.json("regions").isString() ? p.json("regions").string() : std::string();
    std::replace(text.begin(), text.end(), ';', '\n');
    // The domain: the cell, else the atoms' bounding box.
    Vec3 low{kHuge, kHuge, kHuge}, high{-kHuge, -kHuge, -kHuge};
    if (p.has("cell")) {
        const Mat3 cell = cellMatrix(p.array("cell"));
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const Vec3 corner = rowTimes({static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)}, cell);
                    for (std::size_t d = 0; d < 3; ++d) { low[d] = std::min(low[d], corner[d]); high[d] = std::max(high[d], corner[d]); }
                }
    } else {
        for (const auto& x : positions)
            for (std::size_t d = 0; d < 3; ++d) { low[d] = std::min(low[d], x[d] - 1); high[d] = std::max(high[d], x[d] + 1); }
        if (positions.empty()) throw std::runtime_error("No atoms and no cell to define the domain");
    }
    const auto regions = parseRegions(text, positive(p.number("lattice_scale", 1.0), "lattice_scale"), low, high);
    Json perRegion = Json::object(), counts = Json::object();
    std::vector<bool> any(positions.size(), false);
    for (const auto& r : regions) {
        std::vector<bool> inside(positions.size());
        long long count = 0;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            inside[i] = r.distance(positions[i]) >= 0;
            count += inside[i];
            if (inside[i]) any[i] = true;
        }
        perRegion[r.id] = toJson(inside);
        counts[r.id] = static_cast<double>(count);
    }
    Json result = Json::object();
    result["inside"] = perRegion;
    result["atoms_inside"] = counts;
    result["inside_any"] = toJson(any);
    // Each region's surface within the domain, from its signed distance.
    if (p.boolean("export_mesh", true)) {
        const long long resolution = integer(p.number("resolution", 64), "resolution", 8);
        double longest = 0;
        for (std::size_t d = 0; d < 3; ++d) longest = std::max(longest, high[d] - low[d]);
        const double step = longest / static_cast<double>(resolution);
        std::ostringstream obj;
        obj.precision(8);
        obj << "# AtomForge LAMMPS regions\n";
        std::size_t vertices = 0;
        for (const auto& r : regions) {
            taskCheckpoint();
            atomforge::electronic::Grid field;
            field.periodic = false;
            field.origin = glm::dvec3(low[0] - step, low[1] - step, low[2] - step);
            for (int d = 0; d < 3; ++d) {
                field.shape[static_cast<std::size_t>(d)] = static_cast<int>(std::ceil((high[static_cast<std::size_t>(d)] - low[static_cast<std::size_t>(d)]) / step)) + 3;
                field.cell[d] = glm::dvec3(0.0);
                field.cell[d][d] = step * (field.shape[static_cast<std::size_t>(d)] - 1);
            }
            field.values.resize(static_cast<std::size_t>(field.shape[0]) * field.shape[1] * field.shape[2]);
            for (int z = 0; z < field.shape[2]; ++z)
                for (int y = 0; y < field.shape[1]; ++y)
                    for (int x = 0; x < field.shape[0]; ++x) {
                        const glm::dvec3 q = field.position(x, y, z);
                        // Outside the domain counts as outside every region, closing the surface there.
                        double d = r.distance({q.x, q.y, q.z});
                        for (std::size_t k = 0; k < 3; ++k) d = std::min({d, q[static_cast<glm::length_t>(k)] - low[k], high[k] - q[static_cast<glm::length_t>(k)]});
                        field.values[(static_cast<std::size_t>(z) * field.shape[1] + y) * field.shape[0] + x] = d;
                    }
            const auto mesh = atomforge::electronic::isosurface(field, 0.0);
            obj << "o " << r.id << '\n';
            for (const auto& v : mesh.vertices) obj << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';
            for (std::size_t t = 0; t + 2 < mesh.vertices.size(); t += 3)
                obj << "f " << vertices + t + 1 << ' ' << vertices + t + 2 << ' ' << vertices + t + 3 << '\n';
            vertices += mesh.vertices.size();
        }
        Json files = Json::object();
        files["regions.obj"] = obj.str();
        result["files"] = files;
    }
    return result;
}
}
