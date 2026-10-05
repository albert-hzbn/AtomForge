// Physical-invariant checks for the vacancy, strain, surface, SQS, nanowire
// and core-shell builders (src/algorithms). Neighbour counting uses the
// independent native neighbour list from src/science.
#include "algorithms/NanostructureTools.h"
#include "algorithms/SQSBuilder.h"
#include "algorithms/StrainTool.h"
#include "algorithms/SurfaceBuilder.h"
#include "algorithms/VacancyBuilder.h"
#include "science/ScienceCore.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace atomforge;
using science::Vec3;
using science::Mat3;

namespace
{
int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) throw std::runtime_error(what);
}

void close(double actual, double expected, double tolerance, const std::string& what)
{
    if (!(std::abs(actual - expected) <= tolerance))
        throw std::runtime_error(what + ": expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}

void test(const std::string& name, const std::function<void()>& body)
{
    try {
        body();
        std::cout << "  ok   " << name << '\n';
    } catch (const std::exception& error) {
        ++failures;
        std::cout << "  FAIL " << name << ": " << error.what() << '\n';
    }
}

AtomSite atom(const std::string& symbol, double x, double y, double z)
{
    AtomSite site;
    site.symbol = symbol;
    for (int zNumber = 1; zNumber <= 118; ++zNumber)
        if (symbol == elementSymbol(zNumber)) site.atomicNumber = zNumber;
    site.x = x; site.y = y; site.z = z;
    return site;
}

Structure fcc(double a, int n, const std::string& symbol = "Cu")
{
    Structure s;
    s.hasUnitCell = true;
    s.cellVectors = {{{a * n, 0, 0}, {0, a * n, 0}, {0, 0, a * n}}};
    const double basis[4][3] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (const auto& b : basis) s.atoms.push_back(atom(symbol, (i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a));
    return s;
}

Mat3 cellOf(const Structure& s)
{
    Mat3 cell{};
    for (int r = 0; r < 3; ++r) cell[r] = s.cellVectors[r];
    return cell;
}

std::vector<Vec3> positionsOf(const Structure& s)
{
    std::vector<Vec3> result;
    for (const auto& a : s.atoms) result.push_back({a.x, a.y, a.z});
    return result;
}

std::vector<int> coordination(const Structure& s, double cutoff, science::Pbc pbc = {true, true, true})
{
    std::vector<int> count(s.atoms.size(), 0);
    for (const auto& n : science::neighborList(positionsOf(s), cellOf(s), pbc, cutoff)) ++count[static_cast<std::size_t>(n.i)];
    return count;
}
}

int main()
{
    std::cout << "Builder tool regressions\n";

    test("vacancies: counts, element filter, separation and reproducibility", [] {
        const Structure host = fcc(3.6, 3);  // 108 atoms
        VacancyParams p;
        p.targetPercentage = 5;
        auto r = buildVacancies(host, p);
        check(r.success, r.message);
        check(r.eligibleCount == 108 && r.removedCount == r.requestedCount, "requested vacancies removed");
        close(r.removedCount, std::round(108 * 0.05), 1, "5% of 108 sites");
        check(r.structure.atoms.size() == 108 - static_cast<std::size_t>(r.removedCount), "atom count");
        check(r.structure.hasUnitCell, "cell retained");
        // Survivors are untouched lattice sites.
        for (const auto& a : r.structure.atoms) {
            bool onSite = false;
            for (const auto& b : host.atoms) onSite = onSite || (std::abs(a.x - b.x) + std::abs(a.y - b.y) + std::abs(a.z - b.z) < 1e-12);
            check(onSite, "remaining atoms stay on lattice sites");
        }
        const auto again = buildVacancies(host, p);
        check(again.structure.atoms.size() == r.structure.atoms.size(), "seeded reproducibility");
        bool identical = true;
        for (std::size_t i = 0; i < again.structure.atoms.size(); ++i) identical = identical && again.structure.atoms[i].x == r.structure.atoms[i].x;
        check(identical, "same seed, same vacancies");

        p = {};
        p.targetCount = 6;
        p.minSeparation = 5.0;
        r = buildVacancies(host, p);
        check(r.success && r.removedCount == 6, "six separated vacancies");
        std::vector<Vec3> removed;
        for (const auto& b : host.atoms) {
            bool present = false;
            for (const auto& a : r.structure.atoms) present = present || (a.x == b.x && a.y == b.y && a.z == b.z);
            if (!present) removed.push_back({b.x, b.y, b.z});
        }
        check(removed.size() == 6, "six sites removed");
        const science::MinimumImage mic(cellOf(host), {true, true, true});
        for (std::size_t i = 0; i < removed.size(); ++i)
            for (std::size_t j = i + 1; j < removed.size(); ++j)
                check(science::norm(mic(science::sub(removed[i], removed[j]))) >= 5.0 - 1e-9, "minimum vacancy separation (periodic)");

        Structure alloy = host;
        for (std::size_t i = 0; i < alloy.atoms.size(); i += 2) alloy.atoms[i] = atom("Ni", alloy.atoms[i].x, alloy.atoms[i].y, alloy.atoms[i].z);
        p = {};
        p.element = "Ni";
        p.targetCount = 10;
        r = buildVacancies(alloy, p);
        int ni = 0;
        for (const auto& a : r.structure.atoms) ni += a.symbol == "Ni";
        check(r.success && r.eligibleCount == 54 && ni == 44, "only the selected element is removed");
        p.targetCount = 1000;
        check(!buildVacancies(alloy, p).success, "cannot remove more atoms than are eligible");
    });

    test("strain: deformation gradient, volume and fractional coordinates", [] {
        const Structure host = fcc(3.6, 2);
        const auto params = engineeringStrain(.02, -.01, .03, .01, 0, .005);
        check(params.f[0][1] == params.f[1][0] && params.f[1][2] == params.f[2][1], "symmetric shear");
        const auto r = applyStrain(host, params);
        check(r.success, r.message);
        const Mat3 before = cellOf(host), after = cellOf(r.structure);
        Mat3 f{};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) f[i][j] = params.f[i][j];
        close(science::cellVolume(after), science::cellVolume(before) * science::determinant(f), 1e-9, "volume scales by det F");
        for (std::size_t i = 0; i < host.atoms.size(); ++i) {
            const Vec3 fa = science::fractional({host.atoms[i].x, host.atoms[i].y, host.atoms[i].z}, before);
            const Vec3 fb = science::fractional({r.structure.atoms[i].x, r.structure.atoms[i].y, r.structure.atoms[i].z}, after);
            for (int k = 0; k < 3; ++k) close(fb[k], fa[k], 1e-12, "fractional coordinates preserved");
        }
        Structure cluster = host;
        cluster.hasUnitCell = false;
        check(!applyStrain(cluster, params).success, "strain needs a periodic cell");
    });

#ifdef ATOMS_ENABLE_SPGLIB
    test("surfaces: layer count, orientation and surface coordination", [] {
        const Structure copper = fcc(3.6, 1);
        const double a = 3.6;
        struct Case { std::array<int, 3> hkl; int surfaceCn; double spacing; };
        for (const Case& c : {Case{{1, 1, 1}, 9, a / std::sqrt(3.0)}, Case{{1, 0, 0}, 8, a / 2}, Case{{1, 1, 0}, 7, a / (2 * std::sqrt(2.0))}}) {
            const auto& hkl = c.hkl;
            const int surfaceCn = c.surfaceCn;
            SurfaceParams p;
            p.h = hkl[0]; p.k = hkl[1]; p.l = hkl[2];
            p.layers = 6;
            p.vacuum = 12;
            const auto r = buildSurface(copper, p);
            check(r.success, r.message);
            check(r.atomsPerLayer > 0 && static_cast<int>(r.structure.atoms.size()) == 6 * r.atomsPerLayer, "atoms = layers x atoms per layer");
            const Mat3 cell = cellOf(r.structure);
            // Slabs are oriented with the surface normal along the third cell vector (z).
            close(science::norm(science::cross(cell[0], cell[1])), std::abs(science::cross(cell[0], cell[1])[2]), 1e-9, "normal along z");
            std::vector<double> heights;
            for (const auto& atomSite : r.structure.atoms) heights.push_back(atomSite.z);
            std::sort(heights.begin(), heights.end());
            std::vector<double> planes = {heights.front()};
            for (double z : heights) if (z - planes.back() > 1e-4) planes.push_back(z);
            check(planes.size() == 6, "six atomic planes");
            for (std::size_t i = 1; i < planes.size(); ++i) close(planes[i] - planes[i - 1], c.spacing, 1e-6, "fcc d(hkl) layer spacing");
            close(cell[2][2] - (planes.back() - planes.front()), 12 + c.spacing, 1e-6, "vacuum plus one interlayer gap");
            const auto cn = coordination(r.structure, 2.8);
            const int low = *std::min_element(cn.begin(), cn.end()), high = *std::max_element(cn.begin(), cn.end());
            check(high == 12, "bulk-like interior (CN range " + std::to_string(low) + "-" + std::to_string(high) + ", " + std::to_string(r.structure.atoms.size()) + " atoms)");
            check(low == surfaceCn, "surface coordination for (" + std::to_string(hkl[0]) + std::to_string(hkl[1]) + std::to_string(hkl[2]) +
                  ") is " + std::to_string(surfaceCn) + ", got " + std::to_string(low));
            int surface = 0;
            for (int c : cn) surface += c == surfaceCn;
            check(surface >= 2 * r.atomsPerLayer, "two exposed faces");
        }
    });
#endif

    test("SQS: exact composition and vanishing short-range order", [] {
        const Structure host = fcc(3.6, 4);  // 256 sites
        SQSParams p;
        p.composition = {{"Cu", 0.5}, {"Ni", 0.5}};
        p.steps = 20000;
        p.shells = 2;
        const auto r = buildSQS(host, p);
        check(r.success, r.message);
        int cu = 0, ni = 0;
        for (const auto& a : r.structure.atoms) { cu += a.symbol == "Cu"; ni += a.symbol == "Ni"; }
        check(cu == 128 && ni == 128, "exact 50/50 composition");
        check(r.finalObjective <= r.initialObjective, "annealing does not worsen the objective");
        // Independent first-shell Warren-Cowley parameter: alpha = 1 - P(Ni|Cu)/c_Ni.
        double unlike = 0, pairs = 0;
        for (const auto& n : science::neighborList(positionsOf(r.structure), cellOf(r.structure), {true, true, true}, 2.8)) {
            if (r.structure.atoms[static_cast<std::size_t>(n.i)].symbol != "Cu") continue;
            pairs += 1;
            unlike += r.structure.atoms[static_cast<std::size_t>(n.j)].symbol == "Ni";
        }
        const double alpha = 1 - (unlike / pairs) / 0.5;
        check(std::abs(alpha) < 0.05, "first-shell Warren-Cowley |alpha| < 0.05 (got " + std::to_string(alpha) + ")");
        SQSParams bad;
        bad.composition = {{"Cu", 0.7}, {"Ni", 0.7}};
        check(!buildSQS(host, bad).success, "composition must sum to one");
    });

    test("nanowire: periodic axis and cross-section", [] {
        const Structure copper = fcc(3.6, 1);
        for (int sides : {0, 6}) {
            NanowireParams p;
            p.axis = 2; p.radius = 8; p.sides = sides; p.vacuum = 10; p.axisRepeats = 3;
            const auto r = buildNanowire(copper, p);
            check(r.success, r.message);
            check(r.atomCount > 50 && static_cast<int>(r.structure.atoms.size()) == r.atomCount, "atoms in the wire");
            const Mat3 cell = cellOf(r.structure);
            close(cell[2][2], 3 * 3.6, 1e-9, "axis length = repeats x lattice");
            double cx = 0, cy = 0;
            for (const auto& a : r.structure.atoms) { cx += a.x; cy += a.y; }
            cx /= r.atomCount; cy /= r.atomCount;
            const double limit = sides >= 3 ? 8 / std::cos(std::acos(-1.0) / sides) : 8;
            double extent = 0;
            for (const auto& a : r.structure.atoms) extent = std::max(extent, std::hypot(a.x - cx, a.y - cy));
            check(extent <= limit + 1e-6 && extent > 0.7 * limit, "atoms fill the cross-section without exceeding it");
            // Cross-section box: 2 (radius + vacuum) along both in-plane cell vectors.
            close(science::norm(cell[0]), 2 * (8 + 10), 1e-9, "in-plane box a");
            close(science::norm(cell[1]), 2 * (8 + 10), 1e-9, "in-plane box b");
            check(science::norm(cell[0]) - 2 * extent >= 2 * 10 * std::cos(std::acos(-1.0) / 6) - 1e-6, "periodic images separated by the vacuum");
            // Periodic along the axis: interior atoms keep full fcc coordination across the boundary.
            const auto cn = coordination(r.structure, 2.8);
            check(*std::max_element(cn.begin(), cn.end()) == 12, "bulk coordination inside a periodic wire");
        }
    });

    test("core-shell: radial assignment", [] {
        Structure cluster = fcc(3.6, 4);
        cluster.hasUnitCell = false;
        CoreShellParams p;
        p.coreRadius = 5;
        p.coreElement = "Au";
        p.shellElement = "Ag";
        const auto r = applyCoreShell(cluster, p);
        check(r.success, r.message);
        check(r.coreCount + r.shellCount == static_cast<int>(cluster.atoms.size()), "every atom assigned");
        double cx = 0, cy = 0, cz = 0;
        for (const auto& a : cluster.atoms) { cx += a.x; cy += a.y; cz += a.z; }
        const double n = static_cast<double>(cluster.atoms.size());
        cx /= n; cy /= n; cz /= n;
        int core = 0;
        for (std::size_t i = 0; i < cluster.atoms.size(); ++i) {
            const auto& a = r.structure.atoms[i];
            const double d = std::sqrt((a.x - cx) * (a.x - cx) + (a.y - cy) * (a.y - cy) + (a.z - cz) * (a.z - cz));
            const bool inside = d < 5.0;
            check(a.symbol == (inside ? "Au" : "Ag"), "species follows the radius");
            check(a.atomicNumber == (inside ? 79 : 47), "atomic number updated with the symbol");
            core += inside;
        }
        check(core == r.coreCount && core > 0, "core count");
    });

    if (failures) {
        std::cout << failures << " builder regression(s) failed\n";
        return 1;
    }
    std::cout << "All builder tool regressions passed\n";
    return 0;
}
