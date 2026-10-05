// Analytical and independent-reference checks for the native scientific tools
// (src/science). The cases mirror python/tests/test_condensed_matter.py so the
// desktop/CLI implementation is held to the same contracts as the Python API.
#include "science/Potentials.h"
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "science/Simulation.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace atomforge::science;

namespace
{
const double kPi = std::acos(-1.0);
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

template<class Work>
void expectError(Work work, const std::string& what)
{
    try { work(); } catch (const std::exception&) { return; }
    throw std::runtime_error(what + ": expected an error");
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

// Deterministic standard normals for synthetic data.
class Normal
{
public:
    explicit Normal(std::uint64_t seed) : m_state(seed) {}
    double operator()()
    {
        auto uniform = [&] {
            m_state = m_state * 6364136223846793005ULL + 1442695040888963407ULL;
            return (static_cast<double>(m_state >> 11) + 0.5) / 9007199254740992.0;
        };
        return std::sqrt(-2 * std::log(uniform())) * std::cos(2 * kPi * uniform());
    }
private:
    std::uint64_t m_state;
};

Json run(const std::string& tool, const Json& request) { return runTool(tool, request).result; }
std::vector<double> values(const Json& value) { return toArray(value, "value").values; }
double number(const Json& value) { return value.number(); }

Json rows(const std::vector<Vec3>& points)
{
    Json result = Json::array();
    for (const auto& p : points) result.push(toJson(p));
    return result;
}

Json matrixJson(const Mat3& m) { return toJson(m); }

Json object(std::initializer_list<std::pair<const char*, Json>> members)
{
    Json result = Json::object();
    for (const auto& [key, value] : members) result[key] = value;
    return result;
}

struct Crystal { std::vector<std::string> symbols; std::vector<Vec3> positions; Mat3 cell{}; };

Crystal fccCubic(double a, int repeat, const std::string& symbol = "Cu")
{
    Crystal crystal;
    const Vec3 basis[4] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
    for (int i = 0; i < repeat; ++i)
        for (int j = 0; j < repeat; ++j)
            for (int k = 0; k < repeat; ++k)
                for (const auto& b : basis) {
                    crystal.symbols.push_back(symbol);
                    crystal.positions.push_back({(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a});
                }
    crystal.cell = {{{a * repeat, 0, 0}, {0, a * repeat, 0}, {0, 0, a * repeat}}};
    return crystal;
}

Mat3 rotation(const Vec3& axis, double angle)
{
    const Vec3 u = scale(axis, 1 / norm(axis));
    const double c = std::cos(angle), s = std::sin(angle), t = 1 - c;
    return {{{t * u[0] * u[0] + c, t * u[0] * u[1] - s * u[2], t * u[0] * u[2] + s * u[1]},
             {t * u[0] * u[1] + s * u[2], t * u[1] * u[1] + c, t * u[1] * u[2] - s * u[0]},
             {t * u[0] * u[2] - s * u[1], t * u[1] * u[2] + s * u[0], t * u[2] * u[2] + c}}};
}

// Applies x -> M x to row vectors: x M^T.
Vec3 transformVector(const Mat3& m, const Vec3& x) { return rowTimes(x, transpose(m)); }

Json structureJsonOf(const Crystal& crystal, bool periodic = true)
{
    Json symbols = Json::array();
    for (const auto& s : crystal.symbols) symbols.push(s);
    Json result = object({{"symbols", symbols}, {"positions", rows(crystal.positions)}});
    if (periodic) result["cell"] = matrixJson(crystal.cell);
    return result;
}

Configuration configuration(const Crystal& crystal, bool periodic = true)
{
    Structure structure;
    for (std::size_t i = 0; i < crystal.symbols.size(); ++i) {
        AtomSite atom;
        atom.symbol = crystal.symbols[i];
        atom.atomicNumber = atomicNumber(atom.symbol);
        atom.x = crystal.positions[i][0]; atom.y = crystal.positions[i][1]; atom.z = crystal.positions[i][2];
        structure.atoms.push_back(atom);
    }
    structure.hasUnitCell = periodic;
    for (int r = 0; r < 3; ++r) structure.cellVectors[r] = crystal.cell[r];
    return configurationFrom(structure, {periodic, periodic, periodic});
}

class CurvedDoubleWell final : public Potential
{
public:
    PotentialResult compute(const Configuration& c, bool) const override
    {
        PotentialResult result;
        for (const auto& p : c.positions) {
            const double x = p[0], y = p[1], z = p[2], residual = y - 0.2 * (1 - x * x);
            result.energy += (x * x - 1) * (x * x - 1) + 5 * residual * residual + z * z;
            result.forces.push_back({-4 * x * (x * x - 1) - 4 * x * residual, -10 * residual, -2 * z});
        }
        return result;
    }
    std::string description() const override { return "curved double well"; }
};

class FreeParticles final : public Potential
{
public:
    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        PotentialResult result;
        result.forces.assign(c.size(), {0, 0, 0});
        result.hasStress = stress;
        return result;
    }
    std::string description() const override { return "free particles"; }
};

std::filesystem::path scratch(const std::string& name)
{
    const auto folder = std::filesystem::temp_directory_path() / "atomforge-science-tests";
    std::filesystem::create_directories(folder);
    return folder / name;
}
}

int main()
{
    std::cout << "Native scientific tool regressions\n";

    test("json round trip", [] {
        const Json value = Json::parse(R"({"a": [1, 2.5, -3e-2, true, null], "b": {"c": "x\nyé"}})");
        check(value.at("a").size() == 5 && value.at("a").items()[3].boolean(), "parsed array");
        check(value.at("b").at("c").string() == "x\ny\xC3\xA9", "unicode escape");
        const Json again = Json::parse(value.dump(2));
        check(again.dump() == value.dump(), "dump/parse round trip");
        check(Json(std::nan("")).dump() == "null", "non-finite as null");
        expectError([] { Json::parse("[1, 2"); }, "unterminated array");
    });

    test("ballistic msd and drift removal", [] {
        const Vec3 velocity[2] = {{.1, .2, .3}, {-.2, .1, 0}};
        Json positions = Json::array(), translated = Json::array();
        for (int t = 0; t < 12; ++t) {
            const double time = 2.0 * t;
            positions.push(rows({scale(velocity[0], time), scale(velocity[1], time)}));
            translated.push(rows({scale(velocity[0], time), scale(velocity[0], time), scale(velocity[0], time)}));
        }
        const auto msd = values(run("msd", object({{"positions", positions}, {"timestep_fs", 2.0}})).at("msd_A2"));
        const double meanSquare = (dot(velocity[0], velocity[0]) + dot(velocity[1], velocity[1])) / 2;
        for (int t = 0; t < 12; ++t) close(msd[t], 4.0 * t * t * meanSquare, 1e-13, "ballistic msd");
        const auto drift = values(run("msd", object({{"positions", translated}, {"timestep_fs", 2.0}, {"remove_drift", true}})).at("msd_A2"));
        for (double v : drift) close(v, 0, 1e-25, "drift removed");
    });

    test("msd unwrapping in a skew cell", [] {
        const Mat3 cell = {{{5, 0, 0}, {4, 3, 0}, {0, 0, 6}}};
        Json unwrapped = Json::array(), wrapped = Json::array();
        for (int t = 0; t < 20; ++t) {
            const Vec3 f = {.13 * t, .07 * t, .01 * t};
            const Vec3 w = {f[0] - std::floor(f[0]), f[1] - std::floor(f[1]), f[2] - std::floor(f[2])};
            unwrapped.push(rows({rowTimes(f, cell)}));
            wrapped.push(rows({rowTimes(w, cell)}));
        }
        const auto a = values(run("msd", object({{"positions", unwrapped}, {"timestep_fs", 1.0}})).at("msd_A2"));
        const auto b = values(run("msd", object({{"positions", wrapped}, {"timestep_fs", 1.0}, {"wrapped", true}, {"cell", matrixJson(cell)}})).at("msd_A2"));
        for (std::size_t i = 0; i < a.size(); ++i) close(b[i], a[i], 1e-11, "unwrapped msd");
        expectError([&] { run("msd", object({{"positions", wrapped}, {"timestep_fs", 1.0}, {"wrapped", true}})); }, "wrapped needs cell");
    });

    test("diffusion units and fit", [] {
        std::vector<double> time, msd, falling;
        for (int i = 0; i <= 100; ++i) { time.push_back(10.0 * i); msd.push_back(6 * .003 * 10.0 * i + .4); falling.push_back(30 - 0.1 * i); }
        const Json result = run("diffusion", object({{"lag_fs", toJson(time)}, {"msd_A2", toJson(msd)}, {"fit_range_fs", Json::array({100.0, 800.0})}}));
        close(number(result.at("D_A2_per_fs")), .003, 1e-12, "D");
        close(number(result.at("D_m2_per_s")), 3e-8, 1e-17, "D SI");
        close(number(result.at("intercept_A2")), .4, 1e-10, "intercept");
        close(number(result.at("r_squared")), 1, 1e-12, "r squared");
        expectError([&] { run("diffusion", object({{"lag_fs", toJson(time)}, {"msd_A2", toJson(falling)}, {"fit_range_fs", Json::array({0.0, 500.0})}})); }, "negative slope");
    });

    test("brownian diffusion", [] {
        Normal normal(824);
        const double diffusion = .004;
        const int frames = 400, atoms = 300;
        NdArray positions({static_cast<std::size_t>(frames), static_cast<std::size_t>(atoms), 3});
        for (int f = 1; f < frames; ++f)
            for (int a = 0; a < atoms; ++a)
                for (int k = 0; k < 3; ++k) positions(f, a, k) = positions(f - 1, a, k) + normal() * std::sqrt(2 * diffusion);
        const Json msd = run("msd", object({{"positions", toJson(positions)}, {"timestep_fs", 1.0}, {"max_lag", 60}}));
        const Json fit = run("diffusion", object({{"lag_fs", msd.at("lag_fs")}, {"msd_A2", msd.at("msd_A2")}, {"fit_range_fs", Json::array({10.0, 60.0})}}));
        check(std::abs(number(fit.at("D_A2_per_fs")) / diffusion - 1) < .06, "Brownian D within 6%");
    });

    test("fft vacf matches direct correlation", [] {
        Normal normal(10);
        NdArray velocity({35, 4, 3});
        for (double& v : velocity.values) v = normal();
        const auto vacf = values(run("vacf", object({{"velocities", toJson(velocity)}, {"timestep_fs", .5}})).at("vacf"));
        for (int lag = 0; lag < 35; ++lag) {
            double sum = 0;
            for (int t = 0; t + lag < 35; ++t)
                for (int a = 0; a < 4; ++a)
                    for (int k = 0; k < 3; ++k) sum += velocity(t, a, k) * velocity(t + lag, a, k);
            close(vacf[lag], sum / (4.0 * (35 - lag)), 2e-13, "vacf lag " + std::to_string(lag));
        }
        close(values(run("vacf", object({{"velocities", toJson(velocity)}, {"timestep_fs", 1.0}, {"normalize", true}})).at("vacf"))[0], 1, 1e-14, "normalized");
        expectError([] { run("vacf", object({{"velocities", toJson(NdArray({35, 4, 3}))}, {"timestep_fs", 1.0}, {"normalize", true}})); }, "zero vacf");
    });

    test("vibrational spectrum frequency and normalization", [] {
        NdArray velocity({1000, 2, 3});
        for (int t = 0; t < 1000; ++t)
            for (int a = 0; a < 2; ++a) velocity(t, a, 0) = std::cos(2 * kPi * .02 * t);
        const Json result = run("vibrational-spectrum", object({{"velocities", toJson(velocity)}, {"timestep_fs", 1.0}, {"masses", Json::array({1.0, 12.0})}}));
        const auto f = values(result.at("frequency_THz")), d = values(result.at("density_per_THz"));
        std::size_t peak = 0;
        for (std::size_t i = 0; i < d.size(); ++i) if (d[i] > d[peak]) peak = i;
        close(f[peak], 20, 1e-9, "peak frequency");
        double area = 0;
        for (std::size_t i = 1; i < d.size(); ++i) area += .5 * (d[i] + d[i - 1]) * (f[i] - f[i - 1]);
        close(area, 1, 1e-12, "unit area");
        expectError([] { run("vibrational-spectrum", object({{"velocities", toJson(NdArray({1000, 2, 3}))}, {"timestep_fs", 1.0}})); }, "zero spectrum");
    });

    test("affine strain, rotation and D2min", [] {
        const Crystal atoms = fccCubic(3.6, 2);
        const Mat3 deformation = {{{1.03, .08, 0}, {0, .98, .01}, {.02, 0, 1.04}}};
        std::vector<Vec3> deformed;
        for (const auto& p : atoms.positions) deformed.push_back(transformVector(deformation, p));
        Mat3 deformedCell{};
        for (int r = 0; r < 3; ++r) deformedCell[r] = transformVector(deformation, atoms.cell[r]);
        const Json periodic = Json::array({true, true, true});
        const Json result = run("local-strain", object({{"reference", rows(atoms.positions)}, {"current", rows(deformed)}, {"cutoff_A", 2.8},
            {"reference_cell", matrixJson(atoms.cell)}, {"current_cell", matrixJson(deformedCell)}, {"pbc", periodic}}));
        const Mat3 ftf = multiply(transpose(deformation), deformation);
        const NdArray strain = toArray(result.at("green_lagrange_strain"), "strain");
        for (std::size_t a = 0; a < atoms.positions.size(); ++a) {
            check(result.at("valid").items()[a].boolean(), "valid environment");
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j) close(strain(a, i, j), (ftf[i][j] - (i == j)) / 2, 1e-13, "green strain");
        }
        for (double v : values(result.at("d2min_A2"))) close(v, 0, 1e-24, "affine D2min");
        const Mat3 r = rotation({0, 0, 1}, .3);
        std::vector<Vec3> rotated;
        for (const auto& p : atoms.positions) rotated.push_back(transformVector(r, p));
        Mat3 rotatedCell{};
        for (int k = 0; k < 3; ++k) rotatedCell[k] = transformVector(r, atoms.cell[k]);
        const Json rigid = run("local-strain", object({{"reference", rows(atoms.positions)}, {"current", rows(rotated)}, {"cutoff_A", 2.8},
            {"reference_cell", matrixJson(atoms.cell)}, {"current_cell", matrixJson(rotatedCell)}, {"pbc", periodic}}));
        for (double v : values(rigid.at("green_lagrange_strain"))) close(v, 0, 1e-13, "rotation is unstrained");
        auto perturbed = atoms.positions;
        perturbed[0][0] += .15;
        const Json nonaffine = run("local-strain", object({{"reference", rows(atoms.positions)}, {"current", rows(perturbed)}, {"cutoff_A", 2.8},
            {"reference_cell", matrixJson(atoms.cell)}, {"pbc", periodic}}));
        double largest = 0;
        for (double v : values(nonaffine.at("d2min_A2"))) largest = std::max(largest, v);
        check(largest > .01, "nonaffine D2min");
        const Json planar = rows({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}});
        const Json planarResult = run("local-strain", object({{"reference", planar}, {"current", planar}, {"cutoff_A", 2.0}}));
        for (const auto& v : planarResult.at("valid").items())
            check(!v.boolean(), "planar environments are rank deficient");
    });

    test("centrosymmetry of fcc and a distortion", [] {
        Crystal atoms = fccCubic(3.6, 2);
        const Json periodic = Json::array({true, true, true});
        for (double v : values(run("centrosymmetry", object({{"positions", rows(atoms.positions)}, {"cutoff_A", 2.8}, {"cell", matrixJson(atoms.cell)}, {"pbc", periodic}})).at("centrosymmetry_A2")))
            close(v, 0, 1e-24, "perfect fcc CSP");
        atoms.positions[0][0] += .1;
        const auto csp = values(run("centrosymmetry", object({{"positions", rows(atoms.positions)}, {"cutoff_A", 2.8}, {"cell", matrixJson(atoms.cell)}, {"pbc", periodic}})).at("centrosymmetry_A2"));
        check(csp[0] > .01, "distorted CSP");
    });

    test("fcc Steinhardt order and rotation invariance", [] {
        // One-atom primitive cell: periodic self images are the neighbours.
        Mat3 cell = {{{0, 1.8, 1.8}, {1.8, 0, 1.8}, {1.8, 1.8, 0}}};
        const Json periodic = Json::array({true, true, true});
        const Json result = run("bond-order", object({{"positions", rows({{0, 0, 0}})}, {"cutoff_A", 2.8}, {"cell", matrixJson(cell)}, {"pbc", periodic}}));
        close(values(result.at("order").at("q4"))[0], .1909406539564933, 1e-12, "q4");
        close(values(result.at("order").at("q6"))[0], .5745242597140697, 1e-12, "q6");
        const Mat3 r = rotation({1, 2, 3}, 37 * kPi / 180);
        for (auto& row : cell) row = transformVector(r, row);
        const Json rotated = run("bond-order", object({{"positions", rows({{0, 0, 0}})}, {"cutoff_A", 2.8}, {"cell", matrixJson(cell)}, {"pbc", periodic}}));
        close(values(rotated.at("order").at("q6"))[0], .5745242597140697, 1e-12, "rotated q6");
    });

    test("Wigner-Seitz Frenkel pair and box mapping", [] {
        const std::vector<Vec3> sites = {{0, 0, 0}, {2, 0, 0}, {4, 0, 0}};
        auto moved = sites;
        moved[1] = {.2, 0, 0};
        const Json result = run("wigner-seitz", object({{"reference_sites", rows(sites)}, {"positions", rows(moved)}}));
        check(result.at("vacancies").number() == 1 && result.at("interstitial_excess").number() == 1, "Frenkel pair");
        check(values(result.at("occupancy")) == std::vector<double>({2, 0, 1}), "occupancy");
        std::vector<Vec3> scaled;
        for (const auto& s : sites) scaled.push_back(scale(s, 1.2));
        const Mat3 box = {{{6, 0, 0}, {0, 6, 0}, {0, 0, 6}}}, big = {{{7.2, 0, 0}, {0, 7.2, 0}, {0, 0, 7.2}}};
        const Json mapped = run("wigner-seitz", object({{"reference_sites", rows(sites)}, {"positions", rows(scaled)}, {"cell", matrixJson(box)},
            {"current_cell", matrixJson(big)}, {"pbc", Json::array({true, true, true})}}));
        check(values(mapped.at("occupancy")) == std::vector<double>({1, 1, 1}), "mapped occupancy");
        check(run("wigner-seitz", object({{"reference_sites", rows(sites)}, {"positions", rows({{1, 0, 0}})}})).at("ambiguous").items()[0].boolean(), "tie flagged");
    });

    test("static structure factor extinction", [] {
        const Json q = rows({{0, 0, 0}, {2 * kPi, 0, 0}, {2 * kPi, 2 * kPi, 0}});
        const auto s = values(run("structure-factor", object({{"positions", rows({{0, 0, 0}, {.5, .5, .5}})}, {"q_vectors", q}})).at("S_q"));
        close(s[0], 2, 1e-13, "S(0)"); close(s[1], 0, 1e-13, "extinct"); close(s[2], 2, 1e-13, "allowed");
        const auto shifted = values(run("structure-factor", object({{"positions", rows({{7.1, -3.4, 2.6}, {7.6, -2.9, 3.1}})}, {"q_vectors", q}})).at("S_q"));
        for (int i = 0; i < 3; ++i) close(shifted[i], s[i], 1e-12, "translation invariance");
    });

    test("band gaps: direct, indirect, metal, spin", [] {
        const Json gap = run("band-gap", object({{"energies_eV", Json::parse("[[-1, 2], [-0.2, 3]]")}, {"fermi_eV", 0.0}}));
        close(number(gap.at("gap_eV")), 2.2, 1e-12, "indirect gap");
        close(number(gap.at("direct_gap_eV")), 3, 1e-12, "direct gap");
        check(!gap.at("metal_on_sampled_mesh").boolean(), "insulator");
        check(run("band-gap", object({{"energies_eV", Json::parse("[[-1, 0.1], [0.1, 2]]")}, {"fermi_eV", 0.0}})).at("metal_on_sampled_mesh").boolean(), "crossing metal");
        check(run("band-gap", object({{"energies_eV", Json::parse("[[-1, 0, 2]]")}, {"fermi_eV", 0.0}})).at("metal_on_sampled_mesh").boolean(), "touching metal");
        close(number(run("band-gap", object({{"energies_eV", Json::parse("[[[-1, 2], [-0.2, 3]], [[-2, 1], [-1, 1.5]]]")}, {"fermi_eV", 0.0}})).at("gap_eV")), 1.2, 1e-12, "spin gap");
    });

    test("effective-mass full tensor", [] {
        const Mat3 r = rotation({.3, -.1, .2}, std::sqrt(.09 + .01 + .04));
        Mat3 mass{};
        const double principal[3] = {.2, .5, 1.4};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) mass[i][j] += r[i][k] * principal[k] * r[j][k];
        Mat3 hessian = inverse(mass);
        for (auto& row : hessian) row = scale(row, 7.619964231073853);
        Normal normal(44);
        std::vector<Vec3> k;
        std::vector<double> energy;
        for (int n = 0; n < 60; ++n) {
            const Vec3 point = {normal() * .02, normal() * .02, normal() * .02};
            k.push_back(point);
            energy.push_back(2.4 + dot(point, {.01, -.02, .005}) + dot(point, rowTimes(point, hessian)) / 2);
        }
        const Json result = run("effective-mass", object({{"kpoints_inv_A", rows(k)}, {"energies_eV", toJson(energy)}, {"center_inv_A", Json::array({0.0, 0.0, 0.0})}}));
        const auto tensor = values(result.at("mass_tensor_m_e"));
        for (int i = 0; i < 9; ++i) close(tensor[i], mass[i / 3][i % 3], 1e-10, "mass tensor");
        const auto gradient = values(result.at("gradient_eV_A"));
        close(gradient[0], .01, 1e-10, "gx"); close(gradient[1], -.02, 1e-10, "gy"); close(gradient[2], .005, 1e-10, "gz");
        check(number(result.at("rms_fit_error_eV")) < 1e-12, "exact fit");
        std::vector<Vec3> line;
        std::vector<double> lineEnergy;
        for (int n = 0; n < 10; ++n) { line.push_back({double(n), 0, 0}); lineEnergy.push_back(n); }
        expectError([&] { run("effective-mass", object({{"kpoints_inv_A", rows(line)}, {"energies_eV", toJson(lineEnergy)}, {"center_inv_A", Json::array({0.0, 0.0, 0.0})}})); }, "line path");
    });

    test("work function and sloped vacuum", [] {
        std::vector<double> x, flat, sloped;
        for (int i = 0; i <= 200; ++i) { x.push_back(i * .1); flat.push_back(i * .1 < 10 ? 2.0 : 6.0); sloped.push_back(flat.back() + .1 * x.back()); }
        const Json result = run("work-function", object({{"distance_A", toJson(x)}, {"potential_eV", toJson(flat)}, {"fermi_eV", 1.5}, {"vacuum_range_A", Json::array({12.0, 18.0})}}));
        close(number(result.at("work_function_eV")), 4.5, 1e-12, "work function");
        check(result.at("flat_vacuum").boolean(), "flat vacuum");
        const Json bad = run("work-function", object({{"distance_A", toJson(x)}, {"potential_eV", toJson(sloped)}, {"fermi_eV", 1.5}, {"vacuum_range_A", Json::array({12.0, 18.0})}}));
        check(!bad.at("flat_vacuum").boolean(), "sloped vacuum flagged");
        close(number(bad.at("slope_eV_per_A")), .1, 1e-10, "slope");
    });

    test("Birch-Murnaghan recovery", [] {
        std::vector<double> volume, energy;
        const double v0 = 16.2, e0 = -4.1, b0 = .7, derivative = 4.3;
        for (int i = 0; i < 13; ++i) {
            const double v = 14 + 4.0 * i / 12;
            const double eta = std::pow(v0 / v, 2.0 / 3) - 1;
            volume.push_back(v);
            energy.push_back(e0 + 9 * v0 * b0 / 16 * (derivative * eta * eta * eta + (6 - 4 * (eta + 1)) * eta * eta));
        }
        const Json result = run("equation-of-state", object({{"volumes_A3", toJson(volume)}, {"energies_eV", toJson(energy)}}));
        close(number(result.at("volume_A3")), v0, 5e-7, "V0");
        close(number(result.at("bulk_modulus_GPa")), b0 * 160.2176634, 5e-6, "B0");
        close(number(result.at("bulk_derivative")), derivative, 5e-6, "B0'");
        expectError([&] { run("equation-of-state", object({{"volumes_A3", toJson(volume)}, {"energies_eV", toJson(volume)}})); }, "no minimum");
    });

    test("isotropic elasticity and singular inputs", [] {
        const double bulk = 120, shear = 50, lame = bulk - 2 * shear / 3;
        double tensor[6][6] = {};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) tensor[i][j] = lame + (i == j ? 2 * shear : 0);
        for (int i = 3; i < 6; ++i) tensor[i][i] = shear;
        auto build = [&](const double (&c)[6][6], bool prestress) {
            Json strains = Json::array(), stresses = Json::array();
            const double base[6] = {.1, -.2, .05, 0, .03, 0};
            for (int n = 0; n < 13; ++n) {
                double e[6] = {};
                if (n > 0) e[(n - 1) % 6] = n <= 6 ? .005 : -.005;
                Json se = Json::array(), ss = Json::array();
                for (int i = 0; i < 6; ++i) {
                    double s = prestress ? base[i] : 0;
                    for (int j = 0; j < 6; ++j) s += c[i][j] * e[j];
                    se.push(e[i]); ss.push(s);
                }
                strains.push(se); stresses.push(ss);
            }
            return object({{"strains", strains}, {"stresses_GPa", stresses}});
        };
        const Json result = run("elastic-tensor", build(tensor, true));
        const auto stiffness = values(result.at("stiffness_GPa"));
        for (int i = 0; i < 36; ++i) close(stiffness[i], tensor[i / 6][i % 6], 1e-9, "stiffness");
        close(number(result.at("moduli").at("bulk_hill_GPa")), bulk, 1e-9, "bulk");
        close(number(result.at("moduli").at("shear_hill_GPa")), shear, 1e-9, "shear");
        close(number(result.at("moduli").at("universal_anisotropy")), 0, 1e-9, "isotropy");
        check(result.at("stable_zero_prestress").boolean(), "stable");
        expectError([] { run("elastic-tensor", object({{"strains", toJson(NdArray({10, 6}))}, {"stresses_GPa", toJson(NdArray({10, 6}))}})); }, "singular design");
        tensor[5][5] = -1;
        const Json unstable = run("elastic-tensor", build(tensor, false));
        check(!unstable.at("stable_zero_prestress").boolean() && unstable.at("moduli").isNull(), "unstable tensor");
    });

    test("phonon DOS integral and weights", [] {
        std::vector<double> grid;
        for (int i = 0; i <= 3000; ++i) grid.push_back(-.05 + .15 * i / 3000);
        const Json result = run("phonon-dos", object({{"energies_eV", Json::parse("[[0.01, 0.02, 0.03], [-0.01, 0.02, 0.03]]")}, {"energy_grid_eV", toJson(grid)}, {"weights", Json::array({3.0, 1.0})}}));
        close(number(result.at("enclosed_modes")), 3, 1e-10, "enclosed modes");
        close(number(result.at("imaginary_weight")), .25, 1e-14, "imaginary weight");
    });

    test("harmonic limits and thermodynamic derivatives", [] {
        const double kb = 8.617333262145e-5;
        const Json result = run("harmonic-thermodynamics", object({{"energies_eV", Json::parse("[[0.01, 0.02, 0.03]]")}, {"temperatures_K", Json::parse("[0, 299.9, 300, 300.1, 1e7]")}}));
        close(number(result.at("zero_point_energy_eV")), .03, 1e-15, "ZPE");
        const auto cv = values(result.at("heat_capacity_eV_per_K")), f = values(result.at("free_energy_eV"));
        const auto u = values(result.at("internal_energy_eV")), s = values(result.at("entropy_eV_per_K"));
        close(cv[0], 0, 0, "Cv(0)");
        close(cv[4], 3 * kb, 1e-11, "classical Cv");
        close(-(f[3] - f[1]) / .2, s[2], 1e-10, "S = -dF/dT");
        close((u[3] - u[1]) / .2, cv[2], 1e-10, "Cv = dU/dT");
        close(number(run("harmonic-thermodynamics", object({{"energies_eV", Json::parse("[[0, 0.01, 0.02]]")}, {"temperatures_K", Json::parse("[300]")}})).at("omitted_zero_mode_weight")), 1, 0, "omitted weight");
        expectError([] { run("harmonic-thermodynamics", object({{"energies_eV", Json::parse("[[-0.001, 0.01]]")}, {"temperatures_K", Json::parse("[300]")}})); }, "unstable phonons");
    });

    test("EMT copper: equilibrium, forces and virial stress", [] {
        const auto emt = makeEmt();
        double best = 0, bestEnergy = HUGE_VAL;
        for (int i = 0; i <= 40; ++i) {
            const double a = 3.50 + .005 * i;
            const double e = emt->compute(configuration(fccCubic(a, 2)), false).energy / 32;
            if (e < bestEnergy) { bestEnergy = e; best = a; }
        }
        check(best > 3.55 && best < 3.63, "EMT Cu lattice constant near 3.59 A (got " + std::to_string(best) + ")");
        check(std::abs(bestEnergy) < .02, "EMT energies are referenced to the fcc ground state");
        Crystal crystal = fccCubic(3.6, 2);
        crystal.positions[3] = add(crystal.positions[3], {.07, -.04, .05});
        crystal.positions[9] = add(crystal.positions[9], {-.03, .06, .02});
        const auto c = configuration(crystal);
        const auto result = emt->compute(c, true);
        const double h = 1e-5;
        for (int atom : {3, 9, 0})
            for (int k = 0; k < 3; ++k) {
                auto plus = c, minus = c;
                plus.positions[atom][k] += h; minus.positions[atom][k] -= h;
                const double numeric = -(emt->compute(plus, false).energy - emt->compute(minus, false).energy) / (2 * h);
                close(result.forces[atom][k], numeric, 1e-6, "EMT force");
            }
        // sigma_ab = (1/V) dE/d(epsilon_ab) under homogeneous strain of atoms and cell.
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                auto strained = [&](double e) {
                    Mat3 f = identity();
                    f[a][b] += e / 2; f[b][a] += e / 2;
                    auto s = c;
                    for (auto& p : s.positions) p = transformVector(f, p);
                    for (auto& row : s.cell) row = transformVector(f, row);
                    return emt->compute(s, false).energy;
                };
                const double numeric = (strained(h) - strained(-h)) / (2 * h) / cellVolume(c.cell);
                close(result.stress[a][b], numeric, 1e-7, "EMT stress");
            }
    });

    test("Lennard-Jones forces, stress and cutoff shift", [] {
        const auto lj = makeLennardJones(.0104, 3.4, 8.5);
        Crystal crystal = fccCubic(5.26, 2, "Ar");
        crystal.positions[1] = add(crystal.positions[1], {.1, -.05, .08});
        const auto c = configuration(crystal);
        const auto result = lj->compute(c, true);
        const double h = 1e-5;
        for (int k = 0; k < 3; ++k) {
            auto plus = c, minus = c;
            plus.positions[1][k] += h; minus.positions[1][k] -= h;
            close(result.forces[1][k], -(lj->compute(plus, false).energy - lj->compute(minus, false).energy) / (2 * h), 1e-8, "LJ force");
        }
        auto expand = [&](double e) { auto s = c; for (auto& p : s.positions) p = scale(p, 1 + e); for (auto& r : s.cell) r = scale(r, 1 + e); return lj->compute(s, false).energy; };
        const double trace = result.stress[0][0] + result.stress[1][1] + result.stress[2][2];
        close(trace, (expand(h) - expand(-h)) / (2 * h) / cellVolume(c.cell), 1e-8, "LJ virial trace");
        Configuration pair;
        pair.symbols = {"Ar", "Ar"}; pair.numbers = {18, 18}; pair.masses = {39.948, 39.948};
        pair.positions = {{0, 0, 0}, {8.4999999, 0, 0}};
        close(lj->compute(pair, false).energy, 0, 1e-10, "energy vanishes at the cutoff");
    });

    test("NEB on a curved double well", [] {
        Configuration initial;
        initial.symbols = {"H"}; initial.numbers = {1}; initial.masses = {1.008};
        initial.positions = {{-1, 0, 0}};
        Configuration final = initial;
        final.positions = {{1, 0, 0}};
        NebOptions options;
        options.images = 7; options.fmax = .002; options.steps = 300;
        const auto result = migrationPath(initial, final, [] { return std::make_unique<CurvedDoubleWell>(); }, options);
        check(result.result.at("converged").boolean(), "NEB converged");
        close(number(result.result.at("forward_barrier_eV")), 1, 5e-6, "barrier");
        close(result.frames[3].atoms[0].y, .2, 5e-4, "curved saddle");
        for (double f : values(result.result.at("endpoint_max_forces_eV_per_A"))) close(f, 0, 1e-12, "relaxed endpoints");
    });

    test("NVT temperature control and reproducibility", [] {
        Configuration gas;
        Normal normal(1);
        for (int i = 0; i < 100; ++i) {
            gas.symbols.push_back("Ar"); gas.numbers.push_back(18); gas.masses.push_back(39.948);
            gas.positions.push_back({std::abs(normal()) * 5, std::abs(normal()) * 5, std::abs(normal()) * 5});
        }
        FreeParticles free;
        DynamicsOptions options;
        options.steps = 400; options.thermostatFs = 10; options.seed = 27; options.sampleInterval = 10;
        const auto result = nvtDynamics(gas, free, options);
        const auto temperature = values(result.result.at("temperature_K"));
        double mean = 0;
        for (std::size_t i = 10; i < temperature.size(); ++i) mean += temperature[i];
        mean /= static_cast<double>(temperature.size() - 10);
        check(std::abs(mean / 300 - 1) < .1, "mean temperature within 10% (got " + std::to_string(mean) + ")");
        options.steps = 10;
        const auto repeat = nvtDynamics(gas, free, options);
        for (int f = 0; f < 2; ++f)
            for (std::size_t a = 0; a < 100; ++a)
                for (int k = 0; k < 3; ++k) close(repeat.velocities[f][a][k], result.velocities[f][a][k], 0, "seeded reproducibility");
        check(result.frames.size() == 41, "41 frames");
        close(values(result.result.at("time_fs")).back(), 400, 0, "final time");
    });

    test("NPT compression with EMT stress", [] {
        DynamicsOptions options;
        options.steps = 30; options.timestepFs = .5; options.temperatureK = 100; options.pressureGPa = 5;
        options.thermostatFs = 20; options.barostatFs = 100; options.sampleInterval = 5;
        const auto emt = makeEmt();
        const auto result = nptDynamics(configuration(fccCubic(3.6, 2)), *emt, options);
        for (double p : values(result.result.at("pressure_GPa"))) check(std::isfinite(p), "finite pressure");
        for (double t : values(result.result.at("temperature_K"))) check(std::isfinite(t), "finite temperature");
        const auto volume = values(result.result.at("volume_A3"));
        check(volume.back() < volume.front(), "volume decreases under 5 GPa");
        check(result.frames.size() == 7, "7 frames");
    });

    test("NPT barostat reaches the target pressure on average", [] {
        DynamicsOptions options;
        options.steps = 3000; options.timestepFs = 2; options.temperatureK = 300; options.pressureGPa = 2;
        options.thermostatFs = 50; options.barostatFs = 300; options.sampleInterval = 10;
        const auto emt = makeEmt();
        const auto result = nptDynamics(configuration(fccCubic(3.6, 3)), *emt, options);
        const auto pressure = values(result.result.at("pressure_GPa")), temperature = values(result.result.at("temperature_K"));
        double p = 0, t = 0;
        for (std::size_t i = pressure.size() / 2; i < pressure.size(); ++i) { p += pressure[i]; t += temperature[i]; }
        p /= static_cast<double>(pressure.size() - pressure.size() / 2);
        t /= static_cast<double>(pressure.size() - pressure.size() / 2);
        check(std::abs(p - 2) < .5, "mean pressure near 2 GPa (got " + std::to_string(p) + ")");
        check(std::abs(t / 300 - 1) < .15, "mean temperature near 300 K (got " + std::to_string(t) + ")");
    });

#ifdef ATOMS_ENABLE_SPGLIB
    test("HPKOT symmetry paths", [] {
        Crystal silicon;
        silicon.cell = {{{5.43, 0, 0}, {0, 5.43, 0}, {0, 0, 5.43}}};
        const Vec3 fcc[4] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
        for (const auto& b : fcc)
            for (const Vec3& shift : {Vec3{0, 0, 0}, Vec3{.25, .25, .25}}) {
                silicon.symbols.push_back("Si");
                silicon.positions.push_back(rowTimes(add(b, shift), silicon.cell));
            }
        const Json si = run("reciprocal-path", object({{"structure", structureJsonOf(silicon)}, {"spacing_inv_A", .1}}));
        check(si.at("spacegroup_number").number() == 227, "Si space group 227");
        check(si.at("bravais_lattice_extended").string() == "cF2", "cF2 path");
        check(si.at("primitive_structure").at("symbols").size() == 2, "2-atom primitive cell");
        const NdArray frac = toArray(si.at("kpoints_fractional"), "k"), cart = toArray(si.at("kpoints_inv_A"), "k"), rec = toArray(si.at("reciprocal_lattice_inv_A"), "b");
        for (std::size_t n = 0; n < frac.shape[0]; ++n)
            for (std::size_t k = 0; k < 3; ++k)
                close(frac(n, 0) * rec(0, k) + frac(n, 1) * rec(1, k) + frac(n, 2) * rec(2, k), cart(n, k), 1e-12, "basis consistency");
        bool gamma = false;
        for (const auto& label : si.at("labels").items()) gamma = gamma || label.string() == "GAMMA";
        check(gamma && si.at("segments").size() > 1, "labelled multi-segment path");
        // The X point of the fcc path sits at |k| = 2 pi / a.
        const NdArray x = toArray(si.at("special_points_fractional").at("X"), "X");
        Vec3 xk{};
        for (int k = 0; k < 3; ++k) xk[k] = x(0) * rec(0, k) + x(1) * rec(1, k) + x(2) * rec(2, k);
        close(norm(xk), 2 * kPi / 5.43, 1e-9, "X point distance");

        // Zincblende lacks inversion: without time reversal the -k path is added.
        Crystal gaas = silicon;
        for (std::size_t i = 1; i < gaas.symbols.size(); i += 2) gaas.symbols[i] = "As";
        for (std::size_t i = 0; i < gaas.symbols.size(); i += 2) gaas.symbols[i] = "Ga";
        const Json zb = run("reciprocal-path", object({{"structure", structureJsonOf(gaas)}, {"time_reversal", false}}));
        check(zb.at("spacegroup_number").number() == 216 && zb.at("augmented_path").boolean(), "augmented zincblende path");
        check(zb.at("bravais_lattice_extended").string() == "cF2", "space groups 207-230 use the cF2 path");

        Crystal magnesium;
        const double a = 3.21, c = 5.21;
        magnesium.cell = {{{a, 0, 0}, {-a / 2, a * std::sqrt(3.0) / 2, 0}, {0, 0, c}}};
        magnesium.symbols = {"Mg", "Mg"};
        magnesium.positions = {rowTimes({1.0 / 3, 2.0 / 3, .25}, magnesium.cell), rowTimes({2.0 / 3, 1.0 / 3, .75}, magnesium.cell)};
        const Json hcp = run("reciprocal-path", object({{"structure", structureJsonOf(magnesium)}}));
        check(hcp.at("spacegroup_number").number() == 194 && hcp.at("bravais_lattice_extended").string() == "hP2", "hcp hP2 path");

        Crystal triclinic;
        triclinic.cell = {{{4.1, 0, 0}, {.7, 5.3, 0}, {-.9, .6, 6.2}}};
        triclinic.symbols = {"Cu", "Ag"};
        triclinic.positions = {rowTimes({.1, .2, .3}, triclinic.cell), rowTimes({.6, .55, .8}, triclinic.cell)};
        const Json ap = run("reciprocal-path", object({{"structure", structureJsonOf(triclinic)}}));
        const std::string kind = ap.at("bravais_lattice_extended").string();
        check(kind == "aP2" || kind == "aP3", "triclinic path");
        Crystal bcc;
        bcc.cell = {{{3.3, 0, 0}, {0, 3.3, 0}, {0, 0, 3.3}}};
        bcc.symbols = {"Nb", "Nb"};
        bcc.positions = {{0, 0, 0}, {1.65, 1.65, 1.65}};
        const Json body = run("reciprocal-path", object({{"structure", structureJsonOf(bcc)}}));
        check(body.at("spacegroup_number").number() == 229 && body.at("primitive_structure").at("symbols").size() == 1, "bcc primitive");
    });
#endif

    test("file references: extXYZ velocities, CSV columns and validation", [] {
        Normal normal(5);
        std::vector<Structure> frames;
        std::vector<std::vector<Vec3>> velocities;
        std::vector<double> times;
        NdArray inline_({12, 3, 3});
        for (int f = 0; f < 12; ++f) {
            Structure frame;
            frame.hasUnitCell = true;
            frame.cellVectors = {{{10, 0, 0}, {0, 10, 0}, {0, 0, 10}}};
            std::vector<Vec3> v;
            for (int a = 0; a < 3; ++a) {
                AtomSite atom;
                atom.symbol = a == 1 ? "O" : "H";
                atom.atomicNumber = a == 1 ? 8 : 1;
                atom.x = f * .1 + a; atom.y = a; atom.z = 0;
                frame.atoms.push_back(atom);
                v.push_back({normal() * .01, normal() * .01, normal() * .01});
                for (int k = 0; k < 3; ++k) inline_(f, a, k) = v.back()[k];
            }
            frames.push_back(frame);
            velocities.push_back(v);
            times.push_back(f * 2.0);
        }
        const auto path = scratch("velocities.extxyz");
        writeExtxyz(path, frames, velocities, times);
        const auto fromFile = values(run("vacf", object({{"velocities", object({{"file", path.u8string()}})}, {"timestep_fs", 2.0}})).at("vacf"));
        const auto fromInline = values(run("vacf", object({{"velocities", toJson(inline_)}, {"timestep_fs", 2.0}})).at("vacf"));
        for (std::size_t i = 0; i < fromFile.size(); ++i) close(fromFile[i], fromInline[i], 1e-12 + 1e-9 * std::abs(fromInline[i]), "velocity round trip");
        expectError([&] { run("vacf", object({{"velocities", object({{"file", path.u8string()}})}, {"timestep_fs", 1.0}})); }, "timestep mismatch");
        const auto msd = values(run("msd", object({{"positions", object({{"file", path.u8string()}})}, {"timestep_fs", 2.0}})).at("msd_A2"));
        close(msd[3], .09, 1e-12, "positions from extXYZ");

        const auto csv = scratch("table.csv");
        { std::ofstream out(csv); out << "# lag, msd\n0, 1\n10, 7\n20, 13\n30, 19\n"; }
        const Json fit = run("diffusion", object({{"lag_fs", object({{"file", csv.u8string()}, {"column", 0}})},
            {"msd_A2", object({{"file", csv.u8string()}, {"column", 1}})}, {"fit_range_fs", Json::array({0.0, 30.0})}, {"dimensions", 3}}));
        close(number(fit.at("slope_A2_per_fs")), .6, 1e-12, "CSV columns");
        // VASP EIGENVAL (ISPIN=2): two k points, two bands per spin.
        const auto eigenval = scratch("EIGENVAL");
        { std::ofstream out(eigenval); out << "    2    2    1    2\n  0.1E+02  0.4E-09  0.4E-09  0.4E-09  0.5E-15\n  1.0E-004\n  CAR\n test\n     8     2     2\n\n"
              "  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.5000000E+00\n    1   -1.0000   -2.0000   1.0000   1.0000\n    2    2.0000    1.0000   0.0000   0.0000\n\n"
              "  0.5000000E+00  0.0000000E+00  0.0000000E+00  0.5000000E+00\n    1   -0.2000   -1.0000   1.0000   1.0000\n    2    3.0000    1.5000   0.0000   0.0000\n"; }
        close(number(run("band-gap", object({{"energies_eV", object({{"file", eigenval.u8string()}})}, {"fermi_eV", 0.0}})).at("gap_eV")), 1.2, 1e-12, "EIGENVAL spin gap");
        expectError([] { run("msd", object({{"positions", Json::parse("[[[0,0,0]],[[1,0,0]]]")}, {"timestep_fs", 1.0}, {"bogus", 1}})); }, "unknown parameter");
        expectError([] { run("msd", object({{"timestep_fs", 1.0}})); }, "missing required parameter");
        expectError([] { run("nvt", object({{"structure", Json::parse(R"({"symbols":["Cu","Cu"],"positions":[[0,0,0],[2.5,0,0]]})")}, {"calculator", object({{"module", "ase.calculators.vasp"}, {"attribute", "Vasp"}})}})); }, "unsupported calculator");
        const std::string report = resultReport("band-gap", run("band-gap", object({{"energies_eV", Json::parse("[[-1, 2], [-0.2, 3]]")}, {"fermi_eV", 0.0}})));
        check(report.find("gap_eV: 2.2") != std::string::npos, "report summary");
    });

    test("catalog defaults are valid requests", [] {
        check(scienceToolCatalog().size() == 20, "twenty tools");
        for (const auto& tool : scienceToolCatalog())
            for (const auto& parameter : tool.parameters) {
                if (!parameter.value[0]) continue;
                const Json value = Json::parse(parameter.value);
                if (std::string(parameter.kind) == "calculator") potentialFactory(value)();
            }
    });

    test("simulation tools through the request interface", [] {
        const Crystal copper = fccCubic(3.6, 2);
        const Json nvt = runTool("nvt", object({{"structure", structureJsonOf(copper)}, {"calculator", object({{"potential", "EMT"}})},
            {"steps", 20}, {"sample_interval", 5}, {"temperature_K", 50.0}})).result;
        check(nvt.at("frames").size() == 5 && nvt.at("ensemble").string() == "nvt", "NVT request");
        // The ASE-style EMT description from earlier releases maps to the native potential.
        Crystal moved = copper;
        moved.positions[0] = add(moved.positions[0], {.9, .9, 0});
        const Json neb = runTool("neb", object({{"initial", structureJsonOf(copper)}, {"final", structureJsonOf(copper)},
            {"calculator_factory", object({{"module", "ase.calculators.emt"}, {"attribute", "EMT"}})}, {"images", 3}, {"steps", 5}})).result;
        check(neb.at("images").size() == 3, "NEB request");
    });

    if (failures) {
        std::cout << failures << " scientific tool regression(s) failed\n";
        return 1;
    }
    std::cout << "All native scientific tool regressions passed\n";
    return 0;
}
