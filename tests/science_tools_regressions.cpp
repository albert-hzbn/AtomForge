// Analytical and independent-reference checks for the native scientific tools
// (src/science). The cases mirror python/tests/test_condensed_matter.py so the
// desktop/CLI implementation is held to the same contracts as the Python API.
#include "science/AtomProperties.h"
#include "science/LammpsExport.h"
#include "science/Batch.h"
#include "science/GifWriter.h"
#include "science/Phonons.h"
#include "science/Potentials.h"
#include "science/ResultPlots.h"
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "science/Simulation.h"
#include "science/VaspElectronic.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
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

    test("relaxation: cell, pressure, shape, positions and clusters", [] {
        const auto emt = makeEmt();
        // Reference lattice constant from a fine energy scan of the perfect crystal.
        double reference = 0, lowest = HUGE_VAL;
        for (int i = 0; i <= 2000; ++i) {
            const double a = 3.55 + 0.00005 * i;
            const double e = emt->compute(configuration(fccCubic(a, 2)), false).energy;
            if (e < lowest) { lowest = e; reference = a; }
        }
        RelaxOptions options;
        options.fmax = 1e-4; options.steps = 3000; options.relaxCell = true;
        const auto relaxed = relaxConfiguration(configuration(fccCubic(3.70, 2)), *emt, options);
        check(relaxed.converged, "cell relaxation converged");
        const double a = std::cbrt(cellVolume(relaxed.configuration.cell) / 8);
        close(a, reference, 2e-4, "relaxed lattice constant matches the energy scan");
        close(relaxed.forces.energy, lowest, 1e-6, "relaxed energy matches the scan minimum");
        const auto& s = relaxed.forces.stress;
        close(-(s[0][0] + s[1][1] + s[2][2]) / 3 * 160.2176634, 0, 0.02, "zero residual pressure");

        options.pressureGPa = 5;
        const auto compressed = relaxConfiguration(configuration(fccCubic(3.60, 2)), *emt, options);
        check(compressed.converged, "pressurised relaxation converged");
        const auto& sp = compressed.forces.stress;
        close(-(sp[0][0] + sp[1][1] + sp[2][2]) / 3 * 160.2176634, 5, 0.02, "target pressure reached");
        check(cellVolume(compressed.configuration.cell) < cellVolume(relaxed.configuration.cell), "compression reduces volume");

        options.pressureGPa = 0;
        Crystal distorted = fccCubic(3.6, 2);
        const Mat3 shear = {{{1.03, .02, 0}, {0, .97, .01}, {0, 0, 1.01}}};
        for (auto& p : distorted.positions) p = transformVector(shear, p);
        for (auto& r : distorted.cell) r = transformVector(shear, r);
        const auto shaped = relaxConfiguration(configuration(distorted), *emt, options);
        check(shaped.converged, "shape relaxation converged");
        const Mat3& c = shaped.configuration.cell;
        close(norm(c[0]), norm(c[1]), 1e-3, "|a| = |b|");
        close(norm(c[1]), norm(c[2]), 1e-3, "|b| = |c|");
        close(dot(c[0], c[1]) / norm(c[0]) / norm(c[1]), 0, 1e-4, "gamma = 90");
        close(dot(c[0], c[2]) / norm(c[0]) / norm(c[2]), 0, 1e-4, "beta = 90");
        close(shaped.forces.energy, lowest, 1e-5, "recovers the cubic minimum");

        Crystal rattled = fccCubic(reference, 2);
        Normal normal(3);
        for (auto& p : rattled.positions) p = add(p, {normal() * .08, normal() * .08, normal() * .08});
        RelaxOptions positions;
        positions.fmax = 1e-4; positions.steps = 3000;
        const auto settled = relaxConfiguration(configuration(rattled), *emt, positions);
        check(settled.converged, "position relaxation converged");
        check(settled.enthalpies.back() < settled.enthalpies.front(), "energy decreases");
        close(settled.forces.energy, lowest, 1e-6, "rattled crystal returns to the perfect energy");

        Configuration trimer;
        trimer.symbols = {"Ar", "Ar", "Ar"}; trimer.numbers = {18, 18, 18}; trimer.masses = {39.948, 39.948, 39.948};
        trimer.positions = {{0, 0, 0}, {3.6, 0.2, 0}, {1.5, 3.3, 0.3}};
        const auto lj = makeLennardJones(.0104, 3.4, 20);
        RelaxOptions tight = positions;
        tight.fmax = 1e-8; tight.steps = 20000;
        const auto cluster = relaxConfiguration(trimer, *lj, tight);
        check(cluster.converged, "trimer converged");
        const double r0 = std::pow(2.0, 1.0 / 6) * 3.4;
        for (int i = 0; i < 3; ++i)
            close(norm(sub(cluster.configuration.positions[i], cluster.configuration.positions[(i + 1) % 3])), r0, 1e-4, "equilateral LJ trimer");
        RelaxOptions bad = positions;
        bad.relaxCell = true;
        expectError([&] { relaxConfiguration(trimer, *lj, bad); }, "cell relaxation needs a cell");

        const Json request = runTool("relax", object({{"structure", structureJsonOf(fccCubic(3.7, 2))}, {"calculator", object({{"potential", "EMT"}})},
            {"relax_cell", true}, {"fmax", 1e-3}})).result;
        check(request.at("converged").boolean() && std::abs(request.at("pressure_GPa").number()) < .2, "relax request");
        check(request.at("energy_change_eV").number() < 0, "relax lowers the energy");
    });

    test("result plots", [] {
        NdArray positions({6, 2, 3});
        for (std::size_t f = 0; f < 6; ++f) positions(f, 0, 0) = positions(f, 1, 1) = 0.5 * f;
        const auto msd = resultPlots("msd", run("msd", object({{"positions", toJson(positions)}, {"timestep_fs", 2.0}})));
        check(msd.size() == 1 && msd[0].series.size() == 4 && msd[0].series[0].x.size() == 6, "MSD plot with total and components");
        close(msd[0].series[0].x[5], 10, 0, "MSD x axis is lag time");

        std::vector<double> volume, energy;
        for (int i = 0; i < 9; ++i) { volume.push_back(14 + .5 * i); energy.push_back(.02 * std::pow(volume.back() - 16.1, 2) - 3); }
        const auto eos = resultPlots("equation-of-state", run("equation-of-state", object({{"volumes_A3", toJson(volume)}, {"energies_eV", toJson(energy)}})));
        check(eos.size() == 1 && eos[0].series.size() == 2 && eos[0].series[0].points && eos[0].series[0].x.size() == 9, "EOS data and fit");
        const auto& fit = eos[0].series[1];
        const std::size_t best = static_cast<std::size_t>(std::min_element(fit.y.begin(), fit.y.end()) - fit.y.begin());
        close(fit.x[best], 16.1, .03, "fit curve minimum at V0");

        std::vector<double> lag, line;
        for (int i = 0; i <= 50; ++i) { lag.push_back(i); line.push_back(.6 * i + 2); }
        const auto diffusion = resultPlots("diffusion", run("diffusion", object({{"lag_fs", toJson(lag)}, {"msd_A2", toJson(line)}, {"fit_range_fs", Json::array({10.0, 40.0})}})));
        check(diffusion.size() == 1 && diffusion[0].series[1].x.size() == 2, "diffusion fit line");
        close(diffusion[0].series[1].y[1], .6 * 40 + 2, 1e-9, "fit line endpoint");

        RelaxOptions options;
        options.fmax = 1e-3;
        Crystal rattled = fccCubic(3.6, 2);
        rattled.positions[0][0] += .1;
        const auto emt = makeEmt();
        const auto relax = resultPlots("relax", relaxStructure(configuration(rattled), *emt, options).result);
        check(relax.size() == 2 && relax[1].logY, "relaxation energy and log force plots");

        Configuration initial;
        initial.symbols = {"H"}; initial.numbers = {1}; initial.masses = {1.008};
        initial.positions = {{-1, 0, 0}};
        Configuration final = initial;
        final.positions = {{1, 0, 0}};
        NebOptions neb;
        neb.images = 5; neb.fmax = .01;
        const auto path = resultPlots("neb", migrationPath(initial, final, [] { return std::make_unique<CurvedDoubleWell>(); }, neb).result);
        check(path.size() == 1, "NEB plot");
        close(path[0].series[0].y.front(), 0, 0, "energies relative to the initial image");
        for (std::size_t i = 1; i < path[0].series[0].x.size(); ++i) check(path[0].series[0].x[i] > path[0].series[0].x[i - 1], "increasing reaction coordinate");
        check(path[0].series[0].x.back() > 2.0, "path length exceeds the straight-line distance on a curved path");

        DynamicsOptions dynamics;
        dynamics.steps = 20; dynamics.sampleInterval = 5;
        check(resultPlots("nvt", nvtDynamics(configuration(fccCubic(3.6, 2)), *emt, dynamics).result).size() == 2, "NVT plots");
        check(resultPlots("npt", nptDynamics(configuration(fccCubic(3.6, 2)), *emt, dynamics).result).size() == 4, "NPT plots");
        check(resultPlots("band-gap", run("band-gap", object({{"energies_eV", Json::parse("[[-1, 2]]")}, {"fermi_eV", 0.0}}))).empty(), "no curve for scalar results");

        const auto h = histogram("values", {1, 2, 2, 3, std::nan(""), 3, 3}, 3);
        double total = 0;
        for (double c : h.y) total += c;
        close(total, 6, 0, "histogram counts every finite value once");
        close(h.y[2], 3, 0, "upper bin includes the maximum");
        const std::string csv = plotCsv(msd[0]);
        check(std::count(csv.begin(), csv.end(), '\n') == 7, "CSV header plus one row per sample");
    });

    test("per-atom properties and viewport colouring", [] {
        const auto low = viridis(0), high = viridis(1);
        // Polynomial fit to matplotlib viridis: anchors within its ~0.015 fit error.
        close(low[0], .267, .015, "viridis(0) r"); close(low[2], .329, .015, "viridis(0) b");
        close(high[0], .993, .015, "viridis(1) r"); close(high[1], .906, .015, "viridis(1) g");
        double previous = -1;
        for (int i = 0; i <= 20; ++i) {
            const auto c = viridis(i / 20.0);
            const double luminance = .2126 * c[0] + .7152 * c[1] + .0722 * c[2];
            check(luminance > previous, "viridis luminance increases monotonically");
            previous = luminance;
        }
        // Simple shear gamma: E = (F^T F - I)/2 gives a known von Mises invariant.
        const Crystal atoms = fccCubic(3.6, 2);
        const double gamma = .04;
        const Mat3 shear = {{{1, gamma, 0}, {0, 1, 0}, {0, 0, 1}}};
        std::vector<Vec3> sheared;
        for (const auto& p : atoms.positions) sheared.push_back(transformVector(shear, p));
        Mat3 shearedCell{};
        for (int r = 0; r < 3; ++r) shearedCell[r] = transformVector(shear, atoms.cell[r]);
        const Json strain = run("local-strain", object({{"reference", rows(atoms.positions)}, {"current", rows(sheared)}, {"cutoff_A", 2.8},
            {"reference_cell", matrixJson(atoms.cell)}, {"current_cell", matrixJson(shearedCell)}, {"pbc", Json::array({true, true, true})}}));
        const auto properties = perAtomProperties("local-strain", strain);
        check(properties.size() == 4 && properties[1].name == "Von Mises shear strain", "local-strain properties");
        const double exy = gamma / 2, eyy = gamma * gamma / 2;
        const double expected = std::sqrt(exy * exy + (eyy * eyy + eyy * eyy) / 6);
        for (double v : properties[1].values) close(v, expected, 1e-12, "von Mises shear strain");
        check(properties[0].values.size() == atoms.positions.size(), "aligned with atoms");
        const Json csp = run("centrosymmetry", object({{"positions", rows(atoms.positions)}, {"cutoff_A", 2.8}, {"cell", matrixJson(atoms.cell)}, {"pbc", Json::array({true, true, true})}}));
        check(perAtomProperties("centrosymmetry", csp).size() == 1, "CSP property");
        check(perAtomProperties("band-gap", Json::object()).empty(), "no per-atom data for scalar tools");

        const std::vector<double> values = {0, 1, 2, 3, std::nan("")};
        std::vector<std::array<float, 3>> colours;
        std::vector<bool> visible;
        PropertyDisplay display;
        colourByProperty(values, display, colours, visible);
        check(colours[0] == viridis(0) && colours[3] == viridis(1), "automatic range spans the data");
        check(colours[4][0] == colours[4][1] && visible[4], "invalid atoms grey and shown by default");
        display.autoRange = false;
        display.range = {1, 2};
        display.hideOutside = true;
        display.hideInvalid = true;
        colourByProperty(values, display, colours, visible);
        check(!visible[0] && visible[1] && visible[2] && !visible[3] && !visible[4], "range and invalid filtering");
        check(colours[2] == viridis(1), "manual range");

        Structure structure;
        for (int i = 0; i < 3; ++i) structure.atoms.push_back({"Cu", 29, double(i), 0, 0});
        structure.atomProperty = {10, 20, 30};
        structure.eraseAtom(1);
        check(structure.atomProperty == std::vector<double>({10, 30}), "property stays aligned after deleting an atom");
    });

    test("VASP and structure-file inputs", [] {
        const auto xdatcar = scratch("XDATCAR");
        { std::ofstream out(xdatcar); out << "Cu\n1.0\n4 0 0\n0 4 0\n0 0 4\nCu\n2\nDirect configuration= 1\n0 0 0\n0.5 0.5 0.5\n"
                                           "Direct configuration= 2\n0.1 0 0\n0.5 0.5 0.6\n"; }
        auto frames = readFrames(xdatcar);
        check(frames.size() == 2 && frames[1].structure.atoms.size() == 2, "constant-cell XDATCAR frames");
        close(frames[1].structure.atoms[1].z, 2.4, 1e-12, "direct to Cartesian");
        const auto variable = scratch("XDATCAR_variable");
        { std::ofstream out(variable); out << "Cu\n1.0\n4 0 0\n0 4 0\n0 0 4\nCu\n1\nDirect configuration= 1\n0.5 0 0\n"
                                            "Cu\n1.0\n5 0 0\n0 5 0\n0 0 5\nCu\n1\nDirect configuration= 2\n0.5 0 0\n"; }
        frames = readFrames(variable);
        check(frames.size() == 2, "variable-cell XDATCAR frames");
        close(frames[1].structure.atoms[0].x, 2.5, 1e-12, "repeated header updates the cell");
        expectError([&] { run("msd", object({{"positions", object({{"file", variable.u8string()}})}, {"timestep_fs", 1.0},
            {"wrapped", true}, {"cell", object({{"file", variable.u8string()}})}})); }, "variable cell rejected for fixed-cell input");

        const auto poscar = scratch("POSCAR");
        { std::ofstream out(poscar); out << "NaCl\n-179.406144\n1 0 0\n0 1 0\n0 0 1\nNa Cl\n1 1\nSelective dynamics\nCartesian\n0 0 0 T T T\n2.82 2.82 2.82 F F F\n"; }
        frames = readFrames(poscar);
        check(frames.size() == 1 && frames[0].structure.atoms[1].symbol == "Cl", "POSCAR species");
        close(frames[0].structure.cellVectors[0][0], std::cbrt(179.406144), 1e-9, "negative scale is the cell volume");
        close(frames[0].structure.atoms[1].x, 2.82, 1e-12, "Cartesian coordinates are not rescaled");

        // An active structure saved as JSON supplies its positions to array inputs.
        const Crystal copper = fccCubic(3.6, 2);
        const auto active = scratch("active.json");
        { std::ofstream out(active); out << structureJsonOf(copper).dump(); }
        const Json fromFile = run("centrosymmetry", object({{"positions", object({{"file", active.u8string()}})}, {"cutoff_A", 2.8},
            {"cell", matrixJson(copper.cell)}, {"pbc", Json::array({true, true, true})}}));
        check(fromFile.at("valid").size() == copper.positions.size(), "active-structure JSON positions");
        for (double v : values(fromFile.at("centrosymmetry_A2"))) close(v, 0, 1e-20, "perfect crystal from active JSON");
    });

    test("EAM setfl, Finnis-Sinclair and funcfl tables", [] {
        // Analytic two-element model; the tables must reproduce it exactly up to interpolation error.
        const double rc = 5.5;
        auto cutoff = [&](double r) { return r < rc ? std::pow((rc - r) / rc, 3) : 0.0; };
        auto embed = [](int e, double rho) { return e == 0 ? 0.08 * rho * rho - 1.1 * rho : 0.05 * rho * rho - 1.4 * rho; };
        auto dens = [&](int from, int at, double r) {  // FS: density at an `at` atom from a `from` atom
            const double base = (from == 0 ? 1.0 : 1.3) * std::exp(-1.2 * (r - 2.5)) * cutoff(r);
            return base * (from == at ? 1.0 : 0.8);
        };
        auto pair = [&](int a, int b, double r) {
            const double d = a + b == 0 ? .30 : (a + b == 1 ? .35 : .40);
            const double x = std::exp(-1.5 * (r - 2.55));
            return d * (x * x - 2 * x) * cutoff(r);
        };
        const int nrho = 4001, nr = 6001;
        const double drho = 40.0 / (nrho - 1), dr = rc / (nr - 1);
        auto writeTable = [](std::ofstream& out, const std::function<double(int)>& f, int n) {
            out << std::setprecision(16);
            for (int k = 0; k < n; ++k) out << f(k) << ((k % 5 == 4) ? '\n' : ' ');
            out << '\n';
        };
        auto writeAlloy = [&](const std::filesystem::path& path, bool fs) {
            std::ofstream out(path);
            out << std::setprecision(16) << "analytic test potential\nline 2\nline 3\n2 Cu Ni\n" << nrho << ' ' << drho << ' ' << nr << ' ' << dr << ' ' << rc << '\n';
            for (int e = 0; e < 2; ++e) {
                out << (e == 0 ? "29 63.546 3.6 fcc\n" : "28 58.693 3.52 fcc\n");
                writeTable(out, [&](int k) { return embed(e, k * drho); }, nrho);
                for (int at = 0; at < (fs ? 2 : 1); ++at)
                    writeTable(out, [&](int k) { return fs ? dens(e, at, k * dr) : dens(e, e, k * dr); }, nr);
            }
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b <= a; ++b) writeTable(out, [&](int k) { return k * dr * pair(a, b, k * dr); }, nr);
        };
        const auto setfl = scratch("analytic.eam.alloy"), fs = scratch("analytic.eam.fs");
        writeAlloy(setfl, false);
        writeAlloy(fs, true);

        Crystal alloy = fccCubic(3.58, 3);
        Normal normal(12);
        std::vector<int> type;
        for (std::size_t i = 0; i < alloy.positions.size(); ++i) {
            const bool nickel = normal() > 0;
            alloy.symbols[i] = nickel ? "Ni" : "Cu";
            type.push_back(nickel ? 1 : 0);
            alloy.positions[i] = add(alloy.positions[i], {normal() * .05, normal() * .05, normal() * .05});
        }
        const auto c = configuration(alloy);
        auto analytic = [&](bool finnisSinclair) {
            std::vector<double> rho(c.size(), 0.0);
            double energy = 0;
            for (const auto& nb : neighborList(c.positions, c.cell, c.pbc, rc)) {
                const double r = norm(nb.vector);
                const int i = type[static_cast<std::size_t>(nb.i)], j = type[static_cast<std::size_t>(nb.j)];
                rho[static_cast<std::size_t>(nb.i)] += finnisSinclair ? dens(j, i, r) : dens(j, j, r);
                energy += 0.5 * pair(i, j, r);
            }
            for (std::size_t i = 0; i < c.size(); ++i) energy += embed(type[i], rho[i]);
            return energy;
        };
        for (const bool finnisSinclair : {false, true}) {
            const auto eam = makeEam(finnisSinclair ? fs : setfl);
            const auto result = eam->compute(c, true);
            close(result.energy, analytic(finnisSinclair), 1e-6, std::string(finnisSinclair ? "FS" : "setfl") + " energy matches the analytic model");
            const double h = 1e-5;
            for (int atom : {0, 7, 50})
                for (int k = 0; k < 3; ++k) {
                    auto plus = c, minus = c;
                    plus.positions[atom][k] += h; minus.positions[atom][k] -= h;
                    close(result.forces[atom][k], -(eam->compute(plus, false).energy - eam->compute(minus, false).energy) / (2 * h), 2e-6, "EAM force");
                }
            for (int a = 0; a < 3; ++a)
                for (int b = a; b < 3; ++b) {
                    auto strained = [&](double e) {
                        Mat3 f = identity();
                        f[a][b] += e / 2; f[b][a] += e / 2;
                        auto s = c;
                        for (auto& p : s.positions) p = transformVector(f, p);
                        for (auto& row : s.cell) row = transformVector(f, row);
                        return eam->compute(s, false).energy;
                    };
                    close(result.stress[a][b], (strained(h) - strained(-h)) / (2 * h) / cellVolume(c.cell), 1e-7, "EAM stress");
                }
        }

        // funcfl: single element with a repulsive pair from the effective charge Z(r).
        const auto funcfl = scratch("analytic.eam");
        auto repulsive = [&](double r) { return 3.0 * std::exp(-1.1 * r) * cutoff(r); };
        {
            std::ofstream out(funcfl);
            out << std::setprecision(16) << "single element analytic\n29 63.546 3.6 fcc\n" << nrho << ' ' << drho << ' ' << nr << ' ' << dr << ' ' << rc << '\n';
            writeTable(out, [&](int k) { return embed(0, k * drho); }, nrho);
            writeTable(out, [&](int k) { return std::sqrt(k * dr * repulsive(k * dr) / (27.2 * 0.529)); }, nr);
            writeTable(out, [&](int k) { return dens(0, 0, k * dr); }, nr);
        }
        const auto copper = configuration(fccCubic(3.6, 3));
        double expected = 0;
        std::vector<double> rho(copper.size(), 0.0);
        for (const auto& nb : neighborList(copper.positions, copper.cell, copper.pbc, rc)) {
            const double r = norm(nb.vector);
            rho[static_cast<std::size_t>(nb.i)] += dens(0, 0, r);
            expected += 0.5 * repulsive(r);
        }
        for (double value : rho) expected += embed(0, value);
        close(makeEam(funcfl)->compute(copper, false).energy, expected, 1e-6, "funcfl energy");

        // Request interface: relative file resolved against the request directory.
        const Json relaxed = runTool("relax", object({{"structure", structureJsonOf(fccCubic(3.58, 2))},
            {"calculator", object({{"potential", "EAM"}, {"file", setfl.filename().u8string()}})}, {"relax_cell", true}, {"fmax", 1e-3}}),
            setfl.parent_path()).result;
        check(relaxed.at("converged").boolean() && relaxed.at("potential").string().find("EAM") == 0, "EAM relaxation through a request");
        Crystal gold = fccCubic(4.08, 1, "Au");
        expectError([&] { makeEam(setfl)->compute(configuration(gold), false); }, "element missing from the EAM file");
        expectError([&] { potentialFactory(object({{"potential", "EAM"}}))(); }, "EAM needs a file");
    });

    test("phonons: analytic Lennard-Jones dispersion and EMT copper", [] {
        // Primitive fcc argon with Lennard-Jones; cutoff between the 2nd and 3rd shells.
        const double a = 5.26, epsilon = .0104, sigma = 3.4, rc = 6.0, mass = 39.948;
        Configuration unit;
        unit.symbols = {"Ar"}; unit.numbers = {18}; unit.masses = {mass}; unit.pbc = {true, true, true};
        unit.positions = {{0, 0, 0}};
        unit.cell = {{{0, a / 2, a / 2}, {a / 2, 0, a / 2}, {a / 2, a / 2, 0}}};
        const auto lj = makeLennardJones(epsilon, sigma, rc);
        const auto constants = forceConstants(unit, *lj, supercellFor(unit, 13.0), 0.002);
        // Independent lattice sum: D(q) = (1/m) sum_R k(R) (1 - cos q.R).
        auto analytic = [&](const Vec3& q) {
            Mat3 d{};
            for (int i = -4; i <= 4; ++i)
                for (int j = -4; j <= 4; ++j)
                    for (int k = -4; k <= 4; ++k) {
                        const Vec3 r = rowTimes({double(i), double(j), double(k)}, unit.cell);
                        const double length = norm(r);
                        if (length < 1e-9 || length >= rc) continue;
                        const double s6 = std::pow(sigma / length, 6);
                        const double d1 = 4 * epsilon * (-12 * s6 * s6 + 6 * s6) / length;
                        const double d2 = 4 * epsilon * (156 * s6 * s6 - 42 * s6) / (length * length);
                        const double factor = 1 - std::cos(dot(q, r));
                        for (int x = 0; x < 3; ++x)
                            for (int y = 0; y < 3; ++y)
                                d[x][y] += factor * ((d2 - d1 / length) * r[x] * r[y] / (length * length) + (x == y ? d1 / length : 0)) / mass;
                    }
            Eigen::Matrix3d m;
            for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y) m(x, y) = d[x][y];
            const Eigen::Vector3d lambda = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(m).eigenvalues();
            std::vector<double> f;
            const double acceleration = 1.602176634e-19 / (1e-10 * 1.66053906660e-27) * 1e10 * 1e-30;
            for (int x = 0; x < 3; ++x) f.push_back(std::sqrt(std::max(0.0, lambda(x)) * acceleration) / (2 * kPi) * 1000);
            return f;
        };
        const double g = 2 * kPi / a;
        for (const Vec3& q : {Vec3{g, 0, 0}, Vec3{g / 2, g / 2, g / 2}, Vec3{.37 * g, .21 * g, .05 * g}, Vec3{.8 * g, .4 * g, 0}}) {
            const auto numeric = phononFrequencies(constants, q);
            const auto exact = analytic(q);
            for (int m = 0; m < 3; ++m) close(numeric[static_cast<std::size_t>(m)], exact[static_cast<std::size_t>(m)], 2e-3, "LJ phonon frequency");
        }
        for (double f : phononFrequencies(constants, {0, 0, 0})) close(f, 0, 1e-6, "acoustic sum rule at Gamma");

        const Json copper = runTool("phonons", object({{"structure", structureJsonOf(fccCubic(3.59, 1))}, {"calculator", object({{"potential", "EMT"}})},
            {"mesh", Json::array({6, 6, 6})}})).result;
        check(copper.at("unit_cell_atoms").number() == 1, "primitive cell used for the symmetry path");
        check(copper.at("imaginary_mesh_modes").number() == 0, "stable fcc copper");
        const double top = copper.at("max_frequency_THz").number();
        check(top > 6 && top < 9.5, "EMT Cu maximum frequency near the measured 7.2 THz (got " + std::to_string(top) + ")");
        close(copper.at("dos").at("enclosed_modes").number(), 3, 0.02, "DOS holds three modes per atom");
        const auto cv = values(copper.at("thermodynamics").at("heat_capacity_eV_per_K"));
        close(cv.back() / (3 * 8.617333262145e-5), 1, 0.02, "Dulong-Petit limit at 1000 K");
        const auto plots = resultPlots("phonons", copper);
        check(plots.size() == 4 && plots[0].series.size() == 3 && !plots[0].markers.empty(), "dispersion with labelled points, DOS and thermodynamics plots");
    });

    test("NEB: automatic images, endpoint relaxation and restart", [] {
        auto well = [] { return std::make_unique<CurvedDoubleWell>(); };
        Configuration initial;
        initial.symbols = {"H"}; initial.numbers = {1}; initial.masses = {1.008};
        initial.positions = {{-1, 0, 0}};
        Configuration final = initial;
        final.positions = {{1, 0, 0}};
        NebOptions automatic;
        automatic.images = 0; automatic.imageSpacing = 0.1; automatic.steps = 1;
        check(migrationPath(initial, final, well, automatic).result.at("image_count").number() == 21, "ceil(2/0.1)+1 images");
        automatic.imageSpacing = 1.0;
        check(migrationPath(initial, final, well, automatic).result.at("image_count").number() == 5, "at least five images");

        // Endpoints displaced from the minima are relaxed back before the band forms.
        Configuration rough = initial, roughFinal = final;
        rough.positions = {{-0.85, 0.1, 0.05}};
        roughFinal.positions = {{1.1, -0.08, 0}};
        NebOptions relaxed;
        relaxed.fmax = .002; relaxed.relaxEndpoints = true; relaxed.endpointFmax = 1e-5;
        const auto withRelax = migrationPath(rough, roughFinal, well, relaxed);
        check(withRelax.result.at("endpoint_relaxation").at("initial_converged").boolean(), "endpoint relaxed");
        close(withRelax.frames.front().atoms[0].x, -1, 1e-4, "initial endpoint at its minimum");
        close(number(withRelax.result.at("forward_barrier_eV")), 1, 1e-4, "barrier from relaxed endpoints");

        // A short unconverged run, saved and restarted, continues to the same answer.
        NebOptions partial;
        partial.fmax = .002; partial.steps = 8;
        const auto first = migrationPath(initial, final, well, partial);
        check(!first.result.at("converged").boolean(), "first segment unconverged");
        const auto saved = scratch("neb_images.extxyz");
        writeExtxyz(saved, first.frames, {}, {});
        NebOptions restart;
        restart.fmax = .002; restart.steps = 300;
        for (const auto& frame : readFrames(saved)) restart.restart.push_back(configurationFrom(frame.structure, frame.pbc));
        const auto resumed = migrationPath(initial, final, well, restart);
        check(resumed.result.at("restarted").boolean() && resumed.result.at("converged").boolean(), "restart converged");
        close(number(resumed.result.at("forward_barrier_eV")), 1, 1e-5, "restart barrier");
        close(resumed.frames[3].atoms[0].y, .2, 1e-3, "restart saddle on the curved path");

        // Request-level restart of a copper vacancy hop without endpoint structures.
        Crystal copper = fccCubic(3.61, 2);
        copper.positions.erase(copper.positions.begin());
        copper.symbols.erase(copper.symbols.begin());
        Crystal hopped = copper;
        hopped.positions[0] = {0, 0, 0};  // a neighbour moves into the vacancy
        const Json firstRun = runTool("neb", object({{"initial", structureJsonOf(copper)}, {"final", structureJsonOf(hopped)},
            {"calculator_factory", object({{"potential", "EMT"}})}, {"images", 5}, {"steps", 3}})).result;
        std::vector<Structure> images;
        for (const auto& image : firstRun.at("images").items()) {
            const auto path = scratch("image.json");
            { std::ofstream out(path); out << image.dump(); }
            images.push_back(readFrames(path).front().structure);
        }
        const auto restartFile = scratch("cu_hop_images.extxyz");
        writeExtxyz(restartFile, images, {}, {});
        const Json resumedRun = runTool("neb", object({{"restart_images", object({{"file", restartFile.u8string()}})},
            {"calculator_factory", object({{"potential", "EMT"}})}, {"steps", 2}})).result;
        check(resumedRun.at("restarted").boolean() && resumedRun.at("image_count").number() == 5, "request restart keeps the image count");
        close(resumedRun.at("energies_eV").items()[2].number(), firstRun.at("energies_eV").items()[2].number(), 0.5, "restart starts from the saved band");
        expectError([] { runTool("neb", object({{"calculator_factory", object({{"potential", "EMT"}})}})); }, "NEB needs endpoints or a restart");
    });

    test("MD: NVE production, energy conservation and block statistics", [] {
        auto stats = blockStatistics(std::vector<double>(50, 2.5));
        close(stats.mean, 2.5, 0, "constant mean");
        close(stats.standardDeviation, 0, 0, "constant spread");
        close(stats.standardError, 0, 0, "constant standard error");
        std::vector<double> alternating;
        for (int i = 0; i < 100; ++i) alternating.push_back(i % 2 ? 1.0 : -1.0);
        alternating.push_back(std::nan(""));
        stats = blockStatistics(alternating);
        check(stats.samples == 100, "NaN samples ignored");
        close(stats.mean, 0, 1e-15, "alternating mean");
        close(stats.standardError, 0, 1e-15, "fast fluctuations average out within blocks");
        std::vector<double> ramp;
        for (int i = 0; i < 100; ++i) ramp.push_back(i);
        stats = blockStatistics(ramp);
        close(stats.standardError, 20 * std::sqrt(2.5) / std::sqrt(5.0), 1e-9, "block-mean standard error of a ramp");
        check(std::isnan(blockStatistics({1, 2, 3}).standardError), "too few samples for blocks");

        const auto emt = makeEmt();
        DynamicsOptions options;
        options.steps = 200; options.timestepFs = 1; options.temperatureK = 300; options.thermostatFs = 20;
        options.sampleInterval = 5; options.productionSteps = 500; options.seed = 4;
        const auto run = nvtDynamics(configuration(fccCubic(3.59, 3)), *emt, options);
        const Json& result = run.result;
        check(result.at("production_start_fs").number() == 200, "production follows the thermostatted run");
        check(run.frames.size() == 41 + 100, "frames from both stages");
        close(values(result.at("time_fs")).back(), 700, 0, "production time");
        const double drift = result.at("nve_energy_drift_eV_per_atom_per_ps").number();
        check(std::abs(drift) < 2e-3, "NVE energy conserved (drift " + std::to_string(drift) + " eV/atom/ps)");
        const auto& energy = result.at("statistics").at("total_energy_eV");
        check(energy.at("samples").number() == 101, "statistics use the production samples");
        check(energy.at("standard_deviation").number() / 108 < 1e-3, "NVE total energy fluctuations are small");
        const auto& temperature = result.at("statistics").at("temperature_K");
        check(temperature.at("mean").number() > 100 && temperature.at("mean").number() < 400, "production temperature near the thermostat target");
        check(temperature.at("standard_error").isNumber(), "temperature standard error");
        const auto plots = resultPlots("nvt", result);
        check(!plots[0].markers.empty() && plots[0].markers[0].label == "NVE", "production start marked on plots");

        options.productionSteps = 0;
        options.steps = 100;
        options.equilibrationFs = 50;
        const auto windowed = nvtDynamics(configuration(fccCubic(3.59, 2)), *emt, options);
        check(windowed.result.at("statistics").at("temperature_K").at("samples").number() == 11, "equilibration window excludes early samples");
        check(!windowed.result.contains("production_start_fs"), "no production stage when not requested");

        DynamicsOptions npt;
        npt.steps = 60; npt.timestepFs = 1; npt.temperatureK = 300; npt.pressureGPa = 0; npt.thermostatFs = 20;
        npt.barostatFs = 200; npt.sampleInterval = 5; npt.productionSteps = 40;
        const auto pressured = nptDynamics(configuration(fccCubic(3.59, 2)), *emt, npt);
        check(pressured.result.at("statistics").contains("volume_A3") && pressured.result.at("statistics").contains("pressure_GPa"), "NPT statistics");
        const auto volumes = values(pressured.result.at("volume_A3"));
        close(volumes.back(), volumes[12], 1e-9, "production keeps the final cell fixed");
    });

    test("Ackland-Jones structure types", [] {
        auto classify = [](const Crystal& crystal, bool periodic, double cutoff) {
            Json request = object({{"positions", rows(crystal.positions)}, {"cutoff_A", cutoff}});
            if (periodic) { request["cell"] = matrixJson(crystal.cell); request["pbc"] = Json::array({true, true, true}); }
            return run("structure-type", request);
        };
        auto fraction = [](const Json& result, const char* type) { return result.at("fractions").at(type).number(); };
        close(fraction(classify(fccCubic(3.61, 3), true, 0.0), "fcc"), 1, 0, "perfect fcc");
        Crystal bcc;
        const double a = 2.87;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k)
            for (const Vec3& b : {Vec3{0, 0, 0}, Vec3{.5, .5, .5}}) { bcc.symbols.push_back("Fe"); bcc.positions.push_back({(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a}); }
        bcc.cell = {{{4 * a, 0, 0}, {0, 4 * a, 0}, {0, 0, 4 * a}}};
        close(fraction(classify(bcc, true, 0.0), "bcc"), 1, 0, "perfect bcc");
        Crystal hcp;
        const double ah = 3.21, ch = ah * std::sqrt(8.0 / 3);
        const Mat3 hcell = {{{ah, 0, 0}, {-ah / 2, ah * std::sqrt(3.0) / 2, 0}, {0, 0, ch}}};
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 3; ++k)
            for (const Vec3& b : {Vec3{1.0 / 3, 2.0 / 3, .25}, Vec3{2.0 / 3, 1.0 / 3, .75}}) {
                hcp.symbols.push_back("Mg");
                hcp.positions.push_back(rowTimes({i + b[0], j + b[1], k + b[2]}, hcell));
            }
        hcp.cell = {{scale(hcell[0], 4), scale(hcell[1], 4), scale(hcell[2], 3)}};
        close(fraction(classify(hcp, true, 0.0), "hcp"), 1, 0, "ideal hcp");
        // Icosahedral 13-atom cluster: the centre is icosahedral, surface atoms are not.
        Crystal ico;
        const double phi = (1 + std::sqrt(5.0)) / 2, radius = 2.5;
        ico.symbols.push_back("Cu");
        ico.positions.push_back({0, 0, 0});
        for (int s1 : {-1, 1}) for (int s2 : {-1, 1})
            for (const Vec3& v : {Vec3{0, double(s1), s2 * phi}, Vec3{double(s1), s2 * phi, 0}, Vec3{s2 * phi, 0, double(s1)}}) {
                ico.symbols.push_back("Cu");
                ico.positions.push_back(scale(v, radius / std::sqrt(1 + phi * phi)));
            }
        const Json icoResult = classify(ico, false, 3.5);
        check(icoResult.at("structure_type").items()[0].number() == 4, "icosahedral centre");
        check(icoResult.at("counts").at("icosahedral").number() == 1, "only the centre is icosahedral");
        // Thermal noise (sigma 0.07 A, ~2.7% of the bond) keeps fcc identifiable.
        Crystal noisy = fccCubic(3.61, 3);
        Normal normal(31);
        for (auto& p : noisy.positions) p = add(p, {normal() * .07, normal() * .07, normal() * .07});
        check(fraction(classify(noisy, true, 0.0), "fcc") > 0.95, "noisy fcc stays fcc");
        Crystal gas;
        for (int i = 0; i < 200; ++i) {
            gas.symbols.push_back("Ar");
            gas.positions.push_back({std::abs(normal()) * 8, std::abs(normal()) * 8, std::abs(normal()) * 8});
        }
        check(fraction(classify(gas, false, 4.0), "other") > 0.9, "random gas is unclassified");
        check(perAtomProperties("structure-type", icoResult).size() == 1, "per-atom structure type");
        expectError([&] { classify(ico, false, 0.0); }, "automatic cutoff needs a cell");
    });

    test("cluster analysis: separation, periodic wrap, mask and percolation", [] {
        // Two dimers and a trimer; the trimer straddles the periodic x boundary.
        const std::vector<Vec3> atoms = {{1, 1, 1}, {2, 1, 1}, {5, 5, 5}, {5, 6, 5}, {9.6, 3, 3}, {0.4, 3, 3}, {1.4, 3, 3}};
        const Mat3 box = {{{10, 0, 0}, {0, 10, 0}, {0, 0, 10}}};
        const Json periodic = run("cluster-analysis", object({{"positions", rows(atoms)}, {"cutoff_A", 1.5}, {"cell", matrixJson(box)}, {"pbc", Json::array({true, true, true})}}));
        check(periodic.at("cluster_count").number() == 3, "three clusters with periodic bonds");
        const Json& trimer = periodic.at("clusters").items()[0];
        check(trimer.at("size").number() == 3, "largest is the wrapped trimer");
        const double cx = trimer.at("centroid_A").items()[0].number();
        const double wrapped = std::fmod(cx + 100.0, 10.0);
        close(wrapped, (9.6 + 10.4 + 11.4) / 3 - 10, 1e-9, "unwrapped centroid across the boundary");
        const double mean = (9.6 + 10.4 + 11.4) / 3;
        const double gyration = std::sqrt((std::pow(9.6 - mean, 2) + std::pow(10.4 - mean, 2) + std::pow(11.4 - mean, 2)) / 3);
        close(trimer.at("radius_of_gyration_A").number(), gyration, 1e-9, "trimer radius of gyration");
        check(run("cluster-analysis", object({{"positions", rows(atoms)}, {"cutoff_A", 1.5}})).at("cluster_count").number() == 4, "open boundaries split the trimer");
        const Json masked = run("cluster-analysis", object({{"positions", rows(atoms)}, {"cutoff_A", 1.5}, {"mask", Json::array({1, 1, 0, 0, 0, 0, 0})}}));
        check(masked.at("cluster_count").number() == 1 && masked.at("cluster_id").items()[2].number() == -1, "mask selects atoms");
        const Crystal copper = fccCubic(3.61, 2);
        const Json bulk = run("cluster-analysis", object({{"positions", rows(copper.positions)}, {"cutoff_A", 2.8}, {"cell", matrixJson(copper.cell)}, {"pbc", Json::array({true, true, true})}}));
        check(bulk.at("cluster_count").number() == 1 && bulk.at("clusters").items()[0].at("percolating").boolean(), "bulk crystal percolates");
    });

    test("void analysis: vacancies, divacancy merge and spherical cavity", [] {
        const Crystal perfect = fccCubic(3.61, 3);
        auto voids = [](const Crystal& c, double grid = 0.25) {
            return run("void-analysis", object({{"positions", rows(c.positions)}, {"cell", matrixJson(c.cell)}, {"pbc", Json::array({true, true, true})},
                                               {"grid_spacing_A", grid}}));
        };
        check(voids(perfect).at("void_count").number() == 0, "no voids in a perfect crystal with a 1 A probe");
        Crystal vacancy = perfect;
        const Vec3 site = vacancy.positions[13];
        vacancy.positions.erase(vacancy.positions.begin() + 13);
        vacancy.symbols.pop_back();
        const Json one = voids(vacancy);
        check(one.at("void_count").number() == 1, "one vacancy void");
        const Vec3 centre = {one.at("voids").items()[0].at("centroid_A").items()[0].number(), one.at("voids").items()[0].at("centroid_A").items()[1].number(),
                             one.at("voids").items()[0].at("centroid_A").items()[2].number()};
        check(norm(sub(centre, site)) < 0.2, "void centred on the vacant site");
        int lining = 0;
        for (const auto& v : one.at("lining_void_id").items()) lining += v.number() == 0;
        check(lining == 12, "twelve atoms line an fcc vacancy (got " + std::to_string(lining) + ")");
        // Neighbouring vacancies merge into one void; distant ones stay separate.
        Crystal pair = perfect, apart = perfect;
        std::size_t neighbour = 0, far = 0;
        for (std::size_t i = 0; i < perfect.positions.size(); ++i) {
            const double d = norm(sub(perfect.positions[i], site));
            if (std::abs(d - 3.61 / std::sqrt(2.0)) < 1e-6 && !neighbour) neighbour = i;
            if (d > 7.0 && !far) far = i;
        }
        for (Crystal* c : {&pair, &apart}) {
            const std::size_t other = c == &pair ? neighbour : far;
            std::vector<Vec3> kept;
            for (std::size_t i = 0; i < c->positions.size(); ++i) if (i != 13 && i != other) kept.push_back(c->positions[i]);
            c->positions = kept;
            c->symbols.resize(kept.size());
        }
        // The channel between adjacent fcc vacancies has a free radius of
        // a/2 - r (as narrow as an octahedral hole), so a 1 A probe sees two
        // regions whose total equals two isolated vacancies.
        const Json divacancy = voids(pair);
        check(divacancy.at("void_count").number() == 2, "adjacent vacancies: two probe-accessible regions");
        close(divacancy.at("accessible_volume_A3").number(), 2 * one.at("accessible_volume_A3").number(), 0.2 * one.at("accessible_volume_A3").number(), "divacancy volume");
        const Json smallProbe = run("void-analysis", object({{"positions", rows(pair.positions)}, {"cell", matrixJson(pair.cell)},
            {"pbc", Json::array({true, true, true})}, {"probe_radius_A", 0.0}, {"grid_spacing_A", 0.2}}));
        check(smallProbe.at("void_count").number() >= 1, "a point probe sees the connected empty space");
        check(voids(apart).at("void_count").number() == 2, "distant vacancies are two voids");
        // Spherical cavity of radius 5 A carved from a larger crystal.
        Crystal cavity = fccCubic(3.61, 5);
        const Vec3 middle = scale(add(add(cavity.cell[0], cavity.cell[1]), cavity.cell[2]), 0.5);
        std::vector<Vec3> kept;
        for (const auto& x : cavity.positions) if (norm(sub(x, middle)) > 5.0) kept.push_back(x);
        cavity.positions = kept;
        cavity.symbols.resize(kept.size());
        const Json pore = voids(cavity, 0.3);
        check(pore.at("void_count").number() == 1, "one cavity");
        const double radius = pore.at("voids").items()[0].at("pore_radius_estimate_A").number();
        check(radius > 4.5 && radius < 6.2, "pore radius near the carved 5 A (got " + std::to_string(radius) + ")");
        close(pore.at("void_fraction").number(), pore.at("accessible_volume_A3").number() / cellVolume(cavity.cell), 1e-12, "void fraction");
    });

    test("powder XRD and electron diffraction", [] {
        const double a = 3.615, lambda = 1.54184;
        const Json window = Json::array({10.0, 100.0});  // Cu (311) lies at 90.02 and (222) at 95.2 degrees
        const Json copper = run("powder-xrd", object({{"structure", structureJsonOf(fccCubic(a, 1))}, {"two_theta_range_deg", window}}));
        const auto& peaks = copper.at("peaks").items();
        const int expected[5][3] = {{1, 1, 1}, {2, 0, 0}, {2, 2, 0}, {3, 1, 1}, {2, 2, 2}};
        const int multiplicity[5] = {8, 6, 12, 24, 8};
        check(peaks.size() == 5, "five fcc reflections below 100 degrees (got " + std::to_string(peaks.size()) + ")");
        for (int i = 0; i < 5; ++i) {
            const double q = std::sqrt(double(expected[i][0] * expected[i][0] + expected[i][1] * expected[i][1] + expected[i][2] * expected[i][2]));
            close(peaks[static_cast<std::size_t>(i)].at("two_theta_deg").number(), 2 * std::asin(lambda * q / (2 * a)) * 180 / kPi, 1e-6, "Bragg angle");
            check(peaks[static_cast<std::size_t>(i)].at("multiplicity").number() == multiplicity[i], "multiplicity");
            std::array<int, 3> indices{};
            for (int k = 0; k < 3; ++k) indices[static_cast<std::size_t>(k)] = std::abs(static_cast<int>(peaks[static_cast<std::size_t>(i)].at("hkl").items()[static_cast<std::size_t>(k)].number()));
            std::sort(indices.rbegin(), indices.rend());
            check(indices[0] == expected[i][0] && indices[1] == expected[i][1] && indices[2] == expected[i][2], "reflection family");
        }
        close(peaks[0].at("intensity").number(), 100, 1e-9, "(111) strongest");
        const double i200 = peaks[1].at("intensity").number(), i220 = peaks[2].at("intensity").number(), i311 = peaks[3].at("intensity").number();
        check(i200 > 35 && i200 < 60 && i220 > 12 && i220 < 35 && i311 > 12 && i311 < 35, "Cu relative intensities (200 " + std::to_string(i200) + ")");
        // Rock salt: (111) is the weak difference reflection, (200) the strongest.
        Crystal salt;
        salt.cell = {{{5.64, 0, 0}, {0, 5.64, 0}, {0, 0, 5.64}}};
        const Vec3 fcc[4] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
        for (const auto& b : fcc) {
            salt.symbols.push_back("Na"); salt.positions.push_back(rowTimes(b, salt.cell));
            salt.symbols.push_back("Cl"); salt.positions.push_back(rowTimes(add(b, {.5, 0, 0}), salt.cell));
        }
        const Json saltResult = run("powder-xrd", object({{"structure", structureJsonOf(salt)}}));
        const auto& saltPeaks = saltResult.at("peaks").items();
        check(saltPeaks[0].at("intensity").number() < 20 && saltPeaks[1].at("intensity").number() == 100, "NaCl (111) weak, (200) strongest");
        // bcc tungsten: (110) first and strongest, no (100).
        Crystal tungsten;
        tungsten.cell = {{{3.165, 0, 0}, {0, 3.165, 0}, {0, 0, 3.165}}};
        tungsten.symbols = {"W", "W"};
        tungsten.positions = {{0, 0, 0}, {1.5825, 1.5825, 1.5825}};
        const Json wResult = run("powder-xrd", object({{"structure", structureJsonOf(tungsten)}}));
        const auto& wPeaks = wResult.at("peaks").items();
        close(wPeaks[0].at("two_theta_deg").number(), 2 * std::asin(lambda * std::sqrt(2.0) / (2 * 3.165)) * 180 / kPi, 1e-6, "bcc (110) first");
        close(wPeaks[0].at("intensity").number(), 100, 1e-9, "bcc (110) strongest");
        const Json warmResult = run("powder-xrd", object({{"structure", structureJsonOf(fccCubic(a, 1))}, {"two_theta_range_deg", window}, {"debye_waller_A2", 1.0}}));
        const auto& warm = warmResult.at("peaks").items();
        check(warm[4].at("intensity").number() < peaks[4].at("intensity").number(), "Debye-Waller damps high angles");

        const Json zone = run("electron-diffraction", object({{"structure", structureJsonOf(fccCubic(a, 1))}, {"zone_axis", Json::array({0, 0, 1})}, {"g_max_inv_A", 0.8}}));
        close(zone.at("wavelength_A").number(), 0.025079, 2e-6, "200 kV relativistic wavelength");
        int first = 0;
        for (const auto& spot : zone.at("spots").items()) {
            const auto& h = spot.at("hkl").items();
            check(h[2].number() == 0, "zero-order Laue zone");
            const int hh = static_cast<int>(h[0].number()), kk = static_cast<int>(h[1].number());
            check(hh % 2 == 0 && kk % 2 == 0, "fcc [001]: only all-even reflections");
            if (std::abs(spot.at("g_inv_A").number() - 2 / a) < 1e-9) ++first;
        }
        check(first == 4, "four {200} spots around [001]");
        const Json zone110 = run("electron-diffraction", object({{"structure", structureJsonOf(fccCubic(a, 1))}, {"zone_axis", Json::array({1, -1, 0})}, {"g_max_inv_A", 0.6}}));
        int odd = 0;
        for (const auto& spot : zone110.at("spots").items()) odd += std::abs(static_cast<int>(spot.at("hkl").items()[0].number())) % 2;
        check(odd == 4, "four {111} spots in the [1-10] zone");
        check(resultPlots("powder-xrd", copper).at(0).markers.size() == 5, "pattern plot labels the peaks");
    });

    test("VASP band structure, DOS and fat bands", [] {
        const auto poscar = scratch("vasp_POSCAR"), kpoints = scratch("vasp_KPOINTS"), eigenval = scratch("vasp_EIGENVAL");
        const auto doscar = scratch("vasp_DOSCAR"), procar = scratch("vasp_PROCAR");
        { std::ofstream out(poscar); out << "CuO test\n1.0\n4 0 0\n0 4 0\n0 0 4\nCu O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n"; }
        { std::ofstream out(kpoints); out << "path\n5\nLine-mode\nReciprocal\n0 0 0 ! \\Gamma\n0.5 0 0 ! X\n\n0.5 0 0 ! X\n0.5 0.5 0 ! M\n\n0.5 0.5 0.5 ! R\n0 0 0 ! \\Gamma\n"; }
        // Segment end points (fractional); five points per segment.
        const Vec3 ends[3][2] = {{{0, 0, 0}, {.5, 0, 0}}, {{.5, 0, 0}, {.5, .5, 0}}, {{.5, .5, .5}, {0, 0, 0}}};
        std::vector<Vec3> ks;
        for (const auto& seg : ends)
            for (int i = 0; i < 5; ++i) ks.push_back(add(seg[0], scale(sub(seg[1], seg[0]), i / 4.0)));
        {
            std::ofstream out(eigenval);
            out << "    2    2    1    1\n  0.1E+02\n  1E-4\n  CAR\n test\n     9    15     3\n";
            for (const auto& k : ks)
                out << "\n  " << k[0] << ' ' << k[1] << ' ' << k[2] << "  0.0666\n    1  " << -1 - k[0] << "  1.0\n    2  " << 1 + k[1] << "  0.0\n    3  2.0  0.0\n";
        }
        {
            std::ofstream out(doscar);
            out << "    2    2    1    0\n  x\n  x\n  CAR\n test\n  2.0 -2.0 5 0.5 1.0\n";
            for (int i = 0; i < 5; ++i) out << -2 + i << ' ' << 1.0 + i << ' ' << i << '\n';
            for (int a = 0; a < 2; ++a) {
                out << "  2.0 -2.0 5 0.5 1.0\n";
                for (int i = 0; i < 5; ++i) out << -2 + i << ' ' << 0.1 * (a + 1) << ' ' << 0.2 * (a + 1) << ' ' << 0.3 * (a + 1) << '\n';
            }
        }
        {
            std::ofstream out(procar);
            out << "PROCAR lm decomposed\n# of k-points:   15         # of bands:    3         # of ions:     2\n";
            for (std::size_t k = 0; k < ks.size(); ++k) {
                out << "\n k-point    " << k + 1 << " :    " << ks[k][0] << ' ' << ks[k][1] << ' ' << ks[k][2] << "     weight = 0.0666\n";
                for (int b = 0; b < 3; ++b) {
                    out << "\nband     " << b + 1 << " # energy   0.0 # occ.  1.0\n\nion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot\n";
                    out << "    1  0.000  0.000  0.000  0.000  0.600  0.000  0.000  0.000  0.000  0.600\n";
                    out << "    2  0.000  0.200  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.200\n";
                    out << "tot    0.000  0.200  0.000  0.000  0.600  0.000  0.000  0.000  0.000  0.800\n";
                }
            }
        }
        auto request = [&](const char* orbital) {
            return object({{"eigenval_file", eigenval.u8string()}, {"kpoints_file", kpoints.u8string()}, {"structure", object({{"file", poscar.u8string()}})},
                           {"doscar_file", doscar.u8string()}, {"procar_file", procar.u8string()}, {"projection_orbital", orbital}});
        };
        const Json r = run("vasp-electronic", request("all"));
        check(r.at("fermi_source").string() == "DOSCAR", "Fermi level from DOSCAR");
        const auto distance = values(r.at("distance"));
        const double unit = 2 * kPi / 4;
        close(distance[4], 0.5 * unit, 1e-9, "Gamma-X length");
        close(distance[9], unit, 1e-9, "X-M length");
        close(distance[10], distance[9], 1e-12, "path break M|R has no length");
        close(distance.back(), unit + 0.5 * std::sqrt(3.0) * unit, 1e-9, "R-Gamma length");
        const auto& labels = r.at("labels").items();
        check(labels.size() == 4 && labels[0].at("label").string() == "G" && labels[2].at("label").string() == "M|R", "labels with merged break: " + r.at("labels").dump());
        close(r.at("energies_minus_fermi_eV").items()[0].items()[0].items()[0].number(), -1.5, 1e-12, "energies relative to E_F");
        close(r.at("gap").at("gap_eV").number(), 2.0, 1e-9, "band gap from the sampled bands");
        const auto cu = values(r.at("dos").at("projected_by_element").at("Cu").items()[0]);
        close(cu[0], 0.6, 1e-12, "Cu projected DOS sums its orbitals");
        close(r.at("fat_band_weights").at("Cu").items()[0].items()[3].items()[1].number(), 0.75, 1e-12, "Cu fat-band weight");
        close(run("vasp-electronic", request("d")).at("fat_band_weights").at("Cu").items()[0].items()[0].items()[0].number(), 0.75, 1e-12, "Cu d weight");
        close(run("vasp-electronic", request("s")).at("fat_band_weights").at("Cu").items()[0].items()[0].items()[0].number(), 0, 1e-12, "Cu s weight");
        const auto plots = resultPlots("vasp-electronic", r);
        check(plots.size() == 3 && plots[0].markers.size() == 4 && plots[0].series.size() == 3, "band, fat-band and DOS plots");
    });

#ifdef ATOMS_ENABLE_SPGLIB
    test("DFT band-structure inputs round trip", [] {
        Crystal silicon;
        silicon.cell = {{{5.43, 0, 0}, {0, 5.43, 0}, {0, 0, 5.43}}};
        const Vec3 fcc[4] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
        for (const auto& b : fcc)
            for (const Vec3& shift : {Vec3{0, 0, 0}, Vec3{.25, .25, .25}}) {
                silicon.symbols.push_back("Si");
                silicon.positions.push_back(rowTimes(add(b, shift), silicon.cell));
            }
        const Json result = run("dft-inputs", object({{"structure", structureJsonOf(silicon)}, {"points_per_segment", 40}}));
        const auto& files = result.at("files");
        const auto poscar = scratch("generated_POSCAR"), kpoints = scratch("generated_KPOINTS");
        { std::ofstream out(poscar); out << files.at("POSCAR").string(); }
        { std::ofstream out(kpoints); out << files.at("KPOINTS").string(); }
        const auto frames = readFrames(poscar);
        check(frames.size() == 1 && frames[0].structure.atoms.size() == 2, "primitive POSCAR with two atoms");
        const auto& primitive = result.at("primitive_structure").at("cell").items();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                close(frames[0].structure.cellVectors[r][c], primitive[static_cast<std::size_t>(r)].items()[static_cast<std::size_t>(c)].number(), 1e-9, "POSCAR cell");
        int perSegment = 0;
        const auto labels = readLineModeLabels(kpoints, perSegment);
        const auto& path = result.at("path").items();
        check(perSegment == 40 && labels.size() == 2 * path.size(), "KPOINTS segments");
        for (std::size_t s = 0; s < path.size(); ++s) {
            const std::string from = path[s].items()[0].string() == "GAMMA" ? "G" : path[s].items()[0].string();
            const std::string to = path[s].items()[1].string() == "GAMMA" ? "G" : path[s].items()[1].string();
            check(labels[2 * s] == from && labels[2 * s + 1] == to, "KPOINTS labels follow the path");
        }
        // QE crystal_b: 40 steps along joined segments, 1 for jumps and for the last point.
        std::istringstream qe(files.at("qe_band_cards.in").string());
        std::string line;
        while (std::getline(qe, line) && line.rfind("K_POINTS crystal_b", 0) != 0) {}
        std::getline(qe, line);
        const int count = std::stoi(line);
        std::vector<int> weights;
        for (int i = 0; i < count; ++i) {
            std::getline(qe, line);
            std::istringstream row(line);
            double x, y, z;
            int w;
            row >> x >> y >> z >> w;
            weights.push_back(w);
        }
        int jumps = 0;
        for (std::size_t s = 1; s < path.size(); ++s) jumps += path[s].items()[0].string() != path[s - 1].items()[1].string();
        check(count == static_cast<int>(path.size()) + 1 + jumps, "one QE point per path vertex");
        check(weights.back() == 1 && std::count(weights.begin(), weights.end(), 40) == static_cast<long>(path.size()), "40 steps per segment, 1 for jumps");
        check(files.at("qe_band_cards.in").string().find("Si 28.") != std::string::npos, "QE species mass");
    });
#endif

    test("batch runs: sweeps, file globs, grids, failures and CSV", [] {
        const Json pressures = Json::parse(R"({"tool": "relax", "base": {"calculator": {"potential": "EMT"}, "relax_cell": true, "fmax": 0.001},
            "sweep": {"pressure_GPa": [0, 5, 10]}, "collect": ["volume_A3", "pressure_GPa", "converged"]})");
        Json batch = pressures;
        batch["base"]["structure"] = structureJsonOf(fccCubic(3.6, 2));
        const Json result = runBatch(batch, ".");
        check(result.at("run_count").number() == 3 && result.at("failures").number() == 0, "three relaxations");
        const auto& rows = result.at("table").at("rows").items();
        for (std::size_t i = 0; i < 3; ++i) {
            close(rows[i].items()[2].number(), 5.0 * i, 0.05, "reached each target pressure");
            check(rows[i].items()[3].boolean(), "converged");
            if (i) check(rows[i].items()[1].number() < rows[i - 1].items()[1].number(), "volume falls with pressure");
        }
        // One run per matching file, in sorted order.
        const auto folder = scratch("batch_gaps");
        std::filesystem::create_directories(folder);
        const double gaps[3] = {1.0, 2.5, 0.4};
        for (int i = 0; i < 3; ++i) {
            std::ofstream out(folder / ("gap_" + std::to_string(i) + ".json"));
            out << "[[-0.05, " << gaps[i] - 0.05 << "]]";  // edges straddle E_F = 0
        }
        { std::ofstream out(folder / "other.json"); out << "[[0]]"; }
        const Json files = runBatch(Json::parse(R"({"tool": "band-gap", "base": {"fermi_eV": 0},
            "files": {"energies_eV": "gap_*.json"}, "collect": ["gap_eV"]})"), folder);
        check(files.at("run_count").number() == 3, "glob matches three files");
        for (int i = 0; i < 3; ++i) {
            check(files.at("table").at("rows").items()[static_cast<std::size_t>(i)].items()[0].string() == "gap_" + std::to_string(i) + ".json", "sorted file order");
            close(files.at("table").at("rows").items()[static_cast<std::size_t>(i)].items()[1].number(), gaps[i], 1e-12, "collected gap");
        }
        // Cartesian product and recorded failures.
        const Json grid = runBatch(Json::parse(R"({"tool": "harmonic-thermodynamics", "base": {"energies_eV": [[0.01, 0.02]]},
            "sweep": {"temperatures_K": [[100], [300]], "zero_tolerance_eV": [1e-8, -1, 1e-6]}, "collect": ["heat_capacity_eV_per_K.0"]})"), ".");
        check(grid.at("run_count").number() == 6 && grid.at("failures").number() == 2, "2 x 3 grid with two invalid tolerances");
        check(!grid.at("table").at("rows").items()[1].items()[3].isNull(), "failure message recorded");
        expectError([] { runBatch(Json::parse(R"({"tool": "band-gap", "base": {"fermi_eV": 0}, "sweep": {"energies_eV": [[[-1, 1]], [[0]]]},
            "continue_on_error": false})"), "."); }, "stop on first failure when requested");
        check(valueAt(Json::parse(R"({"a": {"b": [1, {"c": 7}]}})"), "a.b.1.c").number() == 7, "dotted result path");
        check(valueAt(Json::parse(R"({"a": 1})"), "a.missing").isNull(), "missing path is null");
        const std::string csv = batchCsv(files);
        check(csv.rfind("energies_eV,gap_eV,error\n", 0) == 0 && std::count(csv.begin(), csv.end(), '\n') == 4, "CSV table");
    });

    test("LAMMPS export: triclinic box, coordinates, potentials and tasks", [] {
        // Parse a data file back: box rows and atom positions.
        auto parse = [](const std::string& text, Mat3& box, std::vector<Vec3>& positions, std::vector<int>& types) {
            std::istringstream in(text);
            std::string line;
            box = Mat3{};
            positions.clear();
            types.clear();
            bool atoms = false;
            while (std::getline(in, line)) {
                std::istringstream row(line);
                double a = 0, b = 0, c = 0;
                if (line.find("xlo xhi") != std::string::npos) { row >> a >> b; box[0][0] = b - a; }
                else if (line.find("ylo yhi") != std::string::npos) { row >> a >> b; box[1][1] = b - a; }
                else if (line.find("zlo zhi") != std::string::npos) { row >> a >> b; box[2][2] = b - a; }
                else if (line.find("xy xz yz") != std::string::npos) { row >> a >> b >> c; box[1][0] = a; box[2][0] = b; box[2][1] = c; }
                else if (line.rfind("Atoms", 0) == 0) atoms = true;
                else if (atoms) {
                    int id = 0, type = 0;
                    if (row >> id >> type >> a >> b >> c) { positions.push_back({a, b, c}); types.push_back(type); }
                }
            }
        };
        Crystal hexagonal;
        const double a = 3.21, c = 5.21;
        hexagonal.cell = {{{a, 0, 0}, {-a / 2, a * std::sqrt(3.0) / 2, 0}, {0, 0, c}}};
        hexagonal.symbols = {"Mg", "Mg", "Al"};
        hexagonal.positions = {rowTimes({1.0 / 3, 2.0 / 3, .25}, hexagonal.cell), rowTimes({2.0 / 3, 1.0 / 3, .75}, hexagonal.cell), rowTimes({.1, .2, .5}, hexagonal.cell)};
        const Json result = run("lammps-export", object({{"structure", structureJsonOf(hexagonal)}}));
        check(result.at("triclinic").boolean(), "hexagonal cell is triclinic");
        Mat3 box;
        std::vector<Vec3> positions;
        std::vector<int> types;
        parse(result.at("files").at("data.lammps").string(), box, positions, types);
        close(cellVolume(box), cellVolume(hexagonal.cell), 1e-6, "box volume");
        check(std::abs(box[1][0]) <= box[0][0] / 2 + 1e-9, "reduced xy tilt");
        check(types == std::vector<int>({1, 1, 2}), "atom types in order of first appearance");
        for (std::size_t i = 0; i < positions.size(); ++i) {
            Vec3 f0 = fractional(hexagonal.positions[i], hexagonal.cell), f1 = fractional(positions[i], box);
            for (int k = 0; k < 3; ++k) {
                const double d = f1[k] - f0[k];
                close(d - std::round(d), 0, 1e-9, "fractional coordinates preserved");
            }
        }
        // A strongly sheared cell is reduced to legal tilts without changing its lattice.
        const Mat3 sheared = {{{4, 0, 0}, {9, 4, 0}, {7, -11, 4}}};
        const Mat3 lmp = lammpsCell(sheared);
        check(std::abs(lmp[1][0]) <= lmp[0][0] / 2 + 1e-9 && std::abs(lmp[2][0]) <= lmp[0][0] / 2 + 1e-9 && std::abs(lmp[2][1]) <= lmp[1][1] / 2 + 1e-9, "tilts reduced");
        close(cellVolume(lmp), cellVolume(sheared), 1e-9, "reduction keeps the volume");
        expectError([] { lammpsCell({{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}}}); }, "left-handed cell rejected");

        const Crystal copper = fccCubic(3.61, 2);
        const Json eam = run("lammps-export", object({{"structure", structureJsonOf(copper)}, {"calculator", object({{"potential", "EAM"}, {"file", "Cu_u3.eam"}})},
                                                      {"task", "npt"}, {"pressure_GPa", 2.0}}));
        const std::string input = eam.at("files").at("in.lammps").string();
        check(input.find("pair_style eam\npair_coeff 1 1 Cu_u3.eam") != std::string::npos, "funcfl EAM pair style");
        check(input.find("iso 20000 20000") != std::string::npos, "2 GPa = 20000 bar");
        check(input.find("dump_modify traj element Cu") != std::string::npos, "dump with element names");
        const Json lj = run("lammps-export", object({{"structure", structureJsonOf(copper)}, {"task", "nvt"},
            {"calculator", object({{"potential", "LennardJones"}, {"epsilon", 0.0104}, {"sigma", 3.4}, {"cutoff", 8.5}})}}));
        check(lj.at("files").at("in.lammps").string().find("pair_style lj/cut 8.5\npair_coeff * * 0.0104 3.4\npair_modify shift yes") != std::string::npos, "LJ pair style");
        check(run("lammps-export", object({{"structure", structureJsonOf(copper)}, {"calculator", object({{"potential", "EAM"}, {"file", "CuNi.eam.alloy"}})}}))
                  .at("files").at("in.lammps").string().find("pair_style eam/alloy\npair_coeff * * CuNi.eam.alloy Cu") != std::string::npos, "setfl pair style");
        Crystal moved = copper;
        moved.positions[0] = add(moved.positions[0], {0.3, 0.1, 0});
        const Json neb = run("lammps-export", object({{"structure", structureJsonOf(copper)}, {"final", structureJsonOf(moved)}, {"task", "neb"}}));
        std::istringstream finalFile(neb.at("files").at("final.neb").string());
        std::string header;
        std::getline(finalFile, header);
        check(std::stoi(header) == static_cast<int>(copper.positions.size()), "final.neb atom count");
        int id;
        Vec3 first;
        finalFile >> id >> first[0] >> first[1] >> first[2];
        Mat3 cubeBox;
        std::vector<Vec3> cubePositions;
        std::vector<int> cubeTypes;
        parse(neb.at("files").at("data.lammps").string(), cubeBox, cubePositions, cubeTypes);
        close(norm(sub(first, cubePositions[0])), norm(Vec3{0.3, 0.1, 0}), 1e-9, "final.neb carries the displacement");
        expectError([&] { run("lammps-export", object({{"structure", structureJsonOf(copper)}, {"task", "neb"}})); }, "NEB needs a final structure");
        const Json open = run("lammps-export", object({{"structure", Json::parse(R"({"symbols":["Ar","Ar"],"positions":[[0,0,0],[3.8,0,0]]})")}}));
        check(open.at("files").at("in.lammps").string().find("boundary f f f") != std::string::npos, "open structure uses fixed boundaries");
    });

    test("animated GIF encoding round trip", [] {
        const auto path = std::filesystem::temp_directory_path() / "atomforge_gif_test.gif";
        const int w = 97, h = 61;
        std::vector<std::vector<unsigned char>> frames;
        for (int f = 0; f < 3; ++f) {
            std::vector<unsigned char> rgba(static_cast<std::size_t>(w * h * 4));
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x) {
                    unsigned char* p = &rgba[static_cast<std::size_t>((y * w + x) * 4)];
                    // Few flat colours (exact after quantisation) plus a moving stripe.
                    const bool stripe = (x + 7 * f) % 20 < 4;
                    p[0] = stripe ? 200 : static_cast<unsigned char>(x < w / 2 ? 16 : 240);
                    p[1] = stripe ? 40 : static_cast<unsigned char>(y < h / 2 ? 80 : 160);
                    p[2] = static_cast<unsigned char>(8 * f);
                    p[3] = 255;
                }
            frames.push_back(rgba);
        }
        {
            GifWriter writer(path, w, h, 7);
            for (const auto& frame : frames) writer.addFrame(frame);
            check(writer.frames() == 3, "frame count while writing");
        }
        const GifImage image = readGif(path);
        check(image.width == w && image.height == h && image.looping, "GIF header and loop extension");
        check(image.frames.size() == 3 && image.delays == std::vector<int>({7, 7, 7}), "frames and delays decoded");
        int mismatches = 0;
        for (std::size_t f = 0; f < 3; ++f)
            for (std::size_t i = 0; i < static_cast<std::size_t>(w * h); ++i)
                for (std::size_t k = 0; k < 3; ++k)
                    mismatches += image.frames[f][3 * i + k] != frames[f][4 * i + k];
        check(mismatches == 0, "flat colours reproduced exactly");

        // A smooth gradient with far more than 256 colours exercises the LZW table
        // reset and nearest-colour mapping.
        const int gw = 300, gh = 200;
        std::vector<unsigned char> rgb(static_cast<std::size_t>(gw * gh * 3));
        for (int y = 0; y < gh; ++y)
            for (int x = 0; x < gw; ++x) {
                unsigned char* p = &rgb[static_cast<std::size_t>((y * gw + x) * 3)];
                p[0] = static_cast<unsigned char>(x * 255 / (gw - 1));
                p[1] = static_cast<unsigned char>(y * 255 / (gh - 1));
                p[2] = static_cast<unsigned char>((x * y) % 256);
            }
        {
            GifWriter writer(path, gw, gh, 4);
            writer.addFrame(rgb, 3);
        }
        const GifImage gradient = readGif(path);
        check(gradient.frames.size() == 1, "gradient frame decoded");
        double error = 0;
        for (std::size_t i = 0; i < rgb.size(); ++i) error += std::abs(int(gradient.frames[0][i]) - int(rgb[i]));
        error /= static_cast<double>(rgb.size());
        check(error < 40, "gradient mean error " + std::to_string(error));
        expectError([&] { GifWriter bad(path, 4, 4, 5); bad.addFrame(std::vector<unsigned char>(10)); }, "short frame rejected");
        std::filesystem::remove(path);
    });

    test("catalog defaults are valid requests", [] {
        check(scienceToolCatalog().size() >= 20, "catalog tools");
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
