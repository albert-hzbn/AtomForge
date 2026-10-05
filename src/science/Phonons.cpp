#include "science/Phonons.h"
#include "science/Analysis.h"
#include "science/ReciprocalPath.h"
#include "util/TaskControl.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);
// eV/(A^2 amu) -> 1/fs^2.
const double kAcceleration = 1.602176634e-19 / (1e-10 * 1.66053906660e-27) * 1e10 * 1e-30;

Mat3 reciprocal(const Mat3& cell)
{
    const Mat3 inv = inverse(cell);
    Mat3 result{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) result[i][j] = 2 * kPi * inv[j][i];
    return result;
}

double frequencyFromEigenvalue(double lambda)
{
    const double omega = std::sqrt(std::abs(lambda) * kAcceleration);  // rad/fs
    return (lambda < 0 ? -1.0 : 1.0) * omega / (2 * kPi) * 1000.0;     // THz
}
}

std::array<int, 3> supercellFor(const Configuration& unit, double minimumLength)
{
    positive(minimumLength, "supercell_min_A");
    const Mat3 inv = inverse(unit.cell);
    std::array<int, 3> result{};
    for (int k = 0; k < 3; ++k) {
        const double spacing = 1.0 / norm({inv[0][k], inv[1][k], inv[2][k]});
        result[static_cast<std::size_t>(k)] = std::max(1, static_cast<int>(std::ceil(minimumLength / spacing - 1e-9)));
    }
    return result;
}

ForceConstants forceConstants(const Configuration& unit, const Potential& potential,
                              const std::array<int, 3>& supercell, double displacement)
{
    if (!unit.fullyPeriodic() || std::abs(determinant(unit.cell)) < 1e-12) throw std::runtime_error("Phonons require a full periodic cell");
    positive(displacement, "displacement_A");
    for (int n : supercell)
        if (n < 1 || n > 50) throw std::runtime_error("Supercell repetitions must be between 1 and 50");
    ForceConstants result;
    result.unit = unit;
    result.supercell = supercell;
    Configuration& super = result.super;
    super.pbc = {true, true, true};
    for (int r = 0; r < 3; ++r) super.cell[r] = scale(unit.cell[r], supercell[static_cast<std::size_t>(r)]);
    for (int i = 0; i < supercell[0]; ++i)
        for (int j = 0; j < supercell[1]; ++j)
            for (int k = 0; k < supercell[2]; ++k) {
                const Vec3 shift = rowTimes({double(i), double(j), double(k)}, unit.cell);
                for (std::size_t a = 0; a < unit.size(); ++a) {
                    super.symbols.push_back(unit.symbols[a]);
                    super.numbers.push_back(unit.numbers[a]);
                    super.masses.push_back(unit.masses[a]);
                    super.positions.push_back(add(unit.positions[a], shift));
                }
            }
    const std::size_t n = unit.size(), total = super.size();
    result.phi.assign(n, std::vector<Mat3>(total, Mat3{}));
    for (std::size_t a = 0; a < n; ++a)
        for (int alpha = 0; alpha < 3; ++alpha) {
            taskProgress(static_cast<double>(3 * a + static_cast<std::size_t>(alpha)) / static_cast<double>(3 * n));
            Configuration plus = super, minus = super;
            plus.positions[a][alpha] += displacement;
            minus.positions[a][alpha] -= displacement;
            const auto fPlus = potential.compute(plus, false).forces, fMinus = potential.compute(minus, false).forces;
            for (std::size_t b = 0; b < total; ++b)
                for (int beta = 0; beta < 3; ++beta)
                    result.phi[a][b][alpha][beta] = -(fPlus[b][beta] - fMinus[b][beta]) / (2 * displacement);
        }
    // Acoustic sum rule: a rigid translation produces no force.
    for (std::size_t a = 0; a < n; ++a)
        for (int alpha = 0; alpha < 3; ++alpha)
            for (int beta = 0; beta < 3; ++beta) {
                double sum = 0;
                for (std::size_t b = 0; b < total; ++b) sum += result.phi[a][b][alpha][beta];
                result.phi[a][a][alpha][beta] -= sum;
            }
    return result;
}

std::vector<double> phononFrequencies(const ForceConstants& constants, const Vec3& q)
{
    const std::size_t n = constants.unit.size(), total = constants.super.size();
    const Mat3& superCell = constants.super.cell;
    Eigen::MatrixXcd dynamical = Eigen::MatrixXcd::Zero(static_cast<int>(3 * n), static_cast<int>(3 * n));
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t bs = 0; bs < total; ++bs) {
            const std::size_t b = bs % n;
            const Vec3 base = sub(constants.super.positions[bs], constants.unit.positions[a]);
            // Minimum-image vectors (shared equally between equidistant images).
            std::vector<Vec3> images;
            double best = HUGE_VAL;
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j)
                    for (int k = -1; k <= 1; ++k) {
                        const Vec3 d = add(base, rowTimes({double(i), double(j), double(k)}, superCell));
                        const double length = norm(d);
                        if (length < best - 1e-5) { best = length; images.assign(1, d); }
                        else if (std::abs(length - best) <= 1e-5) images.push_back(d);
                    }
            std::complex<double> phase = 0;
            for (const auto& d : images) phase += std::polar(1.0, dot(q, d));
            phase /= static_cast<double>(images.size());
            const double massFactor = 1.0 / std::sqrt(constants.unit.masses[a] * constants.unit.masses[b]);
            const Mat3& f = constants.phi[a][bs];
            for (int alpha = 0; alpha < 3; ++alpha)
                for (int beta = 0; beta < 3; ++beta)
                    dynamical(static_cast<int>(3 * a) + alpha, static_cast<int>(3 * b) + beta) += f[alpha][beta] * massFactor * phase;
        }
    const Eigen::MatrixXcd hermitian = (dynamical + dynamical.adjoint()) / 2.0;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hermitian, Eigen::EigenvaluesOnly);
    std::vector<double> frequencies;
    for (int i = 0; i < solver.eigenvalues().size(); ++i) frequencies.push_back(frequencyFromEigenvalue(solver.eigenvalues()(i)));
    return frequencies;
}

ToolOutput phononCalculation(const StructureInput& structure, const Potential& potential, const Parameters& p)
{
    Json result = Json::object();
    Configuration unit = configurationFrom(structure.structure, structure.pbc);
    std::vector<Vec3> pathFractional;
    std::vector<double> pathDistance;
    Json labels = Json::array();
    ToolOutput output;
    if (p.boolean("symmetry_path", true)) {
        // Phonons are computed in the standardized primitive cell so the
        // HPKOT k-points refer to the same reciprocal basis.
        const ToolOutput path = reciprocalPath(structure, p.number("spacing_inv_A", 0.05), p.number("symprec_A", 1e-5), true);
        unit = configurationFrom(path.frames.front(), {true, true, true});
        const NdArray k = toArray(path.result.at("kpoints_fractional"), "kpoints");
        for (std::size_t i = 0; i < k.shape[0]; ++i) pathFractional.push_back({k(i, 0), k(i, 1), k(i, 2)});
        pathDistance = toArray(path.result.at("distance_inv_A"), "distance").values;
        labels = path.result.at("labels");
        result["spacegroup_symbol"] = path.result.at("spacegroup_symbol");
        result["path"] = path.result.at("path");
        result["primitive_structure"] = path.result.at("primitive_structure");
    } else if (p.has("q_points_fractional")) {
        for (const auto& q : points(p.array("q_points_fractional"), "q_points_fractional")) {
            if (!pathFractional.empty()) {
                const Mat3 b = reciprocal(unit.cell);
                pathDistance.push_back(pathDistance.back() + norm(rowTimes(sub(q, pathFractional.back()), b)));
            } else pathDistance.push_back(0.0);
            pathFractional.push_back(q);
            labels.push("");
        }
    }
    std::array<int, 3> supercell{};
    if (p.has("supercell")) {
        const auto values = finiteArray(p.array("supercell"), "supercell", 1).values;
        if (values.size() != 3) throw std::runtime_error("supercell needs three repetitions");
        for (int k = 0; k < 3; ++k) supercell[static_cast<std::size_t>(k)] = static_cast<int>(integer(values[static_cast<std::size_t>(k)], "supercell"));
    } else supercell = supercellFor(unit, p.number("supercell_min_A", 12.0));
    const ForceConstants constants = forceConstants(unit, potential, supercell, p.number("displacement_A", 0.01));
    const Mat3 b = reciprocal(unit.cell);

    Json bands = Json::array();
    for (const auto& q : pathFractional) {
        taskCheckpoint();
        bands.push(toJson(phononFrequencies(constants, rowTimes(q, b))));
    }
    std::array<int, 3> mesh = {8, 8, 8};
    if (p.has("mesh")) {
        const auto values = finiteArray(p.array("mesh"), "mesh", 1).values;
        if (values.size() != 3) throw std::runtime_error("mesh needs three sizes");
        for (int k = 0; k < 3; ++k) mesh[static_cast<std::size_t>(k)] = static_cast<int>(integer(values[static_cast<std::size_t>(k)], "mesh"));
    }
    Json meshEnergies = Json::array();
    double lowest = HUGE_VAL, highest = -HUGE_VAL;
    int imaginary = 0;
    for (int i = 0; i < mesh[0]; ++i)
        for (int j = 0; j < mesh[1]; ++j)
            for (int k = 0; k < mesh[2]; ++k) {
                taskCheckpoint();
                const Vec3 q = {double(i) / mesh[0], double(j) / mesh[1], double(k) / mesh[2]};
                Json row = Json::array();
                for (double f : phononFrequencies(constants, rowTimes(q, b))) {
                    lowest = std::min(lowest, f);
                    highest = std::max(highest, f);
                    imaginary += f < -0.05;
                    row.push(f * kThzToEv);
                }
                meshEnergies.push(row);
            }
    result["supercell"] = Json::array({supercell[0], supercell[1], supercell[2]});
    result["unit_cell_atoms"] = unit.size();
    result["displacement_A"] = p.number("displacement_A", 0.01);
    if (!pathFractional.empty()) {
        Json q = Json::array();
        for (const auto& k : pathFractional) q.push(toJson(k));
        result["q_points_fractional"] = q;
        result["distance_inv_A"] = toJson(pathDistance);
        result["labels"] = labels;
        result["frequencies_THz"] = bands;
    }
    result["mesh"] = Json::array({mesh[0], mesh[1], mesh[2]});
    result["mesh_energies_eV"] = meshEnergies;
    result["min_frequency_THz"] = lowest;
    result["max_frequency_THz"] = highest;
    result["imaginary_mesh_modes"] = imaginary;

    // Density of states and thermodynamics reuse the tabulated-mode tools.
    const double sigma = positive(p.number("sigma_eV", 0.001), "sigma_eV");
    Json grid = Json::array();
    const double from = std::min(0.0, lowest * kThzToEv) - 5 * sigma, to = highest * kThzToEv + 5 * sigma;
    for (int i = 0; i <= 1000; ++i) grid.push(from + (to - from) * i / 1000.0);
    Json dosRequest = Json::object();
    dosRequest["energies_eV"] = meshEnergies;
    dosRequest["energy_grid_eV"] = grid;
    dosRequest["sigma_eV"] = sigma;
    result["dos"] = phononDos(Parameters("phonon-dos", dosRequest, "."));
    if (imaginary == 0) {
        Json thermoRequest = Json::object();
        thermoRequest["energies_eV"] = meshEnergies;
        if (p.has("temperatures_K")) thermoRequest["temperatures_K"] = toJson(p.array("temperatures_K"));
        else {
            Json temperatures = Json::array();
            for (int t = 0; t <= 1000; t += 25) temperatures.push(t);
            thermoRequest["temperatures_K"] = temperatures;
        }
        // Acoustic modes at Gamma are zero only to finite-difference accuracy.
        thermoRequest["zero_tolerance_eV"] = 2e-4;
        result["thermodynamics"] = harmonicThermodynamics(Parameters("harmonic-thermodynamics", thermoRequest, "."));
    } else {
        result["thermodynamics"] = Json();
        result["thermodynamics_note"] = "Imaginary modes on the mesh: the structure is dynamically unstable with this potential; harmonic thermodynamics is not defined.";
    }
    result["potential"] = potential.description();
    output.frames.push_back(structureFrom(unit));
    output.result = result;
    return output;
}
}
