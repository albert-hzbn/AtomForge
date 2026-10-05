#include "science/Analysis.h"
#include "util/TaskControl.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kBoltzmann = 8.617333262145e-5;  // eV/K
const double kEvPerA3ToGPa = 160.2176634;

NdArray trajectory(const NdArray& values, const std::string& name)
{
    const NdArray result = finiteArray(values, name, 3);
    if (result.shape[0] < 2 || result.shape[2] != 3)
        throw std::runtime_error(name + " must have shape (at least two frames, atoms, 3)");
    return result;
}

std::vector<double> vector1d(const NdArray& values, const std::string& name)
{
    return finiteArray(values, name, 1).values;
}

std::optional<long long> optionalInteger(const Parameters& p, const std::string& name, long long minimum)
{
    if (!p.has(name)) return std::nullopt;
    return integer(p.number(name), name, minimum);
}

std::vector<double> lagTimes(long long limit, double dt)
{
    std::vector<double> lag(static_cast<std::size_t>(limit + 1));
    for (std::size_t i = 0; i < lag.size(); ++i) lag[i] = static_cast<double>(i) * dt;
    return lag;
}

std::vector<int> origins(std::size_t count, long long limit)
{
    std::vector<int> result(static_cast<std::size_t>(limit + 1));
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = static_cast<int>(count - i);
    return result;
}

double trapezoid(const std::vector<double>& y, const std::vector<double>& x)
{
    double area = 0;
    for (std::size_t i = 1; i < y.size(); ++i) area += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    return area;
}

// NumPy-compatible numerical rank (SVD with max(M,N)*eps*sigma_max tolerance).
int matrixRank(const Eigen::MatrixXd& matrix)
{
    if (!matrix.size()) return 0;
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    const auto singular = svd.singularValues();
    if (!singular.size()) return 0;
    const double tolerance = singular(0) * static_cast<double>(std::max(matrix.rows(), matrix.cols())) * std::numeric_limits<double>::epsilon();
    int rank = 0;
    for (int i = 0; i < singular.size(); ++i) rank += singular(i) > tolerance;
    return rank;
}

Eigen::MatrixXd leastSquares(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
{
    return a.completeOrthogonalDecomposition().solve(b);
}

Json matrixJson(const Eigen::MatrixXd& matrix)
{
    Json result = Json::array();
    for (int i = 0; i < matrix.rows(); ++i) {
        Json row = Json::array();
        for (int j = 0; j < matrix.cols(); ++j) row.push(matrix(i, j));
        result.push(row);
    }
    return result;
}

Json vectorJson(const Eigen::VectorXd& vector)
{
    Json result = Json::array();
    for (int i = 0; i < vector.size(); ++i) result.push(vector(i));
    return result;
}

struct Linear { double slope, intercept, r, standardError; };

Linear linearRegression(const std::vector<double>& x, const std::vector<double>& y)
{
    const double n = static_cast<double>(x.size());
    const double xm = std::accumulate(x.begin(), x.end(), 0.0) / n;
    const double ym = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double ssxm = 0, ssym = 0, ssxym = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        ssxm += (x[i] - xm) * (x[i] - xm);
        ssym += (y[i] - ym) * (y[i] - ym);
        ssxym += (x[i] - xm) * (y[i] - ym);
    }
    ssxm /= n; ssym /= n; ssxym /= n;
    double r = (ssxm == 0 || ssym == 0) ? 0.0 : ssxym / std::sqrt(ssxm * ssym);
    r = std::clamp(r, -1.0, 1.0);
    const double slope = ssxym / ssxm;
    const double df = n - 2;
    const double standardError = std::sqrt(std::max(0.0, (1 - r * r) * ssym / ssxm / df));
    return {slope, ym - slope * xm, r, standardError};
}

std::vector<Neighbor> localNeighbors(const std::vector<Vec3>& positions, double cutoff, const Parameters& p,
                                     const std::string& cellName, const Pbc& pbc, Mat3& cell)
{
    cell = p.has(cellName) ? cellMatrix(p.array(cellName)) : identity();
    if ((pbc[0] || pbc[1] || pbc[2]) && !p.has(cellName)) throw std::runtime_error("Periodic analysis requires a cell");
    auto neighbors = neighborList(positions, cell, pbc, cutoff);
    for (const auto& neighbor : neighbors)
        if (norm(neighbor.vector) < 1e-10) throw std::runtime_error("Coincident atoms make local environment analysis undefined");
    return neighbors;
}

std::vector<std::vector<std::size_t>> byCenter(const std::vector<Neighbor>& neighbors, std::size_t count)
{
    std::vector<std::vector<std::size_t>> result(count);
    for (std::size_t n = 0; n < neighbors.size(); ++n) result[static_cast<std::size_t>(neighbors[n].i)].push_back(n);
    return result;
}

// Orthonormal spherical harmonic magnitude basis: Ylm(theta, 0) with Condon-Shortley phase.
double normalizedLegendre(int l, int m, double x)
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

Json meanSquareDisplacement(const Parameters& p)
{
    NdArray coordinates = trajectory(p.array("positions"), "positions");
    const bool wrapped = p.boolean("wrapped", false);
    const bool removeDrift = p.boolean("remove_drift", false);
    const double dt = positive(p.number("timestep_fs"), "timestep_fs");
    const std::size_t count = coordinates.shape[0], atoms = coordinates.shape[1];
    if (wrapped) {
        if (!p.has("cell")) throw std::runtime_error("cell must be a nonempty finite array of the required dimension");
        const MinimumImage mic(cellMatrix(p.array("cell")), p.pbc("pbc", {true, true, true}));
        const NdArray raw = coordinates;
        for (std::size_t f = 1; f < count; ++f)
            for (std::size_t a = 0; a < atoms; ++a) {
                const Vec3 step = mic({raw(f, a, 0) - raw(f - 1, a, 0), raw(f, a, 1) - raw(f - 1, a, 1), raw(f, a, 2) - raw(f - 1, a, 2)});
                for (std::size_t k = 0; k < 3; ++k) coordinates(f, a, k) = coordinates(f - 1, a, k) + step[k];
            }
    }
    if (removeDrift)
        for (std::size_t f = 0; f < count; ++f)
            for (std::size_t k = 0; k < 3; ++k) {
                double mean = 0;
                for (std::size_t a = 0; a < atoms; ++a) mean += coordinates(f, a, k);
                mean /= static_cast<double>(atoms);
                for (std::size_t a = 0; a < atoms; ++a) coordinates(f, a, k) -= mean;
            }
    const long long limit = optionalInteger(p, "max_lag", 1).value_or(static_cast<long long>(count) - 1);
    if (limit >= static_cast<long long>(count)) throw std::runtime_error("max_lag must be smaller than frame count");
    NdArray covariance({static_cast<std::size_t>(limit + 1), 3, 3});
    for (long long lag = 1; lag <= limit; ++lag) {
        taskProgress(static_cast<double>(lag) / static_cast<double>(limit));
        double sum[3][3] = {};
        const std::size_t samples = count - static_cast<std::size_t>(lag);
        for (std::size_t t = 0; t < samples; ++t)
            for (std::size_t a = 0; a < atoms; ++a) {
                double d[3];
                for (std::size_t k = 0; k < 3; ++k) d[k] = coordinates(t + static_cast<std::size_t>(lag), a, k) - coordinates(t, a, k);
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) sum[i][j] += d[i] * d[j];
            }
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                covariance(static_cast<std::size_t>(lag), i, j) = sum[i][j] / static_cast<double>(samples * atoms);
    }
    std::vector<double> msd(static_cast<std::size_t>(limit + 1));
    NdArray components({static_cast<std::size_t>(limit + 1), 3});
    for (std::size_t lag = 0; lag < msd.size(); ++lag) {
        for (std::size_t k = 0; k < 3; ++k) components(lag, k) = covariance(lag, k, k);
        msd[lag] = components(lag, 0) + components(lag, 1) + components(lag, 2);
    }
    Json result = Json::object();
    result["lag_fs"] = toJson(lagTimes(limit, dt));
    result["msd_A2"] = toJson(msd);
    result["components_A2"] = toJson(components);
    result["displacement_tensor_A2"] = toJson(covariance);
    result["origins"] = toJson(origins(count, limit));
    return result;
}

Json diffusionCoefficient(const Parameters& p)
{
    const auto times = vector1d(p.array("lag_fs"), "lag_fs");
    const auto values = vector1d(p.array("msd_A2"), "msd_A2");
    const double dimensionValue = p.number("dimensions", 3);
    if (dimensionValue != 1 && dimensionValue != 2 && dimensionValue != 3) throw std::runtime_error("dimensions must be 1, 2 or 3");
    const int dimensions = static_cast<int>(dimensionValue);
    bool valid = times.size() == values.size() && times[0] >= 0;
    for (std::size_t i = 0; valid && i < times.size(); ++i) valid = values[i] >= 0 && (i == 0 || times[i] > times[i - 1]);
    if (!valid) throw std::runtime_error("MSD must be nonnegative with matching, increasing nonnegative times");
    const auto interval = vector1d(p.array("fit_range_fs"), "fit_range_fs");
    if (interval.size() != 2 || interval[0] >= interval[1] || interval[0] < times.front() || interval[1] > times.back())
        throw std::runtime_error("fit_range_fs must lie inside the sampled time interval");
    std::vector<double> x, y;
    for (std::size_t i = 0; i < times.size(); ++i)
        if (times[i] >= interval[0] && times[i] <= interval[1]) { x.push_back(times[i]); y.push_back(values[i]); }
    if (x.size() < 3) throw std::runtime_error("At least three fit points are required");
    const Linear fit = linearRegression(x, y);
    if (fit.slope < 0) throw std::runtime_error("Negative MSD slope is not a physical diffusion estimate");
    const double diffusion = fit.slope / (2 * dimensions);
    Json result = Json::object();
    result["D_A2_per_fs"] = diffusion;
    result["D_m2_per_s"] = diffusion * 1e-5;
    result["slope_A2_per_fs"] = fit.slope;
    result["intercept_A2"] = fit.intercept;
    result["r_squared"] = std::isfinite(fit.r) ? fit.r * fit.r : 0.0;
    result["ols_slope_stderr"] = fit.standardError;
    result["fit_points"] = x.size();
    result["fit_range_fs"] = toJson(interval);
    return result;
}

Json velocityAutocorrelation(const Parameters& p)
{
    NdArray velocity = trajectory(p.array("velocities"), "velocities");
    const bool removeDrift = p.boolean("remove_drift", false);
    const bool normalize = p.boolean("normalize", false);
    const double dt = positive(p.number("timestep_fs"), "timestep_fs");
    const std::size_t count = velocity.shape[0], atoms = velocity.shape[1];
    if (removeDrift)
        for (std::size_t f = 0; f < count; ++f)
            for (std::size_t k = 0; k < 3; ++k) {
                double mean = 0;
                for (std::size_t a = 0; a < atoms; ++a) mean += velocity(f, a, k);
                mean /= static_cast<double>(atoms);
                for (std::size_t a = 0; a < atoms; ++a) velocity(f, a, k) -= mean;
            }
    const long long limit = optionalInteger(p, "max_lag", 1).value_or(static_cast<long long>(count) - 1);
    if (limit >= static_cast<long long>(count)) throw std::runtime_error("max_lag must be smaller than frame count");
    // Zero-padded FFT correlation avoids circular wrap-around.
    std::size_t size = 1;
    while (size < 2 * count - 1) size <<= 1;
    std::vector<double> correlation(static_cast<std::size_t>(limit + 1), 0.0);
    std::vector<std::complex<double>> series(size);
    for (std::size_t a = 0; a < atoms; ++a) {
        taskProgress(static_cast<double>(a) / static_cast<double>(atoms));
        // Two real series per complex transform: x + i y.
        for (std::size_t k = 0; k < 3; k += 2) {
            std::fill(series.begin(), series.end(), std::complex<double>());
            for (std::size_t f = 0; f < count; ++f)
                series[f] = {velocity(f, a, k), k + 1 < 3 ? velocity(f, a, k + 1) : 0.0};
            auto spectrum = fft(series);
            std::vector<std::complex<double>> power(size);
            for (std::size_t q = 0; q < size; ++q) {
                const auto x = 0.5 * (spectrum[q] + std::conj(spectrum[(size - q) % size]));
                const auto y = std::complex<double>(0, -0.5) * (spectrum[q] - std::conj(spectrum[(size - q) % size]));
                power[q] = std::norm(x) + std::norm(y);
            }
            const auto back = fft(power, true);
            for (std::size_t lag = 0; lag < correlation.size(); ++lag) correlation[lag] += back[lag].real();
        }
    }
    for (std::size_t lag = 0; lag < correlation.size(); ++lag)
        correlation[lag] /= static_cast<double>(atoms) * static_cast<double>(count - lag);
    if (normalize) {
        if (correlation[0] <= 0) throw std::runtime_error("Cannot normalize zero velocity autocorrelation");
        const double zero = correlation[0];
        for (double& value : correlation) value /= zero;
    }
    Json result = Json::object();
    result["lag_fs"] = toJson(lagTimes(limit, dt));
    result["vacf"] = toJson(correlation);
    result["unit"] = normalize ? "dimensionless" : "A^2/fs^2";
    result["origins"] = toJson(origins(count, limit));
    return result;
}

Json vibrationalSpectrum(const Parameters& p)
{
    const bool removeMean = p.boolean("remove_mean", true);
    NdArray velocity = trajectory(p.array("velocities"), "velocities");
    const double dt = positive(p.number("timestep_fs"), "timestep_fs");
    const std::size_t count = velocity.shape[0], atoms = velocity.shape[1];
    if (count < 4) throw std::runtime_error("At least four velocity frames are required");
    if (removeMean)
        for (std::size_t a = 0; a < atoms; ++a)
            for (std::size_t k = 0; k < 3; ++k) {
                double mean = 0;
                for (std::size_t f = 0; f < count; ++f) mean += velocity(f, a, k);
                mean /= static_cast<double>(count);
                for (std::size_t f = 0; f < count; ++f) velocity(f, a, k) -= mean;
            }
    if (p.has("masses")) {
        const auto weights = vector1d(p.array("masses"), "masses");
        bool valid = weights.size() == atoms;
        for (double weight : weights) valid = valid && weight > 0;
        if (!valid) throw std::runtime_error("One positive mass is required per atom");
        for (std::size_t f = 0; f < count; ++f)
            for (std::size_t a = 0; a < atoms; ++a)
                for (std::size_t k = 0; k < 3; ++k) velocity(f, a, k) *= std::sqrt(weights[a]);
    }
    const std::string window = p.has("window") ? p.json("window").string() : "hann";
    if (window != "hann" && window != "none") throw std::runtime_error("window must be hann or none");
    std::vector<double> taper(count, 1.0);
    if (window == "hann")
        for (std::size_t f = 0; f < count; ++f) taper[f] = 0.5 - 0.5 * std::cos(2 * kPi * static_cast<double>(f) / static_cast<double>(count - 1));
    const std::size_t bins = count / 2 + 1;
    std::vector<double> power(bins, 0.0);
    std::vector<std::complex<double>> series(count);
    for (std::size_t a = 0; a < atoms; ++a) {
        taskProgress(static_cast<double>(a) / static_cast<double>(atoms));
        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t f = 0; f < count; ++f) series[f] = velocity(f, a, k) * taper[f];
            const auto spectrum = fft(series);
            for (std::size_t q = 0; q < bins; ++q) power[q] += std::norm(spectrum[q]);
        }
    }
    const std::size_t end = count % 2 == 0 ? bins - 1 : bins;
    for (std::size_t q = 1; q < end; ++q) power[q] *= 2;
    std::vector<double> frequency(bins);
    for (std::size_t q = 0; q < bins; ++q) frequency[q] = static_cast<double>(q) / (static_cast<double>(count) * dt) * 1000;
    const double area = trapezoid(power, frequency);
    if (!(area > 0)) throw std::runtime_error("No nonzero vibrational spectral weight");
    for (double& value : power) value /= area;
    Json result = Json::object();
    result["frequency_THz"] = toJson(frequency);
    result["density_per_THz"] = toJson(power);
    result["resolution_THz"] = 1000 / (static_cast<double>(count) * dt);
    result["nyquist_THz"] = 500 / dt;
    return result;
}

Json localStrain(const Parameters& p)
{
    const auto reference = points(p.array("reference"), "positions");
    const double cutoff = positive(p.number("cutoff_A"), "cutoff_A");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    Mat3 referenceCell{};
    const auto neighbors = localNeighbors(reference, cutoff, p, "reference_cell", pbc, referenceCell);
    const auto now = points(p.array("current"), "current");
    if (now.size() != reference.size()) throw std::runtime_error("Reference/current atom ordering and counts must correspond");
    const Mat3 cell = p.has("current_cell") ? cellMatrix(p.array("current_cell")) : referenceCell;
    const MinimumImage mic(cell, pbc);
    const Mat3 mapping = multiply(inverse(referenceCell), cell);
    const std::size_t count = reference.size();
    NdArray deformation({count, 3, 3}, kNaN), strain({count, 3, 3}, kNaN);
    std::vector<double> residual(count, kNaN);
    std::vector<bool> valid(count, false);
    std::vector<int> coordination(count, 0);
    const auto groups = byCenter(neighbors, count);
    for (std::size_t index = 0; index < count; ++index) {
        if (index % 64 == 0) taskProgress(static_cast<double>(index) / static_cast<double>(count));
        const auto& group = groups[index];
        coordination[index] = static_cast<int>(group.size());
        if (group.size() < 3) continue;
        Eigen::MatrixXd original(static_cast<int>(group.size()), 3), displaced(static_cast<int>(group.size()), 3);
        for (std::size_t n = 0; n < group.size(); ++n) {
            const auto& neighbor = neighbors[group[n]];
            Vec3 d = add(sub(now[static_cast<std::size_t>(neighbor.j)], now[index]),
                         rowTimes({double(neighbor.shift[0]), double(neighbor.shift[1]), double(neighbor.shift[2])}, cell));
            if (periodic) {
                const Vec3 prediction = rowTimes(neighbor.vector, mapping);
                d = add(prediction, mic(sub(d, prediction)));
            }
            for (int k = 0; k < 3; ++k) {
                original(static_cast<int>(n), k) = neighbor.vector[k];
                displaced(static_cast<int>(n), k) = d[k];
            }
        }
        if (matrixRank(original) < 3) continue;
        // Row-vector displacements Y = X F^T.
        const Eigen::MatrixXd transposed = leastSquares(original, displaced);
        const Eigen::Matrix3d green = (transposed * transposed.transpose() - Eigen::Matrix3d::Identity()) / 2;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                deformation(index, i, j) = transposed(static_cast<int>(j), static_cast<int>(i));
                strain(index, i, j) = green(static_cast<int>(i), static_cast<int>(j));
            }
        residual[index] = (displaced - original * transposed).squaredNorm();
        valid[index] = true;
    }
    Json result = Json::object();
    result["deformation_gradient"] = toJson(deformation);
    result["green_lagrange_strain"] = toJson(strain);
    result["d2min_A2"] = toJson(residual);
    result["valid"] = toJson(valid);
    result["coordination"] = toJson(coordination);
    return result;
}

Json centrosymmetry(const Parameters& p)
{
    const long long neighborCount = integer(p.number("neighbors", 12), "neighbors", 2);
    if (neighborCount % 2) throw std::runtime_error("neighbors must be even");
    const auto positions = points(p.array("positions"), "positions");
    const double cutoff = positive(p.number("cutoff_A"), "cutoff_A");
    Mat3 cell{};
    const auto neighbors = localNeighbors(positions, cutoff, p, "cell", p.pbc("pbc", {false, false, false}), cell);
    const auto groups = byCenter(neighbors, positions.size());
    const std::size_t n = static_cast<std::size_t>(neighborCount);
    std::vector<double> values(positions.size(), kNaN);
    std::vector<bool> valid(positions.size(), false);
    for (std::size_t index = 0; index < positions.size(); ++index) {
        if (index % 256 == 0) taskProgress(static_cast<double>(index) / static_cast<double>(positions.size()));
        if (groups[index].size() < n) continue;
        std::vector<Vec3> bonds;
        for (std::size_t g : groups[index]) bonds.push_back(neighbors[g].vector);
        std::stable_sort(bonds.begin(), bonds.end(), [](const Vec3& a, const Vec3& b) { return norm(a) < norm(b); });
        bonds.resize(n);
        std::vector<double> costs;
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = a + 1; b < n; ++b) {
                const Vec3 sum = add(bonds[a], bonds[b]);
                costs.push_back(dot(sum, sum));
            }
        std::sort(costs.begin(), costs.end());
        values[index] = std::accumulate(costs.begin(), costs.begin() + static_cast<std::ptrdiff_t>(n / 2), 0.0);
        valid[index] = true;
    }
    Json result = Json::object();
    result["centrosymmetry_A2"] = toJson(values);
    result["valid"] = toJson(valid);
    return result;
}

Json bondOrder(const Parameters& p)
{
    std::vector<int> orders;
    if (p.has("degrees")) {
        for (double degree : p.array("degrees").values) orders.push_back(static_cast<int>(integer(degree, "degree", 0)));
    } else orders = {4, 6};
    std::vector<int> sorted = orders;
    std::sort(sorted.begin(), sorted.end());
    if (orders.empty() || std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) throw std::runtime_error("Supply distinct degrees");
    for (int order : orders)
        if (order > 32) throw std::runtime_error("degree must be <= 32");
    const auto positions = points(p.array("positions"), "positions");
    const double cutoff = positive(p.number("cutoff_A"), "cutoff_A");
    Mat3 cell{};
    const auto neighbors = localNeighbors(positions, cutoff, p, "cell", p.pbc("pbc", {false, false, false}), cell);
    const auto groups = byCenter(neighbors, positions.size());
    std::vector<int> coordination(positions.size());
    std::vector<bool> valid(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        coordination[i] = static_cast<int>(groups[i].size());
        valid[i] = !groups[i].empty();
    }
    Json orderResult = Json::object();
    for (int l : orders) {
        std::vector<double> q(positions.size(), kNaN);
        for (std::size_t index = 0; index < positions.size(); ++index) {
            if (index % 256 == 0) taskCheckpoint();
            const auto& group = groups[index];
            if (group.empty()) continue;
            double total = 0;
            for (int m = 0; m <= l; ++m) {
                std::complex<double> mean = 0;
                for (std::size_t g : group) {
                    const Vec3& bond = neighbors[g].vector;
                    const double cosine = std::clamp(bond[2] / norm(bond), -1.0, 1.0);
                    const double azimuth = std::atan2(bond[1], bond[0]);
                    mean += normalizedLegendre(l, m, cosine) * std::polar(1.0, m * azimuth);
                }
                mean /= static_cast<double>(group.size());
                total += (m ? 2.0 : 1.0) * std::norm(mean);
            }
            q[index] = std::sqrt(4 * kPi / (2 * l + 1) * total);
        }
        orderResult["q" + std::to_string(l)] = toJson(q);
    }
    Json result = Json::object();
    result["order"] = orderResult;
    result["coordination"] = toJson(coordination);
    result["valid"] = toJson(valid);
    return result;
}

Json wignerSeitz(const Parameters& p)
{
    const auto sites = points(p.array("reference_sites"), "reference_sites");
    std::vector<Vec3> atoms;
    const NdArray raw = p.array("positions");
    if (!(raw.size() == 0 && (raw.shape == std::vector<std::size_t>{0} || raw.shape == std::vector<std::size_t>{0, 3})))
        atoms = points(raw, "positions");
    const Pbc pbc = p.pbc("pbc", {false, false, false});
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    Mat3 cell = identity();
    if (p.has("cell")) cell = cellMatrix(p.array("cell"));
    else if (periodic || p.has("current_cell")) throw std::runtime_error("cell must be a nonempty finite array of the required dimension");
    if (p.has("current_cell")) {
        const Mat3 mapping = multiply(inverse(cellMatrix(p.array("current_cell"))), cell);
        for (auto& atom : atoms) atom = rowTimes(atom, mapping);
    }
    const MinimumImage mic(cell, pbc);
    std::vector<int> assignments, occupancy(sites.size(), 0);
    std::vector<double> distances;
    std::vector<bool> ties;
    for (std::size_t a = 0; a < atoms.size(); ++a) {
        if (a % 64 == 0) taskProgress(static_cast<double>(a) / static_cast<double>(atoms.size()));
        std::vector<double> norms(sites.size());
        for (std::size_t s = 0; s < sites.size(); ++s) norms[s] = norm(mic(sub(sites[s], atoms[a])));
        const std::size_t selected = static_cast<std::size_t>(std::min_element(norms.begin(), norms.end()) - norms.begin());
        int close = 0;
        for (double value : norms) close += std::abs(value - norms[selected]) <= 1e-10 + 1e-10 * std::abs(norms[selected]);
        assignments.push_back(static_cast<int>(selected));
        distances.push_back(norms[selected]);
        ties.push_back(close > 1);
        ++occupancy[selected];
    }
    std::vector<int> vacancySites, interstitialSites;
    int vacancies = 0, excess = 0;
    for (std::size_t s = 0; s < sites.size(); ++s) {
        if (occupancy[s] == 0) { vacancySites.push_back(static_cast<int>(s)); ++vacancies; }
        if (occupancy[s] > 1) { interstitialSites.push_back(static_cast<int>(s)); excess += occupancy[s] - 1; }
    }
    Json result = Json::object();
    result["occupancy"] = toJson(occupancy);
    result["site_index"] = toJson(assignments);
    result["distance_A"] = toJson(distances);
    result["ambiguous"] = toJson(ties);
    result["vacancy_sites"] = toJson(vacancySites);
    result["interstitial_sites"] = toJson(interstitialSites);
    result["vacancies"] = vacancies;
    result["interstitial_excess"] = excess;
    return result;
}

Json staticStructureFactor(const Parameters& p)
{
    NdArray coordinates = finiteArray(p.array("positions"), "positions");
    if (coordinates.ndim() == 2) coordinates.shape.insert(coordinates.shape.begin(), 1);
    if (coordinates.ndim() != 3 || coordinates.shape[2] != 3)
        throw std::runtime_error("positions must have shape (atoms,3) or (frames,atoms,3)");
    const auto wavevectors = points(p.array("q_vectors"), "q_vectors");
    const std::size_t frames = coordinates.shape[0], atoms = coordinates.shape[1];
    std::vector<double> factors(atoms, 1.0);
    if (p.has("weights")) factors = p.array("weights").values;
    double total = 0;
    bool finite = factors.size() == atoms;
    for (double factor : factors) { finite = finite && std::isfinite(factor); total += factor * factor; }
    if (!finite || total == 0) throw std::runtime_error("weights must contain one finite scattering amplitude per atom, not all zero");
    std::vector<double> intensity(wavevectors.size(), 0.0);
    for (std::size_t f = 0; f < frames; ++f)
        for (std::size_t q = 0; q < wavevectors.size(); ++q) {
            if (q % 64 == 0) taskProgress((static_cast<double>(f) + static_cast<double>(q) / static_cast<double>(wavevectors.size())) / static_cast<double>(frames));
            std::complex<double> amplitude = 0;
            for (std::size_t a = 0; a < atoms; ++a) {
                const double phase = wavevectors[q][0] * coordinates(f, a, 0) + wavevectors[q][1] * coordinates(f, a, 1) + wavevectors[q][2] * coordinates(f, a, 2);
                amplitude += factors[a] * std::polar(1.0, phase);
            }
            intensity[q] += std::norm(amplitude);
        }
    for (double& value : intensity) value /= static_cast<double>(frames) * total;
    Json qJson = Json::array();
    for (const auto& q : wavevectors) qJson.push(toJson(q));
    Json result = Json::object();
    result["q_vectors_rad_per_A"] = qJson;
    result["S_q"] = toJson(intensity);
    return result;
}

Json bandGap(const Parameters& p)
{
    NdArray energy = finiteArray(p.array("energies_eV"), "energies_eV");
    const double fermi = p.number("fermi_eV");
    const double tolerance = positive(p.number("tolerance_eV", 1e-6), "tolerance_eV");
    if (energy.ndim() == 2) energy.shape.insert(energy.shape.begin(), 1);
    if (energy.ndim() != 3) throw std::runtime_error("energies_eV must have shape (k,band) or (spin,k,band)");
    const std::size_t spins = energy.shape[0], kpoints = energy.shape[1], bands = energy.shape[2];
    bool anyBelow = false, anyAbove = false, touching = false, crossing = false;
    double valence = -HUGE_VAL, conduction = HUGE_VAL;
    for (double value : energy.values) {
        const bool below = value < fermi - tolerance, above = value > fermi + tolerance;
        anyBelow = anyBelow || below; anyAbove = anyAbove || above;
        touching = touching || std::abs(value - fermi) <= tolerance;
        if (below) valence = std::max(valence, value);
        if (above) conduction = std::min(conduction, value);
    }
    if (!anyBelow || !anyAbove) throw std::runtime_error("Include both occupied and unoccupied bands");
    for (std::size_t s = 0; s < spins; ++s)
        for (std::size_t b = 0; b < bands; ++b) {
            bool below = false, above = false;
            for (std::size_t k = 0; k < kpoints; ++k) {
                below = below || energy(s, k, b) < fermi - tolerance;
                above = above || energy(s, k, b) > fermi + tolerance;
            }
            crossing = crossing || (below && above);
        }
    const bool metal = crossing || touching;
    double direct = HUGE_VAL;
    for (std::size_t s = 0; s < spins; ++s)
        for (std::size_t k = 0; k < kpoints; ++k) {
            double occupied = -HUGE_VAL, empty = HUGE_VAL;
            for (std::size_t b = 0; b < bands; ++b) {
                const double value = energy(s, k, b);
                if (value < fermi - tolerance) occupied = std::max(occupied, value);
                if (value > fermi + tolerance) empty = std::min(empty, value);
            }
            if (std::isfinite(occupied) && std::isfinite(empty)) direct = std::min(direct, empty - occupied);
        }
    Json vbm = Json::array(), cbm = Json::array();
    for (std::size_t s = 0; s < spins; ++s)
        for (std::size_t k = 0; k < kpoints; ++k)
            for (std::size_t b = 0; b < bands; ++b) {
                const double value = energy(s, k, b);
                const Json index = Json::array({Json(s), Json(k), Json(b)});
                if (value < fermi - tolerance && std::abs(value - valence) <= tolerance) vbm.push(index);
                if (value > fermi + tolerance && std::abs(value - conduction) <= tolerance) cbm.push(index);
            }
    Json result = Json::object();
    result["metal_on_sampled_mesh"] = metal;
    result["gap_eV"] = metal ? 0.0 : conduction - valence;
    result["direct_gap_eV"] = metal ? Json(0.0) : (std::isfinite(direct) ? Json(direct) : Json());
    result["vbm_eV"] = valence;
    result["cbm_eV"] = conduction;
    result["vbm_indices_spin_k_band"] = vbm;
    result["cbm_indices_spin_k_band"] = cbm;
    return result;
}

Json effectiveMass(const Parameters& p)
{
    const auto k = points(p.array("kpoints_inv_A"), "kpoints_inv_A");
    const auto energy = vector1d(p.array("energies_eV"), "energies_eV");
    const auto center = vector1d(p.array("center_inv_A"), "center_inv_A");
    if (center.size() != 3 || energy.size() != k.size() || k.size() < 10)
        throw std::runtime_error("At least ten 3D k points, matching energies and one center are required");
    const int rows = static_cast<int>(k.size());
    Eigen::MatrixXd design(rows, 10);
    Eigen::VectorXd values(rows);
    for (int i = 0; i < rows; ++i) {
        const double x = k[static_cast<std::size_t>(i)][0] - center[0], y = k[static_cast<std::size_t>(i)][1] - center[1], z = k[static_cast<std::size_t>(i)][2] - center[2];
        design.row(i) << 1, x, y, z, x * x / 2, y * y / 2, z * z / 2, x * y, x * z, y * z;
        values(i) = energy[static_cast<std::size_t>(i)];
    }
    // Scaling protects tiny neighborhoods against misleading rank decisions.
    const Eigen::RowVectorXd scale = design.colwise().norm();
    if ((scale.array() == 0).any()) throw std::runtime_error("k points must span a full 3D quadratic fit; a line path is insufficient");
    const Eigen::MatrixXd scaled = design * scale.cwiseInverse().asDiagonal();
    if (matrixRank(scaled) != 10) throw std::runtime_error("k points must span a full 3D quadratic fit; a line path is insufficient");
    const Eigen::VectorXd coefficients = scale.cwiseInverse().transpose().asDiagonal() * Eigen::VectorXd(leastSquares(scaled, values));
    Eigen::Matrix3d hessian;
    hessian << coefficients(4), coefficients(7), coefficients(8),
               coefficients(7), coefficients(5), coefficients(9),
               coefficients(8), coefficients(9), coefficients(6);
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(hessian);
    const Eigen::Vector3d curvature = solver.eigenvalues();
    if ((curvature.array().abs() < 1e-10).any()) throw std::runtime_error("Singular band curvature: finite effective mass cannot be determined");
    // hbar^2 / m_e in eV A^2 (CODATA conversion).
    const double constant = 7.619964231073853;
    const Eigen::VectorXd fitted = design * coefficients;
    Json result = Json::object();
    result["principal_masses_m_e"] = vectorJson((constant / curvature.array()).matrix());
    result["principal_axes_columns"] = matrixJson(solver.eigenvectors());
    result["mass_tensor_m_e"] = matrixJson(constant * hessian.inverse());
    result["hessian_eV_A2"] = matrixJson(hessian);
    result["gradient_eV_A"] = vectorJson(coefficients.segment(1, 3));
    result["energy_at_center_eV"] = coefficients(0);
    result["rms_fit_error_eV"] = std::sqrt((fitted - values).squaredNorm() / rows);
    return result;
}

Json workFunction(const Parameters& p)
{
    const auto position = vector1d(p.array("distance_A"), "distance_A");
    const auto potential = vector1d(p.array("potential_eV"), "potential_eV");
    const auto region = vector1d(p.array("vacuum_range_A"), "vacuum_range_A");
    const double fermi = p.number("fermi_eV");
    const double limit = positive(p.number("max_slope_eV_per_A", 0.01), "max_slope_eV_per_A");
    bool increasing = position.size() == potential.size();
    for (std::size_t i = 1; increasing && i < position.size(); ++i) increasing = position[i] > position[i - 1];
    if (!increasing) throw std::runtime_error("Distances must increase and match the potential");
    if (region.size() != 2 || region[0] >= region[1] || region[0] < position.front() || region[1] > position.back())
        throw std::runtime_error("Vacuum interval must lie within the supplied profile");
    std::vector<double> x, y;
    for (std::size_t i = 0; i < position.size(); ++i)
        if (position[i] >= region[0] && position[i] <= region[1]) { x.push_back(position[i]); y.push_back(potential[i]); }
    if (x.size() < 3) throw std::runtime_error("At least three vacuum samples are required");
    const double slope = linearRegression(x, y).slope;
    const double vacuum = std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(y.size());
    double variance = 0;
    for (double value : y) variance += (value - vacuum) * (value - vacuum);
    Json result = Json::object();
    result["work_function_eV"] = vacuum - fermi;
    result["vacuum_level_eV"] = vacuum;
    result["vacuum_std_eV"] = std::sqrt(variance / static_cast<double>(y.size()));
    result["slope_eV_per_A"] = slope;
    result["flat_vacuum"] = std::abs(slope) <= limit;
    result["samples"] = x.size();
    return result;
}

namespace
{
// Third-order Birch-Murnaghan E(V) with parameters (E0, B0, B0', V0), as in ASE.
double birchMurnaghan(double volume, const Eigen::Vector4d& q, Eigen::Vector4d* gradient = nullptr)
{
    const double x = std::pow(q(3) / volume, 2.0 / 3.0) - 1.0;
    const double shape = 2 * x * x + (q(2) - 4) * x * x * x;
    const double prefactor = 9.0 / 16.0;
    if (gradient) {
        (*gradient)(0) = 1;
        (*gradient)(1) = prefactor * q(3) * shape;
        (*gradient)(2) = prefactor * q(1) * q(3) * x * x * x;
        (*gradient)(3) = prefactor * q(1) * (shape + 2.0 / 3.0 * (x + 1) * (4 * x + 3 * (q(2) - 4) * x * x));
    }
    return q(0) + prefactor * q(1) * q(3) * shape;
}
}

Json equationOfState(const Parameters& p)
{
    const auto volumes = vector1d(p.array("volumes_A3"), "volumes_A3");
    const auto energies = vector1d(p.array("energies_eV"), "energies_eV");
    std::vector<double> unique = volumes;
    std::sort(unique.begin(), unique.end());
    const bool distinct = std::adjacent_find(unique.begin(), unique.end()) == unique.end();
    if (volumes.size() != energies.size() || volumes.size() < 5 || unique.front() <= 0 || !distinct)
        throw std::runtime_error("Provide at least five distinct positive volumes and matching energies");
    const int n = static_cast<int>(volumes.size());
    Eigen::MatrixXd quadraticDesign(n, 3);
    Eigen::VectorXd e(n);
    for (int i = 0; i < n; ++i) {
        const double v = volumes[static_cast<std::size_t>(i)];
        quadraticDesign.row(i) << v * v, v, 1;
        e(i) = energies[static_cast<std::size_t>(i)];
    }
    const Eigen::VectorXd quadratic = leastSquares(quadraticDesign, e);
    if (quadratic(0) <= 0) throw std::runtime_error("Energy-volume samples must contain a convex minimum");
    const double vmin = unique.front(), vmax = unique.back();
    const double v0 = -quadratic(1) / (2 * quadratic(0));
    if (!(vmin < v0 && v0 < vmax)) throw std::runtime_error("Volume samples must bracket the equilibrium minimum");
    Eigen::Vector4d q(quadratic(0) * v0 * v0 + quadratic(1) * v0 + quadratic(2), 2 * quadratic(0) * v0, 4.0, v0);
    const Eigen::Vector4d lower(-HUGE_VAL, 1e-12, -20, vmin), upper(HUGE_VAL, HUGE_VAL, 30, vmax);
    auto evaluate = [&](const Eigen::Vector4d& parameters, Eigen::VectorXd& residual, Eigen::MatrixXd& jacobian) {
        residual.resize(n);
        jacobian.resize(n, 4);
        for (int i = 0; i < n; ++i) {
            Eigen::Vector4d gradient;
            residual(i) = e(i) - birchMurnaghan(volumes[static_cast<std::size_t>(i)], parameters, &gradient);
            jacobian.row(i) = gradient.transpose();
        }
        return residual.squaredNorm();
    };
    // Bound-projected Levenberg-Marquardt with Marquardt diagonal scaling.
    Eigen::VectorXd residual;
    Eigen::MatrixXd jacobian;
    double cost = evaluate(q, residual, jacobian);
    double lambda = 1e-3;
    for (int iteration = 0; iteration < 20000; ++iteration) {
        taskCheckpoint();
        const Eigen::Matrix4d normal = jacobian.transpose() * jacobian;
        const Eigen::Vector4d gradient = jacobian.transpose() * residual;
        bool improved = false;
        while (lambda < 1e20) {
            Eigen::Matrix4d damped = normal;
            for (int k = 0; k < 4; ++k) damped(k, k) += lambda * std::max(normal(k, k), 1e-30);
            Eigen::Vector4d candidate = q + damped.ldlt().solve(gradient);
            for (int k = 0; k < 4; ++k) candidate(k) = std::clamp(candidate(k), lower(k), upper(k));
            Eigen::VectorXd candidateResidual;
            Eigen::MatrixXd candidateJacobian;
            const double candidateCost = evaluate(candidate, candidateResidual, candidateJacobian);
            if (std::isfinite(candidateCost) && candidateCost < cost) {
                const double change = (candidate - q).cwiseAbs().maxCoeff();
                const double relative = (cost - candidateCost) / std::max(cost, 1e-300);
                q = candidate; residual = candidateResidual; jacobian = candidateJacobian; cost = candidateCost;
                lambda = std::max(lambda / 10, 1e-15);
                improved = true;
                if (change < 1e-15 * (1 + q.cwiseAbs().maxCoeff()) || relative < 1e-15) iteration = 1 << 30;
                break;
            }
            lambda *= 10;
        }
        if (!improved || cost == 0.0) break;
    }
    const double e0 = q(0), modulus = q(1), derivative = q(2), volume = q(3);
    if (!(vmin + 1e-8 < volume && volume < vmax - 1e-8))
        throw std::runtime_error("Fitted equilibrium lies on the volume boundary; sample a wider range");
    std::vector<double> residuals(volumes.size());
    for (std::size_t i = 0; i < volumes.size(); ++i) residuals[i] = energies[i] - birchMurnaghan(volumes[i], q);
    double squared = 0;
    for (double value : residuals) squared += value * value;
    Eigen::Matrix4d covariance = Eigen::Matrix4d::Constant(HUGE_VAL);
    const Eigen::Matrix4d normal = jacobian.transpose() * jacobian;
    if (n > 4 && matrixRank(normal) == 4) covariance = normal.inverse() * (squared / (n - 4));
    Json result = Json::object();
    result["volume_A3"] = volume;
    result["energy_eV"] = e0;
    result["bulk_modulus_GPa"] = modulus * kEvPerA3ToGPa;
    result["bulk_derivative"] = derivative;
    result["residuals_eV"] = toJson(residuals);
    result["rms_error_eV"] = std::sqrt(squared / static_cast<double>(residuals.size()));
    result["parameter_covariance"] = matrixJson(covariance);
    return result;
}

Json elasticTensor(const Parameters& p)
{
    const NdArray strain = finiteArray(p.array("strains"), "strains", 2);
    const NdArray stress = finiteArray(p.array("stresses_GPa"), "stresses_GPa", 2);
    if (strain.shape != stress.shape || strain.shape[1] != 6 || strain.shape[0] < 7)
        throw std::runtime_error("Provide matching (samples>=7,6) engineering strains and stresses");
    const int n = static_cast<int>(strain.shape[0]);
    Eigen::MatrixXd design(n, 7), target(n, 6);
    for (int i = 0; i < n; ++i) {
        design(i, 0) = 1;
        for (int j = 0; j < 6; ++j) {
            design(i, j + 1) = strain(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
            target(i, j) = stress(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
        }
    }
    if (matrixRank(design) != 7) throw std::runtime_error("Strains must independently sample all six strain components");
    const Eigen::MatrixXd coefficients = leastSquares(design, target);
    const Eigen::MatrixXd raw = coefficients.bottomRows(6).transpose();
    const Eigen::MatrixXd stiffness = (raw + raw.transpose()) / 2;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(stiffness);
    const Eigen::VectorXd eigenvalues = solver.eigenvalues();
    const bool stable = (eigenvalues.array() > 0).all();
    Json moduli;
    if (stable) {
        const Eigen::MatrixXd compliance = stiffness.inverse();
        auto diagonal = [](const Eigen::MatrixXd& m) {
            return std::array<double, 3>{m(0, 0) + m(1, 1) + m(2, 2), m(0, 1) + m(0, 2) + m(1, 2), m(3, 3) + m(4, 4) + m(5, 5)};
        };
        auto [a, b, c] = diagonal(stiffness);
        const double bulkV = (a + 2 * b) / 9, shearV = (a - b + 3 * c) / 15;
        const auto s = diagonal(compliance);
        const double bulkR = 1 / (s[0] + 2 * s[1]), shearR = 15 / (4 * s[0] - 4 * s[1] + 3 * s[2]);
        const double bulk = (bulkV + bulkR) / 2, shear = (shearV + shearR) / 2;
        moduli = Json::object();
        moduli["bulk_voigt_GPa"] = bulkV;
        moduli["bulk_reuss_GPa"] = bulkR;
        moduli["shear_voigt_GPa"] = shearV;
        moduli["shear_reuss_GPa"] = shearR;
        moduli["bulk_hill_GPa"] = bulk;
        moduli["shear_hill_GPa"] = shear;
        moduli["young_hill_GPa"] = 9 * bulk * shear / (3 * bulk + shear);
        moduli["poisson_hill"] = (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear));
        moduli["universal_anisotropy"] = 5 * shearV / shearR + bulkV / bulkR - 6;
    }
    const Eigen::MatrixXd prediction = design.col(0) * coefficients.row(0) + design.rightCols(6) * stiffness.transpose();
    const Eigen::MatrixXd residual = target - prediction;
    Json result = Json::object();
    result["stiffness_GPa"] = matrixJson(stiffness);
    result["prestress_GPa"] = vectorJson(coefficients.row(0).transpose());
    result["stable_zero_prestress"] = stable;
    result["stiffness_eigenvalues_GPa"] = vectorJson(eigenvalues);
    result["antisymmetric_norm_GPa"] = (raw - raw.transpose()).norm();
    result["rms_stress_error_GPa"] = std::sqrt(residual.squaredNorm() / static_cast<double>(residual.size()));
    result["moduli"] = moduli;
    return result;
}

namespace
{
struct Modes { NdArray energies; std::vector<double> weights; };

Modes modes(const Parameters& p)
{
    Modes result;
    result.energies = finiteArray(p.array("energies_eV"), "energies_eV", 2);
    const std::size_t rows = result.energies.shape[0];
    if (p.has("weights")) result.weights = vector1d(p.array("weights"), "weights");
    else result.weights.assign(rows, 1.0);
    double total = 0;
    bool valid = result.weights.size() == rows;
    for (double weight : result.weights) { valid = valid && weight >= 0; total += weight; }
    if (!valid || total <= 0) throw std::runtime_error("Supply one nonnegative q-point weight per row, with positive total");
    for (double& weight : result.weights) weight /= total;
    return result;
}
}

Json phononDos(const Parameters& p)
{
    const Modes m = modes(p);
    const auto grid = vector1d(p.array("energy_grid_eV"), "energy_grid_eV");
    bool increasing = grid.size() >= 2;
    for (std::size_t i = 1; increasing && i < grid.size(); ++i) increasing = grid[i] > grid[i - 1];
    if (!increasing) throw std::runtime_error("Energy grid must contain at least two increasing samples");
    const double sigma = positive(p.number("sigma_eV", 0.001), "sigma_eV");
    std::vector<double> density(grid.size(), 0.0);
    double imaginary = 0;
    for (std::size_t q = 0; q < m.energies.shape[0]; ++q) {
        taskProgress(static_cast<double>(q) / static_cast<double>(m.energies.shape[0]));
        for (std::size_t b = 0; b < m.energies.shape[1]; ++b) {
            const double mode = m.energies(q, b);
            if (mode < 0) imaginary += m.weights[q];
            for (std::size_t g = 0; g < grid.size(); ++g)
                density[g] += m.weights[q] * std::exp(-0.5 * std::pow((grid[g] - mode) / sigma, 2)) / (sigma * std::sqrt(2 * kPi));
        }
    }
    Json result = Json::object();
    result["energy_eV"] = toJson(grid);
    result["dos_per_eV"] = toJson(density);
    result["total_modes"] = m.energies.shape[1];
    result["enclosed_modes"] = trapezoid(density, grid);
    result["imaginary_weight"] = imaginary;
    return result;
}

Json harmonicThermodynamics(const Parameters& p)
{
    const Modes m = modes(p);
    const auto temperatures = vector1d(p.array("temperatures_K"), "temperatures_K");
    const double tolerance = positive(p.number("zero_tolerance_eV", 1e-8), "zero_tolerance_eV");
    bool valid = true;
    for (double t : temperatures) valid = valid && t >= 0;
    for (double e : m.energies.values) valid = valid && e >= -tolerance;
    if (!valid) throw std::runtime_error("Nonnegative temperatures and dynamically stable phonons are required");
    std::vector<double> modeEnergies, modeWeights;
    double omitted = 0;
    for (std::size_t q = 0; q < m.energies.shape[0]; ++q)
        for (std::size_t b = 0; b < m.energies.shape[1]; ++b) {
            if (m.energies(q, b) > tolerance) { modeEnergies.push_back(m.energies(q, b)); modeWeights.push_back(m.weights[q]); }
            else omitted += m.weights[q];
        }
    if (modeEnergies.empty()) throw std::runtime_error("No positive-frequency modes remain");
    double zpe = 0;
    for (std::size_t i = 0; i < modeEnergies.size(); ++i) zpe += modeWeights[i] * modeEnergies[i] / 2;
    std::vector<double> free, internal, entropy, capacity;
    for (double temperature : temperatures) {
        double f = zpe, u = zpe, s = 0, cv = 0;
        if (temperature > 0) {
            double logSum = 0, energySum = 0, heat = 0;
            for (std::size_t i = 0; i < modeEnergies.size(); ++i) {
                const double x = modeEnergies[i] / (kBoltzmann * temperature);
                const double decay = std::exp(-x);
                const double denominator = -std::expm1(-x);
                logSum += modeWeights[i] * std::log(denominator);
                energySum += modeWeights[i] * modeEnergies[i] * decay / denominator;
                // Avoid inf*0 for very low T by omitting frozen modes (x>700).
                if (x < 700) heat += modeWeights[i] * x * x * decay / (denominator * denominator);
            }
            f = zpe + kBoltzmann * temperature * logSum;
            u = zpe + energySum;
            s = (u - f) / temperature;
            cv = kBoltzmann * heat;
        }
        free.push_back(f); internal.push_back(u); entropy.push_back(s); capacity.push_back(cv);
    }
    Json result = Json::object();
    result["temperature_K"] = toJson(temperatures);
    result["free_energy_eV"] = toJson(free);
    result["internal_energy_eV"] = toJson(internal);
    result["entropy_eV_per_K"] = toJson(entropy);
    result["heat_capacity_eV_per_K"] = toJson(capacity);
    result["zero_point_energy_eV"] = zpe;
    result["omitted_zero_mode_weight"] = omitted;
    return result;
}
}
