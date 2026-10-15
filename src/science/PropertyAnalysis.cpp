#include "science/PropertyAnalysis.h"
#include "util/TaskControl.h"

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
const double kNaN = std::numeric_limits<double>::quiet_NaN();

// Per-atom values (1 for every atom when absent); NaN marks atoms without a value.
std::vector<double> perAtom(const Parameters& p, const std::string& name, std::size_t count)
{
    if (!p.has(name)) return std::vector<double>(count, 1.0);
    const NdArray array = p.array(name);
    if (array.ndim() != 1) throw std::runtime_error(name + " must be a list of numbers");
    if (array.values.size() != count) throw std::runtime_error(name + " needs one value per atom");
    return array.values;
}

std::array<long long, 3> binCounts(const Parameters& p, std::array<long long, 3> fallback)
{
    if (!p.has("bins")) return fallback;
    const auto values = p.array("bins").values;
    if (values.size() != 3) throw std::runtime_error("bins needs three counts (1 = not binned along that cell vector)");
    std::array<long long, 3> out{};
    for (std::size_t d = 0; d < 3; ++d) out[d] = integer(values[d], "bins", 1);
    if (out[0] * out[1] * out[2] > 50'000'000) throw std::runtime_error("Too many bins");
    return out;
}

// Discrete Fourier transform of a 3D grid (row-major, axes n0 n1 n2) in place.
void fft3(std::vector<std::complex<double>>& grid, const std::array<std::size_t, 3>& n, bool inverse)
{
    const std::array<std::size_t, 3> stride = {n[1] * n[2], n[2], 1};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (n[axis] == 1) continue;
        const std::size_t a = (axis + 1) % 3, b = (axis + 2) % 3;
        std::vector<std::complex<double>> line(n[axis]);
        for (std::size_t i = 0; i < n[a]; ++i)
            for (std::size_t j = 0; j < n[b]; ++j) {
                const std::size_t base = i * stride[a] + j * stride[b];
                for (std::size_t k = 0; k < n[axis]; ++k) line[k] = grid[base + k * stride[axis]];
                line = fft(line, inverse);
                for (std::size_t k = 0; k < n[axis]; ++k) grid[base + k * stride[axis]] = line[k];
            }
        taskCheckpoint();
    }
}
}

Json spatialBinning(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    if (!p.has("cell")) throw std::runtime_error("Spatial binning needs the cell");
    const Mat3 cell = cellMatrix(p.array("cell"));
    const Pbc pbc = p.pbc("pbc", {true, true, true});
    const auto values = perAtom(p, "values", positions.size());
    const auto bins = binCounts(p, {20, 1, 1});
    const std::string operation = p.has("operation") ? p.json("operation").string() : std::string("mean");
    if (operation != "mean" && operation != "sum" && operation != "density" && operation != "min" && operation != "max" && operation != "count")
        throw std::runtime_error("operation must be mean, sum, density, min, max or count");
    const std::size_t total = static_cast<std::size_t>(bins[0] * bins[1] * bins[2]);
    std::vector<double> sum(total, 0.0), low(total, std::numeric_limits<double>::infinity()), high(total, -std::numeric_limits<double>::infinity());
    std::vector<long long> count(total, 0);
    std::size_t outside = 0;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        Vec3 f = fractional(positions[i], cell);
        std::array<long long, 3> index{};
        bool inside = true;
        for (std::size_t d = 0; d < 3; ++d) {
            if (pbc[d]) {
                f[d] -= std::floor(f[d]);
                if (f[d] >= 1) f[d] = 0;  // a tiny negative coordinate wraps to exactly 1
            }
            if (f[d] < 0 || f[d] >= 1) { inside = false; break; }
            index[d] = std::min(bins[d] - 1, static_cast<long long>(f[d] * static_cast<double>(bins[d])));
        }
        if (!inside) { ++outside; continue; }
        if (std::isnan(values[i])) continue;
        const std::size_t k = static_cast<std::size_t>((index[0] * bins[1] + index[1]) * bins[2] + index[2]);
        sum[k] += values[i];
        low[k] = std::min(low[k], values[i]);
        high[k] = std::max(high[k], values[i]);
        ++count[k];
    }
    const double binVolume = cellVolume(cell) / static_cast<double>(total);
    std::vector<double> out(total);
    for (std::size_t k = 0; k < total; ++k) {
        if (operation == "mean") out[k] = count[k] ? sum[k] / static_cast<double>(count[k]) : kNaN;
        else if (operation == "sum") out[k] = sum[k];
        else if (operation == "density") out[k] = sum[k] / binVolume;
        else if (operation == "min") out[k] = count[k] ? low[k] : kNaN;
        else if (operation == "max") out[k] = count[k] ? high[k] : kNaN;
        else out[k] = static_cast<double>(count[k]);
    }
    Json result = Json::object();
    result["bins"] = Json::array({static_cast<double>(bins[0]), static_cast<double>(bins[1]), static_cast<double>(bins[2])});
    result["values"] = toJson(out);
    std::vector<double> counts(count.begin(), count.end());
    result["atom_counts"] = toJson(counts);
    result["operation"] = operation;
    result["bin_volume_A3"] = binVolume;
    result["atoms_outside"] = static_cast<double>(outside);
    // Bin centres along each binned cell vector (fractional and Angstrom).
    Json centers = Json::object();
    const char* axes[3] = {"a", "b", "c"};
    for (std::size_t d = 0; d < 3; ++d) {
        if (bins[d] == 1) continue;
        const double length = norm(cell[d]);
        std::vector<double> fractionalCenters, distance;
        for (long long b = 0; b < bins[d]; ++b) {
            const double f = (static_cast<double>(b) + 0.5) / static_cast<double>(bins[d]);
            fractionalCenters.push_back(f);
            distance.push_back(f * length);
        }
        Json axis = Json::object();
        axis["fractional"] = toJson(fractionalCenters);
        axis["distance_A"] = toJson(distance);
        centers[axes[d]] = axis;
    }
    result["centers"] = centers;
    return result;
}

Json propertyHistogram(const Parameters& p)
{
    const NdArray input = p.array("values");
    if (input.ndim() != 1) throw std::runtime_error("values must be a list of numbers");
    std::vector<double> values = input.values;
    if (p.has("mask")) {
        const auto mask = finiteArray(p.array("mask"), "mask", 1).values;
        if (mask.size() != values.size()) throw std::runtime_error("mask needs one value per atom");
        std::vector<double> kept;
        for (std::size_t i = 0; i < values.size(); ++i)
            if (mask[i] != 0) kept.push_back(values[i]);
        values = std::move(kept);
    }
    values.erase(std::remove_if(values.begin(), values.end(), [](double v) { return std::isnan(v); }), values.end());
    if (values.empty()) throw std::runtime_error("No values to analyse");
    const long long bins = integer(p.number("bins", 50), "bins");
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    double low = sorted.front(), high = sorted.back();
    if (p.has("range")) {
        const auto range = p.array("range").values;
        if (range.size() != 2 || !(range[1] > range[0])) throw std::runtime_error("range must be [low, high] with high > low");
        low = range[0];
        high = range[1];
    }
    // Equal values (up to rounding) get one unit-wide range around them.
    if (high - low <= 1e-9 * std::max(1.0, std::abs(high))) { low -= 0.5; high += 0.5; }
    const double width = (high - low) / static_cast<double>(bins);
    std::vector<double> counts(static_cast<std::size_t>(bins), 0.0), centers(static_cast<std::size_t>(bins));
    for (long long b = 0; b < bins; ++b) centers[static_cast<std::size_t>(b)] = low + (static_cast<double>(b) + 0.5) * width;
    std::size_t outside = 0;
    for (double v : values) {
        if (v < low || v > high) { ++outside; continue; }
        const long long b = std::min(bins - 1, static_cast<long long>((v - low) / width));
        counts[static_cast<std::size_t>(b)] += 1;
    }
    const double n = static_cast<double>(values.size());
    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    const double mean = sum / n;
    double squares = 0;
    for (double v : values) squares += (v - mean) * (v - mean);
    const std::size_t half = sorted.size() / 2;
    Json result = Json::object();
    result["centers"] = toJson(centers);
    result["counts"] = toJson(counts);
    result["bin_width"] = width;
    result["values_outside_range"] = static_cast<double>(outside);
    result["count"] = n;
    result["sum"] = sum;
    result["mean"] = mean;
    result["std"] = std::sqrt(squares / n);
    result["min"] = sorted.front();
    result["max"] = sorted.back();
    result["median"] = sorted.size() % 2 ? sorted[half] : 0.5 * (sorted[half - 1] + sorted[half]);
    return result;
}

Json spatialCorrelation(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    if (!p.has("cell")) throw std::runtime_error("The spatial correlation function needs a periodic cell");
    const Mat3 cell = cellMatrix(p.array("cell"));
    const std::size_t atoms = positions.size();
    const auto a = perAtom(p, "values_a", atoms);
    const auto b = p.has("values_b") ? perAtom(p, "values_b", atoms) : a;
    const double spacing = positive(p.number("grid_spacing_A", 1.0), "grid_spacing_A");
    const long long radialBins = integer(p.number("bins", 50), "bins");
    // Largest unambiguous distance: half the smallest distance between opposite cell faces.
    const double volume = cellVolume(cell);
    double height = std::numeric_limits<double>::infinity();
    for (std::size_t d = 0; d < 3; ++d) height = std::min(height, volume / norm(cross(cell[(d + 1) % 3], cell[(d + 2) % 3])));
    const double rMax = p.number("r_max_A", height / 2);
    if (!(rMax > 0)) throw std::runtime_error("r_max_A must be positive");
    std::array<std::size_t, 3> n{};
    for (std::size_t d = 0; d < 3; ++d) n[d] = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(norm(cell[d]) / spacing)));
    const std::size_t cells = n[0] * n[1] * n[2];
    if (cells > 32'000'000) throw std::runtime_error("FFT grid too fine; increase grid_spacing_A");
    std::vector<std::complex<double>> gridA(cells), gridB(cells), gridN(cells);
    std::size_t used = 0;
    for (std::size_t i = 0; i < atoms; ++i) {
        if (std::isnan(a[i]) || std::isnan(b[i])) continue;  // atoms without a value take no part
        ++used;
        Vec3 f = fractional(positions[i], cell);
        std::array<std::size_t, 3> g{};
        for (std::size_t d = 0; d < 3; ++d) {
            f[d] -= std::floor(f[d]);
            g[d] = static_cast<std::size_t>(std::llround(f[d] * static_cast<double>(n[d]))) % n[d];
        }
        const std::size_t k = (g[0] * n[1] + g[1]) * n[2] + g[2];
        gridA[k] += a[i];
        gridB[k] += b[i];
        gridN[k] += 1.0;
    }
    // Cross-correlations: sum_g A(g) B(g + delta) = IFFT(conj(FFT A) FFT B).
    fft3(gridA, n, false);
    fft3(gridB, n, false);
    fft3(gridN, n, false);
    std::vector<std::complex<double>> ab(cells), nn(cells);
    for (std::size_t k = 0; k < cells; ++k) {
        ab[k] = std::conj(gridA[k]) * gridB[k];
        nn[k] = std::conj(gridN[k]) * gridN[k];
    }
    fft3(ab, n, true);
    fft3(nn, n, true);
    std::vector<double> numerator(static_cast<std::size_t>(radialBins), 0.0), pairs(static_cast<std::size_t>(radialBins), 0.0);
    const double width = rMax / static_cast<double>(radialBins);
    for (std::size_t i = 0; i < n[0]; ++i)
        for (std::size_t j = 0; j < n[1]; ++j)
            for (std::size_t k = 0; k < n[2]; ++k) {
                const auto wrapped = [](std::size_t index, std::size_t size) {
                    const double v = static_cast<double>(index) / static_cast<double>(size);
                    return v >= 0.5 ? v - 1 : v;
                };
                const Vec3 f = {wrapped(i, n[0]), wrapped(j, n[1]), wrapped(k, n[2])};
                const double r = norm(rowTimes(f, cell));
                if (r >= rMax) continue;
                const std::size_t bin = static_cast<std::size_t>(r / width);
                const std::size_t index = (i * n[1] + j) * n[2] + k;
                numerator[bin] += ab[index].real();
                pairs[bin] += nn[index].real();
            }
    if (!used) throw std::runtime_error("No atom has values");
    double meanA = 0, meanB = 0, meanAB = 0;
    for (std::size_t i = 0; i < atoms; ++i)
        if (!std::isnan(a[i]) && !std::isnan(b[i])) { meanA += a[i]; meanB += b[i]; meanAB += a[i] * b[i]; }
    meanA /= static_cast<double>(used);
    meanB /= static_cast<double>(used);
    meanAB /= static_cast<double>(used);
    const double covariance = meanAB - meanA * meanB;
    std::vector<double> radius(static_cast<std::size_t>(radialBins)), correlation(static_cast<std::size_t>(radialBins), kNaN),
        normalized(static_cast<std::size_t>(radialBins), kNaN);
    for (std::size_t bin = 0; bin < radius.size(); ++bin) {
        radius[bin] = (static_cast<double>(bin) + 0.5) * width;
        if (pairs[bin] > 0.5) {
            correlation[bin] = numerator[bin] / pairs[bin];
            if (std::abs(covariance) > 0) normalized[bin] = (correlation[bin] - meanA * meanB) / covariance;
        }
    }
    Json result = Json::object();
    result["r_A"] = toJson(radius);
    result["correlation"] = toJson(correlation);
    result["normalized_correlation"] = toJson(normalized);
    result["pair_counts"] = toJson(pairs);
    result["mean_a"] = meanA;
    result["mean_b"] = meanB;
    result["covariance"] = covariance;
    result["grid"] = Json::array({static_cast<double>(n[0]), static_cast<double>(n[1]), static_cast<double>(n[2])});
    result["r_max_A"] = rMax;
    return result;
}

Json scatterPlot(const Parameters& p)
{
    const NdArray xs = p.array("values_x"), ys = p.array("values_y");
    if (xs.ndim() != 1 || ys.ndim() != 1 || xs.values.size() != ys.values.size())
        throw std::runtime_error("values_x and values_y must be lists of the same length");
    std::vector<double> x, y;
    for (std::size_t i = 0; i < xs.values.size(); ++i)
        if (std::isfinite(xs.values[i]) && std::isfinite(ys.values[i])) { x.push_back(xs.values[i]); y.push_back(ys.values[i]); }
    if (x.size() < 2) throw std::runtime_error("At least two atoms need both values");
    const double n = static_cast<double>(x.size());
    const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n, my = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double sxx = 0, syy = 0, sxy = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
        sxy += (x[i] - mx) * (y[i] - my);
    }
    Json result = Json::object();
    result["x"] = toJson(x);
    result["y"] = toJson(y);
    result["points"] = n;
    result["pearson_r"] = sxx > 0 && syy > 0 ? sxy / std::sqrt(sxx * syy) : kNaN;
    result["slope"] = sxx > 0 ? sxy / sxx : kNaN;
    result["intercept"] = sxx > 0 ? my - sxy / sxx * mx : kNaN;
    return result;
}

Json ambientOcclusion(const Parameters& p)
{
    const auto positions = points(p.array("positions"), "positions");
    const std::size_t n = positions.size();
    if (n == 0) throw std::runtime_error("No atoms");
    std::vector<double> radius(n, 0.0);
    if (p.has("radii")) {
        const NdArray radii = p.array("radii");
        if (radii.values.size() != n) throw std::runtime_error("radii needs one value per atom");
        radius = radii.values;
    } else {
        // Half the typical nearest-neighbour spacing.
        Vec3 low = positions[0], high = positions[0];
        for (const auto& x : positions)
            for (std::size_t d = 0; d < 3; ++d) { low[d] = std::min(low[d], x[d]); high[d] = std::max(high[d], x[d]); }
        double volume = 1;
        for (std::size_t d = 0; d < 3; ++d) volume *= std::max(high[d] - low[d], 1.0);
        std::fill(radius.begin(), radius.end(), 0.55 * std::cbrt(volume / static_cast<double>(n)));
    }
    for (double r : radius)
        if (!(r > 0)) throw std::runtime_error("radii must be positive");
    const long long directions = integer(p.number("directions", 40), "directions");
    const long long resolution = integer(p.number("resolution", 256), "resolution", 16);
    // Directions spread evenly over the sphere (Fibonacci lattice).
    const double golden = std::acos(-1.0) * (3 - std::sqrt(5.0));
    std::vector<double> light(n, 0.0);
    double maxRadius = *std::max_element(radius.begin(), radius.end());
    for (long long k = 0; k < directions; ++k) {
        taskProgress(static_cast<double>(k) / static_cast<double>(directions));
        const double z = 1 - 2 * (static_cast<double>(k) + 0.5) / static_cast<double>(directions);
        const double ring = std::sqrt(std::max(0.0, 1 - z * z));
        const Vec3 dir = {ring * std::cos(golden * static_cast<double>(k)), ring * std::sin(golden * static_cast<double>(k)), z};
        const Vec3 helper = std::abs(dir[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        const Vec3 u = scale(cross(dir, helper), 1 / norm(cross(dir, helper)));
        const Vec3 v = cross(dir, u);
        // Orthographic projection along -dir: the front-most atom at each pixel receives light.
        double lowU = std::numeric_limits<double>::infinity(), lowV = lowU, highU = -lowU, highV = -lowU;
        for (const auto& x : positions) {
            lowU = std::min(lowU, dot(x, u)); highU = std::max(highU, dot(x, u));
            lowV = std::min(lowV, dot(x, v)); highV = std::max(highV, dot(x, v));
        }
        const double extent = std::max(highU - lowU, highV - lowV) + 2 * maxRadius;
        const double pixel = extent / static_cast<double>(resolution);
        const std::size_t side = static_cast<std::size_t>(resolution);
        std::vector<double> depth(side * side, -std::numeric_limits<double>::infinity());
        std::vector<int> owner(side * side, -1);
        for (std::size_t i = 0; i < n; ++i) {
            const double cu = (dot(positions[i], u) - lowU + maxRadius) / pixel, cv = (dot(positions[i], v) - lowV + maxRadius) / pixel;
            const double r = radius[i] / pixel, h = dot(positions[i], dir);
            const long long u0 = std::max(0LL, static_cast<long long>(std::floor(cu - r))), u1 = std::min(static_cast<long long>(side) - 1, static_cast<long long>(std::ceil(cu + r)));
            const long long v0 = std::max(0LL, static_cast<long long>(std::floor(cv - r))), v1 = std::min(static_cast<long long>(side) - 1, static_cast<long long>(std::ceil(cv + r)));
            for (long long a = u0; a <= u1; ++a)
                for (long long b = v0; b <= v1; ++b) {
                    const double du = static_cast<double>(a) + 0.5 - cu, dv = static_cast<double>(b) + 0.5 - cv;
                    const double inside = r * r - du * du - dv * dv;
                    if (inside < 0) continue;
                    // Height of the sphere's front surface at this pixel.
                    const double front = h + std::sqrt(inside) * pixel;
                    const std::size_t k2 = static_cast<std::size_t>(a) * side + static_cast<std::size_t>(b);
                    if (front > depth[k2]) { depth[k2] = front; owner[k2] = static_cast<int>(i); }
                }
        }
        for (int o : owner)
            if (o >= 0) light[static_cast<std::size_t>(o)] += pixel * pixel;
    }
    // Brightness relative to the best-lit atom per unit of its cross-section.
    std::vector<double> brightness(n);
    double best = 0;
    for (std::size_t i = 0; i < n; ++i) {
        brightness[i] = light[i] / (radius[i] * radius[i]);
        best = std::max(best, brightness[i]);
    }
    for (double& b : brightness) b = best > 0 ? b / best : 0;
    Json result = Json::object();
    result["brightness"] = toJson(brightness);
    result["directions"] = static_cast<double>(directions);
    return result;
}
}
