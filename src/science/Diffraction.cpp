#include "science/Diffraction.h"
#include "science/ScatteringFactors.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <array>
#include <map>
#include <tuple>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);

const scattering::Coefficients& coefficients(const std::string& symbol)
{
    for (const auto& entry : scattering::table)
        if (symbol == entry.symbol) return entry;
    throw std::runtime_error("No scattering factor coefficients for " + symbol);
}

struct Crystal
{
    Mat3 cell{};
    Mat3 reciprocal{};  // rows b_i with a_i . b_j = delta_ij (no 2 pi)
    std::vector<std::string> symbols;
    std::vector<Vec3> fractional;
};

Crystal crystalFrom(const StructureInput& input)
{
    const Structure& s = input.structure;
    if (!s.hasUnitCell || s.atoms.empty()) throw std::runtime_error("Diffraction requires a nonempty periodic structure");
    Crystal c;
    for (int r = 0; r < 3; ++r) c.cell[r] = s.cellVectors[r];
    c.reciprocal = transpose(inverse(c.cell));
    for (const auto& atom : s.atoms) {
        c.symbols.push_back(atom.symbol);
        c.fractional.push_back(fractional({atom.x, atom.y, atom.z}, c.cell));
        coefficients(atom.symbol);
    }
    return c;
}

std::complex<double> structureFactor(const Crystal& c, const std::array<int, 3>& hkl, double s,
                                     double (*factor)(const std::string&, double), double debyeWaller)
{
    std::complex<double> total = 0;
    std::map<std::string, double> cache;
    for (std::size_t i = 0; i < c.symbols.size(); ++i) {
        auto found = cache.find(c.symbols[i]);
        if (found == cache.end()) found = cache.emplace(c.symbols[i], factor(c.symbols[i], s) * std::exp(-debyeWaller * s * s)).first;
        const double phase = 2 * kPi * (hkl[0] * c.fractional[i][0] + hkl[1] * c.fractional[i][1] + hkl[2] * c.fractional[i][2]);
        total += found->second * std::polar(1.0, phase);
    }
    return total;
}

// All (h,k,l) with |g| <= gMax (g = h b1 + k b2 + l b3, 1/Angstrom).
std::vector<std::array<int, 3>> reflections(const Crystal& c, double gMax)
{
    int range[3];
    for (int k = 0; k < 3; ++k) range[k] = static_cast<int>(std::ceil(gMax * norm(c.cell[k]))) + 1;
    std::vector<std::array<int, 3>> result;
    for (int h = -range[0]; h <= range[0]; ++h)
        for (int k = -range[1]; k <= range[1]; ++k)
            for (int l = -range[2]; l <= range[2]; ++l) {
                if (!h && !k && !l) continue;
                if (norm(rowTimes({double(h), double(k), double(l)}, c.reciprocal)) <= gMax) result.push_back({h, k, l});
            }
    return result;
}

Json hklJson(const std::array<int, 3>& hkl) { return Json::array({hkl[0], hkl[1], hkl[2]}); }
}

double xrayScatteringFactor(const std::string& symbol, double s)
{
    const auto& c = coefficients(symbol);
    double sum = 0;
    for (int i = 0; i < 4; ++i) sum += c.a[i] * std::exp(-c.b[i] * s * s);
    int z = atomicNumber(symbol);
    return z - 41.78214 * s * s * sum;
}

double electronScatteringFactor(const std::string& symbol, double s)
{
    // Mott-Bethe: f_e = 0.023934 (Z - f_x) / s^2 = 0.023934 * 41.78214 * sum a_i exp(-b_i s^2).
    const auto& c = coefficients(symbol);
    double sum = 0;
    for (int i = 0; i < 4; ++i) sum += c.a[i] * std::exp(-c.b[i] * s * s);
    return 0.023934 * 41.78214 * sum;
}

ToolOutput powderXrd(const StructureInput& structure, const Parameters& p)
{
    const Crystal c = crystalFrom(structure);
    const double wavelength = positive(p.number("wavelength_A", 1.54184), "wavelength_A");
    NdArray range;
    double low = 10, high = 90;
    if (p.has("two_theta_range_deg")) {
        const auto values = finiteArray(p.array("two_theta_range_deg"), "two_theta_range_deg", 1).values;
        if (values.size() != 2 || values[0] < 0 || values[1] > 180 || values[0] >= values[1]) throw std::runtime_error("two_theta_range_deg must be [low, high] within 0-180");
        low = values[0]; high = values[1];
    }
    const double debyeWaller = p.number("debye_waller_A2", 0.0);
    const double fwhm = positive(p.number("fwhm_deg", 0.1), "fwhm_deg");
    const double gMax = 2 * std::sin(high / 2 * kPi / 180) / wavelength;
    // Merge reflections sharing a Bragg angle (within 1e-5 degrees).
    struct Peak { double twoTheta, d, intensity; std::vector<std::array<int, 3>> hkls; };
    std::map<long long, Peak> peaks;
    const auto all = reflections(c, gMax);
    for (std::size_t n = 0; n < all.size(); ++n) {
        if (n % 512 == 0) taskProgress(static_cast<double>(n) / static_cast<double>(all.size()));
        const auto& hkl = all[n];
        const double g = norm(rowTimes({double(hkl[0]), double(hkl[1]), double(hkl[2])}, c.reciprocal));
        const double theta = std::asin(std::min(1.0, wavelength * g / 2));
        const double twoTheta = 2 * theta * 180 / kPi;
        if (twoTheta < low || twoTheta > high) continue;
        const auto f = structureFactor(c, hkl, g / 2, xrayScatteringFactor, debyeWaller);
        const double lorentz = (1 + std::pow(std::cos(2 * theta), 2)) / (std::pow(std::sin(theta), 2) * std::cos(theta));
        const long long key = std::llround(twoTheta / 1e-5);
        auto& peak = peaks[key];
        if (peak.hkls.empty()) { peak.twoTheta = twoTheta; peak.d = 1 / g; }
        peak.intensity += std::norm(f) * lorentz;
        peak.hkls.push_back(hkl);
    }
    double maximum = 0;
    for (const auto& [key, peak] : peaks) maximum = std::max(maximum, peak.intensity);
    Json list = Json::array();
    std::vector<std::pair<double, double>> kept;
    for (const auto& [key, peak] : peaks) {
        const double scaled = maximum > 0 ? 100 * peak.intensity / maximum : 0;
        if (scaled < 1e-3) continue;  // systematic absences and negligible peaks
        // Representative index: the reflection with the largest (h, k, l) in sorted order.
        auto representative = *std::max_element(peak.hkls.begin(), peak.hkls.end(), [](const auto& a, const auto& b) {
            std::array<int, 3> x = {std::abs(a[0]), std::abs(a[1]), std::abs(a[2])}, y = {std::abs(b[0]), std::abs(b[1]), std::abs(b[2])};
            std::sort(x.rbegin(), x.rend()); std::sort(y.rbegin(), y.rend());
            return x < y || (x == y && a < b);
        });
        Json entry = Json::object();
        entry["two_theta_deg"] = peak.twoTheta;
        entry["d_A"] = peak.d;
        entry["intensity"] = scaled;
        entry["hkl"] = hklJson(representative);
        entry["multiplicity"] = peak.hkls.size();
        list.push(entry);
        kept.push_back({peak.twoTheta, scaled});
    }
    // Gaussian-broadened profile on a 0.02 degree grid.
    const double sigma = fwhm / (2 * std::sqrt(2 * std::log(2.0)));
    std::vector<double> x, y;
    for (double t = low; t <= high + 1e-9; t += 0.02) {
        double value = 0;
        for (const auto& [position, intensity] : kept)
            if (std::abs(t - position) < 6 * sigma) value += intensity * std::exp(-0.5 * std::pow((t - position) / sigma, 2));
        x.push_back(t);
        y.push_back(value);
    }
    Json result = Json::object();
    result["wavelength_A"] = wavelength;
    result["peaks"] = list;
    result["peak_count"] = list.size();
    result["two_theta_deg"] = toJson(x);
    result["profile"] = toJson(y);
    ToolOutput output;
    output.result = result;
    return output;
}

ToolOutput electronDiffraction(const StructureInput& structure, const Parameters& p)
{
    const Crystal c = crystalFrom(structure);
    std::vector<double> axis = {0, 0, 1};
    if (p.has("zone_axis")) axis = finiteArray(p.array("zone_axis"), "zone_axis", 1).values;
    if (axis.size() != 3) throw std::runtime_error("zone_axis needs three indices [u, v, w]");
    const Vec3 beam = rowTimes({axis[0], axis[1], axis[2]}, c.cell);
    if (norm(beam) < 1e-9) throw std::runtime_error("zone_axis cannot be zero");
    const double voltage = positive(p.number("voltage_kV", 200.0), "voltage_kV") * 1000;
    // Relativistic electron wavelength (Angstrom).
    const double wavelength = 12.2642598 / std::sqrt(voltage * (1 + 0.978475e-6 * voltage));
    const double gMax = positive(p.number("g_max_inv_A", 1.0), "g_max_inv_A");
    const Vec3 direction = scale(beam, 1 / norm(beam));
    // Detector axes: e1 along the shortest allowed reflection, e2 = beam x e1.
    std::vector<std::pair<std::array<int, 3>, Vec3>> zone;
    for (const auto& hkl : reflections(c, gMax)) {
        if (std::abs(hkl[0] * axis[0] + hkl[1] * axis[1] + hkl[2] * axis[2]) > 1e-9) continue;
        zone.push_back({hkl, rowTimes({double(hkl[0]), double(hkl[1]), double(hkl[2])}, c.reciprocal)});
    }
    std::vector<std::tuple<std::array<int, 3>, Vec3, double>> spots;
    double maximum = 0;
    for (const auto& [hkl, g] : zone) {
        const double intensity = std::norm(structureFactor(c, hkl, norm(g) / 2, electronScatteringFactor, 0.0));
        maximum = std::max(maximum, intensity);
        spots.push_back({hkl, g, intensity});
    }
    std::sort(spots.begin(), spots.end(), [](const auto& a, const auto& b) { return norm(std::get<1>(a)) < norm(std::get<1>(b)); });
    Vec3 e1{0, 0, 0};
    for (const auto& spot : spots)
        if (std::get<2>(spot) > 1e-6 * maximum) { e1 = scale(std::get<1>(spot), 1 / norm(std::get<1>(spot))); break; }
    if (norm(e1) == 0) throw std::runtime_error("No allowed reflections in this zone within g_max");
    const Vec3 e2 = cross(direction, e1);
    Json list = Json::array();
    std::vector<double> xs, ys;
    for (const auto& [hkl, g, intensity] : spots) {
        const double scaled = 100 * intensity / maximum;
        if (scaled < 1e-3) continue;
        Json entry = Json::object();
        entry["hkl"] = hklJson(hkl);
        entry["g_inv_A"] = norm(g);
        entry["d_A"] = 1 / norm(g);
        entry["x_inv_A"] = dot(g, e1);
        entry["y_inv_A"] = dot(g, e2);
        entry["intensity"] = scaled;
        entry["angle_mrad"] = 2 * std::asin(wavelength * norm(g) / 2) * 1000;
        list.push(entry);
        xs.push_back(dot(g, e1));
        ys.push_back(dot(g, e2));
    }
    Json result = Json::object();
    result["wavelength_A"] = wavelength;
    result["zone_axis"] = toJson(axis);
    result["spots"] = list;
    result["spot_count"] = list.size();
    result["spot_x_inv_A"] = toJson(xs);
    result["spot_y_inv_A"] = toJson(ys);
    ToolOutput output;
    output.result = result;
    return output;
}
}
