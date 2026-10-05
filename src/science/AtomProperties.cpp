#include "science/AtomProperties.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace atomforge::science
{
namespace
{
const double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<double> numbers(const Json* value)
{
    std::vector<double> result;
    if (!value || !value->isArray()) return result;
    for (const auto& item : value->items())
        result.push_back(item.isNumber() ? item.number() : kNaN);
    return result;
}

double element(const Json& tensor, int i, int j)
{
    if (!tensor.isArray() || tensor.size() != 3) return kNaN;
    const Json& row = tensor.items()[static_cast<std::size_t>(i)];
    if (!row.isArray() || row.size() != 3) return kNaN;
    const Json& value = row.items()[static_cast<std::size_t>(j)];
    return value.isNumber() ? value.number() : kNaN;
}

void add(std::vector<AtomProperty>& properties, const std::string& name, std::vector<double> values)
{
    if (!values.empty()) properties.push_back({name, std::move(values)});
}
}

std::vector<AtomProperty> perAtomProperties(const std::string& tool, const Json& result)
{
    std::vector<AtomProperty> properties;
    if (!result.isObject()) return properties;
    if (tool == "cluster-analysis") {
        add(properties, "Cluster id (-1 unselected)", numbers(result.find("cluster_id")));
    } else if (tool == "void-analysis") {
        add(properties, "Lining void id (-1 none)", numbers(result.find("lining_void_id")));
    } else if (tool == "dislocation-lines") {
        add(properties, "Dislocation line id (-1 none)", numbers(result.find("line_id")));
        add(properties, "Structure type (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)", numbers(result.find("structure_type")));
    } else if (tool == "structure-type") {
        add(properties, "Structure type (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)", numbers(result.find("structure_type")));
    } else if (tool == "centrosymmetry") {
        add(properties, "Centrosymmetry (A^2)", numbers(result.find("centrosymmetry_A2")));
    } else if (tool == "local-strain") {
        add(properties, "D2min (A^2)", numbers(result.find("d2min_A2")));
        std::vector<double> volumetric, shear;
        if (const Json* strain = result.find("green_lagrange_strain"); strain && strain->isArray())
            for (const auto& e : strain->items()) {
                const double xx = element(e, 0, 0), yy = element(e, 1, 1), zz = element(e, 2, 2);
                const double xy = element(e, 0, 1), xz = element(e, 0, 2), yz = element(e, 1, 2);
                volumetric.push_back(xx + yy + zz);
                // Von Mises local shear invariant (Shimizu, Ogata, Li 2007).
                shear.push_back(std::sqrt(xy * xy + xz * xz + yz * yz +
                    ((xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) + (xx - zz) * (xx - zz)) / 6));
            }
        add(properties, "Von Mises shear strain", shear);
        add(properties, "Volumetric strain trace(E)", volumetric);
        add(properties, "Coordination", numbers(result.find("coordination")));
    } else if (tool == "bond-order") {
        if (const Json* order = result.find("order"); order && order->isObject())
            for (const auto& [name, values] : order->members()) add(properties, name, numbers(&values));
        add(properties, "Coordination", numbers(result.find("coordination")));
    } else if (tool == "wigner-seitz") {
        add(properties, "Distance to assigned site (A)", numbers(result.find("distance_A")));
        add(properties, "Assigned site index", numbers(result.find("site_index")));
    }
    return properties;
}

std::array<float, 3> viridis(double t)
{
    if (!std::isfinite(t)) t = 0;
    t = std::clamp(t, 0.0, 1.0);
    const double r = 0.2777273272234177 + t * (0.1050930431085774 + t * (-0.3308618287255563 + t * (-4.634230498983486 + t * (6.228269936347081 + t * (4.776384997670288 + t * -5.435455855934631)))));
    const double g = 0.005407344544966578 + t * (1.404613529898575 + t * (0.214847559468213 + t * (-5.799100973351585 + t * (14.17993336680509 + t * (-13.74514537774601 + t * 4.645852612178535)))));
    const double b = 0.3340998053353061 + t * (1.384590162594685 + t * (0.09509516302823659 + t * (-19.33244095627987 + t * (56.69055260068105 + t * (-65.35303263337234 + t * 26.3124352495832)))));
    return {static_cast<float>(std::clamp(r, 0.0, 1.0)), static_cast<float>(std::clamp(g, 0.0, 1.0)), static_cast<float>(std::clamp(b, 0.0, 1.0))};
}

PropertyRange finiteRange(const std::vector<double>& values)
{
    PropertyRange range{HUGE_VAL, -HUGE_VAL};
    for (double v : values)
        if (std::isfinite(v)) { range.low = std::min(range.low, v); range.high = std::max(range.high, v); }
    if (!(range.low <= range.high)) return {0.0, 1.0};
    return range;
}

void colourByProperty(const std::vector<double>& values, const PropertyDisplay& display,
                      std::vector<std::array<float, 3>>& colours, std::vector<bool>& visible)
{
    const PropertyRange range = display.autoRange ? finiteRange(values) : display.range;
    const double span = range.high - range.low;
    colours.assign(values.size(), {0.55f, 0.55f, 0.55f});
    visible.assign(values.size(), true);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double v = values[i];
        if (!std::isfinite(v)) { visible[i] = !display.hideInvalid; continue; }
        colours[i] = viridis(span > 0 ? (v - range.low) / span : 0.5);
        if (display.hideOutside && (v < range.low || v > range.high)) visible[i] = false;
    }
}
}
