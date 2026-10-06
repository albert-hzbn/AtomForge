#include "science/ResultPlots.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>

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

std::vector<double> column(const Json* rows, std::size_t k)
{
    std::vector<double> result;
    if (!rows || !rows->isArray()) return result;
    for (const auto& row : rows->items())
        result.push_back(row.isArray() && row.size() > k && row.items()[k].isNumber() ? row.items()[k].number() : kNaN);
    return result;
}

std::vector<double> indices(std::size_t count)
{
    std::vector<double> result(count);
    std::iota(result.begin(), result.end(), 0.0);
    return result;
}

PlotSeries line(const std::string& name, std::vector<double> x, std::vector<double> y, bool points = false)
{
    PlotSeries series;
    series.name = name;
    series.x = std::move(x);
    series.y = std::move(y);
    series.points = points;
    return series;
}

PlotSpec plot(const std::string& title, const std::string& xLabel, const std::string& yLabel, std::vector<PlotSeries> series)
{
    PlotSpec spec;
    spec.title = title;
    spec.xLabel = xLabel;
    spec.yLabel = yLabel;
    spec.series = std::move(series);
    return spec;
}

double number(const Json& result, const char* key, double fallback = kNaN)
{
    const Json* value = result.find(key);
    return value && value->isNumber() ? value->number() : fallback;
}

void addHistogram(std::vector<PlotSpec>& plots, const std::string& title, const std::string& xLabel, const std::vector<double>& values)
{
    bool any = false;
    for (double v : values) any = any || std::isfinite(v);
    if (any) plots.push_back(plot(title, xLabel, "Atoms", {histogram(xLabel, values, 40)}));
}

using PlotBuilder = void (*)(const Json& result, std::vector<PlotSpec>& plots);

void msdPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto lag = numbers(result.find("lag_fs"));
    const Json* components = result.find("components_A2");
    plots.push_back(plot("Mean-square displacement", "Lag time (fs)", "MSD (A^2)",
        {line("total", lag, numbers(result.find("msd_A2"))), line("x", lag, column(components, 0)),
         line("y", lag, column(components, 1)), line("z", lag, column(components, 2))}));
}

void diffusionPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto lag = numbers(result.find("lag_fs"));
    const auto range = numbers(result.find("fit_range_fs"));
    if (!lag.empty() && range.size() == 2) {
        const double slope = number(result, "slope_A2_per_fs"), intercept = number(result, "intercept_A2");
        plots.push_back(plot("Einstein diffusion fit", "Lag time (fs)", "MSD (A^2)",
            {line("MSD", lag, numbers(result.find("msd_A2"))),
             line("linear fit", range, {intercept + slope * range[0], intercept + slope * range[1]})}));
    }
}

void vacfPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    auto spec = plot("Velocity autocorrelation", "Lag time (fs)", "VACF",
        {line("VACF", numbers(result.find("lag_fs")), numbers(result.find("vacf")))});
    spec.horizontal = {0.0};
    plots.push_back(spec);
}

void vibrationalSpectrumPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    plots.push_back(plot("Vibrational spectrum", "Frequency (THz)", "Spectral density (1/THz)",
        {line("density", numbers(result.find("frequency_THz")), numbers(result.find("density_per_THz")))}));
}

void structureFactorPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const Json* q = result.find("q_vectors_rad_per_A");
    const auto s = numbers(result.find("S_q"));
    std::vector<std::pair<double, double>> pairs;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const double qx = column(q, 0)[i], qy = column(q, 1)[i], qz = column(q, 2)[i];
        pairs.push_back({std::sqrt(qx * qx + qy * qy + qz * qz), s[i]});
    }
    std::sort(pairs.begin(), pairs.end());
    std::vector<double> x, y;
    for (const auto& [qq, ss] : pairs) { x.push_back(qq); y.push_back(ss); }
    plots.push_back(plot("Static structure factor", "|q| (rad/A)", "S(q)", {line("S(q)", x, y, true)}));
}

void localStrainPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    addHistogram(plots, "Non-affine displacement", "D2min (A^2)", numbers(result.find("d2min_A2")));
    std::vector<double> volumetric;
    if (const Json* strain = result.find("green_lagrange_strain"); strain && strain->isArray())
        for (const auto& tensor : strain->items()) {
            double trace = 0;
            for (int k = 0; k < 3; ++k) {
                const auto row = tensor.isArray() && tensor.size() == 3 ? tensor.items()[static_cast<std::size_t>(k)] : Json();
                trace += row.isArray() && row.items()[static_cast<std::size_t>(k)].isNumber() ? row.items()[static_cast<std::size_t>(k)].number() : kNaN;
            }
            volumetric.push_back(trace);
        }
    addHistogram(plots, "Volumetric strain", "trace(E)", volumetric);
}

void vaspElectronicPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto distance = numbers(result.find("distance"));
    const Json* energies = result.find("energies_minus_fermi_eV");
    if (energies && energies->isArray() && distance.size() > 1) {
        std::vector<PlotSeries> series;
        for (std::size_t s = 0; s < energies->size(); ++s) {
            const Json& spin = energies->items()[s];
            const std::size_t count = spin.items()[0].size();
            for (std::size_t b = 0; b < count; ++b)
                series.push_back(line((energies->size() > 1 ? (s == 0 ? "up " : "down ") : "band ") + std::to_string(b + 1), distance, column(&spin, b)));
        }
        auto spec = plot("Band structure", "Wavevector distance", "E - E_F (eV)", series);
        if (const Json* labels = result.find("labels"); labels && labels->isArray())
            for (const auto& marker : labels->items()) spec.markers.push_back({marker.at("distance").number(), marker.at("label").string()});
        spec.horizontal = {0.0};
        plots.push_back(spec);
        if (const Json* fat = result.find("fat_band_weights"); fat && fat->isObject() && fat->size() > 0) {
            const auto& [element, weights] = fat->members()[0];
            std::vector<double> x, y;
            const Json& spin = energies->items()[0];
            for (std::size_t k = 0; k < distance.size(); ++k)
                for (std::size_t b = 0; b < spin.items()[k].size(); ++b)
                    if (weights.items()[0].items()[k].items()[b].number() > 0.3) {
                        x.push_back(distance[k]);
                        y.push_back(spin.items()[k].items()[b].number());
                    }
            auto fatSpec = spec;
            fatSpec.title = "Fat bands: states with >30% " + element + " character";
            fatSpec.series.push_back(line(element + " > 30%", x, y, true));
            plots.push_back(fatSpec);
        }
    }
    if (const Json* dos = result.find("dos"); dos && dos->isObject()) {
        const auto energy = numbers(dos->find("energy_minus_fermi_eV"));
        std::vector<PlotSeries> series;
        const Json* total = dos->find("total");
        for (std::size_t s = 0; total && s < total->size(); ++s) series.push_back(line(s ? "total (down)" : "total", energy, numbers(&total->items()[s])));
        if (const Json* projected = dos->find("projected_by_element"); projected && projected->isObject())
            for (const auto& [element, data] : projected->members()) series.push_back(line(element, energy, numbers(&data.items()[0])));
        auto spec = plot("Density of states", "E - E_F (eV)", "DOS (states/eV)", series);
        spec.markers.push_back({0.0, "E_F"});
        plots.push_back(spec);
    }
}

void powderXrdPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    auto spec = plot("Powder diffraction pattern", "2-theta (degrees)", "Intensity", {line("profile", numbers(result.find("two_theta_deg")), numbers(result.find("profile")))});
    if (const Json* peaks = result.find("peaks"); peaks && peaks->isArray())
        for (const auto& peak : peaks->items())
            if (peak.at("intensity").number() >= 5) {
                const auto& h = peak.at("hkl").items();
                spec.markers.push_back({peak.at("two_theta_deg").number(),
                    std::to_string(static_cast<int>(h[0].number())) + std::to_string(static_cast<int>(h[1].number())) + std::to_string(static_cast<int>(h[2].number()))});
            }
    plots.push_back(spec);
}

void electronDiffractionPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    plots.push_back(plot("Diffraction spots", "x (1/A)", "y (1/A)", {line("spots", numbers(result.find("spot_x_inv_A")), numbers(result.find("spot_y_inv_A")), true)}));
}

void clusterAnalysisPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const Json* histogram = result.find("size_histogram");
    if (histogram && histogram->isArray() && histogram->size() > 0) {
        std::vector<double> x, y;
        for (const auto& row : histogram->items()) { x.push_back(row.items()[0].number()); y.push_back(row.items()[1].number()); }
        if (x.size() == 1) { x.push_back(x[0] + 1); y.push_back(0); }
        plots.push_back(plot("Cluster size distribution", "Cluster size (atoms)", "Clusters", {line("clusters", x, y, true)}));
    }
}

void voidAnalysisPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    std::vector<double> rank, volume;
    if (const Json* voids = result.find("voids"); voids && voids->isArray())
        for (const auto& v : voids->items()) { rank.push_back(static_cast<double>(rank.size() + 1)); volume.push_back(v.at("accessible_volume_A3").number()); }
    if (rank.size() >= 2) plots.push_back(plot("Void volumes", "Void (largest first)", "Accessible volume (A^3)", {line("voids", rank, volume, true)}));
}

void structureTypePlots(const Json& result, std::vector<PlotSpec>& plots)
{
    if (const Json* counts = result.find("counts"); counts && counts->isObject()) {
        std::vector<double> x, y;
        for (std::size_t t = 0; t < counts->size(); ++t) {
            x.push_back(static_cast<double>(t));
            y.push_back(counts->members()[t].second.number());
        }
        auto spec = plot("Structure types (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)", "Type", "Atoms", {line("atoms", x, y, true)});
        plots.push_back(spec);
    }
}

void centrosymmetryPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    addHistogram(plots, "Centrosymmetry distribution", "CSP (A^2)", numbers(result.find("centrosymmetry_A2")));
}

void bondOrderPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    std::vector<PlotSeries> series;
    if (const Json* order = result.find("order"); order && order->isObject())
        for (const auto& [name, values] : order->members()) {
            auto h = histogram(name, numbers(&values), 40);
            if (!h.x.empty()) series.push_back(h);
        }
    if (!series.empty()) plots.push_back(plot("Steinhardt order distribution", "q_l", "Atoms", series));
}

void wignerSeitzPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    addHistogram(plots, "Distance to assigned site", "Distance (A)", numbers(result.find("distance_A")));
}

void workFunctionPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto distance = numbers(result.find("distance_A"));
    if (!distance.empty()) {
        auto spec = plot("Planar potential", "Distance (A)", "Potential energy (eV)",
            {line("potential", distance, numbers(result.find("potential_eV")))});
        spec.horizontal = {number(result, "vacuum_level_eV"), number(result, "fermi_eV")};
        plots.push_back(spec);
    }
}

void equationOfStatePlots(const Json& result, std::vector<PlotSpec>& plots)
{
    plots.push_back(plot("Birch-Murnaghan equation of state", "Volume (A^3)", "Energy (eV)",
        {line("data", numbers(result.find("volumes_A3")), numbers(result.find("energies_eV")), true),
         line("fit", numbers(result.find("fit_volumes_A3")), numbers(result.find("fit_energies_eV")))}));
}

void phononDosPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    plots.push_back(plot("Phonon density of states", "Energy (eV)", "DOS (1/eV)",
        {line("DOS", numbers(result.find("energy_eV")), numbers(result.find("dos_per_eV")))}));
}

void harmonicThermodynamicsPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto t = numbers(result.find("temperature_K"));
    plots.push_back(plot("Harmonic free and internal energy", "Temperature (K)", "Energy (eV/cell)",
        {line("F", t, numbers(result.find("free_energy_eV"))), line("U", t, numbers(result.find("internal_energy_eV")))}));
    plots.push_back(plot("Entropy and heat capacity", "Temperature (K)", "eV/(cell K)",
        {line("S", t, numbers(result.find("entropy_eV_per_K"))), line("Cv", t, numbers(result.find("heat_capacity_eV_per_K")))}));
}

void phononsPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const Json* bands = result.find("frequencies_THz");
    const auto distance = numbers(result.find("distance_inv_A"));
    if (bands && bands->isArray() && bands->size() == distance.size() && bands->size() > 1) {
        std::vector<PlotSeries> branches;
        const std::size_t count = bands->items()[0].size();
        for (std::size_t b = 0; b < count; ++b)
            branches.push_back(line("branch " + std::to_string(b + 1), distance, column(bands, b)));
        auto spec = plot("Phonon dispersion", "Wavevector distance (1/A)", "Frequency (THz)", branches);
        const Json* labels = result.find("labels");
        for (std::size_t i = 0; labels && labels->isArray() && i < labels->size() && i < distance.size(); ++i)
            if (labels->items()[i].isString() && !labels->items()[i].string().empty()) {
                std::string label = labels->items()[i].string();
                if (label == "GAMMA") label = "G";
                spec.markers.push_back({distance[i], label});
            }
        spec.horizontal = {0.0};
        plots.push_back(spec);
    }
    if (const Json* dos = result.find("dos"); dos && dos->isObject()) {
        auto energy = numbers(dos->find("energy_eV"));
        auto density = numbers(dos->find("dos_per_eV"));
        for (double& e : energy) e /= 4.135667696e-3;
        for (double& d : density) d *= 4.135667696e-3;
        plots.push_back(plot("Phonon density of states", "Frequency (THz)", "DOS (states/THz)", {line("DOS", energy, density)}));
    }
    if (const Json* thermo = result.find("thermodynamics"); thermo && thermo->isObject()) {
        const auto t = numbers(thermo->find("temperature_K"));
        plots.push_back(plot("Harmonic free and internal energy", "Temperature (K)", "Energy (eV/cell)",
            {line("F", t, numbers(thermo->find("free_energy_eV"))), line("U", t, numbers(thermo->find("internal_energy_eV")))}));
        plots.push_back(plot("Heat capacity", "Temperature (K)", "Cv (eV/(cell K))", {line("Cv", t, numbers(thermo->find("heat_capacity_eV_per_K")))}));
    }
}

void relaxPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto enthalpy = numbers(result.find("enthalpy_history_eV"));
    const auto force = numbers(result.find("max_force_history"));
    plots.push_back(plot("Relaxation energy", "Step", "E (+PV) (eV)", {line("energy", indices(enthalpy.size()), enthalpy)}));
    auto spec = plot("Convergence", "Step", "Max force (eV/A)", {line("max force", indices(force.size()), force)});
    spec.logY = true;
    plots.push_back(spec);
}

void nebPlots(const Json& result, std::vector<PlotSpec>& plots)
{
    auto energies = numbers(result.find("energies_eV"));
    if (!energies.empty()) {
        const double start = energies.front();
        for (double& e : energies) e -= start;
        auto x = numbers(result.find("reaction_coordinate_A"));
        if (x.size() != energies.size()) x = indices(energies.size());
        auto spec = plot("Minimum-energy path", "Reaction coordinate (A)", "E - E(initial) (eV)",
            {line("path", x, energies), line("images", x, energies, true)});
        spec.horizontal = {0.0};
        plots.push_back(spec);
    }
}

void dynamicsPlots(const Json& result, std::vector<PlotSpec>& plots, bool npt)
{
    const auto time = numbers(result.find("time_fs"));
    auto temperature = plot("Temperature", "Time (fs)", "T (K)", {line("T", time, numbers(result.find("temperature_K")))});
    temperature.horizontal = {number(result, "target_temperature_K")};
    auto energy = plot("Total energy", "Time (fs)", "E (eV)", {line("E", time, numbers(result.find("total_energy_eV")))});
    if (const Json* start = result.find("production_start_fs"); start && start->isNumber()) {
        temperature.markers.push_back({start->number(), "NVE"});
        energy.markers.push_back({start->number(), "NVE"});
    }
    plots.push_back(temperature);
    plots.push_back(energy);
    if (npt) {
        plots.push_back(plot("Volume", "Time (fs)", "V (A^3)", {line("V", time, numbers(result.find("volume_A3")))}));
        plots.push_back(plot("Pressure", "Time (fs)", "P (GPa)", {line("P", time, numbers(result.find("pressure_GPa")))}));
    }
}

void trajectoryStructurePlots(const Json& result, std::vector<PlotSpec>& plots)
{
    const auto r = numbers(result.find("r_A"));
    std::vector<PlotSeries> rdf = {line("total", r, numbers(result.find("g_total")))};
    if (const Json* partial = result.find("g_partial"); partial && partial->isObject() && partial->size() > 1)
        for (const auto& [pair, values] : partial->members()) rdf.push_back(line(pair, r, numbers(&values)));
    auto spec = plot("Radial distribution function", "r (A)", "g(r)", rdf);
    spec.horizontal = {1.0};
    plots.push_back(spec);
    if (const Json* distribution = result.find("coordination_distribution"); distribution && distribution->isObject()) {
        std::vector<PlotSeries> series;
        for (const auto& [element, rows] : distribution->members()) {
            std::vector<double> x, y;
            for (const auto& row : rows.items()) { x.push_back(row.items()[0].number()); y.push_back(row.items()[1].number()); }
            if (x.size() == 1) { x.push_back(x[0] + 1); y.push_back(0); }
            series.push_back(line(element, x, y, true));
        }
        plots.push_back(plot("Coordination distribution", "Neighbours within the cutoff", "Fraction of atoms", series));
    }
    plots.push_back(plot("Bond-angle distribution", "Angle (degrees)", "Density (1/degree)",
        {line("angles", numbers(result.find("angle_deg")), numbers(result.find("angle_density_per_deg")))}));
}

const std::map<std::string, PlotBuilder>& plotBuilders()
{
    static const std::map<std::string, PlotBuilder> builders = {
        {"msd", msdPlots},
        {"diffusion", diffusionPlots},
        {"vacf", vacfPlots},
        {"vibrational-spectrum", vibrationalSpectrumPlots},
        {"structure-factor", structureFactorPlots},
        {"local-strain", localStrainPlots},
        {"vasp-electronic", vaspElectronicPlots},
        {"powder-xrd", powderXrdPlots},
        {"electron-diffraction", electronDiffractionPlots},
        {"cluster-analysis", clusterAnalysisPlots},
        {"void-analysis", voidAnalysisPlots},
        {"structure-type", structureTypePlots},
        {"centrosymmetry", centrosymmetryPlots},
        {"bond-order", bondOrderPlots},
        {"wigner-seitz", wignerSeitzPlots},
        {"work-function", workFunctionPlots},
        {"equation-of-state", equationOfStatePlots},
        {"eos-scan", equationOfStatePlots},
        {"trajectory-structure", trajectoryStructurePlots},
        {"phonon-dos", phononDosPlots},
        {"harmonic-thermodynamics", harmonicThermodynamicsPlots},
        {"phonons", phononsPlots},
        {"relax", relaxPlots},
        {"neb", nebPlots},
        {"nvt", [](const Json& result, std::vector<PlotSpec>& plots) { dynamicsPlots(result, plots, false); }},
        {"npt", [](const Json& result, std::vector<PlotSpec>& plots) { dynamicsPlots(result, plots, true); }},
    };
    return builders;
}
}

PlotSeries histogram(const std::string& name, const std::vector<double>& values, int bins)
{
    PlotSeries series;
    series.name = name;
    double low = HUGE_VAL, high = -HUGE_VAL;
    for (double v : values)
        if (std::isfinite(v)) { low = std::min(low, v); high = std::max(high, v); }
    if (!(low <= high) || bins < 1) return series;
    if (high == low) { high = low + 1e-12 + std::abs(low) * 1e-9; }
    const double width = (high - low) / bins;
    std::vector<double> counts(static_cast<std::size_t>(bins), 0.0);
    for (double v : values) {
        if (!std::isfinite(v)) continue;
        const int b = std::min(bins - 1, static_cast<int>((v - low) / width));
        counts[static_cast<std::size_t>(b)] += 1;
    }
    for (int b = 0; b < bins; ++b) {
        series.x.push_back(low + (b + 0.5) * width);
        series.y.push_back(counts[static_cast<std::size_t>(b)]);
    }
    return series;
}

std::vector<PlotSpec> resultPlots(const std::string& tool, const Json& result)
{
    std::vector<PlotSpec> plots;
    if (!result.isObject()) return plots;
    if (const auto builder = plotBuilders().find(tool); builder != plotBuilders().end()) builder->second(result, plots);
    // Drop plots without at least two finite points.
    plots.erase(std::remove_if(plots.begin(), plots.end(), [](const PlotSpec& spec) {
        std::size_t finite = 0;
        for (const auto& series : spec.series)
            for (std::size_t i = 0; i < std::min(series.x.size(), series.y.size()); ++i)
                finite += std::isfinite(series.x[i]) && std::isfinite(series.y[i]);
        return finite < 2;
    }), plots.end());
    return plots;
}

PlotSpec overlayPlots(const PlotSpec& current, const std::string& currentLabel,
                      const std::vector<std::pair<std::string, std::vector<PlotSpec>>>& others)
{
    PlotSpec combined = current;
    combined.series.clear();
    const auto addRun = [&](const std::string& label, const PlotSpec& spec) {
        for (auto series : spec.series) {
            series.name = label + ": " + series.name;
            combined.series.push_back(std::move(series));
        }
    };
    for (const auto& [label, plots] : others)
        for (const auto& spec : plots)
            if (spec.title == current.title) { addRun(label, spec); break; }
    addRun(currentLabel, current);
    return combined;
}

std::string plotCsv(const PlotSpec& plot)
{
    std::string text;
    std::size_t rows = 0;
    for (std::size_t s = 0; s < plot.series.size(); ++s) {
        const std::string name = plot.series[s].name;
        text += (s ? "," : "") + plot.xLabel + " [" + name + "]," + plot.yLabel + " [" + name + "]";
        rows = std::max(rows, plot.series[s].x.size());
    }
    text += "\n";
    char buffer[64];
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t s = 0; s < plot.series.size(); ++s) {
            const auto& series = plot.series[s];
            if (s) text += ",";
            if (r < series.x.size() && r < series.y.size()) {
                std::snprintf(buffer, sizeof(buffer), "%.10g,%.10g", series.x[r], series.y[r]);
                text += buffer;
            } else text += ",";
        }
        text += "\n";
    }
    return text;
}

std::vector<std::string> plottedTools()
{
    std::vector<std::string> tools;
    for (const auto& entry : plotBuilders()) tools.push_back(entry.first);
    return tools;
}
}
