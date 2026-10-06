#include "science/TrajectoryStructure.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace atomforge::science
{
Json trajectoryStructure(const Parameters& p)
{
    if (!p.has("trajectory_file")) throw std::runtime_error("trajectory-structure needs trajectory_file");
    const TrajectoryStream stream(std::filesystem::u8path(p.json("trajectory_file").string()));
    const long long first = integer(p.number("first_frame", 0), "first_frame", 0);
    const long long every = integer(p.number("every", 1), "every");
    const double rMax = positive(p.number("r_max_A", 8.0), "r_max_A");
    const long long bins = integer(p.number("bins", 200), "bins", 10);
    const double cutoff = positive(p.number("cutoff_A", 3.0), "cutoff_A");
    const long long angleBins = integer(p.number("angle_bins", 180), "angle_bins", 10);
    if (static_cast<std::size_t>(first) >= stream.size()) throw std::runtime_error("first_frame is beyond the trajectory");
    const double dr = rMax / static_cast<double>(bins);
    const double pi = std::acos(-1.0);

    std::vector<std::string> elements;
    // pairCounts[a][b][bin]: b atoms in each shell around a atoms, summed over frames.
    std::map<std::string, std::map<std::string, std::vector<double>>> pairCounts;
    std::map<std::string, std::map<std::string, double>> pairNormalization;  // sum over frames of N_a (N_b - delta_ab) / V
    std::map<std::string, std::map<std::string, double>> coordinationSum;     // b neighbours around a, summed over a atoms and frames
    std::map<std::string, double> centreCount;                                // a atoms, summed over frames
    std::map<std::string, std::map<long long, double>> coordinationHistogram; // per element: coordination -> count
    std::vector<double> angles(static_cast<std::size_t>(angleBins), 0.0);
    double angleTotal = 0;
    long long used = 0;
    for (std::size_t index = static_cast<std::size_t>(first); index < stream.size(); index += static_cast<std::size_t>(every)) {
        taskProgress(static_cast<double>(index) / static_cast<double>(stream.size()));
        const FrameData frame = stream.frame(index);
        const Structure& s = frame.structure;
        if (!s.hasUnitCell || !(frame.pbc[0] && frame.pbc[1] && frame.pbc[2]))
            throw std::runtime_error("Trajectory-averaged g(r) needs periodic frames (frame " + std::to_string(index) + ")");
        Mat3 cell{};
        for (int r = 0; r < 3; ++r) cell[r] = s.cellVectors[r];
        const double volume = cellVolume(cell);
        std::vector<Vec3> positions;
        std::map<std::string, double> count;
        for (const auto& atom : s.atoms) {
            positions.push_back({atom.x, atom.y, atom.z});
            count[atom.symbol] += 1;
            if (std::find(elements.begin(), elements.end(), atom.symbol) == elements.end()) elements.push_back(atom.symbol);
        }
        for (const auto& [a, na] : count)
            for (const auto& [b, nb] : count) pairNormalization[a][b] += na * (nb - (a == b ? 1.0 : 0.0)) / volume;
        for (const auto& [a, na] : count) centreCount[a] += na;
        // Neighbour list out to the larger of r_max and the coordination cutoff.
        std::vector<std::vector<Vec3>> bonds(s.atoms.size());
        std::vector<std::map<std::string, long long>> coordination(s.atoms.size());
        for (const auto& n : neighborList(positions, cell, frame.pbc, std::max(rMax, cutoff))) {
            const double r = norm(n.vector);
            const std::string& a = s.atoms[static_cast<std::size_t>(n.i)].symbol;
            const std::string& b = s.atoms[static_cast<std::size_t>(n.j)].symbol;
            if (r < rMax) {
                auto& histogram = pairCounts[a][b];
                if (histogram.empty()) histogram.assign(static_cast<std::size_t>(bins), 0.0);
                histogram[std::min(static_cast<std::size_t>(r / dr), static_cast<std::size_t>(bins) - 1)] += 1;
            }
            if (r < cutoff) {
                bonds[static_cast<std::size_t>(n.i)].push_back(n.vector);
                ++coordination[static_cast<std::size_t>(n.i)][b];
            }
        }
        for (std::size_t i = 0; i < s.atoms.size(); ++i) {
            const std::string& a = s.atoms[i].symbol;
            long long total = 0;
            for (const auto& [b, c] : coordination[i]) { coordinationSum[a][b] += static_cast<double>(c); total += c; }
            coordinationHistogram[a][total] += 1;
            const auto& v = bonds[i];
            for (std::size_t j = 0; j < v.size(); ++j)
                for (std::size_t k = j + 1; k < v.size(); ++k) {
                    const double cosine = std::clamp(dot(v[j], v[k]) / (norm(v[j]) * norm(v[k])), -1.0, 1.0);
                    const double degrees = std::acos(cosine) * 180.0 / pi;
                    angles[std::min(static_cast<std::size_t>(degrees / 180.0 * static_cast<double>(angleBins)), static_cast<std::size_t>(angleBins) - 1)] += 1;
                    angleTotal += 1;
                }
        }
        ++used;
    }
    std::sort(elements.begin(), elements.end());
    std::vector<double> r(static_cast<std::size_t>(bins)), shell(static_cast<std::size_t>(bins));
    for (long long k = 0; k < bins; ++k) {
        const double lo = static_cast<double>(k) * dr, hi = lo + dr;
        r[static_cast<std::size_t>(k)] = lo + 0.5 * dr;
        shell[static_cast<std::size_t>(k)] = 4.0 / 3.0 * pi * (hi * hi * hi - lo * lo * lo);
    }
    // g_ab(r) = <counts of b around a> / (N_a rho_b shell volume), averaged over frames;
    // the total g(r) weights pairs by their number.
    Json partial = Json::object();
    std::vector<double> totalCounts(static_cast<std::size_t>(bins), 0.0);
    double totalNormalization = 0;
    for (const auto& a : elements)
        for (const auto& b : elements) {
            const double normalization = pairNormalization[a][b];
            std::vector<double> g(static_cast<std::size_t>(bins), 0.0);
            const auto found = pairCounts[a].find(b);
            for (long long k = 0; k < bins; ++k) {
                const double counts = found == pairCounts[a].end() ? 0.0 : found->second[static_cast<std::size_t>(k)];
                if (normalization > 0) g[static_cast<std::size_t>(k)] = counts / (normalization * shell[static_cast<std::size_t>(k)]);
                totalCounts[static_cast<std::size_t>(k)] += counts;
            }
            totalNormalization += normalization;
            partial[a + "-" + b] = toJson(g);
        }
    std::vector<double> total(static_cast<std::size_t>(bins), 0.0);
    for (long long k = 0; k < bins; ++k)
        if (totalNormalization > 0) total[static_cast<std::size_t>(k)] = totalCounts[static_cast<std::size_t>(k)] / (totalNormalization * shell[static_cast<std::size_t>(k)]);
    Json coordinationJson = Json::object(), pairCoordination = Json::object(), distribution = Json::object();
    double allNeighbours = 0, allCentres = 0;
    for (const auto& a : elements) {
        double sum = 0;
        for (const auto& b : elements) {
            const double value = coordinationSum[a][b] / centreCount[a];
            pairCoordination[b + " around " + a] = value;
            sum += coordinationSum[a][b];
        }
        coordinationJson[a] = sum / centreCount[a];
        allNeighbours += sum;
        allCentres += centreCount[a];
        Json rows = Json::array();
        for (const auto& [n, c] : coordinationHistogram[a]) rows.push(Json::array({static_cast<double>(n), c / centreCount[a]}));
        distribution[a] = rows;
    }
    std::vector<double> angleGrid(static_cast<std::size_t>(angleBins)), angleDensity(static_cast<std::size_t>(angleBins), 0.0);
    const double width = 180.0 / static_cast<double>(angleBins);
    for (long long k = 0; k < angleBins; ++k) {
        angleGrid[static_cast<std::size_t>(k)] = (static_cast<double>(k) + 0.5) * width;
        if (angleTotal > 0) angleDensity[static_cast<std::size_t>(k)] = angles[static_cast<std::size_t>(k)] / (angleTotal * width);
    }
    Json result = Json::object();
    result["frames_used"] = used;
    result["frame_count"] = static_cast<long long>(stream.size());
    result["elements"] = [&] { Json list = Json::array(); for (const auto& e : elements) list.push(e); return list; }();
    result["r_A"] = toJson(r);
    result["g_total"] = toJson(total);
    result["g_partial"] = partial;
    result["cutoff_A"] = cutoff;
    result["mean_coordination"] = allCentres > 0 ? allNeighbours / allCentres : 0.0;
    result["coordination"] = coordinationJson;
    result["pair_coordination"] = pairCoordination;
    result["coordination_distribution"] = distribution;
    result["angle_deg"] = toJson(angleGrid);
    result["angle_density_per_deg"] = toJson(angleDensity);
    result["bond_angles_counted"] = angleTotal;
    return result;
}
}
