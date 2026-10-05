#include "science/DislocationLines.h"
#include "science/Analysis.h"
#include "util/TaskControl.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);

std::size_t find(std::vector<std::size_t>& parent, std::size_t i)
{
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}

std::vector<Vec3> positionsOf(const Structure& s)
{
    std::vector<Vec3> result;
    for (const auto& a : s.atoms) result.push_back({a.x, a.y, a.z});
    return result;
}

std::vector<int> structureTypes(const std::vector<Vec3>& positions, const Mat3& cell, const Pbc& pbc, double cutoff)
{
    const auto neighbors = neighborList(positions, cell, pbc, cutoff);
    std::vector<std::vector<Vec3>> bonds(positions.size());
    for (const auto& n : neighbors) bonds[static_cast<std::size_t>(n.i)].push_back(n.vector);
    std::vector<int> types(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) types[i] = acklandJonesType(bonds[i]);
    return types;
}

// Removes the lattice translation closest to d (greedy descent over the
// nearest-neighbour vectors, which generate fcc/bcc/hcp-basal lattices).
Vec3 reduceToLattice(Vec3 d, const std::vector<Vec3>& lattice)
{
    for (int iteration = 0; iteration < 64; ++iteration) {
        const Vec3* best = nullptr;
        double bestNorm = norm(d) - 1e-9;
        for (const auto& v : lattice) {
            const double candidate = norm(sub(d, v));
            if (candidate < bestNorm) { bestNorm = candidate; best = &v; }
        }
        if (!best) break;
        d = sub(d, *best);
    }
    return d;
}
}

ToolOutput dislocationLines(const StructureInput& reference, const StructureInput& deformed, const Parameters& p)
{
    const Structure& ref = reference.structure;
    const Structure& def = deformed.structure;
    if (ref.atoms.size() != def.atoms.size() || ref.atoms.empty())
        throw std::runtime_error("Reference and deformed structures need the same atoms in the same order");
    const bool periodic = ref.hasUnitCell && reference.pbc[0] && reference.pbc[1] && reference.pbc[2];
    Mat3 cell = identity();
    if (ref.hasUnitCell) for (int r = 0; r < 3; ++r) cell[r] = ref.cellVectors[r];
    const Pbc pbc = {periodic, periodic, periodic};
    const MinimumImage mic = periodic ? MinimumImage(cell, pbc) : MinimumImage();
    const auto refPositions = positionsOf(ref), positions = positionsOf(def);
    const std::size_t n = positions.size();

    // Nearest-neighbour spacing and lattice vectors from the reference crystal.
    std::vector<double> nearest(n, HUGE_VAL);
    const auto probe = neighborList(refPositions, cell, pbc, 6.0);
    for (const auto& nb : probe) nearest[static_cast<std::size_t>(nb.i)] = std::min(nearest[static_cast<std::size_t>(nb.i)], norm(nb.vector));
    std::vector<double> sorted;
    for (double d : nearest) if (std::isfinite(d)) sorted.push_back(d);
    if (sorted.empty()) throw std::runtime_error("No neighbours within 6 Angstrom");
    std::sort(sorted.begin(), sorted.end());
    const double spacing = sorted[sorted.size() / 2];
    double cutoff = p.number("cutoff_A", 0.0);
    if (cutoff <= 0) cutoff = 1.6 * spacing;  // beyond the second bcc shell (1.155 d)
    const auto referenceTypes = structureTypes(refPositions, cell, pbc, cutoff);
    taskProgress(0.3);
    const auto deformedTypes = structureTypes(positions, def.hasUnitCell ? cell : identity(), pbc, cutoff);
    taskProgress(0.6);
    std::vector<Vec3> lattice;
    for (std::size_t i = 0; i < n && lattice.empty(); ++i) {
        if (referenceTypes[i] == 0) continue;
        for (const auto& nb : probe)
            if (static_cast<std::size_t>(nb.i) == i && norm(nb.vector) < 1.2 * spacing) lattice.push_back(nb.vector);
    }
    if (lattice.empty()) throw std::runtime_error("The reference has no crystalline (fcc, hcp, bcc) atoms");

    // Defect atoms: crystalline in the reference, changed structure type now.
    std::vector<std::size_t> defect;
    for (std::size_t i = 0; i < n; ++i)
        if (referenceTypes[i] != 0 && deformedTypes[i] != referenceTypes[i]) defect.push_back(i);
    if (defect.empty()) throw std::runtime_error("No atoms changed crystal structure: no dislocation cores found");
    const double clusterCutoff = p.has("cluster_cutoff_A") && p.number("cluster_cutoff_A") > 0 ? p.number("cluster_cutoff_A") : 1.5 * spacing;
    const long long minimumAtoms = integer(p.number("min_atoms", 4), "min_atoms", 1);
    std::vector<Vec3> defectPositions;
    for (std::size_t i : defect) defectPositions.push_back(refPositions[i]);
    const auto links = neighborList(defectPositions, cell, pbc, clusterCutoff);
    std::vector<std::size_t> parent(defect.size());
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& link : links) parent[find(parent, static_cast<std::size_t>(link.i))] = find(parent, static_cast<std::size_t>(link.j));
    std::map<std::size_t, std::vector<std::size_t>> clusters;
    for (std::size_t k = 0; k < defect.size(); ++k) clusters[find(parent, k)].push_back(k);
    std::vector<std::vector<std::pair<std::size_t, Vec3>>> adjacency(defect.size());
    for (const auto& link : links) adjacency[static_cast<std::size_t>(link.i)].push_back({static_cast<std::size_t>(link.j), link.vector});
    std::vector<std::vector<std::size_t>> ordered;
    for (auto& [root, members] : clusters)
        if (static_cast<long long>(members.size()) >= minimumAtoms) ordered.push_back(members);
    std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });

    struct Line { Vec3 centroid, direction; double length, radius; bool periodicLine; std::string source; std::vector<Vec3> points; std::size_t atoms; };
    std::vector<Line> found;
    std::vector<int> lineId(n, -1);
    for (std::size_t c = 0; c < ordered.size(); ++c) {
        const auto& members = ordered[c];
        std::map<std::size_t, Vec3> unwrapped;
        unwrapped[members.front()] = defectPositions[members.front()];
        std::queue<std::size_t> queue;
        queue.push(members.front());
        Vec3 period{0, 0, 0};
        double periodLength = HUGE_VAL;
        while (!queue.empty()) {
            const std::size_t k = queue.front();
            queue.pop();
            for (const auto& [j, vector] : adjacency[k]) {
                const Vec3 reached = add(unwrapped[k], vector);
                if (!unwrapped.count(j)) { unwrapped[j] = reached; queue.push(j); continue; }
                const Vec3 mismatch = sub(reached, unwrapped[j]);
                if (norm(mismatch) > 0.5 * spacing && norm(mismatch) < periodLength) { period = mismatch; periodLength = norm(mismatch); }
            }
        }
        Line line;
        line.atoms = members.size();
        Vec3 mean{0, 0, 0};
        for (std::size_t k : members) { mean = add(mean, unwrapped[k]); lineId[defect[k]] = static_cast<int>(c); }
        mean = scale(mean, 1.0 / static_cast<double>(members.size()));
        Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
        for (std::size_t k : members) {
            const Vec3 d = sub(unwrapped[k], mean);
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) covariance(i, j) += d[i] * d[j];
        }
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> shape(covariance);
        Eigen::Vector3d t = shape.eigenvectors().col(2);
        line.source = "shape";
        line.periodicLine = std::isfinite(periodLength);
        if (line.periodicLine) { t = Eigen::Vector3d(period[0], period[1], period[2]).normalized(); line.source = "periodic"; }
        int largest = 0;
        for (int k = 1; k < 3; ++k) if (std::abs(t(k)) > std::abs(t(largest))) largest = k;
        if (t(largest) < 0) t = -t;
        line.direction = {t(0), t(1), t(2)};
        line.centroid = mean;
        double low = HUGE_VAL, high = -HUGE_VAL, radius = 0;
        for (std::size_t k : members) {
            const Vec3 d = sub(unwrapped[k], mean);
            const double s = dot(d, line.direction);
            low = std::min(low, s);
            high = std::max(high, s);
            radius = std::max(radius, norm(sub(d, scale(line.direction, s))));
        }
        line.length = line.periodicLine ? periodLength : high - low + spacing;
        line.radius = radius;
        const double bin = std::max(1.5 * spacing, (high - low) / 20.0);
        std::map<long long, std::pair<Vec3, int>> bins;
        for (std::size_t k : members) {
            auto& entry = bins[static_cast<long long>(std::floor((dot(sub(unwrapped[k], mean), line.direction) - low) / bin))];
            entry.first = add(entry.first, unwrapped[k]);
            ++entry.second;
        }
        for (const auto& [index, entry] : bins) line.points.push_back(scale(entry.first, 1.0 / entry.second));
        found.push_back(line);
    }

    // Burgers circuits: right-handed loops about each line, through atoms
    // near a circle in the plane normal to the line, in the reference frame.
    // Displacement steps are reduced modulo lattice vectors, so the slip
    // jump across the cut drops out and the closure failure is b.
    std::vector<Vec3> displacement(n);
    for (std::size_t i = 0; i < n; ++i) displacement[i] = sub(positions[i], refPositions[i]);
    const double requestedRadius = p.number("circuit_radius_A", 0.0);
    Json lines = Json::array();
    std::size_t otherClusters = 0;
    for (std::size_t c = 0; c < found.size(); ++c) {
        taskCheckpoint();
        const Line& line = found[c];
        double radius = requestedRadius > 0 ? requestedRadius : line.radius + 2.0 * spacing;
        // Keep the circuit clear of neighbouring lines.
        for (std::size_t o = 0; o < found.size(); ++o) {
            if (o == c) continue;
            Vec3 d = mic(sub(found[o].centroid, line.centroid));
            d = sub(d, scale(line.direction, dot(d, line.direction)));
            radius = std::min(radius, std::max(spacing, 0.5 * norm(d)));
        }
        const Vec3 helper = std::abs(line.direction[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        Vec3 e1 = cross(line.direction, helper);
        e1 = scale(e1, 1 / norm(e1));
        const Vec3 e2 = cross(line.direction, e1);
        std::vector<std::size_t> circuit;
        const int samples = std::max(16, static_cast<int>(2 * kPi * radius / (0.7 * spacing)));
        for (int k = 0; k < samples; ++k) {
            const double angle = 2 * kPi * k / samples;
            const Vec3 target = add(line.centroid, add(scale(e1, radius * std::cos(angle)), scale(e2, radius * std::sin(angle))));
            std::size_t best = n;
            double bestDistance = HUGE_VAL;
            for (std::size_t i = 0; i < n; ++i) {
                const Vec3 d = mic(sub(refPositions[i], target));
                if (std::abs(dot(d, line.direction)) > 0.75 * spacing) continue;
                const double distance = norm(d);
                if (distance < bestDistance) { bestDistance = distance; best = i; }
            }
            if (best < n && bestDistance < spacing && (circuit.empty() || circuit.back() != best)) circuit.push_back(best);
        }
        if (circuit.size() > 1 && circuit.front() == circuit.back()) circuit.pop_back();
        Vec3 burgers{0, 0, 0};
        const bool closed = circuit.size() >= 6;
        if (closed)
            for (std::size_t k = 0; k < circuit.size(); ++k) {
                const std::size_t a = circuit[k], b = circuit[(k + 1) % circuit.size()];
                burgers = add(burgers, reduceToLattice(sub(displacement[b], displacement[a]), lattice));
            }
        const double magnitude = norm(burgers);
        // A closed circuit without closure failure encloses no dislocation
        // (e.g. a surface or point-defect cluster); report it separately.
        if (closed && magnitude < 0.1 * spacing) {
            ++otherClusters;
            for (std::size_t i = 0; i < n; ++i) if (lineId[i] == static_cast<int>(c)) lineId[i] = -2;
            continue;
        }
        Json entry = Json::object();
        entry["atoms"] = line.atoms;
        entry["direction"] = toJson(line.direction);
        entry["direction_source"] = line.source;
        entry["burgers_vector_A"] = toJson(burgers);
        entry["burgers_magnitude_A"] = magnitude;
        entry["character_angle_deg"] = magnitude > 1e-6 ? std::acos(std::min(1.0, std::abs(dot(burgers, line.direction)) / magnitude)) * 180 / kPi : Json();
        entry["centroid_A"] = toJson(line.centroid);
        entry["length_A"] = line.length;
        entry["periodic_line"] = line.periodicLine;
        entry["core_radius_A"] = line.radius;
        entry["circuit_radius_A"] = radius;
        entry["circuit_atoms"] = circuit.size();
        entry["circuit_closed"] = closed;
        Json polyline = Json::array();
        for (const auto& point : line.points) polyline.push(toJson(point));
        entry["line_points_A"] = polyline;
        lines.push(entry);
    }
    Json result = Json::object();
    result["lines"] = lines;
    result["line_count"] = lines.size();
    result["non_dislocation_clusters"] = otherClusters;
    result["defect_atoms"] = defect.size();
    result["nearest_neighbour_A"] = spacing;
    result["cluster_cutoff_A"] = clusterCutoff;
    std::map<int, int> renumber;
    for (std::size_t c = 0, kept = 0; c < found.size(); ++c) {
        bool used = false;
        for (int id : lineId) used = used || id == static_cast<int>(c);
        if (used) renumber[static_cast<int>(c)] = static_cast<int>(kept++);
    }
    for (int& id : lineId) id = id >= 0 ? renumber[id] : -1;
    result["structure_type"] = toJson(deformedTypes);
    result["line_id"] = toJson(lineId);
    ToolOutput output;
    output.result = result;
    return output;
}
}
