#include "science/TrajectoryTools.h"
#include "science/AtomProperties.h"
#include "science/ScienceTools.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Selection
{
    std::size_t first = 0;
    std::size_t every = 1;
};

Selection frameSelection(const Parameters& p, const TrajectoryStream& stream)
{
    Selection s;
    s.first = static_cast<std::size_t>(integer(p.number("first_frame", 0), "first_frame", 0));
    s.every = static_cast<std::size_t>(integer(p.number("every", 1), "every"));
    if (s.first >= stream.size()) throw std::runtime_error("first_frame is beyond the trajectory (" + std::to_string(stream.size()) + " frames)");
    return s;
}

TrajectoryStream openTrajectory(const Parameters& p)
{
    if (!p.has("trajectory_file")) throw std::runtime_error("Choose a trajectory_file");
    return TrajectoryStream(std::filesystem::u8path(p.json("trajectory_file").string()));
}

Mat3 cellOf(const Structure& s)
{
    Mat3 m{};
    for (std::size_t r = 0; r < 3; ++r) m[r] = s.cellVectors[r];
    return m;
}

std::vector<Vec3> positionsOf(const Structure& s)
{
    std::vector<Vec3> out;
    out.reserve(s.atoms.size());
    for (const auto& a : s.atoms) out.push_back({a.x, a.y, a.z});
    return out;
}

// Continues `previous` (unwrapped) to the new wrapped positions by the
// minimum-image step of each atom in the frame's own cell.
void unwrapStep(std::vector<Vec3>& unwrapped, const std::vector<Vec3>& previousWrapped, const std::vector<Vec3>& current, const FrameData& frame)
{
    const bool periodic = frame.structure.hasUnitCell && (frame.pbc[0] || frame.pbc[1] || frame.pbc[2]);
    const MinimumImage mic = periodic ? MinimumImage(cellOf(frame.structure), frame.pbc) : MinimumImage();
    for (std::size_t i = 0; i < current.size(); ++i) {
        Vec3 step = sub(current[i], previousWrapped[i]);
        if (periodic) step = mic(step);
        unwrapped[i] = add(unwrapped[i], step);
    }
}

void requireSameAtoms(std::size_t expected, const FrameData& frame, std::size_t index)
{
    if (frame.structure.atoms.size() != expected)
        throw std::runtime_error("Frame " + std::to_string(index) + " has " + std::to_string(frame.structure.atoms.size()) +
                                 " atoms; the first has " + std::to_string(expected) + " (the atoms must stay the same)");
}
}

ToolOutput trajectoryAverage(const Parameters& p)
{
    const TrajectoryStream stream = openTrajectory(p);
    const Selection sel = frameSelection(p, stream);
    const bool unwrap = p.boolean("unwrap", true);
    FrameData first = stream.frame(sel.first);
    const std::size_t n = first.structure.atoms.size();
    std::vector<Vec3> wrapped = positionsOf(first.structure), unwrapped = wrapped;
    std::vector<Vec3> sum(n, Vec3{}), sumSquares(n, Vec3{});
    Mat3 cellSum{};
    std::size_t used = 0;
    for (std::size_t index = sel.first; index < stream.size(); index += sel.every) {
        taskProgress(static_cast<double>(index) / static_cast<double>(stream.size()));
        const FrameData frame = index == sel.first ? first : stream.frame(index);
        requireSameAtoms(n, frame, index);
        const auto current = positionsOf(frame.structure);
        if (index != sel.first) {
            if (unwrap) unwrapStep(unwrapped, wrapped, current, frame);
            else unwrapped = current;
        }
        wrapped = current;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t d = 0; d < 3; ++d) {
                sum[i][d] += unwrapped[i][d];
                sumSquares[i][d] += unwrapped[i][d] * unwrapped[i][d];
            }
        if (frame.structure.hasUnitCell)
            for (std::size_t r = 0; r < 3; ++r) cellSum[r] = add(cellSum[r], frame.structure.cellVectors[r]);
        ++used;
    }
    Structure average = first.structure;
    std::vector<double> fluctuation(n);
    double meanFluctuation = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 mean = scale(sum[i], 1.0 / static_cast<double>(used));
        double variance = 0;
        for (std::size_t d = 0; d < 3; ++d) variance += std::max(0.0, sumSquares[i][d] / static_cast<double>(used) - mean[d] * mean[d]);
        fluctuation[i] = std::sqrt(variance);
        meanFluctuation += fluctuation[i] / static_cast<double>(n);
        average.atoms[i].x = mean[0];
        average.atoms[i].y = mean[1];
        average.atoms[i].z = mean[2];
    }
    if (average.hasUnitCell)
        for (std::size_t r = 0; r < 3; ++r) average.cellVectors[r] = scale(cellSum[r], 1.0 / static_cast<double>(used));
    ToolOutput output;
    output.result = Json::object();
    output.result["frames_averaged"] = static_cast<double>(used);
    output.result["rms_fluctuation_A"] = toJson(fluctuation);
    output.result["mean_rms_fluctuation_A"] = meanFluctuation;
    output.result["unwrapped"] = unwrap;
    output.frames.push_back(std::move(average));
    return output;
}

ToolOutput unwrapTrajectory(const Parameters& p)
{
    const TrajectoryStream stream = openTrajectory(p);
    const Selection sel = frameSelection(p, stream);
    const long long window = integer(p.number("smooth_window", 1), "smooth_window");
    FrameData first = stream.frame(sel.first);
    const std::size_t n = first.structure.atoms.size();
    std::vector<Vec3> wrapped = positionsOf(first.structure), unwrapped = wrapped;
    std::vector<Structure> frames;
    std::vector<double> times;
    for (std::size_t index = sel.first; index < stream.size(); index += sel.every) {
        taskProgress(static_cast<double>(index) / static_cast<double>(stream.size()));
        FrameData frame = index == sel.first ? first : stream.frame(index);
        requireSameAtoms(n, frame, index);
        const auto current = positionsOf(frame.structure);
        if (index != sel.first) unwrapStep(unwrapped, wrapped, current, frame);
        wrapped = current;
        for (std::size_t i = 0; i < n; ++i) {
            frame.structure.atoms[i].x = unwrapped[i][0];
            frame.structure.atoms[i].y = unwrapped[i][1];
            frame.structure.atoms[i].z = unwrapped[i][2];
        }
        frames.push_back(std::move(frame.structure));
        times.push_back(frame.timeFs ? *frame.timeFs : static_cast<double>(index));
    }
    ToolOutput output;
    if (window > 1) {
        // Centred moving average of the unwrapped positions (shorter at the ends).
        const std::size_t half = static_cast<std::size_t>(window / 2);
        std::vector<Structure> smoothed = frames;
        for (std::size_t f = 0; f < frames.size(); ++f) {
            const std::size_t low = f >= half ? f - half : 0, high = std::min(frames.size() - 1, f + half);
            for (std::size_t i = 0; i < n; ++i) {
                Vec3 mean{};
                for (std::size_t g = low; g <= high; ++g)
                    mean = add(mean, Vec3{frames[g].atoms[i].x, frames[g].atoms[i].y, frames[g].atoms[i].z});
                mean = scale(mean, 1.0 / static_cast<double>(high - low + 1));
                smoothed[f].atoms[i].x = mean[0];
                smoothed[f].atoms[i].y = mean[1];
                smoothed[f].atoms[i].z = mean[2];
            }
        }
        frames = std::move(smoothed);
    }
    // Net displacement of each atom over the trajectory (unwrapped).
    std::vector<double> travelled(n, 0.0);
    if (!frames.empty())
        for (std::size_t i = 0; i < n; ++i)
            travelled[i] = norm(sub(Vec3{frames.back().atoms[i].x, frames.back().atoms[i].y, frames.back().atoms[i].z},
                                    Vec3{frames.front().atoms[i].x, frames.front().atoms[i].y, frames.front().atoms[i].z}));
    output.result = Json::object();
    output.result["frames"] = static_cast<double>(frames.size());
    output.result["smooth_window"] = static_cast<double>(window);
    output.result["net_displacement_A"] = toJson(travelled);
    output.frames = std::move(frames);
    output.times = std::move(times);
    return output;
}

Json timeSeries(const Parameters& p)
{
    const TrajectoryStream stream = openTrajectory(p);
    const Selection sel = frameSelection(p, stream);
    const std::string tool = p.has("analysis") ? p.json("analysis").string() : std::string("cna");
    const std::string wanted = p.has("property") ? p.json("property").string() : std::string();
    const double cutoff = p.number("cutoff_A", 3.0);
    std::vector<double> frameIndex, time, volume, density, meanProperty;
    std::map<std::string, std::vector<double>> counts;
    std::vector<std::string> countOrder;
    std::string propertyName;
    std::size_t step = 0;
    for (std::size_t index = sel.first; index < stream.size(); index += sel.every, ++step) {
        taskProgress(static_cast<double>(index) / static_cast<double>(stream.size()));
        const FrameData frame = stream.frame(index);
        const Structure& s = frame.structure;
        frameIndex.push_back(static_cast<double>(index));
        time.push_back(frame.timeFs ? *frame.timeFs : kNaN);
        const double v = s.hasUnitCell ? cellVolume(cellOf(s)) : kNaN;
        volume.push_back(v);
        density.push_back(s.hasUnitCell ? static_cast<double>(s.atoms.size()) / v : kNaN);
        if (tool == "none") continue;
        // The same request as the per-atom analysis of the view.
        Json request = Json::object();
        Json positions = Json::array();
        for (const auto& a : s.atoms) positions.push(Json::array({a.x, a.y, a.z}));
        request["positions"] = positions;
        if (s.hasUnitCell) {
            Json cell = Json::array();
            for (const auto& row : s.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
            request["cell"] = cell;
            request["pbc"] = Json::array({frame.pbc[0], frame.pbc[1], frame.pbc[2]});
        }
        const bool ownRange = tool == "ptm" || tool == "cna" || tool == "voronoi" || tool == "structure-type" || tool == "diamond-structure" ||
                              tool == "chill-plus" || tool == "local-entropy" || tool == "fcc-planar-faults";
        if (!ownRange) request["cutoff_A"] = cutoff;
        if (tool == "local-entropy" && !(frame.pbc[0] && frame.pbc[1] && frame.pbc[2])) request["local_density"] = true;
        const Json result = runTool(tool, request).result;
        if (const Json* c = result.find("counts"); c && c->isObject()) {
            for (const auto& [name, value] : c->members()) {
                if (!counts.count(name)) { counts[name] = std::vector<double>(step, kNaN); countOrder.push_back(name); }
                counts[name].push_back(value.number());
            }
        }
        const auto properties = perAtomProperties(tool, result);
        double mean = kNaN;
        for (const auto& property : properties) {
            if (!wanted.empty() && property.name != wanted) continue;
            propertyName = property.name;
            double sum = 0;
            std::size_t valid = 0;
            for (double x : property.values)
                if (std::isfinite(x)) { sum += x; ++valid; }
            mean = valid ? sum / static_cast<double>(valid) : kNaN;
            break;
        }
        meanProperty.push_back(mean);
    }
    Json result = Json::object();
    result["frame"] = toJson(frameIndex);
    result["time_fs"] = toJson(time);
    result["volume_A3"] = toJson(volume);
    result["number_density_per_A3"] = toJson(density);
    result["analysis"] = tool;
    if (!meanProperty.empty()) {
        result["mean_property"] = toJson(meanProperty);
        result["property"] = propertyName;
    }
    Json countJson = Json::object();
    for (const auto& name : countOrder) countJson[name] = toJson(counts[name]);
    result["counts"] = countJson;
    return result;
}

namespace
{
// The first (or named) per-atom property of an analysis of one frame.
std::vector<double> analyseFrameProperty(const FrameData& frame, const std::string& tool, const std::string& wanted, double cutoff, std::string& name)
{
    const Structure& s = frame.structure;
    Json request = Json::object();
    Json positions = Json::array();
    for (const auto& a : s.atoms) positions.push(Json::array({a.x, a.y, a.z}));
    request["positions"] = positions;
    if (s.hasUnitCell) {
        Json cell = Json::array();
        for (const auto& row : s.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
        request["cell"] = cell;
        request["pbc"] = Json::array({frame.pbc[0], frame.pbc[1], frame.pbc[2]});
    }
    const bool ownRange = tool == "ptm" || tool == "cna" || tool == "voronoi" || tool == "vorotop" || tool == "structure-type" ||
                          tool == "diamond-structure" || tool == "chill-plus" || tool == "local-entropy" || tool == "fcc-planar-faults";
    if (!ownRange) request["cutoff_A"] = cutoff;
    if (tool == "local-entropy" && !(frame.pbc[0] && frame.pbc[1] && frame.pbc[2])) request["local_density"] = true;
    for (const auto& property : perAtomProperties(tool, runTool(tool, request).result))
        if (wanted.empty() || property.name == wanted) { name = property.name; return property.values; }
    throw std::runtime_error(tool + (wanted.empty() ? " gives no per-atom property" : " has no property " + wanted));
}
}

Json frameDifference(const Parameters& p)
{
    const TrajectoryStream stream = openTrajectory(p);
    const long long count = static_cast<long long>(stream.size());
    long long reference = static_cast<long long>(p.number("reference_frame", 0));
    long long current = static_cast<long long>(p.number("frame", -1));
    if (reference < 0) reference += count;
    if (current < 0) current += count;
    if (reference < 0 || reference >= count || current < 0 || current >= count)
        throw std::runtime_error("Frames must lie in the trajectory (" + std::to_string(count) + " frames; negative counts from the end)");
    const FrameData a = stream.frame(static_cast<std::size_t>(reference)), b = stream.frame(static_cast<std::size_t>(current));
    if (a.structure.atoms.size() != b.structure.atoms.size()) throw std::runtime_error("The two frames have different numbers of atoms");
    const std::string tool = p.has("analysis") ? p.json("analysis").string() : std::string("positions");
    const std::size_t n = a.structure.atoms.size();
    std::vector<double> difference(n), before, after;
    std::string name = "displacement (A)";
    if (tool == "positions") {
        const bool periodic = b.structure.hasUnitCell && (b.pbc[0] || b.pbc[1] || b.pbc[2]);
        const MinimumImage mic = periodic ? MinimumImage(cellOf(b.structure), b.pbc) : MinimumImage();
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 d = sub(Vec3{b.structure.atoms[i].x, b.structure.atoms[i].y, b.structure.atoms[i].z},
                         Vec3{a.structure.atoms[i].x, a.structure.atoms[i].y, a.structure.atoms[i].z});
            if (periodic) d = mic(d);
            difference[i] = norm(d);
        }
    } else {
        const std::string wanted = p.has("property") ? p.json("property").string() : std::string();
        const double cutoff = p.number("cutoff_A", 3.0);
        before = analyseFrameProperty(a, tool, wanted, cutoff, name);
        after = analyseFrameProperty(b, tool, wanted, cutoff, name);
        for (std::size_t i = 0; i < n; ++i) difference[i] = after[i] - before[i];
    }
    double sum = 0, squares = 0, largest = 0;
    std::size_t changed = 0, valid = 0;
    for (double d : difference) {
        if (!std::isfinite(d)) continue;
        ++valid;
        sum += d;
        squares += d * d;
        largest = std::max(largest, std::abs(d));
        changed += std::abs(d) > 1e-12;
    }
    Json result = Json::object();
    result["difference"] = toJson(difference);
    if (!before.empty()) { result["reference_values"] = toJson(before); result["frame_values"] = toJson(after); }
    result["quantity"] = name;
    result["reference_frame"] = static_cast<double>(reference);
    result["frame"] = static_cast<double>(current);
    result["mean_difference"] = valid ? sum / static_cast<double>(valid) : std::numeric_limits<double>::quiet_NaN();
    result["rms_difference"] = valid ? std::sqrt(squares / static_cast<double>(valid)) : std::numeric_limits<double>::quiet_NaN();
    result["max_abs_difference"] = largest;
    result["atoms_changed"] = static_cast<double>(changed);
    return result;
}
}
