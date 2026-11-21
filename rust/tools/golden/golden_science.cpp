// Reference outputs of the C++ science core and data readers (science.json)
// for the Rust af-science tests. Usage: golden_science <fixtures> <output dir>.

#include "science/ScienceCore.h"
#include "science/ScienceData.h"

#include <algorithm>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace atomforge::science;
namespace fs = std::filesystem;

namespace
{
template <typename F>
Json attempt(F&& f)
{
    try {
        return f();
    } catch (const std::exception& e) {
        Json error = Json::object();
        error["error"] = std::string(e.what());
        return error;
    }
}

Json frameJson(const FrameData& frame)
{
    Json result = Json::object();
    result["structure"] = structureJson(frame.structure);
    result["has_cell"] = frame.structure.hasUnitCell;
    result["pbc"] = Json::array({frame.pbc[0], frame.pbc[1], frame.pbc[2]});
    Json velocities = Json::array();
    for (const auto& v : frame.velocities) velocities.push(toJson(v));
    result["velocities"] = velocities;
    result["masses"] = toJson(frame.masses);
    result["time_fs"] = frame.timeFs ? Json(*frame.timeFs) : Json();
    Json property = Json::array();
    for (double p : frame.structure.atomProperty) property.push(p);
    result["atom_property"] = property;
    Json colours = Json::array();
    for (const auto& a : frame.structure.atoms) colours.push(Json::array({double(a.r), double(a.g), double(a.b)}));
    result["colours"] = colours;
    return result;
}

Json framesJson(const std::vector<FrameData>& frames)
{
    Json list = Json::array();
    for (const auto& f : frames) list.push(frameJson(f));
    return list;
}

Json mat(const Mat3& m) { return toJson(m); }

Mat3 randomCell(std::mt19937& rng)
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Mat3 cell{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) cell[r][c] = (r == c ? 3.0 + 2.0 * (u(rng) + 1.0) : 1.5 * u(rng));
    return cell;
}

Json describeParameters(const Parameters& p, const std::vector<std::string>& names)
{
    Json r = Json::object();
    for (const auto& n : names) {
        Json item = Json::object();
        item["has"] = p.has(n);
        item["array"] = attempt([&] { return toJson(p.array(n)); });
        item["json"] = attempt([&] { return p.json(n); });
        item["structure"] = attempt([&] {
            const StructureInput& s = p.structure(n);
            Json o = Json::object();
            o["structure"] = structureJson(s.structure);
            o["pbc"] = Json::array({s.pbc[0], s.pbc[1], s.pbc[2]});
            return o;
        });
        item["frames"] = attempt([&] { return framesJson(p.frames(n)); });
        r[n] = item;
    }
    const Pbc pbc = p.pbc("pbc", {false, true, false});
    r["pbc_fallback"] = Json::array({pbc[0], pbc[1], pbc[2]});
    return r;
}

Json parameterCases(const fs::path& base, const std::string& name)
{
    Json cases = Json::array();
    const std::vector<std::pair<std::string, std::string>> requests = {
        {"x", "{\"data\": {\"file\": \"@\"}}"},
        {"x", "{\"data\": {\"file\": \"@\", \"column\": 1}}"},
        {"x", "{\"data\": {\"file\": \"@\", \"column\": 0}}"},
        {"x", "{\"data\": {\"file\": \"@\", \"field\": \"positions\"}}"},
        {"x", "{\"data\": {\"file\": \"@\", \"field\": \"cell.1\"}}"},
        {"x", "{\"energies_eV\": {\"file\": \"@\"}}"},
        {"x", "{\"velocities\": {\"file\": \"@\"}, \"masses\": {\"file\": \"@\"}, \"cell\": {\"file\": \"@\"}}"},
        {"msd", "{\"positions\": {\"file\": \"@\"}, \"timestep_fs\": 2.5}"},
        {"x", "{\"structure\": {\"file\": \"@\"}, \"initial\": {\"file\": \"./@\"}}"},
        {"x", "{\"restart_images\": {\"file\": \"@\"}, \"out_file\": \"sub/../@\", \"other\": {\"file\": \"@\", \"bogus\": 1}}"},
    };
    for (const auto& [tool, pattern] : requests) {
        std::string text = pattern;
        for (std::size_t at; (at = text.find('@')) != std::string::npos;) text.replace(at, 1, name);
        Json item = Json::object();
        item["tool"] = tool;
        item["request"] = text;
        item["result"] = attempt([&] {
            const Parameters p(tool, Json::parse(text), base);
            return describeParameters(p, {"data", "energies_eV", "velocities", "masses", "cell", "positions", "structure", "initial", "restart_images", "out_file"});
        });
        cases.push(item);
    }
    return cases;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: golden_science <fixtures> <output>\n");
        return 1;
    }
    const fs::path fixtures = argv[1];
    Json root = Json::object();

    // Matrix algebra and minimum image on random cells.
    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Json algebra = Json::array();
    for (int c = 0; c < 12; ++c) {
        const Mat3 cell = randomCell(rng);
        const Pbc pbc{c % 4 != 1, c % 4 != 2, c % 3 != 0};
        Json item = Json::object();
        item["cell"] = mat(cell);
        item["pbc"] = Json::array({pbc[0], pbc[1], pbc[2]});
        item["det"] = determinant(cell);
        item["inverse"] = mat(inverse(cell));
        item["product"] = mat(multiply(cell, transpose(cell)));
        Json vectors = Json::array(), images = Json::array(), fracs = Json::array();
        const MinimumImage image(cell, pbc);
        for (int v = 0; v < 20; ++v) {
            const Vec3 vec{8 * u(rng), 8 * u(rng), 8 * u(rng)};
            vectors.push(toJson(vec));
            images.push(toJson(image(vec)));
            fracs.push(toJson(fractional(vec, cell)));
        }
        item["vectors"] = vectors;
        item["images"] = images;
        item["fractional"] = fracs;
        // Neighbour list of random points (some outside the cell).
        Json positions = Json::array();
        std::vector<Vec3> points;
        for (int a = 0; a < 25; ++a) {
            const Vec3 p = rowTimes({1.4 * u(rng) + 0.3, 1.4 * u(rng) + 0.3, 1.4 * u(rng) + 0.3}, cell);
            points.push_back(p);
            positions.push(toJson(p));
        }
        item["positions"] = positions;
        const double cutoff = 2.0 + c * 0.35;
        item["cutoff"] = cutoff;
        item["neighbors"] = attempt([&] {
            Json list = Json::array();
            for (const auto& n : neighborList(points, cell, pbc, cutoff))
                list.push(Json::array({n.i, n.j, n.shift[0], n.shift[1], n.shift[2], n.vector[0], n.vector[1], n.vector[2]}));
            return list;
        });
        algebra.push(item);
    }
    root["algebra"] = algebra;

    // FFTs of several lengths.
    Json ffts = Json::array();
    for (int n : {1, 2, 3, 5, 8, 12, 17, 64, 100}) {
        std::vector<std::complex<double>> data(static_cast<std::size_t>(n));
        Json input = Json::array();
        for (auto& d : data) {
            d = {u(rng), u(rng)};
            input.push(Json::array({d.real(), d.imag()}));
        }
        Json item = Json::object();
        item["input"] = input;
        for (bool inv : {false, true}) {
            Json out = Json::array();
            for (const auto& d : fft(data, inv)) out.push(Json::array({d.real(), d.imag()}));
            item[inv ? "inverse" : "forward"] = out;
        }
        ffts.push(item);
    }
    root["fft"] = ffts;

    // Array conversion and validation helpers.
    Json arrays = Json::array();
    for (const char* text : {"3.5", "[1,2,3]", "[[1,2],[3,4]]", "[[1,2],[3]]", "[]", "[[],[]]", "[true,false]", "[\"a\"]",
                             "[[[1],[2]],[[3],[4]]]", "[[1,[2]],[3,4]]"}) {
        Json item = Json::object();
        item["text"] = std::string(text);
        item["array"] = attempt([&] {
            const NdArray a = toArray(Json::parse(text), "x");
            Json r = Json::object();
            Json shape = Json::array();
            for (auto s : a.shape) shape.push(s);
            r["shape"] = shape;
            r["values"] = toJson(a.values);
            r["json"] = toJson(a);
            return r;
        });
        item["cell"] = attempt([&] { return mat(cellMatrix(toArray(Json::parse(text), "cell"))); });
        item["pbc"] = attempt([&] { const Pbc p = parsePbc(Json::parse(text)); return Json::array({p[0], p[1], p[2]}); });
        arrays.push(item);
    }
    for (const char* text : {"[[1,0,0],[0,1,0],[0,0,1]]", "[[1,0,0],[2,0,0],[0,0,1]]", "[[1e6,0,0],[0,1,0],[0,0,1e-5]]", "[[1,0],[0,1]]"}) {
        Json item = Json::object();
        item["text"] = std::string(text);
        item["cell"] = attempt([&] { return mat(cellMatrix(toArray(Json::parse(text), "cell"))); });
        arrays.push(item);
    }
    root["arrays"] = arrays;
    Json checks = Json::array();
    for (double v : {1.0, 0.0, -2.0, 2.5, 1e16, 3.0}) {
        Json item = Json::object();
        item["value"] = v;
        item["positive"] = attempt([&] { return Json(positive(v, "x")); });
        item["integer"] = attempt([&] { return Json(static_cast<double>(integer(v, "n", 1))); });
        item["integer0"] = attempt([&] { return Json(static_cast<double>(integer(v, "n", 0))); });
        checks.push(item);
    }
    root["checks"] = checks;

    // Frame readers, streams and data files.
    Json files = Json::object();
    std::vector<fs::path> paths;
    for (const auto& entry : fs::directory_iterator(fixtures)) paths.push_back(entry.path());
    std::sort(paths.begin(), paths.end());
    for (const auto& path : paths) {
        const std::string name = path.filename().u8string();
        Json item = Json::object();
        item["frames"] = attempt([&] { return framesJson(readFrames(path)); });
        item["stream"] = attempt([&] {
            TrajectoryStream stream(path);
            Json r = Json::object();
            r["streamed"] = stream.streamed();
            std::vector<FrameData> frames;
            for (std::size_t i = stream.size(); i-- > 0;) frames.insert(frames.begin(), stream.frame(i));
            r["frames"] = framesJson(frames);
            return r;
        });
        item["parameters"] = parameterCases(fixtures, name);
        files[name] = item;
    }
    root["files"] = files;

    std::ofstream(fs::path(argv[2]) / "science.json") << root.dump(1);
    std::printf("wrote %s/science.json\n", argv[2]);
    return 0;
}
