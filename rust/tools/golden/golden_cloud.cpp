// Reference atom clouds of the C++ cloud library (rust/tests/golden/cloud/ and
// cloud.json) for the Rust af-cloud tests: generated crystals, XYZ and LAMMPS
// dumps built with one thread (byte-reproducible files), read back, and the
// error cases. Usage (from the repository root): golden_cloud <fixture dir> <output dir>.

#include "cloud/AtomCloud.h"
#include "science/Json.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>

namespace cloud = atomforge::cloud;
using atomforge::science::Json;

namespace
{
Json vec(const std::array<double, 3>& v) { return Json::array({v[0], v[1], v[2]}); }
Json infoJson(const cloud::CloudInfo& info)
{
    Json o = Json::object();
    o["atoms"] = static_cast<double>(info.atoms);
    o["lower"] = vec(info.lower);
    o["upper"] = vec(info.upper);
    o["has_cell"] = info.hasCell;
    o["cell"] = Json::array({vec(info.cell[0]), vec(info.cell[1]), vec(info.cell[2])});
    o["cell_origin"] = vec(info.cellOrigin);
    Json species = Json::array();
    for (const auto& s : info.species) species.push(Json::array({Json(s.name), Json(s.color[0]), Json(s.color[1]), Json(s.color[2]), Json(s.radius)}));
    o["species"] = species;
    Json chunks = Json::array();
    for (const auto& c : info.chunks)
        chunks.push(Json::array({static_cast<double>(c.offset), static_cast<double>(c.count), c.lower[0], c.lower[1], c.lower[2], c.upper[0], c.upper[1], c.upper[2]}));
    o["chunks"] = chunks;
    o["source"] = info.source;
    return o;
}
Json cases = Json::array();
void attempt(const std::string& name, const std::function<Json()>& work)
{
    Json entry = Json::object();
    entry["name"] = name;
    try { entry["result"] = work(); }
    catch (const std::exception& e) { entry["error"] = std::string(e.what()); }
    cases.push(entry);
}
cloud::BuildOptions options(std::uint64_t chunk = 1024, std::uint64_t seed = 12345)
{
    cloud::BuildOptions o;
    o.chunkAtoms = chunk;
    o.threads = 1;
    o.seed = seed;
    return o;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_cloud <fixture dir> <output dir>\n"); return 1; }
    const std::string in = std::string(argv[1]) + "/", out = std::string(argv[2]) + "/cloud/";
    std::filesystem::create_directories(out);
    const auto built = [&](const std::string& name, const std::function<cloud::CloudInfo(const std::filesystem::path&)>& build) {
        attempt("build " + name, [&] { return infoJson(build(out + name + ".afcloud")); });
    };
    built("fcc", [&](const auto& p) { return cloud::generateCrystal(p, "fcc", 3.615, 0, {10, 10, 10}, "Cu", "", {}, options()); });
    built("hcp", [&](const auto& p) { return cloud::generateCrystal(p, "hcp", 3.21, 0, {6, 4, 5}, "Mg", "Zn", {}, options(2000, 7)); });
    built("diamond", [&](const auto& p) { return cloud::generateCrystal(p, "diamond", 5.43, 0, {5, 5, 5}, "si", "", {}, options()); });
    built("bcc", [&](const auto& p) { return cloud::generateCrystal(p, "bcc", 2.87, 0, {3, 3, 3}, "Fe", "type 9", {}, options()); });
    built("sc", [&](const auto& p) { return cloud::generateCrystal(p, "sc", 2.5, 0, {1, 1, 1}, "Po", "", {}, options()); });
    built("hcp_c", [&](const auto& p) { return cloud::generateCrystal(p, "hcp", 2.95, 4.68, {3, 2, 2}, "Ti", "", {}, options()); });
    built("bad_a", [&](const auto& p) { return cloud::generateCrystal(p, "fcc", 0, 0, {1, 1, 1}, "Cu", "", {}, options()); });
    built("bad_lattice", [&](const auto& p) { return cloud::generateCrystal(p, "a15", 3, 0, {1, 1, 1}, "Cu", "", {}, options()); });
    built("bad_cells", [&](const auto& p) { return cloud::generateCrystal(p, "fcc", 3, 0, {1, 0, 1}, "Cu", "", {}, options()); });
    built("xyz", [&](const auto& p) { return cloud::buildFromXyz(in + "big.extxyz", p, 0, {}, options()); });
    built("xyz_frame1", [&](const auto& p) { return cloud::buildFromXyz(in + "two_frames.xyz", p, 1, {}, options(1500, 99)); });
    built("xyz_bad", [&](const auto& p) { return cloud::buildFromXyz(in + "bad_coordinate.xyz", p, 0, {}, options()); });
    built("xyz_short", [&](const auto& p) { return cloud::buildFromXyz(in + "short.xyz", p, 0, {}, options()); });
    built("dump_elements", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "elements.dump", p, {}, 0, {}, options()); });
    built("dump_triclinic", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "triclinic.dump", p, {"Cu", "Ni"}, 0, {}, options()); });
    built("dump_scaled_frame1", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "scaled_two_frames.dump", p, {}, 1, {}, options()); });
    built("dump_no_positions", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "no_positions.dump", p, {}, 0, {}, options()); });
    built("dump_short_line", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "short_line.dump", p, {}, 0, {}, options()); });
    built("dump_not", [&](const auto& p) { return cloud::buildFromLammpsDump(in + "not_a_dump.dump", p, {}, 0, {}, options()); });

    // Reading back.
    for (const char* name : {"fcc", "hcp", "xyz", "dump_triclinic"})
        attempt(std::string("read ") + name, [&] {
            const cloud::CloudFile file(out + name + ".afcloud");
            Json o = Json::object();
            o["info"] = infoJson(file.info());
            Json atoms = Json::array();
            for (std::size_t c = 0; c < file.info().chunks.size(); ++c) {
                std::vector<cloud::PackedAtom> packed(std::min<std::uint64_t>(6, file.info().chunks[c].count));
                file.read(c, file.info().chunks[c].count - packed.size(), packed.size(), packed.data());
                for (const auto& a : packed) {
                    const auto p = file.position(c, a);
                    atoms.push(Json::array({a.x, a.y, a.z, a.species, p[0], p[1], p[2]}));
                }
            }
            o["atoms"] = atoms;
            return o;
        });
    attempt("read past end", [&] {
        const cloud::CloudFile file(out + "sc.afcloud");
        std::vector<cloud::PackedAtom> packed(5);
        file.read(0, 0, 5, packed.data());
        return Json(true);
    });
    {
        std::ofstream incomplete(out + "incomplete.afcloud", std::ios::binary);
        incomplete.write("AFCLOUD1", 8);
        const char zeros[16] = {};
        incomplete.write(zeros, 16);
    }
    for (const std::string& path : {out + "incomplete.afcloud", in + "text.txt", in + "missing.afcloud"})
        attempt("open " + std::filesystem::path(path).filename().string(), [&] { return Json(static_cast<double>(cloud::CloudFile(path).info().atoms)); });
    for (const std::string& path : {out + "fcc.afcloud", in + "text.txt", in + "missing.afcloud", out + "incomplete.afcloud"})
        attempt("is cloud " + std::filesystem::path(path).filename().string(), [&] { return Json(cloud::isCloudFile(path)); });
    for (const char* name : {"big.extxyz", "two_frames.xyz", "elements.dump", "triclinic.dump", "text.txt", "missing.xyz"})
        attempt(std::string("peek ") + name, [&] { return Json(static_cast<double>(cloud::peekAtomCount(in + name))); });

    std::ofstream(std::string(argv[2]) + "/cloud.json") << cases.dump();
    std::printf("wrote %scloud.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}
