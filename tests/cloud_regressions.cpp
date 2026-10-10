// Out-of-core atom clouds: building from generated crystals, LAMMPS dumps and
// XYZ files, reading back, quantization accuracy, chunking and the random
// order that makes chunk prefixes uniform subsamples.
#include "cloud/AtomCloud.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace cloud = atomforge::cloud;

namespace
{
int failures = 0;
void check(bool condition, const std::string& what)
{
    if (!condition) throw std::runtime_error(what);
}
void test(const std::string& name, const std::function<void()>& body)
{
    try { body(); std::cout << "  ok   " << name << "\n"; }
    catch (const std::exception& error) { ++failures; std::cout << "  FAIL " << name << ": " << error.what() << "\n"; }
}

std::filesystem::path folder()
{
    static const auto path = [] {
        auto p = std::filesystem::temp_directory_path() / ("atomforge-cloud-tests-" + std::to_string(std::random_device{}()));
        std::filesystem::create_directories(p);
        return p;
    }();
    return path;
}

// Every atom of a cloud, decoded.
std::vector<std::array<double, 4>> readAll(const cloud::CloudFile& file)
{
    std::vector<std::array<double, 4>> atoms;
    for (std::size_t c = 0; c < file.info().chunks.size(); ++c) {
        std::vector<cloud::PackedAtom> packed(file.info().chunks[c].count);
        file.read(c, 0, packed.size(), packed.data());
        for (const auto& a : packed) {
            const auto p = file.position(c, a);
            atoms.push_back({p[0], p[1], p[2], static_cast<double>(a.species)});
        }
    }
    return atoms;
}

// Whether two sets of atoms agree within `tolerance` (each atom has a match).
bool sameAtoms(const std::vector<std::array<double, 4>>& a, const std::vector<std::array<double, 4>>& b, double tolerance)
{
    if (a.size() != b.size()) return false;
    std::map<std::array<long long, 3>, std::vector<std::size_t>> grid;
    const double cell = 4 * tolerance;
    for (std::size_t i = 0; i < b.size(); ++i)
        grid[{std::llround(b[i][0] / cell), std::llround(b[i][1] / cell), std::llround(b[i][2] / cell)}].push_back(i);
    std::vector<char> used(b.size(), 0);
    for (const auto& p : a) {
        const long long c[3] = {std::llround(p[0] / cell), std::llround(p[1] / cell), std::llround(p[2] / cell)};
        bool found = false;
        for (long long dx = -1; dx <= 1 && !found; ++dx)
            for (long long dy = -1; dy <= 1 && !found; ++dy)
                for (long long dz = -1; dz <= 1 && !found; ++dz) {
                    const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                    if (it == grid.end()) continue;
                    for (std::size_t j : it->second)
                        if (!used[j] && std::abs(b[j][0] - p[0]) <= tolerance && std::abs(b[j][1] - p[1]) <= tolerance && std::abs(b[j][2] - p[2]) <= tolerance) {
                            used[j] = 1;
                            found = true;
                            break;
                        }
                }
        if (!found) return false;
    }
    return true;
}

// Sorted positions, rounded to a grid, for comparing sets of atoms.
std::vector<std::array<long long, 3>> keys(const std::vector<std::array<double, 4>>& atoms, double grid)
{
    std::vector<std::array<long long, 3>> out;
    for (const auto& a : atoms) out.push_back({std::llround(a[0] / grid), std::llround(a[1] / grid), std::llround(a[2] / grid)});
    std::sort(out.begin(), out.end());
    return out;
}
}

int main()
{
    std::cout << "Atom cloud regressions\n";

    test("generated crystal: every atom once, accurate positions, chunked", [] {
        const auto path = folder() / "fcc.afcloud";
        cloud::BuildOptions options;
        options.chunkAtoms = 4096;  // many chunks from a small crystal
        const auto info = cloud::generateCrystal(path, "fcc", 3.615, 0, {20, 20, 20}, "Cu", "", {}, options);
        check(info.atoms == 32000, "atom count " + std::to_string(info.atoms));
        check(cloud::isCloudFile(path), "magic");
        const cloud::CloudFile file(path);
        check(file.info().atoms == 32000 && file.info().chunks.size() >= 8, "chunks " + std::to_string(file.info().chunks.size()));
        for (const auto& c : file.info().chunks) check(c.count <= 4096, "chunk size limit");
        check(file.info().species.size() == 1 && file.info().species[0].name == "Cu", "species");
        check(file.info().species[0].radius > 1.0f && file.info().species[0].radius < 1.6f, "Cu radius from the element table");
        check(file.info().hasCell && std::abs(file.info().cell[0][0] - 72.3) < 1e-9, "cell");
        // Positions within the quantization error of the lattice sites.
        const auto atoms = readAll(file);
        double worst = 0;
        for (const auto& a : atoms)
            for (int k = 0; k < 3; ++k) {
                const double site = std::round(a[static_cast<std::size_t>(k)] / 1.8075) * 1.8075;
                worst = std::max(worst, std::abs(a[static_cast<std::size_t>(k)] - site));
            }
        check(worst < 2e-3, "quantization error " + std::to_string(worst));
        // No duplicates: every site once.
        auto k = keys(atoms, 0.9);
        check(std::adjacent_find(k.begin(), k.end()) == k.end(), "duplicate atoms");
    });

    test("chunk prefixes are uniform random subsamples", [] {
        const auto path = folder() / "prefix.afcloud";
        cloud::BuildOptions options;
        options.chunkAtoms = 1u << 20;
        cloud::generateCrystal(path, "sc", 2.0, 0, {64, 64, 64}, "Fe", "", {}, options);
        const cloud::CloudFile file(path);
        // The first 10% of the largest chunk spreads evenly over its box (octant counts).
        std::size_t biggest = 0;
        for (std::size_t c = 0; c < file.info().chunks.size(); ++c)
            if (file.info().chunks[c].count > file.info().chunks[biggest].count) biggest = c;
        const auto& chunk = file.info().chunks[biggest];
        std::vector<cloud::PackedAtom> prefix(chunk.count / 10);
        file.read(biggest, 0, prefix.size(), prefix.data());
        int octants[8] = {};
        for (const auto& a : prefix) ++octants[(a.x > 32767 ? 1 : 0) | (a.y > 32767 ? 2 : 0) | (a.z > 32767 ? 4 : 0)];
        const double expected = static_cast<double>(prefix.size()) / 8;
        for (int o : octants) check(std::abs(o - expected) < 0.1 * expected, "octant share " + std::to_string(o) + " vs " + std::to_string(expected));
        // Reads at an offset agree with the whole chunk.
        std::vector<cloud::PackedAtom> all(chunk.count), part(100);
        file.read(biggest, 0, all.size(), all.data());
        file.read(biggest, 1000, 100, part.data());
        check(std::equal(part.begin(), part.end(), all.begin() + 1000, [](const auto& p, const auto& q) { return p.x == q.x && p.y == q.y && p.z == q.z && p.species == q.species; }), "offset read");
    });

    test("LAMMPS dumps: orthogonal and scaled coordinates, types, second frame", [] {
        const auto dump = folder() / "atoms.dump";
        std::ofstream out(dump);
        out.precision(12);
        std::mt19937 random(7);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        std::vector<std::array<double, 4>> first, second;
        for (int frame = 0; frame < 2; ++frame) {
            out << "ITEM: TIMESTEP\n" << frame * 100 << "\nITEM: NUMBER OF ATOMS\n5000\nITEM: BOX BOUNDS pp pp pp\n0 40\n-10 30\n5 25\n";
            out << (frame == 0 ? "ITEM: ATOMS id type x y z\n" : "ITEM: ATOMS id type xs ys zs vx\n");
            for (int i = 0; i < 5000; ++i) {
                const double s[3] = {unit(random), unit(random), unit(random)};
                const int type = 1 + i % 3;
                const double r[3] = {s[0] * 40, -10 + s[1] * 40, 5 + s[2] * 20};
                if (frame == 0) { out << i + 1 << " " << type << " " << r[0] << " " << r[1] << " " << r[2] << "\n"; first.push_back({r[0], r[1], r[2], 0}); }
                else { out << i + 1 << " " << type << " " << s[0] << " " << s[1] << " " << s[2] << " 0.5\n"; second.push_back({r[0], r[1], r[2], 0}); }
            }
        }
        out.close();
        check(cloud::peekAtomCount(dump) == 5000, "peek");
        cloud::BuildOptions options;
        options.chunkAtoms = 1024;
        const auto a = folder() / "dump0.afcloud", b = folder() / "dump1.afcloud";
        cloud::buildFromLammpsDump(dump, a, {"Cu", "Ni"}, 0, {}, options);
        cloud::buildFromLammpsDump(dump, b, {}, 1, {}, options);
        const cloud::CloudFile fa(a), fb(b);
        check(fa.info().atoms == 5000 && fb.info().atoms == 5000, "atom counts");
        check(fa.info().species.size() == 3 && fa.info().species[0].name == "Cu" && fa.info().species[1].name == "Ni" && fa.info().species[2].name == "type 3", "type names");
        check(fa.info().hasCell && std::abs(fa.info().cellOrigin[1] + 10) < 1e-12 && std::abs(fa.info().cell[2][2] - 20) < 1e-12, "dump box");
        check(sameAtoms(readAll(fa), first, 2e-3), "frame 0 positions");
        check(sameAtoms(readAll(fb), second, 2e-3), "frame 1 (scaled) positions");
    });

    test("XYZ and extended XYZ", [] {
        const auto xyz = folder() / "atoms.xyz";
        std::ofstream out(xyz);
        out << "4\nLattice=\"10 0 0 0 10 0 0 0 10\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n"
               "O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\nO 5 5 5\n";
        out.close();
        check(cloud::peekAtomCount(xyz) == 4, "peek");
        const auto path = folder() / "xyz.afcloud";
        cloud::buildFromXyz(xyz, path);
        const cloud::CloudFile file(path);
        check(file.info().atoms == 4 && file.info().species.size() == 2 && file.info().hasCell, "atoms, species and lattice");
        std::map<std::string, int> counts;
        for (const auto& a : readAll(file)) ++counts[file.info().species[static_cast<std::size_t>(a[3])].name];
        check(counts["O"] == 2 && counts["H"] == 2, "species per atom");
    });

    test("oversized buckets split on disk", [] {
        // All atoms in a tiny region of a huge declared box: one bucket gets everything.
        const auto path = folder() / "clustered.afcloud";
        cloud::BuildOptions options;
        options.chunkAtoms = 2048;
        cloud::CloudBuilder builder(path, {0, 0, 0}, {1e6, 1e6, 1e6}, 1u << 22, options);
        std::mt19937 random(3);
        std::normal_distribution<double> spread(0.0, 5.0);
        const auto species = builder.species("Al");
        for (int i = 0; i < 100000; ++i) builder.add({500 + spread(random), 500 + spread(random), 500 + spread(random), species});
        const auto info = builder.finish();
        check(info.atoms == 100000, "every atom");
        for (const auto& c : info.chunks) check(c.count <= 2048, "chunk limit after splitting");
    });

    test("errors", [] {
        bool thrown = false;
        try { cloud::CloudFile file(folder() / "missing.afcloud"); } catch (const std::exception&) { thrown = true; }
        check(thrown, "missing file");
        thrown = false;
        try { cloud::generateCrystal(folder() / "bad.afcloud", "quasicrystal", 3, 0, {1, 1, 1}, "Cu"); } catch (const std::invalid_argument&) { thrown = true; }
        check(thrown, "unknown lattice");
    });

    std::error_code ignored;
    std::filesystem::remove_all(folder(), ignored);
    if (failures) { std::cout << failures << " cloud regression(s) failed\n"; return 1; }
    std::cout << "All atom cloud regressions passed\n";
    return 0;
}
