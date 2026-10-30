#include "pipeline/MoleculeLibrary.h"

#include "science/ScienceData.h"

#include "util/ElementData.h"

#include <filesystem>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>

namespace atomforge::pipeline
{
namespace
{
using science::Vec3;
constexpr double kDegree = 3.14159265358979323846 / 180.0;
// The tetrahedral angle's half-angle helpers.
const double kTetrahedral = std::acos(-1.0 / 3.0);

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

using Atoms = std::vector<MoleculeAtom>;

// Four tetrahedral unit directions (alternate cube corners).
std::vector<Vec3> tetrahedron()
{
    const double s = 1.0 / std::sqrt(3.0);
    return {{s, s, s}, {s, -s, -s}, {-s, s, -s}, {-s, -s, s}};
}

Atoms diatomic(const std::string& a, const std::string& b, double length)
{
    return {{a, {0, 0, 0}}, {b, {0, 0, length}}};
}

// A bent X-H2 or X-Y2 (apex at the origin, molecule in the xz plane).
Atoms bent(const std::string& apex, const std::string& arm, double length, double angleDegrees)
{
    const double half = 0.5 * angleDegrees * kDegree;
    return {{apex, {0, 0, 0}}, {arm, {length * std::sin(half), 0, length * std::cos(half)}},
            {arm, {-length * std::sin(half), 0, length * std::cos(half)}}};
}

// XH3 pyramid with H-X-H angle (C3 axis along -z).
Atoms pyramid(const std::string& apex, const std::string& arm, double length, double angleDegrees)
{
    // Angle between the bonds and the axis from the H-X-H angle.
    const double c = std::cos(angleDegrees * kDegree);
    const double sinBeta = std::sqrt(2.0 * (1.0 - c) / 3.0);
    const double cosBeta = std::sqrt(1.0 - sinBeta * sinBeta);
    Atoms out{{apex, {0, 0, 0}}};
    for (int k = 0; k < 3; ++k) {
        const double phi = 2.0 * 3.14159265358979323846 * k / 3.0;
        out.push_back({arm, {length * sinBeta * std::cos(phi), length * sinBeta * std::sin(phi), -length * cosBeta}});
    }
    return out;
}

Atoms tetrahedralHydride(const std::string& centre, const std::string& arm, double length)
{
    Atoms out{{centre, {0, 0, 0}}};
    for (const Vec3& d : tetrahedron()) out.push_back({arm, science::scale(d, length)});
    return out;
}

// Ethane (staggered): C-C along z, C-H tetrahedral.
Atoms ethane()
{
    const double cc = 1.535, ch = 1.094, hcc = 111.2 * kDegree;
    Atoms out{{"C", {0, 0, -cc / 2}}, {"C", {0, 0, cc / 2}}};
    for (int side = 0; side < 2; ++side) {
        const double sign = side == 0 ? -1.0 : 1.0;
        for (int k = 0; k < 3; ++k) {
            const double phi = 2.0 * 3.14159265358979323846 * k / 3.0 + (side ? 3.14159265358979323846 / 3.0 : 0.0);
            // H on carbon `side`, pointing away from the other carbon.
            const double z = sign * (cc / 2 - ch * std::cos(hcc));
            out.push_back({"H", {ch * std::sin(hcc) * std::cos(phi), ch * std::sin(hcc) * std::sin(phi), z}});
        }
    }
    return out;
}

Atoms ethylene()
{
    const double cc = 1.339, ch = 1.086, hcc = 121.3 * kDegree;
    Atoms out{{"C", {0, 0, -cc / 2}}, {"C", {0, 0, cc / 2}}};
    for (int side = 0; side < 2; ++side) {
        const double sign = side == 0 ? -1.0 : 1.0;
        for (int k = 0; k < 2; ++k) {
            const double x = (k ? 1.0 : -1.0) * ch * std::sin(hcc);
            out.push_back({"H", {x, 0, sign * (cc / 2 - ch * std::cos(hcc))}});
        }
    }
    return out;
}

Atoms acetylene()
{
    const double cc = 1.203, ch = 1.063;
    return {{"C", {0, 0, -cc / 2}}, {"C", {0, 0, cc / 2}}, {"H", {0, 0, -cc / 2 - ch}}, {"H", {0, 0, cc / 2 + ch}}};
}

Atoms benzene()
{
    const double cc = 1.397, ch = 1.084;
    Atoms out;
    for (int k = 0; k < 6; ++k) {
        const double phi = k * 60.0 * kDegree;
        out.push_back({"C", {cc * std::cos(phi), cc * std::sin(phi), 0}});
    }
    for (int k = 0; k < 6; ++k) {
        const double phi = k * 60.0 * kDegree;
        out.push_back({"H", {(cc + ch) * std::cos(phi), (cc + ch) * std::sin(phi), 0}});
    }
    return out;
}

// Methanol: C-O along z, hydroxyl H in the xz plane, methyl H tetrahedral.
Atoms methanol()
{
    const double co = 1.427, oh = 0.956, ch = 1.096, coh = 108.5 * kDegree, hco = 109.5 * kDegree;
    Atoms out{{"C", {0, 0, 0}}, {"O", {0, 0, co}}, {"H", {oh * std::sin(coh), 0, co - oh * std::cos(coh)}}};
    for (int k = 0; k < 3; ++k) {
        const double phi = 3.14159265358979323846 + 2.0 * 3.14159265358979323846 * k / 3.0;
        out.push_back({"H", {ch * std::sin(hco) * std::cos(phi), ch * std::sin(hco) * std::sin(phi), ch * std::cos(hco)}});
    }
    return out;
}

// C60 (truncated icosahedron): the even permutations of (0, +-1, +-3p),
// (+-1, +-(2+p), +-2p), (+-p, +-2, +-p^3), scaled to a mean bond of 1.43 A
// (between the experimental 1.40 and 1.45 A).
Atoms fullerene()
{
    const double p = (1.0 + std::sqrt(5.0)) / 2.0;
    const std::vector<Vec3> seeds = {{0, 1, 3 * p}, {1, 2 + p, 2 * p}, {p, 2, p * p * p}};
    std::vector<Vec3> points;
    for (const Vec3& seed : seeds)
        for (int sx = -1; sx <= 1; sx += 2)
            for (int sy = -1; sy <= 1; sy += 2)
                for (int sz = -1; sz <= 1; sz += 2) {
                    const Vec3 v{seed[0] * sx, seed[1] * sy, seed[2] * sz};
                    // Even (cyclic) permutations.
                    for (const Vec3& q : {v, Vec3{v[1], v[2], v[0]}, Vec3{v[2], v[0], v[1]}}) {
                        bool known = false;
                        for (const Vec3& r : points)
                            if (science::norm(science::sub(q, r)) < 1e-6) known = true;
                        if (!known) points.push_back(q);
                    }
                }
    // Edge length 2 in these coordinates.
    Atoms out;
    for (const Vec3& q : points) out.push_back({"C", science::scale(q, 1.43 / 2.0)});
    return out;
}

const std::map<std::string, Atoms (*)()>& table()
{
    static const std::map<std::string, Atoms (*)()> molecules = {
        {"H2", [] { return diatomic("H", "H", 0.741); }},
        {"N2", [] { return diatomic("N", "N", 1.098); }},
        {"O2", [] { return diatomic("O", "O", 1.208); }},
        {"F2", [] { return diatomic("F", "F", 1.412); }},
        {"Cl2", [] { return diatomic("Cl", "Cl", 1.988); }},
        {"CO", [] { return diatomic("C", "O", 1.128); }},
        {"NO", [] { return diatomic("N", "O", 1.151); }},
        {"HF", [] { return diatomic("H", "F", 0.917); }},
        {"HCl", [] { return diatomic("H", "Cl", 1.275); }},
        {"H2O", [] { return bent("O", "H", 0.9572, 104.52); }},
        {"H2S", [] { return bent("S", "H", 1.336, 92.11); }},
        {"CO2", [] { return Atoms{{"C", {0, 0, 0}}, {"O", {0, 0, 1.160}}, {"O", {0, 0, -1.160}}}; }},
        {"SO2", [] { return bent("S", "O", 1.431, 119.3); }},
        {"NO2", [] { return bent("N", "O", 1.193, 134.1); }},
        {"O3", [] { return bent("O", "O", 1.278, 116.8); }},
        {"NH3", [] { return pyramid("N", "H", 1.012, 106.67); }},
        {"PH3", [] { return pyramid("P", "H", 1.420, 93.3); }},
        {"CH4", [] { return tetrahedralHydride("C", "H", 1.087); }},
        {"SiH4", [] { return tetrahedralHydride("Si", "H", 1.480); }},
        {"CF4", [] { return tetrahedralHydride("C", "F", 1.315); }},
        {"C2H2", acetylene},
        {"C2H4", ethylene},
        {"C2H6", ethane},
        {"CH3OH", methanol},
        {"C6H6", benzene},
        {"HCN", [] { return Atoms{{"H", {0, 0, -1.064}}, {"C", {0, 0, 0}}, {"N", {0, 0, 1.156}}}; }},
        {"H2CO", [] {
             const double co = 1.205, ch = 1.111, half = 0.5 * 116.1 * kDegree;
             return Atoms{{"C", {0, 0, 0}}, {"O", {0, 0, co}}, {"H", {ch * std::sin(half), 0, -ch * std::cos(half)}},
                          {"H", {-ch * std::sin(half), 0, -ch * std::cos(half)}}};
         }},
        {"C60", fullerene},
    };
    return molecules;
}

const std::map<std::string, std::string>& aliases()
{
    static const std::map<std::string, std::string> names = {
        {"water", "H2O"}, {"ammonia", "NH3"}, {"methane", "CH4"}, {"silane", "SiH4"}, {"ethane", "C2H6"},
        {"ethylene", "C2H4"}, {"ethene", "C2H4"}, {"acetylene", "C2H2"}, {"ethyne", "C2H2"}, {"benzene", "C6H6"},
        {"methanol", "CH3OH"}, {"fullerene", "C60"}, {"formaldehyde", "H2CO"}, {"hydrogen cyanide", "HCN"}, {"buckminsterfullerene", "C60"}, {"ozone", "O3"},
        {"hydrogen", "H2"}, {"nitrogen", "N2"}, {"oxygen", "O2"}, {"carbon dioxide", "CO2"}, {"carbon monoxide", "CO"},
    };
    return names;
}

void centre(Atoms& atoms)
{
    double total = 0;
    Vec3 c{0, 0, 0};
    for (const auto& a : atoms) {
        const double m = elementAtomicMass(science::atomicNumber(a.symbol));
        total += m;
        c = science::add(c, science::scale(a.position, m));
    }
    c = science::scale(c, 1.0 / total);
    for (auto& a : atoms) a.position = science::sub(a.position, c);
}
}  // namespace

std::vector<std::string> moleculeNames()
{
    std::vector<std::string> names;
    for (const auto& entry : table()) names.push_back(entry.first);
    return names;
}

std::vector<MoleculeAtom> molecule(const std::string& name)
{
    const std::string key = lower(name);
    for (const auto& entry : table())
        if (lower(entry.first) == key) {
            Atoms atoms = entry.second();
            centre(atoms);
            return atoms;
        }
    if (const auto alias = aliases().find(key); alias != aliases().end()) return molecule(alias->second);
    // A structure file (XYZ, POSCAR, ...): its atoms as the molecule.
    if (name.find_first_of("./\\") != std::string::npos && std::filesystem::exists(std::filesystem::u8path(name))) {
        const auto frames = science::readFrames(std::filesystem::u8path(name));
        Atoms atoms;
        for (const auto& a : frames.front().structure.atoms) atoms.push_back({a.symbol, {a.x, a.y, a.z}});
        if (atoms.empty()) throw std::runtime_error("The molecule file " + name + " has no atoms");
        centre(atoms);
        return atoms;
    }
    // A single element.
    std::string symbol = name;
    if (!symbol.empty()) {
        symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
        for (std::size_t k = 1; k < symbol.size(); ++k) symbol[k] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[k])));
        if (science::atomicNumber(symbol) > 0) return {{symbol, {0, 0, 0}}};
    }
    std::string known;
    for (const auto& n : moleculeNames()) known += (known.empty() ? "" : ", ") + n;
    throw std::runtime_error("Unknown molecule '" + name + "'; give an element symbol or one of " + known);
}

double moleculeMass(const std::vector<MoleculeAtom>& atoms)
{
    double mass = 0;
    for (const auto& a : atoms) mass += elementAtomicMass(science::atomicNumber(a.symbol));
    return mass;
}
}  // namespace atomforge::pipeline

namespace atomforge::pipeline
{
std::array<science::Vec3, 3> principalAxes(const std::vector<MoleculeAtom>& atoms)
{
    // Inertia tensor about the centre of mass, diagonalised by Jacobi rotations.
    double total = 0;
    science::Vec3 c{0, 0, 0};
    for (const auto& a : atoms) {
        const double m = elementAtomicMass(science::atomicNumber(a.symbol));
        total += m;
        c = science::add(c, science::scale(a.position, m));
    }
    c = science::scale(c, 1.0 / total);
    double t[3][3] = {};
    for (const auto& a : atoms) {
        const double m = elementAtomicMass(science::atomicNumber(a.symbol));
        const science::Vec3 r = science::sub(a.position, c);
        const double r2 = science::dot(r, r);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) t[i][j] += m * ((i == j ? r2 : 0.0) - r[i] * r[j]);
    }
    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = std::abs(t[0][1]) + std::abs(t[0][2]) + std::abs(t[1][2]);
        if (off < 1e-14) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(t[p][q]) < 1e-18) continue;
                const double theta = 0.5 * std::atan2(2 * t[p][q], t[q][q] - t[p][p]);
                const double cs = std::cos(theta), sn = std::sin(theta);
                for (int k = 0; k < 3; ++k) {
                    const double tkp = t[k][p], tkq = t[k][q];
                    t[k][p] = cs * tkp - sn * tkq;
                    t[k][q] = sn * tkp + cs * tkq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double tpk = t[p][k], tqk = t[q][k];
                    t[p][k] = cs * tpk - sn * tqk;
                    t[q][k] = sn * tpk + cs * tqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = cs * vkp - sn * vkq;
                    v[k][q] = sn * vkp + cs * vkq;
                }
            }
    }
    std::array<int, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](int a, int b) { return t[a][a] < t[b][b]; });
    std::array<science::Vec3, 3> axes{};
    for (int k = 0; k < 3; ++k) axes[static_cast<std::size_t>(k)] = {v[0][order[k]], v[1][order[k]], v[2][order[k]]};
    // Right-handed.
    if (science::dot(science::cross(axes[0], axes[1]), axes[2]) < 0) axes[2] = science::scale(axes[2], -1);
    return axes;
}
}
