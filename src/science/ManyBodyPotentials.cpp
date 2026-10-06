// Covalent and ionic potentials: Tersoff, Stillinger-Weber and Buckingham with
// Ewald-summed Coulomb interactions. Energies follow the LAMMPS definitions
// (pair_style tersoff, sw and buck/coul/long) so parameter files are shared.
//
// Each potential's energy is a sum of local terms of bond vectors d = r_j - r_i;
// with dE/dd known, the force on j is -dE/dd, the force on i is +dE/dd and the
// virial is the sum of sym(dE/dd (x) d), so forces and stress share one path.
#include "science/Potentials.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);
const double kCoulomb = 14.3996454784255;  // e^2 / (4 pi eps0) in eV Angstrom

// Accumulates forces and the virial from gradients with respect to bond vectors.
struct Accumulator
{
    PotentialResult result;
    Mat3 virial{};

    explicit Accumulator(std::size_t atoms) { result.forces.assign(atoms, {0, 0, 0}); }

    void bond(int i, int j, const Vec3& d, const Vec3& gradient)
    {
        auto& fi = result.forces[static_cast<std::size_t>(i)];
        auto& fj = result.forces[static_cast<std::size_t>(j)];
        fi = add(fi, gradient);
        fj = sub(fj, gradient);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) virial[a][b] += 0.5 * (gradient[a] * d[b] + gradient[b] * d[a]);
    }

    PotentialResult finish(const Configuration& c, bool stress)
    {
        if (stress) {
            if (!c.fullyPeriodic()) throw std::runtime_error("Stress requires a full periodic cell");
            const double volume = cellVolume(c.cell);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) result.stress[a][b] = virial[a][b] / volume;
            result.hasStress = true;
        }
        return result;
    }
};

// Full neighbour lists (both directions, periodic images included), grouped by atom.
std::vector<std::vector<Neighbor>> neighboursByAtom(const Configuration& c, double cutoff)
{
    const bool periodic = c.pbc[0] || c.pbc[1] || c.pbc[2];
    std::vector<std::vector<Neighbor>> lists(c.size());
    for (const auto& n : neighborList(c.positions, periodic ? c.cell : identity(), c.pbc, cutoff))
        lists[static_cast<std::size_t>(n.i)].push_back(n);
    return lists;
}

std::vector<std::string> tokensOf(const std::filesystem::path& file)
{
    std::ifstream in(file);
    if (!in) throw std::runtime_error("Cannot read potential file " + file.u8string());
    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(in, line)) {
        line = line.substr(0, line.find('#'));
        std::istringstream words(line);
        for (std::string word; words >> word;) tokens.push_back(word);
    }
    return tokens;
}

double number(const std::string& text, const std::string& file)
{
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(value)) throw std::invalid_argument(text);
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid number '" + text + "' in " + file);
    }
}

using Triplet = std::array<std::string, 3>;

template<class P>
const P& lookup(const std::map<Triplet, P>& table, const std::string& a, const std::string& b, const std::string& c, const char* kind)
{
    const auto found = table.find({a, b, c});
    if (found == table.end()) throw std::runtime_error(std::string(kind) + " has no parameters for " + a + "-" + b + "-" + c);
    return found->second;
}

template<class P>
void requireElements(const std::map<Triplet, P>& table, const Configuration& c, const char* kind)
{
    for (const auto& a : c.symbols)
        for (const auto& b : c.symbols) lookup(table, a, b, b, kind);
}

// ------------------------------------------------------------------ Tersoff
struct TersoffParameters
{
    double m = 3, gamma = 1, lambda3 = 0, c = 0, d = 1, h = 0, n = 1, beta = 0, lambda2 = 0, B = 0, R = 0, D = 0, lambda1 = 0, A = 0;
};

double cutoffFunction(double r, double R, double D, double& derivative)
{
    if (r < R - D) { derivative = 0; return 1; }
    if (r > R + D) { derivative = 0; return 0; }
    const double x = 0.5 * kPi * (r - R) / D;
    derivative = -0.25 * kPi / D * std::cos(x);
    return 0.5 - 0.5 * std::sin(x);
}

class Tersoff final : public Potential
{
public:
    Tersoff(std::map<Triplet, TersoffParameters> table, std::string label) : m_table(std::move(table)), m_label(std::move(label))
    {
        for (const auto& [key, p] : m_table) m_cutoff = std::max(m_cutoff, p.R + p.D);
    }

    std::string description() const override { return "Tersoff (" + m_label + ")"; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        requireElements(m_table, c, "Tersoff");
        Accumulator acc(c.size());
        const auto lists = neighboursByAtom(c, m_cutoff);
        for (std::size_t i = 0; i < c.size(); ++i) {
            const auto& list = lists[i];
            const std::string& si = c.symbols[i];
            for (std::size_t a = 0; a < list.size(); ++a) {
                const Neighbor& nj = list[a];
                const std::string& sj = c.symbols[static_cast<std::size_t>(nj.j)];
                const TersoffParameters& pij = lookup(m_table, si, sj, sj, "Tersoff");
                const Vec3 dij = nj.vector;
                const double rij = norm(dij);
                if (rij >= pij.R + pij.D) continue;
                double dfc = 0;
                const double fc = cutoffFunction(rij, pij.R, pij.D, dfc);
                // Bond order: zeta over the other neighbours k of i.
                double zeta = 0;
                Vec3 dZetaDij{0, 0, 0};
                std::vector<std::pair<std::size_t, Vec3>> dZetaDik;
                for (std::size_t b = 0; b < list.size(); ++b) {
                    if (b == a) continue;
                    const Neighbor& nk = list[b];
                    const std::string& sk = c.symbols[static_cast<std::size_t>(nk.j)];
                    const TersoffParameters& pijk = lookup(m_table, si, sj, sk, "Tersoff");
                    const Vec3 dik = nk.vector;
                    const double rik = norm(dik);
                    if (rik >= pijk.R + pijk.D) continue;
                    double dfck = 0;
                    const double fck = cutoffFunction(rik, pijk.R, pijk.D, dfck);
                    const double cosine = dot(dij, dik) / (rij * rik);
                    const double h = pijk.h - cosine;
                    const double denominator = pijk.d * pijk.d + h * h;
                    const double g = pijk.gamma * (1 + pijk.c * pijk.c / (pijk.d * pijk.d) - pijk.c * pijk.c / denominator);
                    const double dg = -pijk.gamma * 2 * pijk.c * pijk.c * h / (denominator * denominator);  // dg/dcos
                    const double l3m = std::pow(pijk.lambda3, pijk.m);
                    const double delta = rij - rik;
                    const double ex = std::exp(l3m * std::pow(delta, pijk.m));
                    const double dex = ex * l3m * pijk.m * std::pow(delta, pijk.m - 1);  // d ex / d(rij - rik)
                    zeta += fck * g * ex;
                    const Vec3 uij = scale(dij, 1 / rij), uik = scale(dik, 1 / rik);
                    const Vec3 dcosDij = sub(scale(dik, 1 / (rij * rik)), scale(dij, cosine / (rij * rij)));
                    const Vec3 dcosDik = sub(scale(dij, 1 / (rij * rik)), scale(dik, cosine / (rik * rik)));
                    dZetaDij = add(dZetaDij, add(scale(uij, fck * g * dex), scale(dcosDij, fck * dg * ex)));
                    dZetaDik.push_back({b, add(add(scale(uik, dfck * g * ex - fck * g * dex), scale(dcosDik, fck * dg * ex)), Vec3{0, 0, 0})});
                }
                double bij = 1, dbDzeta = 0;
                if (zeta > 0 && pij.beta > 0) {
                    const double t = std::pow(pij.beta * zeta, pij.n);
                    bij = std::pow(1 + t, -0.5 / pij.n);
                    dbDzeta = -0.5 * std::pow(1 + t, -0.5 / pij.n - 1) * t / zeta;
                }
                const double fr = pij.A * std::exp(-pij.lambda1 * rij), dfr = -pij.lambda1 * fr;
                const double fa = -pij.B * std::exp(-pij.lambda2 * rij), dfa = -pij.lambda2 * fa;
                acc.result.energy += 0.5 * fc * (fr + bij * fa);
                const double radial = 0.5 * (dfc * (fr + bij * fa) + fc * (dfr + bij * dfa));
                const double bondOrder = 0.5 * fc * fa * dbDzeta;
                acc.bond(static_cast<int>(i), nj.j, dij, add(scale(dij, radial / rij), scale(dZetaDij, bondOrder)));
                for (const auto& [b, gradient] : dZetaDik) acc.bond(static_cast<int>(i), list[b].j, list[b].vector, scale(gradient, bondOrder));
            }
        }
        return acc.finish(c, stress);
    }

private:
    std::map<Triplet, TersoffParameters> m_table;
    std::string m_label;
    double m_cutoff = 0;
};

std::map<Triplet, TersoffParameters> builtinTersoffSilicon()
{
    // Tersoff, Phys. Rev. B 38, 9902 (1988), as distributed with LAMMPS (Si.tersoff).
    TersoffParameters p;
    p.m = 3; p.gamma = 1; p.lambda3 = 0; p.c = 1.0039e5; p.d = 16.217; p.h = -0.59825; p.n = 0.78734;
    p.beta = 1.1e-6; p.lambda2 = 1.7322; p.B = 471.18; p.R = 2.85; p.D = 0.15; p.lambda1 = 2.4799; p.A = 1830.8;
    return {{{"Si", "Si", "Si"}, p}};
}

std::map<Triplet, TersoffParameters> readTersoff(const std::filesystem::path& file)
{
    const auto tokens = tokensOf(file);
    if (tokens.empty() || tokens.size() % 17) throw std::runtime_error("A Tersoff file needs 17 entries per line: " + file.u8string());
    std::map<Triplet, TersoffParameters> table;
    for (std::size_t t = 0; t < tokens.size(); t += 17) {
        std::array<double, 14> v{};
        for (std::size_t k = 0; k < 14; ++k) v[k] = number(tokens[t + 3 + k], file.u8string());
        TersoffParameters p{v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13]};
        if (p.m != 1 && p.m != 3) throw std::runtime_error("Tersoff m must be 1 or 3");
        table[{tokens[t], tokens[t + 1], tokens[t + 2]}] = p;
    }
    return table;
}

// --------------------------------------------------------- Stillinger-Weber
struct SwParameters
{
    double epsilon = 1, sigma = 1, a = 1, lambda = 0, gamma = 1, cos0 = 0, A = 0, B = 0, p = 4, q = 0, tol = 0;
};

class StillingerWeber final : public Potential
{
public:
    StillingerWeber(std::map<Triplet, SwParameters> table, std::string label) : m_table(std::move(table)), m_label(std::move(label))
    {
        for (const auto& [key, p] : m_table) m_cutoff = std::max(m_cutoff, p.a * p.sigma);
    }

    std::string description() const override { return "Stillinger-Weber (" + m_label + ")"; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        requireElements(m_table, c, "Stillinger-Weber");
        Accumulator acc(c.size());
        const auto lists = neighboursByAtom(c, m_cutoff);
        for (std::size_t i = 0; i < c.size(); ++i) {
            const auto& list = lists[i];
            const std::string& si = c.symbols[i];
            // Two-body term, half per direction.
            for (const auto& nj : list) {
                const SwParameters& p = lookup(m_table, si, c.symbols[static_cast<std::size_t>(nj.j)], c.symbols[static_cast<std::size_t>(nj.j)], "Stillinger-Weber");
                const double r = norm(nj.vector), cut = p.a * p.sigma;
                if (r >= cut) continue;
                const double s = p.sigma / r, sp = std::pow(s, p.p), sq = std::pow(s, p.q);
                const double ex = std::exp(p.sigma / (r - cut));
                const double phi = p.A * p.epsilon * (p.B * sp - sq) * ex;
                const double dphi = p.A * p.epsilon * ex * ((-p.p * p.B * sp + p.q * sq) / r - (p.B * sp - sq) * p.sigma / ((r - cut) * (r - cut)));
                acc.result.energy += 0.5 * phi;
                acc.bond(static_cast<int>(i), nj.j, nj.vector, scale(nj.vector, 0.5 * dphi / r));
            }
            // Three-body term for each pair of neighbours j < k around i.
            for (std::size_t a = 0; a < list.size(); ++a)
                for (std::size_t b = a + 1; b < list.size(); ++b) {
                    const std::string& sj = c.symbols[static_cast<std::size_t>(list[a].j)];
                    const std::string& sk = c.symbols[static_cast<std::size_t>(list[b].j)];
                    const SwParameters& pij = lookup(m_table, si, sj, sj, "Stillinger-Weber");
                    const SwParameters& pik = lookup(m_table, si, sk, sk, "Stillinger-Weber");
                    const SwParameters& pijk = lookup(m_table, si, sj, sk, "Stillinger-Weber");
                    const Vec3 dij = list[a].vector, dik = list[b].vector;
                    const double rij = norm(dij), rik = norm(dik);
                    const double cutj = pij.a * pij.sigma, cutk = pik.a * pik.sigma;
                    if (rij >= cutj || rik >= cutk) continue;
                    const double ej = std::exp(pij.gamma * pij.sigma / (rij - cutj));
                    const double ek = std::exp(pik.gamma * pik.sigma / (rik - cutk));
                    const double cosine = dot(dij, dik) / (rij * rik);
                    const double delta = cosine - pijk.cos0;
                    const double le = pijk.lambda * pijk.epsilon;
                    acc.result.energy += le * delta * delta * ej * ek;
                    const double dEj = -le * delta * delta * ej * ek * pij.gamma * pij.sigma / ((rij - cutj) * (rij - cutj));
                    const double dEk = -le * delta * delta * ej * ek * pik.gamma * pik.sigma / ((rik - cutk) * (rik - cutk));
                    const double dEcos = 2 * le * delta * ej * ek;
                    const Vec3 dcosDij = sub(scale(dik, 1 / (rij * rik)), scale(dij, cosine / (rij * rij)));
                    const Vec3 dcosDik = sub(scale(dij, 1 / (rij * rik)), scale(dik, cosine / (rik * rik)));
                    acc.bond(static_cast<int>(i), list[a].j, dij, add(scale(dij, dEj / rij), scale(dcosDij, dEcos)));
                    acc.bond(static_cast<int>(i), list[b].j, dik, add(scale(dik, dEk / rik), scale(dcosDik, dEcos)));
                }
        }
        return acc.finish(c, stress);
    }

private:
    std::map<Triplet, SwParameters> m_table;
    std::string m_label;
    double m_cutoff = 0;
};

std::map<Triplet, SwParameters> builtinSwSilicon()
{
    // Stillinger and Weber, Phys. Rev. B 31, 5262 (1985), as in LAMMPS Si.sw.
    SwParameters p{2.1683, 2.0951, 1.80, 21.0, 1.20, -1.0 / 3.0, 7.049556277, 0.6022245584, 4.0, 0.0, 0.0};
    return {{{"Si", "Si", "Si"}, p}};
}

std::map<Triplet, SwParameters> readSw(const std::filesystem::path& file)
{
    const auto tokens = tokensOf(file);
    if (tokens.empty() || tokens.size() % 14) throw std::runtime_error("A Stillinger-Weber file needs 14 entries per line: " + file.u8string());
    std::map<Triplet, SwParameters> table;
    for (std::size_t t = 0; t < tokens.size(); t += 14) {
        std::array<double, 11> v{};
        for (std::size_t k = 0; k < 11; ++k) v[k] = number(tokens[t + 3 + k], file.u8string());
        table[{tokens[t], tokens[t + 1], tokens[t + 2]}] = SwParameters{v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10]};
    }
    return table;
}

// ------------------------------------------------- Buckingham + Coulomb
struct BuckinghamPair { double A = 0, rho = 1, C = 0; };

class Buckingham final : public Potential
{
public:
    Buckingham(std::map<std::pair<std::string, std::string>, BuckinghamPair> pairs, std::map<std::string, double> charges,
               double cutoff, double accuracy)
        : m_pairs(std::move(pairs)), m_charges(std::move(charges)), m_cutoff(cutoff), m_accuracy(accuracy) {}

    std::string description() const override { return m_charges.empty() ? "Buckingham" : "Buckingham + Coulomb (Ewald)"; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        Accumulator acc(c.size());
        std::vector<double> q(c.size(), 0.0);
        double total = 0;
        for (std::size_t i = 0; i < c.size(); ++i) {
            if (!m_charges.empty()) {
                const auto found = m_charges.find(c.symbols[i]);
                if (found == m_charges.end()) throw std::runtime_error("No charge given for " + c.symbols[i]);
                q[i] = found->second;
            }
            total += q[i];
        }
        const bool coulomb = !m_charges.empty();
        const bool periodic = c.fullyPeriodic();
        if (coulomb && periodic && std::abs(total) > 1e-8) throw std::runtime_error("Ewald summation needs a charge-neutral cell");
        if (coulomb && !periodic && (c.pbc[0] || c.pbc[1] || c.pbc[2])) throw std::runtime_error("Coulomb interactions need a fully periodic or an open structure");
        // Ewald splitting: real-space cutoff r_c with erfc(alpha r_c) ~ accuracy.
        const double alpha = coulomb && periodic ? std::sqrt(-std::log(m_accuracy)) / m_cutoff : 0.0;
        const double realCutoff = coulomb && !periodic ? 1e9 : m_cutoff;
        const bool open = !(c.pbc[0] || c.pbc[1] || c.pbc[2]);
        // Real-space pairs (each once).
        const auto pairs = [&] {
            std::vector<Neighbor> list;
            if (open && coulomb) {  // direct Coulomb sum over all pairs of an open structure
                for (std::size_t i = 0; i < c.size(); ++i)
                    for (std::size_t j = i + 1; j < c.size(); ++j) {
                        Neighbor n; n.i = static_cast<int>(i); n.j = static_cast<int>(j); n.vector = sub(c.positions[j], c.positions[i]);
                        list.push_back(n);
                    }
                return list;
            }
            for (const auto& n : neighborList(c.positions, open ? identity() : c.cell, c.pbc, realCutoff)) {
                bool keep = n.i < n.j;
                if (n.i == n.j) for (int k = 0; k < 3; ++k) if (n.shift[k] != 0) { keep = n.shift[k] > 0; break; }
                if (keep) list.push_back(n);
            }
            return list;
        }();
        for (const auto& n : pairs) {
            const double r = norm(n.vector);
            double energy = 0, derivative = 0;
            if (r < m_cutoff) {
                const BuckinghamPair* b = pairFor(c.symbols[static_cast<std::size_t>(n.i)], c.symbols[static_cast<std::size_t>(n.j)]);
                if (b) {
                    const double ex = b->A * std::exp(-r / b->rho), r6 = std::pow(r, 6);
                    energy += ex - b->C / r6;
                    derivative += -ex / b->rho + 6 * b->C / (r6 * r);
                }
            }
            if (coulomb) {
                const double qq = kCoulomb * q[static_cast<std::size_t>(n.i)] * q[static_cast<std::size_t>(n.j)];
                if (periodic) {
                    const double e = std::erfc(alpha * r);
                    energy += qq * e / r;
                    derivative += -qq * (e / (r * r) + 2 * alpha / std::sqrt(kPi) * std::exp(-alpha * alpha * r * r) / r);
                } else {
                    energy += qq / r;
                    derivative += -qq / (r * r);
                }
            }
            acc.result.energy += energy;
            acc.bond(n.i, n.j, n.vector, scale(n.vector, derivative / r));
        }
        if (coulomb && periodic) reciprocal(c, q, alpha, acc);
        return acc.finish(c, stress);
    }

private:
    const BuckinghamPair* pairFor(const std::string& a, const std::string& b) const
    {
        auto found = m_pairs.find({a, b});
        if (found == m_pairs.end()) found = m_pairs.find({b, a});
        return found == m_pairs.end() ? nullptr : &found->second;
    }

    // Reciprocal-space and self terms of the Ewald sum, with forces and virial.
    void reciprocal(const Configuration& c, const std::vector<double>& q, double alpha, Accumulator& acc) const
    {
        const double volume = cellVolume(c.cell);
        const Mat3 inv = inverse(c.cell);
        Mat3 b{};  // reciprocal rows (2 pi included)
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) b[r][k] = 2 * kPi * inv[k][r];
        const double kMax = 2 * alpha * std::sqrt(-std::log(m_accuracy));
        std::array<int, 3> limits{};
        for (int r = 0; r < 3; ++r) limits[static_cast<std::size_t>(r)] = static_cast<int>(std::ceil(kMax * norm(c.cell[r]) / (2 * kPi))) + 1;
        double self = 0;
        for (double qi : q) self += qi * qi;
        acc.result.energy -= kCoulomb * alpha / std::sqrt(kPi) * self;
        for (int h = -limits[0]; h <= limits[0]; ++h)
            for (int k = -limits[1]; k <= limits[1]; ++k)
                for (int l = -limits[2]; l <= limits[2]; ++l) {
                    if (h == 0 && k == 0 && l == 0) continue;
                    const Vec3 g = add(add(scale(b[0], h), scale(b[1], k)), scale(b[2], l));
                    const double g2 = dot(g, g);
                    if (g2 > kMax * kMax) continue;
                    double sc = 0, ss = 0;
                    std::vector<double> phase(c.size());
                    for (std::size_t i = 0; i < c.size(); ++i) {
                        phase[i] = dot(g, c.positions[i]);
                        sc += q[i] * std::cos(phase[i]);
                        ss += q[i] * std::sin(phase[i]);
                    }
                    const double prefactor = kCoulomb * 2 * kPi / volume * std::exp(-g2 / (4 * alpha * alpha)) / g2;
                    const double structure2 = sc * sc + ss * ss;
                    // E_g = (2 pi / V) exp(-g^2 / 4 alpha^2) / g^2 |S(g)|^2, summed over all g != 0.
                    const double eg = prefactor * structure2;
                    acc.result.energy += eg;
                    for (std::size_t i = 0; i < c.size(); ++i) {
                        // F_i = -dE/dr_i = 2 prefactor q_i g (sin(g.r_i) S_c - cos(g.r_i) S_s).
                        const double weight = 2 * prefactor * q[i] * (std::sin(phase[i]) * sc - std::cos(phase[i]) * ss);
                        acc.result.forces[i] = add(acc.result.forces[i], scale(g, weight));
                    }
                    // dE_g/d strain = E_g (-delta_ab + 2 (1 + g^2 / 4 alpha^2) g_a g_b / g^2).
                    const double factor = 2 * (1 + g2 / (4 * alpha * alpha)) / g2;
                    for (int a = 0; a < 3; ++a)
                        for (int d = 0; d < 3; ++d) acc.virial[a][d] += eg * (factor * g[a] * g[d] - (a == d ? 1.0 : 0.0));
                }
    }

    std::map<std::pair<std::string, std::string>, BuckinghamPair> m_pairs;
    std::map<std::string, double> m_charges;
    double m_cutoff, m_accuracy;
};
}

std::unique_ptr<Potential> makeTersoff(const std::filesystem::path& file)
{
    if (file.empty()) return std::make_unique<Tersoff>(builtinTersoffSilicon(), "Si, Tersoff 1988");
    return std::make_unique<Tersoff>(readTersoff(file), file.filename().u8string());
}

std::unique_ptr<Potential> makeStillingerWeber(const std::filesystem::path& file)
{
    if (file.empty()) return std::make_unique<StillingerWeber>(builtinSwSilicon(), "Si, Stillinger-Weber 1985");
    return std::make_unique<StillingerWeber>(readSw(file), file.filename().u8string());
}

std::unique_ptr<Potential> makeBuckingham(const Json& options)
{
    std::map<std::pair<std::string, std::string>, BuckinghamPair> pairs;
    if (const Json* list = options.find("pairs"); list && list->isArray())
        for (const auto& entry : list->items()) {
            const auto& elements = entry.at("elements").items();
            if (elements.size() != 2) throw std::runtime_error("Each Buckingham pair needs two elements");
            BuckinghamPair pair{entry.at("A").number(), entry.at("rho").number(), entry.contains("C") ? entry.at("C").number() : 0.0};
            if (!(pair.rho > 0)) throw std::runtime_error("Buckingham rho must be positive");
            pairs[{elements[0].string(), elements[1].string()}] = pair;
        }
    std::map<std::string, double> charges;
    if (const Json* table = options.find("charges"); table && table->isObject())
        for (const auto& [element, charge] : table->members()) charges[element] = charge.number();
    if (pairs.empty() && charges.empty()) throw std::runtime_error("The Buckingham potential needs \"pairs\" and/or \"charges\"");
    const double cutoff = options.contains("cutoff") ? options.at("cutoff").number() : 10.0;
    const double accuracy = options.contains("ewald_accuracy") ? options.at("ewald_accuracy").number() : 1e-6;
    if (!(cutoff > 0) || !(accuracy > 0 && accuracy < 1)) throw std::runtime_error("Buckingham cutoff must be positive and ewald_accuracy in (0, 1)");
    return std::make_unique<Buckingham>(std::move(pairs), std::move(charges), cutoff, accuracy);
}
}
