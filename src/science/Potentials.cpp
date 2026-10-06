#include "science/Potentials.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <array>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
Configuration configurationFrom(const Structure& structure, const Pbc& pbc)
{
    Configuration result;
    for (const auto& atom : structure.atoms) {
        result.symbols.push_back(atom.symbol);
        result.numbers.push_back(atom.atomicNumber ? atom.atomicNumber : atomicNumber(atom.symbol));
        result.positions.push_back({atom.x, atom.y, atom.z});
        result.masses.push_back(atomicMass(atom.symbol));
    }
    if (structure.hasUnitCell)
        for (int r = 0; r < 3; ++r) result.cell[r] = structure.cellVectors[r];
    result.pbc = structure.hasUnitCell ? pbc : Pbc{false, false, false};
    return result;
}

Structure structureFrom(const Configuration& configuration)
{
    Structure structure;
    for (std::size_t i = 0; i < configuration.size(); ++i) {
        AtomSite atom;
        atom.symbol = configuration.symbols[i];
        atom.atomicNumber = configuration.numbers[i];
        atom.x = configuration.positions[i][0];
        atom.y = configuration.positions[i][1];
        atom.z = configuration.positions[i][2];
        getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
        structure.atoms.push_back(atom);
    }
    if (std::abs(determinant(configuration.cell)) > 1e-12) {
        structure.hasUnitCell = true;
        for (int r = 0; r < 3; ++r) structure.cellVectors[r] = configuration.cell[r];
    }
    return structure;
}

namespace
{
// Each unordered pair once, including periodic images (half neighbour list).
std::vector<Neighbor> halfList(const Configuration& configuration, double cutoff)
{
    const bool periodic = configuration.pbc[0] || configuration.pbc[1] || configuration.pbc[2];
    const Mat3 cell = periodic ? configuration.cell : identity();
    std::vector<Neighbor> result;
    for (const auto& neighbor : neighborList(configuration.positions, cell, configuration.pbc, cutoff)) {
        bool keep = neighbor.i < neighbor.j;
        if (neighbor.i == neighbor.j) {
            for (int k = 0; k < 3; ++k)
                if (neighbor.shift[k] != 0) { keep = neighbor.shift[k] > 0; break; }
        }
        if (keep) result.push_back(neighbor);
    }
    return result;
}

void addVirial(Mat3& virial, const Vec3& forceOnI, const Vec3& d)
{
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) virial[a][b] += 0.5 * (forceOnI[a] * d[b] + forceOnI[b] * d[a]);
}

void finishStress(PotentialResult& result, const Configuration& configuration, const Mat3& virial, bool stress)
{
    if (!stress) return;
    if (!configuration.fullyPeriodic()) throw std::runtime_error("Stress requires a full periodic cell");
    const double volume = cellVolume(configuration.cell);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) result.stress[a][b] = virial[a][b] / volume;
    result.hasStress = true;
}

const double kBohr = 0.52917721067;
// ASE's historical rounding of (16 pi / 3)^(1/3) / sqrt(2).
const double kBeta = 1.809;

struct EmtRaw { const char* symbol; double e0, s0, v0, eta2, kappa, lambda, n0; };
// E0 (eV), s0 (bohr), V0 (eV), eta2, kappa, lambda (1/bohr), n0 (1/bohr^3).
const EmtRaw kEmtParameters[] = {
    {"Al", -3.28, 3.00, 1.493, 1.240, 2.000, 1.169, 0.00700},
    {"Cu", -3.51, 2.67, 2.476, 1.652, 2.740, 1.906, 0.00910},
    {"Ag", -2.96, 3.01, 2.132, 1.652, 2.790, 1.892, 0.00547},
    {"Au", -3.80, 3.00, 2.321, 1.674, 2.873, 2.182, 0.00703},
    {"Ni", -4.44, 2.60, 3.673, 1.669, 2.757, 1.948, 0.01030},
    {"Pd", -3.90, 2.87, 2.773, 1.818, 3.107, 2.155, 0.00688},
    {"Pt", -5.85, 2.90, 4.067, 1.812, 3.145, 2.192, 0.00802},
    {"H", -3.21, 1.31, 0.132, 2.652, 2.790, 3.892, 0.00547},
    {"C", -3.50, 1.81, 0.332, 1.652, 2.790, 1.892, 0.01322},
    {"N", -5.10, 1.88, 0.132, 1.652, 2.790, 1.892, 0.01222},
    {"O", -4.60, 1.95, 0.332, 1.652, 2.790, 1.892, 0.00850},
};

struct EmtElement { double e0, s0, v0, eta2, kappa, lambda, n0, gamma1, gamma2; };

class Emt final : public Potential
{
public:
    Emt()
    {
        double maxS0 = 0;
        for (const auto& raw : kEmtParameters) maxS0 = std::max(maxS0, raw.s0 * kBohr);
        m_rc = kBeta * maxS0 * 0.5 * (std::sqrt(3.0) + std::sqrt(4.0));
        const double rr = m_rc * 2 * std::sqrt(4.0) / (std::sqrt(3.0) + std::sqrt(4.0));
        m_acut = std::log(9999.0) / (rr - m_rc);
        m_listCutoff = m_rc + 0.5;
        for (const auto& raw : kEmtParameters) {
            EmtElement p{raw.e0, raw.s0 * kBohr, raw.v0, raw.eta2 / kBohr, raw.kappa / kBohr, raw.lambda / kBohr,
                         raw.n0 / (kBohr * kBohr * kBohr), 0, 0};
            const int shells[3] = {12, 6, 24};
            for (int i = 0; i < 3; ++i) {
                const double r = p.s0 * kBeta * std::sqrt(i + 1.0);
                const double x = shells[i] / (12 * (1.0 + std::exp(m_acut * (r - m_rc))));
                p.gamma1 += x * std::exp(-p.eta2 * (r - kBeta * p.s0));
                p.gamma2 += x * std::exp(-p.kappa / kBeta * (r - kBeta * p.s0));
            }
            m_elements[raw.symbol] = p;
        }
    }

    std::string description() const override { return "EMT (effective-medium theory)"; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        const std::size_t n = c.size();
        std::vector<const EmtElement*> parameters(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto found = m_elements.find(c.symbols[i]);
            if (found == m_elements.end()) throw std::runtime_error("No EMT potential for " + c.symbols[i] + "; EMT supports Al, Cu, Ag, Au, Ni, Pd, Pt, H, C, N and O");
            parameters[i] = &found->second;
        }
        const auto pairs = halfList(c, m_listCutoff);
        PotentialResult result;
        result.forces.assign(n, {0, 0, 0});
        Mat3 virial{};
        std::vector<double> sigma1(n, 0.0), deds(n, 0.0);
        for (const auto& pair : pairs) {
            const auto a1 = static_cast<std::size_t>(pair.i), a2 = static_cast<std::size_t>(pair.j);
            const EmtElement& p1 = *parameters[a1];
            const EmtElement& p2 = *parameters[a2];
            const double ksi = p2.n0 / p1.n0;
            const Vec3& d = pair.vector;
            const double r = norm(d);
            const double x = std::exp(m_acut * (r - m_rc));
            const double theta = 1.0 / (1.0 + x);
            const double y1 = 0.5 * p1.v0 * std::exp(-p2.kappa * (r / kBeta - p2.s0)) * ksi / p1.gamma2 * theta;
            const double y2 = 0.5 * p2.v0 * std::exp(-p1.kappa * (r / kBeta - p1.s0)) / ksi / p2.gamma2 * theta;
            result.energy -= y1 + y2;
            const Vec3 f = scale(d, ((y1 * p2.kappa + y2 * p1.kappa) / kBeta + (y1 + y2) * m_acut * theta * x) / r);
            result.forces[a1] = add(result.forces[a1], f);
            result.forces[a2] = sub(result.forces[a2], f);
            addVirial(virial, f, d);
            sigma1[a1] += std::exp(-p2.eta2 * (r - kBeta * p2.s0)) * ksi * theta / p1.gamma1;
            sigma1[a2] += std::exp(-p1.eta2 * (r - kBeta * p1.s0)) / ksi * theta / p2.gamma1;
        }
        for (std::size_t a = 0; a < n; ++a) {
            const EmtElement& p = *parameters[a];
            if (!(sigma1[a] > 0)) {
                deds[a] = 0.0;
                result.energy -= p.e0;
                continue;
            }
            const double ds = -std::log(sigma1[a] / 12) / (kBeta * p.eta2);
            const double x = p.lambda * ds;
            const double y = std::exp(-x);
            const double z = 6 * p.v0 * std::exp(-p.kappa * ds);
            deds[a] = (x * y * p.e0 * p.lambda + p.kappa * z) / (sigma1[a] * kBeta * p.eta2);
            result.energy += p.e0 * ((1 + x) * y - 1) + z;
        }
        for (const auto& pair : pairs) {
            const auto a1 = static_cast<std::size_t>(pair.i), a2 = static_cast<std::size_t>(pair.j);
            const EmtElement& p1 = *parameters[a1];
            const EmtElement& p2 = *parameters[a2];
            const double ksi = p2.n0 / p1.n0;
            const Vec3& d = pair.vector;
            const double r = norm(d);
            const double x = std::exp(m_acut * (r - m_rc));
            const double theta = 1.0 / (1.0 + x);
            const double y1 = std::exp(-p2.eta2 * (r - kBeta * p2.s0)) * ksi / p1.gamma1 * theta * deds[a1];
            const double y2 = std::exp(-p1.eta2 * (r - kBeta * p1.s0)) / ksi / p2.gamma1 * theta * deds[a2];
            const Vec3 f = scale(d, ((y1 * p2.eta2 + y2 * p1.eta2) + (y1 + y2) * m_acut * theta * x) / r);
            result.forces[a1] = sub(result.forces[a1], f);
            result.forces[a2] = add(result.forces[a2], f);
            addVirial(virial, scale(f, -1), d);
        }
        finishStress(result, c, virial, stress);
        return result;
    }

private:
    double m_rc = 0, m_acut = 0, m_listCutoff = 0;
    std::map<std::string, EmtElement> m_elements;
};

class LennardJones final : public Potential
{
public:
    LennardJones(double epsilon, double sigma, double cutoff) : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
    {
        const double ratio6 = std::pow(sigma / cutoff, 6);
        m_shift = 4 * epsilon * (ratio6 * ratio6 - ratio6);
    }

    std::string description() const override { return "Lennard-Jones"; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        PotentialResult result;
        result.forces.assign(c.size(), {0, 0, 0});
        Mat3 virial{};
        for (const auto& pair : halfList(c, m_cutoff)) {
            const double r2 = dot(pair.vector, pair.vector);
            const double c6 = std::pow(m_sigma * m_sigma / r2, 3), c12 = c6 * c6;
            result.energy += 4 * m_epsilon * (c12 - c6) - m_shift;
            // Force on j along d = r_j - r_i; the reaction acts on i.
            const Vec3 onJ = scale(pair.vector, 24 * m_epsilon * (2 * c12 - c6) / r2);
            result.forces[static_cast<std::size_t>(pair.j)] = add(result.forces[static_cast<std::size_t>(pair.j)], onJ);
            result.forces[static_cast<std::size_t>(pair.i)] = sub(result.forces[static_cast<std::size_t>(pair.i)], onJ);
            addVirial(virial, scale(onJ, -1), pair.vector);
        }
        finishStress(result, c, virial, stress);
        return result;
    }

private:
    double m_epsilon, m_sigma, m_cutoff, m_shift = 0;
};

// LAMMPS pair_eam cubic interpolation of a uniformly tabulated function
// (PairEAM::interpolate), here with zero-based indices.
class Spline
{
public:
    void build(const std::vector<double>& values, double delta)
    {
        const int n = static_cast<int>(values.size());
        if (n < 5 || !(delta > 0)) throw std::runtime_error("EAM tables need at least five samples and a positive spacing");
        m_delta = delta;
        m_c.assign(static_cast<std::size_t>(n), {});
        auto c = [&](int m) -> std::array<double, 7>& { return m_c[static_cast<std::size_t>(m)]; };
        for (int m = 0; m < n; ++m) c(m)[6] = values[static_cast<std::size_t>(m)];
        c(0)[5] = c(1)[6] - c(0)[6];
        c(1)[5] = 0.5 * (c(2)[6] - c(0)[6]);
        c(n - 2)[5] = 0.5 * (c(n - 1)[6] - c(n - 3)[6]);
        c(n - 1)[5] = c(n - 1)[6] - c(n - 2)[6];
        for (int m = 2; m < n - 2; ++m)
            c(m)[5] = ((c(m - 2)[6] - c(m + 2)[6]) + 8.0 * (c(m + 1)[6] - c(m - 1)[6])) / 12.0;
        for (int m = 0; m < n - 1; ++m) {
            c(m)[4] = 3.0 * (c(m + 1)[6] - c(m)[6]) - 2.0 * c(m)[5] - c(m + 1)[5];
            c(m)[3] = c(m)[5] + c(m + 1)[5] - 2.0 * (c(m + 1)[6] - c(m)[6]);
        }
        c(n - 1)[4] = 0.0;
        c(n - 1)[3] = 0.0;
        for (int m = 0; m < n; ++m) {
            c(m)[2] = c(m)[5] / delta;
            c(m)[1] = 2.0 * c(m)[4] / delta;
            c(m)[0] = 3.0 * c(m)[3] / delta;
        }
    }

    void evaluate(double x, double& value, double& derivative) const
    {
        double p = std::max(0.0, x) / m_delta;
        int m = static_cast<int>(p);
        m = std::min(m, static_cast<int>(m_c.size()) - 2);
        p -= m;
        p = std::min(p, 1.0);
        const auto& k = m_c[static_cast<std::size_t>(m)];
        value = ((k[3] * p + k[4]) * p + k[5]) * p + k[6];
        derivative = (k[0] * p + k[1]) * p + k[2];
    }

    double maximum() const { return m_delta * static_cast<double>(m_c.size() - 1); }

private:
    double m_delta = 1.0;
    std::vector<std::array<double, 7>> m_c;
};

class Eam final : public Potential
{
public:
    Eam(const std::filesystem::path& file, std::string format)
    {
        std::ifstream input(file);
        if (!input) throw std::runtime_error("Cannot open EAM potential file " + file.u8string());
        if (format == "auto") {
            const std::string name = file.filename().u8string();
            auto ends = [&](const std::string& suffix) {
                std::string lowerName = name;
                for (char& ch : lowerName) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
                return lowerName.size() >= suffix.size() && lowerName.compare(lowerName.size() - suffix.size(), suffix.size(), suffix) == 0;
            };
            format = ends(".fs") || ends(".eam.fs") ? "fs" : (ends(".alloy") || ends(".setfl") ? "setfl" : "funcfl");
        }
        if (format != "setfl" && format != "fs" && format != "funcfl")
            throw std::runtime_error("EAM format must be setfl, fs or funcfl");
        std::string line;
        std::vector<std::string> tokens;
        auto readTokens = [&](std::istream& stream) {
            for (std::string token; stream >> token;) tokens.push_back(token);
        };
        std::size_t cursor = 0;
        auto number = [&]() {
            if (cursor >= tokens.size()) throw std::runtime_error("Truncated EAM potential file " + file.u8string());
            return std::stod(tokens[cursor++]);
        };
        auto table = [&](int count) {
            std::vector<double> values(static_cast<std::size_t>(count));
            for (double& v : values) v = number();
            return values;
        };
        if (format == "funcfl") {
            std::getline(input, line);  // comment
            readTokens(input);
            const int z = static_cast<int>(number());
            const double mass = number();
            number();          // lattice constant
            ++cursor;          // lattice name
            const int nrho = static_cast<int>(number());
            const double drho = number();
            const int nr = static_cast<int>(number());
            const double dr = number();
            m_cutoff = number();
            m_elements = {elementSymbol(z)};
            m_masses = {mass};
            m_embedding.resize(1);
            m_embedding[0].build(table(nrho), drho);
            m_rhoMax = (nrho - 1) * drho;
            const auto charge = table(nr);
            // funcfl Z(r) in e: r*phi = 27.2 * 0.529 * Z(r)^2 (Hartree-Bohr in eV-Angstrom), as LAMMPS.
            std::vector<double> rphi(charge.size());
            for (std::size_t i = 0; i < charge.size(); ++i) rphi[i] = 27.2 * 0.529 * charge[i] * charge[i];
            m_density.assign(1, std::vector<Spline>(1));
            m_density[0][0].build(table(nr), dr);
            m_pair.resize(1);
            m_pair[0].build(rphi, dr);
        } else {
            for (int k = 0; k < 3; ++k) std::getline(input, line);
            std::getline(input, line);
            std::istringstream header(line);
            int count = 0;
            header >> count;
            if (count < 1 || count > 20) throw std::runtime_error("Invalid element count in EAM file " + file.u8string());
            m_elements.resize(static_cast<std::size_t>(count));
            for (auto& element : m_elements) header >> element;
            readTokens(input);
            const int nrho = static_cast<int>(number());
            const double drho = number();
            const int nr = static_cast<int>(number());
            const double dr = number();
            m_cutoff = number();
            m_rhoMax = (nrho - 1) * drho;
            m_embedding.resize(static_cast<std::size_t>(count));
            m_density.assign(static_cast<std::size_t>(count), std::vector<Spline>(static_cast<std::size_t>(format == "fs" ? count : 1)));
            for (int e = 0; e < count; ++e) {
                number();                 // atomic number
                m_masses.push_back(number());
                number();                 // lattice constant
                ++cursor;                 // lattice name
                m_embedding[static_cast<std::size_t>(e)].build(table(nrho), drho);
                for (auto& density : m_density[static_cast<std::size_t>(e)]) density.build(table(nr), dr);
            }
            m_pair.resize(static_cast<std::size_t>(count * (count + 1) / 2));
            for (int i = 0; i < count; ++i)
                for (int j = 0; j <= i; ++j) m_pair[static_cast<std::size_t>(i * (i + 1) / 2 + j)].build(table(nr), dr);
        }
        m_finnisSinclair = format == "fs";
        if (!(m_cutoff > 0)) throw std::runtime_error("EAM cutoff must be positive");
        m_description = "EAM (" + format + ", " + file.filename().u8string() + ")";
    }

    std::string description() const override { return m_description; }

    PotentialResult compute(const Configuration& c, bool stress) const override
    {
        const std::size_t n = c.size();
        std::vector<std::size_t> type(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto found = std::find(m_elements.begin(), m_elements.end(), c.symbols[i]);
            if (found == m_elements.end()) throw std::runtime_error("The EAM potential has no parameters for " + c.symbols[i]);
            type[i] = static_cast<std::size_t>(found - m_elements.begin());
        }
        const auto pairs = halfList(c, m_cutoff);
        std::vector<double> rho(n, 0.0), embeddingSlope(n, 0.0);
        double value = 0, slope = 0;
        // Density at i contributed by a neighbour of type t: rho_t (setfl) or rho_{t,type(i)} (fs).
        auto density = [&](std::size_t from, std::size_t at) -> const Spline& {
            return m_finnisSinclair ? m_density[from][at] : m_density[from][0];
        };
        for (const auto& pair : pairs) {
            const double r = norm(pair.vector);
            const auto i = static_cast<std::size_t>(pair.i), j = static_cast<std::size_t>(pair.j);
            density(type[j], type[i]).evaluate(r, value, slope);
            rho[i] += value;
            density(type[i], type[j]).evaluate(r, value, slope);
            rho[j] += value;
        }
        PotentialResult result;
        result.forces.assign(n, {0, 0, 0});
        for (std::size_t i = 0; i < n; ++i) {
            m_embedding[type[i]].evaluate(rho[i], value, slope);
            result.energy += value;
            // Linear continuation beyond the tabulated density, as LAMMPS.
            if (rho[i] > m_rhoMax) result.energy += slope * (rho[i] - m_rhoMax);
            embeddingSlope[i] = slope;
        }
        Mat3 virial{};
        for (const auto& pair : pairs) {
            const double r = norm(pair.vector);
            const auto i = static_cast<std::size_t>(pair.i), j = static_cast<std::size_t>(pair.j);
            double rhoJI = 0, dRhoJI = 0, rhoIJ = 0, dRhoIJ = 0, rphi = 0, dRphi = 0;
            density(type[j], type[i]).evaluate(r, rhoJI, dRhoJI);
            density(type[i], type[j]).evaluate(r, rhoIJ, dRhoIJ);
            const std::size_t a = std::max(type[i], type[j]), b = std::min(type[i], type[j]);
            m_pair[a * (a + 1) / 2 + b].evaluate(r, rphi, dRphi);
            const double phi = rphi / r;
            const double dPhi = (dRphi - phi) / r;
            result.energy += phi;
            const double dEdr = embeddingSlope[i] * dRhoJI + embeddingSlope[j] * dRhoIJ + dPhi;
            const Vec3 onI = scale(pair.vector, dEdr / r);
            result.forces[i] = add(result.forces[i], onI);
            result.forces[j] = sub(result.forces[j], onI);
            addVirial(virial, onI, pair.vector);
        }
        finishStress(result, c, virial, stress);
        return result;
    }

private:
    std::vector<std::string> m_elements;
    std::vector<double> m_masses;
    std::vector<Spline> m_embedding;
    std::vector<std::vector<Spline>> m_density;
    std::vector<Spline> m_pair;
    double m_cutoff = 0.0, m_rhoMax = 0.0;
    bool m_finnisSinclair = false;
    std::string m_description;
};

double kwarg(const Json& object, const std::string& name, double fallback)
{
    const Json* value = object.find(name);
    if (!value || value->isNull()) return fallback;
    return positive(value->number(), name);
}
}

std::unique_ptr<Potential> makeEmt()
{
    return std::make_unique<Emt>();
}

std::unique_ptr<Potential> makeEam(const std::filesystem::path& file, const std::string& format)
{
    return std::make_unique<Eam>(file, format);
}

std::unique_ptr<Potential> makeLennardJones(double epsilon, double sigma, double cutoff)
{
    positive(epsilon, "epsilon");
    positive(sigma, "sigma");
    positive(cutoff, "cutoff");
    return std::make_unique<LennardJones>(epsilon, sigma, cutoff);
}

namespace
{
// Parameter tables are parsed once; every image or run shares them (potentials are stateless).
PotentialFactory sharedPotential(std::shared_ptr<const Potential> parsed)
{
    struct Shared final : Potential {
        std::shared_ptr<const Potential> inner;
        PotentialResult compute(const Configuration& c, bool stress) const override { return inner->compute(c, stress); }
        std::string description() const override { return inner->description(); }
    };
    return [parsed] { auto copy = std::make_unique<Shared>(); copy->inner = parsed; return std::unique_ptr<Potential>(std::move(copy)); };
}
}

PotentialFactory potentialFactory(const Json& specification, const std::filesystem::path& base)
{
    if (!specification.isObject()) throw std::runtime_error("The interatomic potential must be an object such as {\"potential\": \"EMT\"}");
    std::string name;
    Json options = Json::object();
    if (const Json* potential = specification.find("potential")) {
        name = potential->string();
        options = specification;
    } else if (specification.contains("module") && specification.contains("attribute")) {
        // Accept the earlier ASE-style description for the native equivalents.
        name = specification.at("attribute").string();
        if (const Json* kwargs = specification.find("kwargs"); kwargs && !kwargs->isNull()) options = *kwargs;
        if (options.contains("rc")) options["cutoff"] = Json(options.at("rc"));
    } else throw std::runtime_error("The interatomic potential must name \"potential\": EMT or LennardJones");
    std::string key;
    for (char c : name) if (std::isalnum(static_cast<unsigned char>(c))) key += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (key == "emt") return [] { return makeEmt(); };
    if (key == "lennardjones" || key == "lj") {
        const double sigma = kwarg(options, "sigma", 1.0);
        const double epsilon = kwarg(options, "epsilon", 1.0);
        const double cutoff = kwarg(options, "cutoff", 3.0 * sigma);
        makeLennardJones(epsilon, sigma, cutoff);
        return [=] { return makeLennardJones(epsilon, sigma, cutoff); };
    }
    if (key == "eam" || key == "eamalloy" || key == "eamfs" || key == "finnissinclair") {
        const Json* file = options.find("file");
        if (!file || !file->isString() || file->string().empty()) throw std::runtime_error("The EAM potential needs a \"file\"");
        const auto path = (base / std::filesystem::u8path(file->string())).lexically_normal();
        std::string format = options.contains("format") ? options.at("format").string() : (key == "eamfs" || key == "finnissinclair" ? "fs" : "auto");
        // Parse once up front so a bad file fails before a calculation starts;
        // each image or run then shares the parsed tables.
        std::shared_ptr<const Potential> parsed(makeEam(path, format).release());
        return sharedPotential(parsed);
    }
    if (key == "tersoff" || key == "stillingerweber" || key == "sw") {
        std::filesystem::path path;
        if (const Json* file = options.find("file"); file && file->isString() && !file->string().empty())
            path = (base / std::filesystem::u8path(file->string())).lexically_normal();
        const bool tersoff = key == "tersoff";
        std::shared_ptr<const Potential> parsed((tersoff ? makeTersoff(path) : makeStillingerWeber(path)).release());
        return sharedPotential(parsed);
    }
    if (key == "buckingham" || key == "buck") {
        std::shared_ptr<const Potential> parsed(makeBuckingham(options).release());
        return sharedPotential(parsed);
    }
    throw std::runtime_error("Unsupported interatomic potential '" + name + "'; AtomForge provides EMT, LennardJones, EAM, Tersoff, StillingerWeber and Buckingham natively");
}
}
