#include "ShortRangeOrderAnalysis.h"
#include <algorithm>
#include <array>
#include <map>
#include <cmath>
#include <glm/glm.hpp>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Every pair within `cutoff` over all periodic images (a small cell counts
// each image of a neighbour, as a large one would), with its distance.
struct SroDistance {
    int i;
    int j;
    double d;
};

static std::vector<SroDistance> sroNeighbours(const Structure& s, double cutoff) {
    std::vector<SroDistance> pairs;
    const std::size_t n = s.atoms.size();
    if (n == 0 || cutoff <= 0.0) return pairs;
    glm::dmat3 cell(0.0);  // columns: a, b, c
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) cell[r][c] = s.cellVectors[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
    const glm::dmat3 inverse = glm::inverse(cell);
    // Images needed along each axis: the planes of constant fractional
    // coordinate k lie 1/|row k of the inverse| apart.
    int range[3];
    for (int k = 0; k < 3; ++k) {
        const glm::dvec3 row(inverse[0][k], inverse[1][k], inverse[2][k]);
        range[k] = static_cast<int>(std::ceil(cutoff * glm::length(row)));
        if (range[k] > 50) throw std::runtime_error("SRO cutoff spans too many periodic images");
    }
    // Atoms wrapped into the cell, binned on a grid of the cutoff.
    std::vector<glm::dvec3> wrapped(n);
    for (std::size_t i = 0; i < n; ++i) {
        glm::dvec3 f = inverse * glm::dvec3(s.atoms[i].x, s.atoms[i].y, s.atoms[i].z);
        f -= glm::floor(f);
        wrapped[i] = cell * f;
    }
    const auto binOf = [cutoff](const glm::dvec3& p) {
        return std::array<long long, 3>{static_cast<long long>(std::floor(p.x / cutoff)),
                                        static_cast<long long>(std::floor(p.y / cutoff)),
                                        static_cast<long long>(std::floor(p.z / cutoff))};
    };
    std::map<std::array<long long, 3>, std::vector<int>> bins;
    for (std::size_t j = 0; j < n; ++j) bins[binOf(wrapped[j])].push_back(static_cast<int>(j));
    const double cutoff2 = cutoff * cutoff;
    for (std::size_t i = 0; i < n; ++i)
        for (int a = -range[0]; a <= range[0]; ++a)
            for (int b = -range[1]; b <= range[1]; ++b)
                for (int c = -range[2]; c <= range[2]; ++c) {
                    const glm::dvec3 shift = cell * glm::dvec3(a, b, c);
                    // Atom j's image at wrapped[j] + shift is near i when wrapped[j] is near wrapped[i] - shift.
                    const auto centre = binOf(wrapped[i] - shift);
                    for (long long dx = -1; dx <= 1; ++dx)
                        for (long long dy = -1; dy <= 1; ++dy)
                            for (long long dz = -1; dz <= 1; ++dz) {
                                const auto found = bins.find({centre[0] + dx, centre[1] + dy, centre[2] + dz});
                                if (found == bins.end()) continue;
                                for (int j : found->second) {
                                    if (j == static_cast<int>(i) && a == 0 && b == 0 && c == 0) continue;
                                    const glm::dvec3 d = wrapped[static_cast<std::size_t>(j)] + shift - wrapped[i];
                                    const double d2 = glm::dot(d, d);
                                    if (d2 < cutoff2 && d2 > 1e-12) pairs.push_back({static_cast<int>(i), j, std::sqrt(d2)});
                                }
                            }
                }
    return pairs;
}

// Clusters sorted distances into shells: each shell averages the distances
// within `tolerance` of its first one.
static std::vector<double> sroClusterShells(std::vector<double> dists, double tolerance, std::size_t maxShells) {
    std::sort(dists.begin(), dists.end());
    std::vector<double> shells;
    std::size_t i = 0;
    while (i < dists.size() && shells.size() < maxShells) {
        double sum = 0.0;
        std::size_t j = i;
        for (; j < dists.size() && dists[j] - dists[i] < tolerance; ++j) sum += dists[j];
        shells.push_back(sum / static_cast<double>(j - i));
        i = j;
    }
    return shells;
}

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

std::map<std::string, double> computeSroComposition(const Structure& structure) {
    std::map<std::string, double> comp;
    int total = static_cast<int>(structure.atoms.size());
    if (total == 0) return comp;
    std::map<std::string, int> counts;
    for (const auto& a : structure.atoms) counts[a.symbol]++;
    for (const auto& [sym, cnt] : counts)
        comp[sym] = static_cast<double>(cnt) / total;
    return comp;
}

std::vector<double> buildSroShellDistances(const Structure& structure,
                                            int maxShells,
                                            double toleranceAng) {
    if (structure.atoms.empty() || !structure.hasUnitCell || maxShells <= 0) return {};
    const glm::dvec3 a(structure.cellVectors[0][0], structure.cellVectors[0][1], structure.cellVectors[0][2]);
    const glm::dvec3 b(structure.cellVectors[1][0], structure.cellVectors[1][1], structure.cellVectors[1][2]);
    const glm::dvec3 c(structure.cellVectors[2][0], structure.cellVectors[2][1], structure.cellVectors[2][2]);
    const double volume = std::abs(glm::dot(a, glm::cross(b, c)));
    if (volume < 1e-9) return {};
    // Start near the second-neighbour distance of a dense packing and widen
    // the search until the requested shells lie wholly inside it.
    double cutoff = 1.6 * std::cbrt(volume / static_cast<double>(structure.atoms.size()));
    std::vector<double> complete;
    for (int attempt = 0; attempt < 8; ++attempt) {
        std::vector<double> dists;
        try {
            for (const auto& pair : sroNeighbours(structure, cutoff)) dists.push_back(pair.d);
        } catch (const std::runtime_error&) {
            break;  // the cutoff spans too many periodic images of a tiny cell
        }
        complete.clear();
        for (double shell : sroClusterShells(dists, toleranceAng, static_cast<std::size_t>(maxShells) + 1))
            if (shell + toleranceAng < cutoff) complete.push_back(shell);
        if (static_cast<int>(complete.size()) >= maxShells) break;
        cutoff *= 1.5;
    }
    if (static_cast<int>(complete.size()) > maxShells) complete.resize(static_cast<std::size_t>(maxShells));
    return complete;
}

std::vector<SroPair> buildSroPairs(const Structure& structure,
                                   const std::vector<double>& shellDistances,
                                   double toleranceAng) {
    std::vector<SroPair> pairs;
    if (structure.atoms.empty() || !structure.hasUnitCell || shellDistances.empty()) return pairs;
    const double cutoff = *std::max_element(shellDistances.begin(), shellDistances.end()) + toleranceAng;
    for (const SroDistance& pair : sroNeighbours(structure, cutoff)) {
        for (int sh = 0; sh < static_cast<int>(shellDistances.size()); ++sh) {
            if (std::abs(pair.d - shellDistances[static_cast<std::size_t>(sh)]) < toleranceAng) {
                pairs.push_back({pair.i, pair.j, sh});
                break;
            }
        }
    }
    return pairs;
}

SroReport analyzeWarrenCowleySRO(const Structure& structure,
                                  int nShells,
                                  double shellTolerance) {
    SroReport report;
    report.nAtoms = static_cast<int>(structure.atoms.size());

    if (structure.atoms.empty() || !structure.hasUnitCell) {
        report.message = "No structure or unit cell defined.";
        return report;
    }

    report.composition = computeSroComposition(structure);
    report.shellDistances = buildSroShellDistances(structure, nShells, shellTolerance);
    report.nShells = static_cast<int>(report.shellDistances.size());

    if (report.nShells == 0) {
        report.message = "Could not identify neighbor shells.";
        return report;
    }

    // Count directional pairs
    std::map<std::tuple<int, std::string, std::string>, int> pairs;
    std::map<std::pair<int, std::string>, int> coordZ;
    for (const SroPair& p : buildSroPairs(structure, report.shellDistances, shellTolerance)) {
        const std::string& si = structure.atoms[static_cast<std::size_t>(p.i)].symbol;
        const std::string& sj = structure.atoms[static_cast<std::size_t>(p.j)].symbol;
        pairs[{p.shell, si, sj}]++;
        coordZ[{p.shell, si}]++;
    }

    for (int s = 0; s < report.nShells; ++s) {
        for (const auto& [ei, ci] : report.composition) {
            int zi = coordZ[{s, ei}];
            if (zi == 0) continue;
            for (const auto& [ej, cj] : report.composition) {
                int nij = pairs[{s, ei, ej}];
                double pij = static_cast<double>(nij) / zi;
                SroEntry e;
                e.shell = s;
                e.shellDistance = report.shellDistances[s];
                e.centerElement = ei;
                e.neighborElement = ej;
                e.nij = nij;
                e.zi = zi;
                e.pij = pij;
                e.cj = cj;
                e.alpha = (cj > 1e-12) ? (1.0 - pij / cj) : 0.0;
                report.entries.push_back(e);
            }
        }
    }

    report.message = "Warren-Cowley SRO completed.";
    return report;
}

SroReportRC analyzeRaoCurtinSRO(const Structure& structure,
                                 int nShells,
                                 double shellTolerance) {
    SroReportRC report;
    report.nAtoms = static_cast<int>(structure.atoms.size());

    if (structure.atoms.empty() || !structure.hasUnitCell) {
        report.message = "No structure or unit cell defined.";
        return report;
    }

    report.composition = computeSroComposition(structure);
    report.shellDistances = buildSroShellDistances(structure, nShells, shellTolerance);
    report.nShells = static_cast<int>(report.shellDistances.size());

    if (report.nShells == 0) {
        report.message = "Could not identify neighbor shells.";
        return report;
    }

    const int N = report.nAtoms;
    std::map<std::tuple<int, std::string, std::string>, int> pairs;
    std::vector<int> totalPerShell(report.nShells, 0);
    for (const SroPair& p : buildSroPairs(structure, report.shellDistances, shellTolerance)) {
        pairs[{p.shell, structure.atoms[static_cast<std::size_t>(p.i)].symbol, structure.atoms[static_cast<std::size_t>(p.j)].symbol}]++;
        totalPerShell[static_cast<std::size_t>(p.shell)]++;
    }

    report.shellCoordination.resize(report.nShells);
    for (int s = 0; s < report.nShells; ++s)
        report.shellCoordination[s] = (N > 0)
            ? static_cast<double>(totalPerShell[s]) / N : 0.0;

    // Unique pairs (i <= j lexicographically), ordered shell-first so the UI
    // can group by shell with a single pass through entries.
    for (int s = 0; s < report.nShells; ++s) {
        double zs = report.shellCoordination[s];
        for (const auto& [ei, ci] : report.composition) {
            for (const auto& [ej, cj] : report.composition) {
                if (ei > ej) continue;
                if (zs < 1e-12) continue;
                int nij = pairs[{s, ei, ej}];
                int nji = (ei != ej) ? pairs[{s, ej, ei}] : nij;
                double pij = static_cast<double>(nij + nji) / (N * zs);
                double denom = 2.0 * ci * cj;
                SroEntryRC e;
                e.shell = s;
                e.shellDistance = report.shellDistances[s];
                e.elem_i = ei;
                e.elem_j = ej;
                e.n_ij = nij;
                e.n_ji = nji;
                e.z_shell = zs;
                e.p_ij = pij;
                e.ci = ci;
                e.cj = cj;
                e.alpha = (denom > 1e-12) ? (1.0 - pij / denom) : 0.0;
                report.entries.push_back(e);
            }
        }
    }

    report.message = "Rao-Curtin SRO (2022) completed.";
    return report;
}

SroAnalysisResult analyzeSRO(const Structure& structure,
                              int nShells,
                              double shellTolerance) {
    SroAnalysisResult result;

    if (structure.atoms.empty() || !structure.hasUnitCell) {
        result.success = false;
        result.message = "Structure is empty or has no unit cell.";
        return result;
    }

    result.warrenCowley = analyzeWarrenCowleySRO(structure, nShells, shellTolerance);
    result.raoCurtin    = analyzeRaoCurtinSRO(structure, nShells, shellTolerance);

    result.success = result.warrenCowley.nShells > 0 && result.raoCurtin.nShells > 0;
    result.message = result.success
        ? "SRO analysis complete (Warren-Cowley + Rao-Curtin)."
        : "SRO analysis failed.";

    return result;
}
