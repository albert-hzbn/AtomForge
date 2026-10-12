#include "ShortRangeOrderAnalysis.h"
#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static glm::dvec3 sroMinimalImageFrac(glm::dvec3 df) {
    for (int i = 0; i < 3; ++i) {
        if (df[i] >  0.5) df[i] -= 1.0;
        if (df[i] < -0.5) df[i] += 1.0;
    }
    return df;
}

static double sroPeriodicDist(const glm::dmat3& cell,
                              const glm::dvec3& fi,
                              const glm::dvec3& fj) {
    return glm::length(cell * sroMinimalImageFrac(fj - fi));
}

// Build column-major cell matrix from row-major cellVectors
static glm::dmat3 buildCell(const Structure& s) {
    glm::dmat3 cell(0.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            cell[j][i] = s.cellVectors[i][j];
    return cell;
}

// Atom positions are Cartesian; the minimum-image distances need them in
// fractional coordinates of the cell.
static std::vector<glm::dvec3> fractionalPositions(const Structure& s, const glm::dmat3& cell) {
    const glm::dmat3 inverse = glm::inverse(cell);
    std::vector<glm::dvec3> frac;
    frac.reserve(s.atoms.size());
    for (const auto& a : s.atoms)
        frac.push_back(inverse * glm::dvec3(a.x, a.y, a.z));
    return frac;
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
    if (structure.atoms.empty() || !structure.hasUnitCell) return {};

    glm::dmat3 cell = buildCell(structure);
    const std::vector<glm::dvec3> frac = fractionalPositions(structure, cell);
    int N = static_cast<int>(structure.atoms.size());

    std::vector<double> dists;
    dists.reserve(N * N);

    for (int i = 0; i < N; ++i) {
        const glm::dvec3& fi = frac[static_cast<std::size_t>(i)];
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            const glm::dvec3& fj = frac[static_cast<std::size_t>(j)];
            double d = sroPeriodicDist(cell, fi, fj);
            if (d > 1e-6) dists.push_back(d);
        }
    }

    if (dists.empty()) return {};
    std::sort(dists.begin(), dists.end());

    std::vector<double> shells;
    std::vector<bool> used(dists.size(), false);

    for (size_t i = 0; i < dists.size() && (int)shells.size() < maxShells; ++i) {
        if (used[i]) continue;
        double sum = dists[i];
        int cnt = 1;
        used[i] = true;
        for (size_t j = i + 1; j < dists.size(); ++j) {
            if (!used[j] && std::abs(dists[j] - dists[i]) < toleranceAng) {
                sum += dists[j];
                ++cnt;
                used[j] = true;
            }
        }
        shells.push_back(sum / cnt);
    }
    return shells;
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

    glm::dmat3 cell = buildCell(structure);
    const std::vector<glm::dvec3> frac = fractionalPositions(structure, cell);
    int N = report.nAtoms;

    // Count directional pairs
    std::map<std::tuple<int, std::string, std::string>, int> pairs;
    std::map<std::pair<int, std::string>, int> coordZ;

    for (int i = 0; i < N; ++i) {
        const glm::dvec3& fi = frac[static_cast<std::size_t>(i)];
        const std::string& si = structure.atoms[i].symbol;

        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            const glm::dvec3& fj = frac[static_cast<std::size_t>(j)];
            const std::string& sj = structure.atoms[j].symbol;
            double d = sroPeriodicDist(cell, fi, fj);
            if (d < 1e-6) continue;

            for (int s = 0; s < report.nShells; ++s) {
                if (std::abs(d - report.shellDistances[s]) < shellTolerance) {
                    pairs[{s, si, sj}]++;
                    coordZ[{s, si}]++;
                    break;
                }
            }
        }
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

    glm::dmat3 cell = buildCell(structure);
    const std::vector<glm::dvec3> frac = fractionalPositions(structure, cell);
    int N = report.nAtoms;

    std::map<std::tuple<int, std::string, std::string>, int> pairs;
    std::vector<int> totalPerShell(report.nShells, 0);

    for (int i = 0; i < N; ++i) {
        const glm::dvec3& fi = frac[static_cast<std::size_t>(i)];
        const std::string& si = structure.atoms[i].symbol;

        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            const glm::dvec3& fj = frac[static_cast<std::size_t>(j)];
            const std::string& sj = structure.atoms[j].symbol;
            double d = sroPeriodicDist(cell, fi, fj);
            if (d < 1e-6) continue;

            for (int s = 0; s < report.nShells; ++s) {
                if (std::abs(d - report.shellDistances[s]) < shellTolerance) {
                    pairs[{s, si, sj}]++;
                    totalPerShell[s]++;
                    break;
                }
            }
        }
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
