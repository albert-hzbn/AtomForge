#pragma once

#include <string>
#include <vector>
#include <map>
#include "model/Structure.h"

/**
 * @struct SroEntry
 * Warren-Cowley SRO parameter entry: alpha_ij^(s) = 1 - P(j|i,s) / c_j
 * Reference: J.M. Cowley, Phys. Rev. 138 (1965) A1384
 */
struct SroEntry {
    int shell = 0;
    double shellDistance = 0.0;
    std::string centerElement;
    std::string neighborElement;
    int nij = 0;           // Count of j-neighbors around i
    int zi = 0;            // Total coordination around i
    double pij = 0.0;      // P(j|i)
    double cj = 0.0;       // Global concentration of j
    double alpha = 0.0;    // Warren-Cowley parameter
};

/**
 * @struct SroEntryRC
 * Rao-Curtin SRO parameter entry: alpha_ij^(s) = 1 - P_ij^(s) / (2*c_i*c_j)
 * Reference: Rao & Curtin, Acta Materialia 226 (2022) 117621
 */
struct SroEntryRC {
    int shell = 0;
    double shellDistance = 0.0;
    std::string elem_i;
    std::string elem_j;
    int n_ij = 0;          // Count of j around i
    int n_ji = 0;          // Count of i around j
    double z_shell = 0.0;  // Average coordination in this shell
    double p_ij = 0.0;     // Symmetric joint probability
    double ci = 0.0, cj = 0.0;  // Compositions
    double alpha = 0.0;    // Rao-Curtin parameter
};

/**
 * @struct SroReport
 * Warren-Cowley SRO analysis report.
 */
struct SroReport {
    int nAtoms = 0;
    int nShells = 0;
    std::vector<double> shellDistances;
    std::map<std::string, double> composition;
    std::vector<SroEntry> entries;
    std::string message;
};

/**
 * @struct SroReportRC
 * Rao-Curtin SRO analysis report.
 */
struct SroReportRC {
    int nAtoms = 0;
    int nShells = 0;
    std::vector<double> shellDistances;
    std::map<std::string, double> composition;
    std::vector<double> shellCoordination;
    std::vector<SroEntryRC> entries;
    std::string message;
};

/**
 * @struct SroAnalysisResult
 * Complete SRO analysis result containing both methods.
 */
struct SroAnalysisResult {
    bool success = false;
    std::string message;
    SroReport warrenCowley;
    SroReportRC raoCurtin;
};

/**
 * Compute composition map from structure.
 */
std::map<std::string, double> computeSroComposition(const Structure& structure);

/**
 * Build unique neighbor shell distances (clustered by tolerance).
 */
std::vector<double> buildSroShellDistances(const Structure& structure,
                                            int maxShells,
                                            double toleranceAng);

/**
 * Perform Warren-Cowley SRO analysis.
 */
SroReport analyzeWarrenCowleySRO(const Structure& structure,
                                  int nShells,
                                  double shellTolerance);

/**
 * Perform Rao-Curtin (2022) SRO analysis.
 */
SroReportRC analyzeRaoCurtinSRO(const Structure& structure,
                                 int nShells,
                                 double shellTolerance);

/**
 * Complete SRO analysis (both methods).
 */
SroAnalysisResult analyzeSRO(const Structure& structure,
                              int nShells,
                              double shellTolerance);
