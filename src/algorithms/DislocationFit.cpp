#include "algorithms/DislocationFit.h"

namespace atomforge
{

DislocationFitResult fitDislocationFromNye(const std::vector<glm::dvec3>& positions,
                                          const NyeTensorResult& nye,
                                          const glm::dvec3& lineDirectionIn,
                                          double areaPerAtom)
{
    DislocationFitResult result;
    if (positions.size() != nye.alpha.size() || positions.empty())
    {
        result.message = "positions and the Nye tensor result must be the same non-empty size.";
        return result;
    }

    const glm::dvec3 lineDir = glm::normalize(lineDirectionIn);

    double totalWeight = 0.0;
    glm::dvec3 weightedCentroid(0.0);
    glm::dvec3 burgers(0.0);

    for (size_t i = 0; i < positions.size(); ++i)
    {
        const double w = nye.alphaNorm[i];
        totalWeight += w;
        weightedCentroid += w * positions[i];

        // Nye's theorem: b = sum_atoms alpha^T . n * areaPerAtom, for the
        // cross-sectional plane with normal n = lineDirection. Transposed
        // relative to the textbook b_i = alpha_ij n_j because NyeTensor.cpp
        // stores alpha[l][k] indexed by (G-tensor column, row) rather than
        // (row, column); verified empirically against a known edge
        // dislocation's Burgers direction (dislocation_fit_regressions.cpp).
        const auto& a = nye.alpha[i];
        for (int row = 0; row < 3; ++row)
        {
            double contracted = 0.0;
            for (int col = 0; col < 3; ++col)
                contracted += a[col][row] * lineDir[col];
            burgers[row] += contracted * areaPerAtom;
        }
    }

    if (totalWeight <= 1e-12)
    {
        result.message = "Nye tensor is zero everywhere; no dislocation signal to fit.";
        return result;
    }

    result.linePosition = weightedCentroid / totalWeight;
    result.burgersVector = burgers;
    result.totalWeight = totalWeight;
    result.success = true;
    result.message = "Fitted from " + std::to_string(positions.size()) + " atoms, total |alpha| weight = "
        + std::to_string(totalWeight) + ".";
    return result;
}

}
