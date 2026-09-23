// Validates fitDislocationFromNye (src/algorithms/DislocationFit.h), which
// serves the same purpose as BABEL's (self-described as "not documented
// yet") displacementFit/vitekFit programs: recovering a dislocation's
// position and Burgers vector from a measured displacement/Nye-tensor
// field, using the Nye tensor's own defining properties (weighted centroid
// = line position, cross-sectional integral = Burgers vector) rather than
// BABEL's iterative nonlinear least-squares fit.
#include "algorithms/DislocationBuilder.h"
#include "algorithms/DislocationFit.h"
#include "algorithms/NyeTensor.h"
#include "model/Structure.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge;

namespace
{
Structure buildFcc(double a, int n)
{
    const double frac[4][3] = {{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
    Structure s;
    s.hasUnitCell = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            s.cellVectors[i][j] = (i == j) ? a * n : 0.0;
    for (int ix = 0; ix < n; ++ix)
    for (int iy = 0; iy < n; ++iy)
    for (int iz = 0; iz < n; ++iz)
    for (int k = 0; k < 4; ++k)
    {
        AtomSite atom;
        atom.symbol = "Cu"; atom.atomicNumber = 29;
        atom.x = (ix + frac[k][0]) * a;
        atom.y = (iy + frac[k][1]) * a;
        atom.z = (iz + frac[k][2]) * a;
        s.atoms.push_back(atom);
    }
    return s;
}
}

int main()
{
    try
    {
        const double a = 3.61;
        Structure reference = buildFcc(a, 13);

        DislocationParams params;
        params.character = DislocationCharacter::Edge;
        params.shape = DislocationShape::Cylinder;
        params.cylinderRadius = 15.0f;
        params.coreRadius = 1.2f;
        params.cutoffRadius = 0.0f;

        const DislocationResult dislo = buildDislocation(reference, params);
        if (!dislo.success)
            throw std::runtime_error("Failed to build reference edge dislocation: " + dislo.message);

        const NyeTensorResult nye = computeNyeTensor(reference, dislo.output, 3.0f, false);
        if (!nye.success)
            throw std::runtime_error("Nye tensor computation failed: " + nye.message);

        std::vector<glm::dvec3> positions(reference.atoms.size());
        for (size_t i = 0; i < positions.size(); ++i)
            positions[i] = glm::dvec3(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);

        const DislocationFitResult fit = fitDislocationFromNye(positions, nye, glm::dvec3(dislo.lineDirection), 1.0);
        if (!fit.success)
            throw std::runtime_error("Fit failed: " + fit.message);

        // Position: perpendicular-plane distance from the fitted point to
        // the TRUE line must be small (well within the core/cutoff region
        // that actually carries Nye-tensor weight).
        const glm::dvec3 lineDir(dislo.lineDirection);
        const glm::dvec3 linePoint(dislo.linePoint);
        const glm::dvec3 rel = fit.linePosition - linePoint;
        const double perpDist = glm::length(rel - lineDir * glm::dot(rel, lineDir));
        std::cout << "Fitted line position perpendicular distance from truth: " << perpDist << " A\n";
        if (!(perpDist < 8.0))
            throw std::runtime_error("Fitted dislocation position is too far from the true line: " + std::to_string(perpDist) + " A");

        // Burgers vector direction: must align closely with the true
        // Burgers direction (magnitude depends on the areaPerAtom
        // calibration, which this test doesn't attempt to match exactly).
        const glm::dvec3 trueBurgersDir = glm::normalize(dislo.burgersDirection);
        const double fitBurgersMag = glm::length(fit.burgersVector);
        std::cout << "Fitted Burgers vector magnitude (uncalibrated): " << fitBurgersMag << "\n";
        if (!(fitBurgersMag > 1e-6))
            throw std::runtime_error("Fitted Burgers vector is essentially zero");
        const glm::dvec3 fitBurgersDir = fit.burgersVector / fitBurgersMag;
        const double alignment = std::abs(glm::dot(fitBurgersDir, trueBurgersDir));
        std::cout << "Burgers direction alignment |cos theta| = " << alignment << "\n";
        if (!(alignment > 0.9))
            throw std::runtime_error("Fitted Burgers vector direction does not align with the true Burgers direction: |cos theta| = " + std::to_string(alignment));

        std::cout << "Dislocation fit regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
