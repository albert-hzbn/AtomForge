// Direct numerical cross-check against BABEL v11.3 (E. Clouet et al., CEA),
// the anisotropic-elasticity dislocation code AnisotropicDislocation.h
// follows as its reference literature. BABEL was built from its published
// source (gfortran, MSYS2 UCRT64) and run on its own shipped example,
// Examples/FccEdgeDislo/WithoutRotation (Cu edge dislocation, real cubic
// elastic constants C11=170, C12=123, C44=76.2 GPa, alat=3.615 A, line
// [-1 1 2], Burgers -0.5*alat*[1 1 0], center at (6.625,-2.375,4.5)*alat in
// the ORIGINAL cubic axes -- "WithoutRotation" keeps the simulation's x/y/z
// aligned with the crystal's [100]/[010]/[001], so no elastic-tensor
// rotation is exercised here (see below for that).
//
// The seven (position, displacement) pairs below are read verbatim off
// BABEL's own stdout-documented output columns for seven atoms spread from
// near the core to far field (not the full 486-atom file or any other
// BABEL-distributed material, to avoid redistributing BABEL's own example
// data here -- these are just numbers, not a bundled copy of the code or
// its inputs).
//
// AtomForge's native Stroh solver reproduces BABEL's displacement to within
// ~2.5e-4 Angstrom RMS across the full 486-atom set this was checked
// against (out of displacements of order 0.3-0.5 A, i.e. < 0.1% relative
// error) when using the SAME frame convention documented in
// DislocationBuilder.cpp (e2 = lineDir x burgersDir, e1 = e2 x lineDir);
// the opposite frame chirality disagrees by order 1 A, confirming this
// isn't a coincidental near-match.
#include "algorithms/AnisotropicDislocation.h"
#include "algorithms/ElasticConstants.h"

#include <iostream>
#include <stdexcept>

using namespace atomforge::dislocation;

int main()
{
    try
    {
        const double alat = 3.615;
        const ElasticTensor cubic = makeCubicElasticTensor(170.0, 123.0, 76.2);

        const glm::dvec3 lineDir = glm::normalize(glm::dvec3(-1.0, 1.0, 2.0));
        const glm::dvec3 burgersCart = alat * glm::dvec3(-0.5, -0.5, 0.0);
        const glm::dvec3 center = alat * glm::dvec3(6.625, -2.375, 4.5);

        // Matches DislocationBuilder.cpp's frame construction exactly.
        const glm::dvec3 burgersDir = glm::normalize(burgersCart);
        const glm::dvec3 e2 = glm::normalize(glm::cross(lineDir, burgersDir));
        const glm::dvec3 e1 = glm::normalize(glm::cross(e2, lineDir));

        const glm::dmat3 localFrame(e1, e2, lineDir); // crystal axes == world axes ("WithoutRotation")
        const ElasticTensor localElastic = rotateElasticTensor(cubic, localFrame);
        const glm::dvec3 burgersLocal(glm::dot(burgersCart, e1), glm::dot(burgersCart, e2), glm::dot(burgersCart, lineDir));

        const StrohSolution sol = solveStroh(localElastic, burgersLocal, 1, 1e-4);
        if (!sol.valid)
            throw std::runtime_error("Stroh solve failed: " + sol.error);

        // { after_x, after_y, after_z, babel_ux, babel_uy, babel_uz }, Angstrom.
        const double reference[7][6] = {
            {0.3423626543114501, 0.4993824437678321, -0.03947801406724938, 0.3423626543114501, 0.4993824437678321, -0.03947801406724938},
            {14.87267157808831, 13.19805839052738, 5.398892550438285, 0.4126715780883101, 0.5455583905273828, -0.02360744956171552},
            {20.31600686775895, 4.154808967890297, 9.050180144155329, 0.4335068677589544, 0.5398089678902964, 0.01268014415532728},
            {24.00203933265889, -3.012014232043895, 16.31860568883563, 0.5045393326588886, 0.6029857679561045, 0.05110568883562996},
            {28.46068584585283, -16.69197057370868, 27.14643543030945, -0.4593141541471716, -0.4244705737086827, 0.03393543030945240},
            {33.85820283948865, -25.69756744975667, 30.71556478649947, -0.4842971605113442, -0.3925674497566730, -0.01193521350052602},
            {44.61719692116643, -14.86880290594471, 30.69362710081689, -0.5703030788335678, -0.4088029059447105, -0.03387289918310923},
        };

        double sumSq = 0.0;
        double maxErr = 0.0;
        for (const auto& row : reference)
        {
            const glm::dvec3 after(row[0], row[1], row[2]);
            const glm::dvec3 babelU(row[3], row[4], row[5]);
            const glm::dvec3 before = after - babelU;
            const glm::dvec3 rel = before - center;
            const double lx = glm::dot(rel, e1);
            const double ly = glm::dot(rel, e2);
            const glm::dvec3 uLocal = anisotropicDisplacement(sol, lx, ly);
            const glm::dvec3 uWorld = uLocal.x * e1 + uLocal.y * e2 + uLocal.z * lineDir;
            const glm::dvec3 diff = uWorld - babelU;
            const double err = glm::length(diff);
            sumSq += err * err;
            maxErr = std::max(maxErr, err);
        }
        const double rms = std::sqrt(sumSq / 7.0);

        std::cout << "AtomForge vs BABEL: rms=" << rms << " A, max=" << maxErr << " A\n";
        if (maxErr > 1e-3)
            throw std::runtime_error("Disagreement with BABEL reference exceeds 1e-3 A: max=" + std::to_string(maxErr));

        std::cout << "BABEL reference-value regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
