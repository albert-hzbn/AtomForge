// Validates the anisotropic (Stroh sextic) dislocation displacement solver
// against AtomForge's independently-implemented isotropic closed form
// (Hirth & Lothe eq. 3-45/3-51) in the near-isotropic limit, and against a
// direct Burgers-vector circuit integral for a genuinely anisotropic
// material. No BABEL binary was available in this environment to compare
// against directly, so these are the primary correctness gates for the
// Stroh implementation; see AnisotropicDislocation.h for the literature
// this follows.
//
// Displacement fields of this kind (line integrals of a log singularity)
// are only defined up to an arbitrary additive rigid translation -- adding
// a constant changes no strain/stress anywhere, so two independently
// derived formulas for the "same" dislocation can (and, empirically, do
// here) differ by a fixed offset while describing the identical physical
// state. Comparisons below are therefore always taken relative to a shared
// reference point rather than as raw displacement values.
#include "algorithms/AnisotropicDislocation.h"
#include "algorithms/ElasticConstants.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge::dislocation;

namespace
{
constexpr double kPi = 3.14159265358979323846;

// Independent reimplementation of the isotropic closed form (matches the
// formulas in src/algorithms/DislocationBuilder.cpp's computeDisplacement,
// written here from the textbook rather than shared, so a bug in one is
// unlikely to be replicated in the other).
glm::dvec3 isotropicDisplacement(double burgersEdge, double burgersScrew, double nu,
                                 double x, double y)
{
    const double r2 = x * x + y * y;
    double ux = 0.0, uy = 0.0, uz = 0.0;
    if (burgersEdge != 0.0)
    {
        const double pref = burgersEdge / (2.0 * kPi);
        ux = pref * (std::atan2(y, x) + (x * y) / (2.0 * (1.0 - nu) * r2));
        const double q = 4.0 * (1.0 - nu);
        uy = -pref * (((1.0 - 2.0 * nu) / q) * std::log(r2) + (x * x - y * y) / (q * r2));
    }
    if (burgersScrew != 0.0)
        uz = (burgersScrew / (2.0 * kPi)) * std::atan2(y, x);
    return glm::dvec3(ux, uy, uz);
}

// Burgers-vector circuit integral: sum displacement steps around a loop,
// discarding the one large step at each term's branch cut (a real jump of
// order |b|, not a smooth field variation).
glm::dvec3 circuitIntegral(const StrohSolution& sol, double radius, int steps)
{
    glm::dvec3 prev = anisotropicDisplacement(sol, radius, 0.0);
    glm::dvec3 accumulated(0.0);
    for (int i = 1; i <= steps; ++i)
    {
        const double theta = 2.0 * kPi * (double)i / (double)steps;
        const glm::dvec3 cur = anisotropicDisplacement(sol, radius * std::cos(theta), radius * std::sin(theta));
        const glm::dvec3 delta = cur - prev;
        if (glm::length(delta) < 1.0)
            accumulated += delta;
        prev = cur;
    }
    return accumulated;
}

// Compares the SHAPE of the two displacement fields -- both evaluated
// relative to a shared reference point, which cancels the arbitrary
// additive-constant gauge freedom described above.
void expectShapeClose(const char* label, const StrohSolution& sol,
                      double burgersEdge, double burgersScrew, double nu,
                      const double (&points)[][2], int count, double tol)
{
    const double refX = 10.0, refY = 0.0;
    const glm::dvec3 anisoRef = anisotropicDisplacement(sol, refX, refY);
    const glm::dvec3 isoRef = isotropicDisplacement(burgersEdge, burgersScrew, nu, refX, refY);

    for (int i = 0; i < count; ++i)
    {
        const double x = points[i][0], y = points[i][1];
        const glm::dvec3 aniso = anisotropicDisplacement(sol, x, y) - anisoRef;
        const glm::dvec3 iso = isotropicDisplacement(burgersEdge, burgersScrew, nu, x, y) - isoRef;
        const double d = glm::length(aniso - iso);
        const double scale = std::max(1.0, glm::length(iso));
        if (d > tol * scale)
        {
            throw std::runtime_error(std::string(label) + " at (" + std::to_string(x) + "," + std::to_string(y)
                + "): anisotropic solver disagrees with the isotropic closed form "
                "(gauge-corrected delta=" + std::to_string(d) + ", tolerance=" + std::to_string(tol * scale) + ")");
        }
    }
}
}

int main()
{
    try
    {
        // E=100 GPa, nu=0.3 expressed as cubic constants, offset by a small
        // deliberate (not random-noise) 3% change in C44 so the sextic roots
        // stay well separated -- exact isotropy makes the sextic degenerate
        // (BABEL documents the same limitation), and perturbing with tiny
        // random noise instead leaves the eigenvectors numerically
        // ill-conditioned (nearly parallel) rather than well resolved. The
        // resulting field is only approximately isotropic, so the tolerance
        // below allows for genuine (not just numerical) O(3%) anisotropic
        // deviation from the pure-isotropic reference.
        const double mu = 100.0 / (2.0 * 1.3);
        const double lambda = 100.0 * 0.3 / (1.3 * 0.4);
        const double c11 = lambda + 2.0 * mu;
        const double c12 = lambda;
        const double c44 = mu;
        const double nu = 0.3;
        const ElasticTensor nearIsotropic = makeCubicElasticTensor(c11, c12, c44 * 1.03);

        const double points[][2] = {{3.0, 1.5}, {-2.0, 4.0}, {1.0, -3.0}, {-4.0, -1.0}, {6.0, 0.5}};
        const int pointCount = (int)(sizeof(points) / sizeof(points[0]));

        // Edge dislocation, Burgers along x, line along z.
        {
            const StrohSolution sol = solveStroh(nearIsotropic, glm::dvec3(2.55, 0, 0));
            if (!sol.valid)
                throw std::runtime_error("Edge Stroh solve failed: " + sol.error);
            const glm::dvec3 circuit = circuitIntegral(sol, 5.0, 3600);
            if (glm::length(circuit - glm::dvec3(2.55, 0, 0)) > 0.02)
                throw std::runtime_error("Edge Burgers-vector circuit did not close");
            expectShapeClose("edge", sol, 2.55, 0.0, nu, points, pointCount, 0.05);
        }

        // Screw dislocation, Burgers along z, line along z.
        {
            const StrohSolution sol = solveStroh(nearIsotropic, glm::dvec3(0, 0, 2.55));
            if (!sol.valid)
                throw std::runtime_error("Screw Stroh solve failed: " + sol.error);
            const glm::dvec3 circuit = circuitIntegral(sol, 5.0, 3600);
            if (glm::length(circuit - glm::dvec3(0, 0, 2.55)) > 0.02)
                throw std::runtime_error("Screw Burgers-vector circuit did not close");
            expectShapeClose("screw", sol, 0.0, 2.55, nu, points, pointCount, 0.05);
        }

        // Mixed dislocation (45 degrees): isotropic elasticity is linear, so
        // this must equal the edge + screw solutions superposed.
        {
            const StrohSolution sol = solveStroh(nearIsotropic, glm::dvec3(1.8, 0, 1.8));
            if (!sol.valid)
                throw std::runtime_error("Mixed Stroh solve failed: " + sol.error);
            expectShapeClose("mixed", sol, 1.8, 1.8, nu, points, pointCount, 0.05);
        }

        // A genuinely anisotropic material (real Cu single-crystal constants,
        // GPa, A = 2*C44/(C11-C12) = 3.2) should still produce a
        // self-consistent solution: the Burgers vector circuit must close
        // even though there is no isotropic reference to compare the field
        // shape against.
        {
            const ElasticTensor copper = makeCubicElasticTensor(168.4, 121.4, 75.4);
            const StrohSolution sol = solveStroh(copper, glm::dvec3(2.556, 0, 0));
            if (!sol.valid)
                throw std::runtime_error("Anisotropic Cu Stroh solve failed: " + sol.error);
            const glm::dvec3 circuit = circuitIntegral(sol, 5.0, 720);
            const double err = glm::length(circuit - glm::dvec3(2.556, 0, 0)) / 2.556;
            if (err > 0.02)
                throw std::runtime_error("Anisotropic Cu Burgers-vector circuit did not close: "
                    "got (" + std::to_string(circuit.x) + "," + std::to_string(circuit.y)
                    + "," + std::to_string(circuit.z) + "), expected (2.556,0,0)");
        }

        // Degenerate (exactly isotropic) input must be refused outright
        // rather than silently returning a garbage solution, unless the
        // caller explicitly asks for the noise workaround.
        {
            const ElasticTensor exactIsotropic = makeCubicElasticTensor(c11, c12, c44);
            const StrohSolution refused = solveStroh(exactIsotropic, glm::dvec3(2.55, 0, 0));
            if (refused.valid)
                throw std::runtime_error("Exactly isotropic input should have been refused");
            const StrohSolution withNoise = solveStroh(exactIsotropic, glm::dvec3(2.55, 0, 0), 7, 1e-3);
            if (!withNoise.valid)
                throw std::runtime_error("Isotropic input with noise workaround should have succeeded: " + withNoise.error);
        }

        std::cout << "Anisotropic dislocation regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
