#pragma once

#include <array>
#include <complex>
#include <string>
#include <vector>
#include <glm/glm.hpp>

namespace atomforge::electronic
{
constexpr double bohr = 0.529177210903;
constexpr double hartree = 27.211386245988;
constexpr double pi = 3.14159265358979323846;

// Lattice vectors are columns. Values are x-fastest. Coordinates are Angstrom.
// Periodic grids exclude the endpoint; finite grids include both endpoints.
struct Grid
{
    std::array<int, 3> shape{2, 2, 2};
    glm::dmat3 cell{1.0};
    glm::dvec3 origin{0.0};
    bool periodic = false;
    std::string name = "field";
    std::string unit = "raw";
    std::vector<double> values;

    void validate() const;
    std::size_t index(int x, int y, int z) const;
    glm::dmat3 steps() const;
    glm::dvec3 position(int x, int y, int z) const;
    double weight(int x, int y, int z) const;
    double sample(const glm::dvec3& point) const;
};

struct Site { int number = 0; glm::dvec3 position{0.0}; };
struct Volume { std::vector<Site> sites; std::vector<Grid> fields; };
struct Mesh { std::vector<glm::dvec3> vertices; std::vector<double> colors; };
struct Basin { double integral = 0.0; double volume = 0.0; };
struct Peak { glm::dvec3 position; double value; };
struct Reflection { glm::ivec3 hkl; std::complex<double> value; };
struct EwaldResult { std::vector<double> potentials; double energy = 0.0; };

Volume loadVolume(const std::string& path, const std::string& quantity = "auto",
                  const std::string& coordinateUnit = "bohr");
void saveVolume(const Volume& volume, const std::string& path, const std::string& format);
void requireAligned(const Grid& a, const Grid& b);
Grid arithmetic(const Grid& a, const Grid& b, const std::string& operation);
Grid scale(const Grid& a, double factor);
Grid resample(const Grid& source, const Grid& target);
Grid periodicEndpoints(const Grid& source, double tolerance = 1e-8);
Grid smooth(const Grid& grid, double sigma, int radius);
std::array<Grid, 3> gradient(const Grid& grid);
Grid laplacian(const Grid& grid);
std::array<Grid, 3> energyDensity(const Grid& density, double floor = 1e-12);
std::vector<glm::dvec2> lineProfile(const Grid& grid, glm::dvec3 start, glm::dvec3 end, int count);
std::vector<glm::dvec2> planarAverage(const Grid& grid, int axis);
std::vector<glm::dvec2> macroscopicAverage(const Grid& grid, int axis, int window);
Grid section(const Grid& grid, glm::dvec3 origin, glm::dvec3 u, glm::dvec3 v, int nu, int nv);
std::vector<glm::dvec3> contours(const Grid& slice, double level);
double integrate(const Grid& grid);
double integrateSphere(const Grid& grid, glm::dvec3 center, double radius);
std::vector<Basin> voronoiIntegrate(const Grid& grid, const std::vector<glm::dvec3>& sites);
std::vector<Peak> peaks(const Grid& grid, int limit = 0);
glm::dvec3 minimumImage(glm::dvec3 delta, const glm::dmat3& cell);
std::vector<Reflection> structureFactors(const Grid& grid, const std::vector<glm::ivec3>& indices);
Grid fourierSynthesis(const Grid& geometry, const std::vector<Reflection>& reflections);
Grid patterson(const Grid& grid);
// Explicit user-supplied scattering factors, occupancies, isotropic B in A^2.
std::vector<Reflection> atomicStructureFactors(const Grid& geometry,
    const std::vector<glm::dvec3>& positions, const std::vector<double>& occupancy,
    const std::vector<double>& bFactors, const std::vector<glm::ivec3>& indices,
    const std::vector<std::complex<double>>& scattering);
EwaldResult ewald(const glm::dmat3& cell, const std::vector<glm::dvec3>& sites,
    const std::vector<double>& charges, double alpha, double realCutoff, double reciprocalCutoff);
Mesh isosurface(const Grid& grid, double level, const Grid* color = nullptr);
void saveMesh(const Mesh& mesh, const std::string& path);
}
