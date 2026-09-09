#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>
#include <cmath>
#include <cstring>
#include <functional>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// -- 2x2 matrix utilities ----------------------------------------------------

struct Mat2 { double m[2][2]; };

int mat2DetInt(const int m[2][2]);

void applyMat(const int n[2][2], const double basis[2][2], double out[2][2]);

Mat2 rotation2D(double theta);

double angleOf(double x, double y);

double wrapDegPm180(double deg);

double vecLen(double x, double y);

// -- Lattice deduplication ----------------------------------------------------

struct LatticeKey
{
    double l1, l2, area;
    bool operator==(const LatticeKey& o) const
    {
        return std::abs(l1 - o.l1) < 1e-6 &&
               std::abs(l2 - o.l2) < 1e-6 &&
               std::abs(area - o.area) < 1e-6;
    }
};

struct LatticeKeyHash
{
    size_t operator()(const LatticeKey& k) const
    {
        auto h = [](double v) -> size_t {
            long long bits = static_cast<long long>(v * 1e6);
            return std::hash<long long>()(bits);
        };
        size_t seed = h(k.l1);
        seed ^= h(k.l2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h(k.area) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

LatticeKey equivalentLatticeKey(const double v[2][2]);

struct SupercellEntry
{
    int mat[2][2];
    double vecs[2][2];
    int det;
};

// -- Supercell generation and strain matching ---------------------------------

std::vector<SupercellEntry> generateUniqueSupercells(
    const double basis[2][2], int nmax, int maxCells);

bool strainComponents(const double v[2][2], const double u[2][2],
                      double& exx, double& eyy, double& exy);

double meanAbsStrain(double exx, double eyy, double exy);

double cubicElasticDensity(double exx, double eyy, double exy, const float c[6][6]);

// -- Structure manipulation ---------------------------------------------------

void get2DBasis(const Structure& s, double basis[2][2]);

void buildMat3From2x2(const int m2[2][2], int m3[3][3]);

Structure makeSupercell2D(const Structure& base, const int mat2[2][2]);

Structure repeatLayersZ(const Structure& base, int layers);

Structure applyTransform2D(const Structure& b, const double F[2][2]);

Structure assembleInterface(const Structure& aSuper, const Structure& bStrained,
                            double zGap, double vacuum);

Structure repeatInterfaceXY(const Structure& iface, int rx, int ry);

// -- Orientation relationship -------------------------------------------------

float orientationAngleFromPlaneDir(const Structure& s,
                                   const float hkl[3],
                                   const float uvw[3]);
