#pragma once

#include "model/Structure.h"

#include <array>
#include <string>
#include <vector>

// ── 3x3 matrix type ────────────────────────────────────────────

struct M3
{
    double v[3][3];
    M3();
    double* operator[](int i) { return v[i]; }
    const double* operator[](int i) const { return v[i]; }
};

M3 eye3();
M3 transpose3(const M3& a);
M3 mul3(const M3& a, const M3& b);
double det3(const M3& m);
M3 inv3(const M3& m);
M3 roundM3(const M3& m, int decimals = 12);
void roundToInt3x3(const M3& m, int out[3][3]);
M3 fromInt3x3(const int a[3][3]);
double dot3(const double* a, const double* b);
double norm3(const double* v);
void solve3x3(const M3& A, const double b[3], double x[3]);
void toIntegerVec3(const double in[3], int out[3]);
void reduceIntVec(int* v, int n);

// ── Rotation and CSL ────────────────────────────────────────────

M3 getRotateMatrix(const int axis[3], double angleDeg);
M3 oLatticeToCsl(M3 oLattice, double n);
M3 reduceCsl(M3 csl);
M3 orthogonalizeCsl(M3 csl, const int axis[3]);
M3 getCslMatrix(int sigma, const M3& rotateMatrix);

// ── Sigma / GB information ──────────────────────────────────────

struct SigmaCandidate
{
    int sigma;
    int m, n;
    float thetaDeg;
    int csl[3][3];
    std::array<int, 3> plane[3];
};

std::string classifyBoundaryType(const int axis[3], const int plane[3]);
std::vector<SigmaCandidate> computeGBInfo(const int axisIn[3], int maxSigma);

// ── Grain type and operations ───────────────────────────────────

struct Grain
{
    std::vector<AtomSite> atoms;
    double cell[3][3];
};

Grain structureToGrain(const Structure& s);
Structure grainToStructure(const Grain& g);

bool reduceToPrimitive(Structure& s, double symprec = 1e-3);
bool reduceToPrimitiveGB(Structure& s, int stackDir, double tol = 0.25);

double cellVecLen(const double v[3]);
void cartToFrac(const double cart[3], const double cellInv[3][3], double frac[3]);
void fracToCart(const double frac[3], const double cell[3][3], double cart[3]);
void invertCell(const double cell[3][3], double inv[3][3]);
double wrapFrac(double f);

void cellToParameters(const double cell[3][3],
                      double& a, double& b, double& c,
                      double& alpha, double& beta, double& gamma);
void parametersToCell(double a, double b, double c,
                      double alpha, double beta, double gamma,
                      double cell[3][3]);

void canonicalizeGrain(Grain& g);
Grain makeSupercell(const Grain& grain, const int M[3][3]);
Grain makeSupercellDiag(const Grain& grain, int na, int nb, int nc);

bool isOrthogonal(const double cell[3][3], double tol = 0.01);
Grain setOrthogonalGrain(const Grain& g, int direction);
Grain getBFromA(const Grain& grainA);
Grain stackGrains(const Grain& grainA, const Grain& grainB, int dir, double vacuum, double gap);
int removeOverlaps(std::vector<AtomSite>& atoms, double distThreshold);

// ── Build result ────────────────────────────────────────────────

struct CSLBuildResult
{
    bool success = false;
    std::string message;
    std::string boundaryType;
    int sigma = 0;
    float thetaDeg = 0.0f;
    int inputAtoms = 0;
    int generatedAtoms = 0;
    int removedOverlap = 0;

    int axis[3];
    int gbPlane[3];
    int csl[3][3];
    int ucA, ucB;
    float vacuumPadding;
    float gap;
    float overlapDist;
};
