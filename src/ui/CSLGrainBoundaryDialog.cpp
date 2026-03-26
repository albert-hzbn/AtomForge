#include "ui/CSLGrainBoundaryDialog.h"

#include "ElementData.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "imgui.h"

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{

// ── Math utilities ──────────────────────────────────────────────

static const double kPi = 3.14159265358979323846;

int gcd2(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0) { int t = a % b; a = b; b = t; }
    return a == 0 ? 1 : a;
}

int gcdN(const int* v, int n)
{
    int g = 0;
    for (int i = 0; i < n; i++) g = gcd2(g, std::abs(v[i]));
    return g == 0 ? 1 : g;
}

bool coPrime(int a, int b) { return gcd2(a, b) <= 1; }

bool isInteger(double x, double tol = 1e-5)
{
    return std::abs(x - std::round(x)) < tol;
}

bool isIntegerVec(const double* v, int n, double tol = 1e-5)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += (v[i] - std::round(v[i])) * (v[i] - std::round(v[i]));
    return std::sqrt(sum) < tol;
}

int getSmallestMultiplier(const double* v, int n, int maxN = 10000)
{
    for (int m = 1; m <= maxN; m++)
    {
        bool ok = true;
        for (int i = 0; i < n; i++)
        {
            if (!isInteger(m * v[i]))
            {
                ok = false;
                break;
            }
        }
        if (ok) return m;
    }
    return 1;
}

void reduceIntVec(int* v, int n)
{
    int g = gcdN(v, n);
    if (g > 1)
        for (int i = 0; i < n; i++) v[i] /= g;
}

// ── 3x3 row-major matrix ───────────────────────────────────────

struct M3
{
    double v[3][3];
    M3() { std::memset(v, 0, sizeof(v)); }
    double* operator[](int i) { return v[i]; }
    const double* operator[](int i) const { return v[i]; }
};

M3 eye3()
{
    M3 m;
    m[0][0] = m[1][1] = m[2][2] = 1.0;
    return m;
}

M3 transpose3(const M3& a)
{
    M3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r[i][j] = a[j][i];
    return r;
}

M3 mul3(const M3& a, const M3& b)
{
    M3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                r[i][j] += a[i][k] * b[k][j];
    return r;
}

double det3(const M3& m)
{
    return m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
         - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
         + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

M3 inv3(const M3& m)
{
    M3 r;
    double d = det3(m);
    if (std::abs(d) < 1e-15) return eye3();
    double id = 1.0 / d;
    r[0][0] =  (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * id;
    r[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1]) * id;
    r[0][2] =  (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * id;
    r[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0]) * id;
    r[1][1] =  (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * id;
    r[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0]) * id;
    r[2][0] =  (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * id;
    r[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0]) * id;
    r[2][2] =  (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * id;
    return r;
}

M3 roundM3(const M3& m, int decimals = 12)
{
    M3 r;
    double scale = 1.0;
    for (int i = 0; i < decimals; i++) scale *= 10.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r[i][j] = std::round(m[i][j] * scale) / scale;
    return r;
}

void roundToInt3x3(const M3& m, int out[3][3])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            out[i][j] = (int)std::round(m[i][j]);
}

M3 fromInt3x3(const int a[3][3])
{
    M3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r[i][j] = (double)a[i][j];
    return r;
}

double dot3(const double* a, const double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double norm3(const double* v)
{
    return std::sqrt(dot3(v, v));
}

// Solve A @ x = b via Cramer's rule
void solve3x3(const M3& A, const double b[3], double x[3])
{
    double d = det3(A);
    if (std::abs(d) < 1e-15) { x[0] = x[1] = x[2] = 0; return; }
    for (int col = 0; col < 3; col++)
    {
        M3 Ai = A;
        for (int row = 0; row < 3; row++)
            Ai[row][col] = b[row];
        x[col] = det3(Ai) / d;
    }
}

// Convert a double vector to the smallest equivalent integer vector.
void toIntegerVec3(const double in[3], int out[3])
{
    int mult = getSmallestMultiplier(in, 3);
    for (int i = 0; i < 3; i++)
        out[i] = (int)std::round(in[i] * mult);
    reduceIntVec(out, 3);
}

// ── CSL computation (ported from aimsgb) ────────────────────────

// Unimodular matrices from aimsgb/grain_bound.py
static const int UNIMODULAR[5][3][3] = {
    {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
    {{ 1, 0, 1}, { 0, 1, 0}, { 0, 1, 1}},
    {{ 1, 0, 1}, { 0, 1, 0}, { 0, 1,-1}},
    {{ 1, 0, 1}, { 0, 1, 0}, {-1, 1, 0}},
    {{ 1, 0, 1}, { 1, 1, 0}, { 1, 1, 1}},
};

// Rodrigues' rotation formula (angle in degrees)
M3 getRotateMatrix(const int axis[3], double angleDeg)
{
    double ax[3] = {(double)axis[0], (double)axis[1], (double)axis[2]};
    double n = norm3(ax);
    if (n < 1e-12) return eye3();
    ax[0] /= n; ax[1] /= n; ax[2] /= n;

    double angle = angleDeg * kPi / 180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);

    // K matrix (skew-symmetric)
    M3 K;
    K[0][1] = -ax[2]; K[0][2] =  ax[1];
    K[1][0] =  ax[2]; K[1][2] = -ax[0];
    K[2][0] = -ax[1]; K[2][1] =  ax[0];

    // R = I + K*sin + K*K*(1-cos)
    M3 KK = mul3(K, K);
    M3 R = eye3();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] += K[i][j] * s + KK[i][j] * (1.0 - c);
    return R;
}

// aimsgb: o_lattice_to_csl (with @transpose_matrix decorator)
// The function body works on the transposed matrix.
M3 oLatticeToCsl(M3 oLattice, double n)
{
    // Apply decorator: transpose input
    M3 csl = transpose3(oLattice);

    if (n < 0)
    {
        for (int j = 0; j < 3; j++) csl[0][j] *= -1.0;
        n *= -1.0;
    }

    // Main loop: make rows of csl[i] integers with product of multipliers <= n
    for (int iter = 0; iter < 200; iter++)
    {
        int m[3];
        for (int i = 0; i < 3; i++)
            m[i] = getSmallestMultiplier(csl[i], 3);
        int mProd = m[0] * m[1] * m[2];

        if (mProd <= (int)std::round(n))
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    csl[i][j] *= m[i];
            if (mProd < (int)std::round(n))
            {
                int nInt = (int)std::round(n);
                int factor = nInt / mProd;
                for (int j = 0; j < 3; j++)
                    csl[0][j] *= factor;
            }
            break;
        }
        else
        {
            bool changed = false;
            for (int i = 0; i < 3 && !changed; i++)
            {
                for (int j = 0; j < 3 && !changed; j++)
                {
                    if (i == j || m[i] == 1 || m[j] == 1) continue;
                    int a = (m[i] <= m[j]) ? i : j;
                    int b = (m[i] <= m[j]) ? j : i;
                    int mB = m[b];
                    // plus_minus_gen(1, mB)
                    for (int ki = 1; ki < mB && !changed; ki++)
                    {
                        for (int sign = 0; sign < 2 && !changed; sign++)
                        {
                            int kk = (sign == 0) ? ki : -ki;
                            double handle[3];
                            for (int c = 0; c < 3; c++)
                                handle[c] = csl[a][c] + kk * csl[b][c];
                            if (getSmallestMultiplier(handle, 3) < m[a])
                            {
                                for (int c = 0; c < 3; c++)
                                    csl[a][c] += kk * csl[b][c];
                                changed = true;
                            }
                        }
                    }
                }
            }
            if (!changed)
            {
                // Rarely happens - just multiply each row by its multiplier
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        csl[i][j] *= m[i];
                break;
            }
        }
    }

    // Round to integers
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            csl[i][j] = std::round(csl[i][j]);

    // Reshape / simplify CSL
    for (int iter = 0; iter < 200; iter++)
    {
        bool changed = false;
        for (int i = 0; i < 3 && !changed; i++)
        {
            for (int j = 0; j < 3 && !changed; j++)
            {
                if (i == j) continue;
                // Try csl[i] += csl[j]
                {
                    double sumAbs[3], origAbs[3];
                    double sumMax = 0, origMax = 0, sumSum = 0, origSum = 0;
                    for (int c = 0; c < 3; c++)
                    {
                        sumAbs[c] = std::abs(csl[i][c] + csl[j][c]);
                        origAbs[c] = std::abs(csl[i][c]);
                        if (sumAbs[c] > sumMax) sumMax = sumAbs[c];
                        if (origAbs[c] > origMax) origMax = origAbs[c];
                        sumSum += sumAbs[c];
                        origSum += origAbs[c];
                    }
                    bool go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    while (go)
                    {
                        for (int c = 0; c < 3; c++)
                            csl[i][c] += csl[j][c];
                        changed = true;
                        sumMax = 0; origMax = 0; sumSum = 0; origSum = 0;
                        for (int c = 0; c < 3; c++)
                        {
                            sumAbs[c] = std::abs(csl[i][c] + csl[j][c]);
                            origAbs[c] = std::abs(csl[i][c]);
                            if (sumAbs[c] > sumMax) sumMax = sumAbs[c];
                            if (origAbs[c] > origMax) origMax = origAbs[c];
                            sumSum += sumAbs[c];
                            origSum += origAbs[c];
                        }
                        go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    }
                }
                // Try csl[i] -= csl[j]
                if (!changed)
                {
                    double sumAbs[3], origAbs[3];
                    double sumMax = 0, origMax = 0, sumSum = 0, origSum = 0;
                    for (int c = 0; c < 3; c++)
                    {
                        sumAbs[c] = std::abs(csl[i][c] - csl[j][c]);
                        origAbs[c] = std::abs(csl[i][c]);
                        if (sumAbs[c] > sumMax) sumMax = sumAbs[c];
                        if (origAbs[c] > origMax) origMax = origAbs[c];
                        sumSum += sumAbs[c];
                        origSum += origAbs[c];
                    }
                    bool go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    while (go)
                    {
                        for (int c = 0; c < 3; c++)
                            csl[i][c] -= csl[j][c];
                        changed = true;
                        sumMax = 0; origMax = 0; sumSum = 0; origSum = 0;
                        for (int c = 0; c < 3; c++)
                        {
                            sumAbs[c] = std::abs(csl[i][c] - csl[j][c]);
                            origAbs[c] = std::abs(csl[i][c]);
                            if (sumAbs[c] > sumMax) sumMax = sumAbs[c];
                            if (origAbs[c] > origMax) origMax = origAbs[c];
                            sumSum += sumAbs[c];
                            origSum += origAbs[c];
                        }
                        go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    }
                }
            }
        }
        if (!changed) break;
    }

    // Apply decorator: transpose output
    return transpose3(csl);
}

// aimsgb: reduce_csl (with @transpose_matrix)
M3 reduceCsl(M3 csl)
{
    // Transpose for decorator
    csl = transpose3(csl);

    // Round and reduce each row
    for (int i = 0; i < 3; i++)
    {
        int iv[3];
        for (int j = 0; j < 3; j++) iv[j] = (int)std::round(csl[i][j]);
        reduceIntVec(iv, 3);
        for (int j = 0; j < 3; j++) csl[i][j] = (double)iv[j];
    }

    return transpose3(csl);
}

// aimsgb: orthogonalize_csl (with @transpose_matrix)
// Sets column 2 of CSL along the rotation axis and orthogonalizes.
M3 orthogonalizeCsl(M3 csl, const int axis[3])
{
    // Transpose for decorator
    csl = transpose3(csl);

    double axisD[3] = {(double)axis[0], (double)axis[1], (double)axis[2]};

    // Solve csl^T @ c = axis  (but csl here is already transposed, so csl^T is the original)
    M3 cslT = transpose3(csl);
    double c[3];
    solve3x3(cslT, axisD, c);

    if (!isIntegerVec(c, 3))
    {
        int mult = getSmallestMultiplier(c, 3);
        for (int i = 0; i < 3; i++) c[i] *= mult;
    }
    for (int i = 0; i < 3; i++) c[i] = std::round(c[i]);

    // Find index with smallest absolute non-zero value
    int ind = -1;
    double minAbsVal = 1e30;
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(c[i]) > 1e-8 && std::abs(c[i]) < minAbsVal)
        {
            minAbsVal = std::abs(c[i]);
            ind = i;
        }
    }
    if (ind < 0) ind = 2;

    // Swap so the picked index goes to position 2
    if (ind != 2)
    {
        double tmp[3];
        std::memcpy(tmp, csl[2], sizeof(tmp));
        for (int j = 0; j < 3; j++) csl[2][j] = -csl[ind][j];
        std::memcpy(csl[ind], tmp, sizeof(tmp));
        double ct = c[2]; c[2] = -c[ind]; c[ind] = ct;
    }

    // Set row 2 = c . csl (linear combination)
    double newRow2[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            newRow2[j] += c[i] * csl[i][j];
    std::memcpy(csl[2], newRow2, sizeof(newRow2));

    if (c[2] < 0)
        for (int j = 0; j < 3; j++) csl[1][j] *= -1;

    // Gram-Schmidt orthogonalization
    double u1[3];
    double n1 = norm3(csl[2]);
    if (n1 < 1e-12) return transpose3(csl);
    for (int i = 0; i < 3; i++) u1[i] = csl[2][i] / n1;

    // Project csl[1] and csl[0] onto the plane perpendicular to u1
    double y2_1[3], y2_2[3];
    double d1 = dot3(csl[1], u1);
    double d0 = dot3(csl[0], u1);
    for (int i = 0; i < 3; i++)
    {
        y2_1[i] = csl[1][i] - d1 * u1[i];
        y2_2[i] = csl[0][i] - d0 * u1[i];
    }

    int c0_1[3], c0_2[3];
    toIntegerVec3(y2_1, c0_1);
    toIntegerVec3(y2_2, c0_2);

    int sum1 = std::abs(c0_1[0]) + std::abs(c0_1[1]) + std::abs(c0_1[2]);
    int sum2 = std::abs(c0_2[0]) + std::abs(c0_2[1]) + std::abs(c0_2[2]);

    if (sum1 > sum2)
    {
        // Use csl[0]'s projection as the second basis
        double n2 = norm3(y2_2);
        if (n2 < 1e-12) return transpose3(csl);
        double u2[3];
        for (int i = 0; i < 3; i++) u2[i] = y2_2[i] / n2;

        double y3[3];
        double d1u1 = dot3(csl[1], u1);
        double d1u2 = dot3(csl[1], u2);
        for (int i = 0; i < 3; i++)
            y3[i] = csl[1][i] - d1u1 * u1[i] - d1u2 * u2[i];

        int y3i[3];
        toIntegerVec3(y3, y3i);
        for (int j = 0; j < 3; j++) csl[1][j] = (double)y3i[j];
        for (int j = 0; j < 3; j++) csl[0][j] = (double)c0_2[j];
    }
    else
    {
        double n2 = norm3(y2_1);
        if (n2 < 1e-12) return transpose3(csl);
        double u2[3];
        for (int i = 0; i < 3; i++) u2[i] = y2_1[i] / n2;

        double y3[3];
        double d0u1 = dot3(csl[0], u1);
        double d0u2 = dot3(csl[0], u2);
        for (int i = 0; i < 3; i++)
            y3[i] = csl[0][i] - d0u1 * u1[i] - d0u2 * u2[i];

        int y3i[3];
        toIntegerVec3(y3, y3i);
        for (int j = 0; j < 3; j++) csl[1][j] = (double)c0_1[j];
        for (int j = 0; j < 3; j++) csl[0][j] = (double)y3i[j];
    }

    // Round all to integers
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            csl[i][j] = std::round(csl[i][j]);

    return transpose3(csl);
}

// aimsgb: get_csl_matrix
M3 getCslMatrix(int sigma, const M3& rotateMatrix)
{
    M3 rInv = inv3(rotateMatrix);

    M3 t;
    for (int ui = 0; ui < 5; ui++)
    {
        M3 u = fromInt3x3(UNIMODULAR[ui]);
        // t = I - U @ inv(S) @ inv(R) @ S = I - U @ inv(R) (since S = I)
        M3 uR = mul3(u, rInv);
        t = eye3();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                t[i][j] -= uR[i][j];
        if (std::abs(det3(t)) > 1e-6)
            break;
    }

    M3 oLattice = roundM3(inv3(t), 12);
    double nVal = std::round((double)sigma / det3(oLattice) * 1e6) / 1e6;
    return oLatticeToCsl(oLattice, nVal);
}

// ── GB information struct ───────────────────────────────────────

struct SigmaCandidate
{
    int sigma;
    int m, n;
    float thetaDeg;
    int csl[3][3];            // integer CSL matrix
    std::array<int, 3> plane[3]; // columns of CSL = possible GB planes
};

std::string classifyBoundaryType(const int axis[3], const int plane[3])
{
    int dotVal = axis[0]*plane[0] + axis[1]*plane[1] + axis[2]*plane[2];
    int cx = axis[1]*plane[2] - axis[2]*plane[1];
    int cy = axis[2]*plane[0] - axis[0]*plane[2];
    int cz = axis[0]*plane[1] - axis[1]*plane[0];
    if (cx == 0 && cy == 0 && cz == 0) return "Twist";
    if (dotVal == 0) return "Tilt";
    return "Mixed";
}

std::vector<SigmaCandidate> computeGBInfo(const int axisIn[3], int maxSigma)
{
    // Reduce axis
    int axis[3] = {axisIn[0], axisIn[1], axisIn[2]};
    reduceIntVec(axis, 3);

    int axisNorm2 = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    if (axisNorm2 <= 0) return std::vector<SigmaCandidate>();

    float axisNorm = std::sqrt((float)axisNorm2);

    // Collect (theta, m, n) per sigma
    struct TMN { double theta; int m; int n; };
    std::map<int, std::vector<TMN>> sigmaTMN;

    int mMax = (int)std::ceil(std::sqrt(4.0 * (double)maxSigma));
    for (int m = 0; m <= mMax; m++)
    {
        for (int n = 0; n <= mMax; n++)
        {
            if (!coPrime(m, n)) continue;
            int sigmaRaw = m * m + n * n * axisNorm2;
            int sigma = sigmaRaw;
            while (sigma != 0 && sigma % 2 == 0) sigma /= 2;
            if (sigma <= 1 || sigma > maxSigma) continue;

            double theta;
            if (m == 0)
                theta = 180.0;
            else
                theta = 2.0 * std::atan((double)n * (double)axisNorm / (double)m) * 180.0 / kPi;

            TMN tmn = {theta, m, n};
            sigmaTMN[sigma].push_back(tmn);
        }
    }

    std::vector<SigmaCandidate> result;
    for (std::map<int, std::vector<TMN>>::iterator it = sigmaTMN.begin(); it != sigmaTMN.end(); ++it)
    {
        int sigma = it->first;
        std::vector<TMN>& entries = it->second;
        std::sort(entries.begin(), entries.end(), [](const TMN& a, const TMN& b) { return a.theta < b.theta; });

        double minTheta = entries[0].theta;

        // Build rotation matrix and CSL
        M3 rotMatrix = getRotateMatrix(axis, minTheta);
        M3 cslRaw = getCslMatrix(sigma, rotMatrix);
        M3 cslOrtho = orthogonalizeCsl(cslRaw, axis);
        M3 cslReduced = reduceCsl(cslOrtho);

        SigmaCandidate cand;
        cand.sigma = sigma;
        cand.m = entries[0].m;
        cand.n = entries[0].n;
        cand.thetaDeg = (float)minTheta;
        roundToInt3x3(cslReduced, cand.csl);

        // GB planes = columns of CSL
        for (int col = 0; col < 3; col++)
        {
            cand.plane[col][0] = cand.csl[0][col];
            cand.plane[col][1] = cand.csl[1][col];
            cand.plane[col][2] = cand.csl[2][col];
        }

        result.push_back(cand);
    }
    return result;
}

// ── Grain construction (ported from aimsgb) ─────────────────────

struct Grain
{
    std::vector<AtomSite> atoms;
    double cell[3][3]; // row-major: cell[i] = i-th lattice vector
};

Grain structureToGrain(const Structure& s)
{
    Grain g;
    g.atoms = s.atoms;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            g.cell[i][j] = s.cellVectors[i][j];
    return g;
}

Structure grainToStructure(const Grain& g)
{
    Structure s;
    s.atoms = g.atoms;
    s.hasUnitCell = true;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            s.cellVectors[i][j] = g.cell[i][j];
    s.cellOffset = {0.0, 0.0, 0.0};
    return s;
}

// Reduce a Structure to its primitive cell using spglib.
// Returns true on success.
bool reduceToPrimitive(Structure& s, double symprec = 1e-3)
{
#ifdef ATOMS_ENABLE_SPGLIB
    if (!s.hasUnitCell || s.atoms.empty()) return false;

    int natom = (int)s.atoms.size();

    // spglib expects lattice as row-major double[3][3] (a, b, c as rows)
    double lattice[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            lattice[i][j] = s.cellVectors[i][j];

    // Compute cell inverse for Cartesian→fractional conversion
    M3 cellM;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cellM[i][j] = lattice[i][j];
    M3 cellInv = inv3(cellM);

    // Fractional positions
    std::vector<double> positions(natom * 3);
    for (int ai = 0; ai < natom; ai++)
    {
        double cart[3] = {s.atoms[ai].x, s.atoms[ai].y, s.atoms[ai].z};
        double frac[3] = {0, 0, 0};
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
                frac[j] += cart[i] * cellInv[i][j];
        for (int j = 0; j < 3; j++)
        {
            frac[j] -= std::floor(frac[j]);
            if (frac[j] >= 1.0 - 1e-8) frac[j] = 0.0;
        }
        positions[ai * 3 + 0] = frac[0];
        positions[ai * 3 + 1] = frac[1];
        positions[ai * 3 + 2] = frac[2];
    }

    // Atom types
    std::vector<int> types(natom);
    for (int ai = 0; ai < natom; ai++)
        types[ai] = s.atoms[ai].atomicNumber;

    // Call spglib to find primitive cell (to_primitive=1, no_idealize=0)
    int num_prim = spg_standardize_cell(lattice,
                                        reinterpret_cast<double(*)[3]>(positions.data()),
                                        types.data(),
                                        natom,
                                        1, // to_primitive
                                        0, // no_idealize
                                        symprec);
    if (num_prim <= 0 || num_prim > natom) return false;

    // Rebuild structure from spglib output
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            s.cellVectors[i][j] = lattice[i][j];

    // Build a map from atomic number to element properties (from original atoms)
    struct ElemInfo { std::string sym; float r, g, b; };
    std::map<int, ElemInfo> zInfo;
    for (const auto& a : s.atoms)
        zInfo[a.atomicNumber] = {a.symbol, a.r, a.g, a.b};

    std::vector<AtomSite> newAtoms(num_prim);
    for (int ai = 0; ai < num_prim; ai++)
    {
        double frac[3] = {positions[ai*3], positions[ai*3+1], positions[ai*3+2]};
        // Convert fractional to Cartesian in new lattice
        double cart[3] = {0, 0, 0};
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
                cart[j] += frac[i] * lattice[i][j];
        newAtoms[ai].x = cart[0];
        newAtoms[ai].y = cart[1];
        newAtoms[ai].z = cart[2];
        newAtoms[ai].atomicNumber = types[ai];
        auto it = zInfo.find(types[ai]);
        if (it != zInfo.end())
        {
            newAtoms[ai].symbol = it->second.sym;
            newAtoms[ai].r = it->second.r;
            newAtoms[ai].g = it->second.g;
            newAtoms[ai].b = it->second.b;
        }
        else
        {
            newAtoms[ai].symbol = elementSymbol(types[ai]);
            newAtoms[ai].r = 0.5f;
            newAtoms[ai].g = 0.5f;
            newAtoms[ai].b = 0.5f;
        }
    }
    s.atoms = newAtoms;
    return true;
#else
    (void)s;
    (void)symprec;
    return false;
#endif
}

double cellVecLen(const double v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void cartToFrac(const double cart[3], const double cellInv[3][3], double frac[3])
{
    // frac = cart @ cellInv (pymatgen convention: frac[j] = sum_i cart[i] * cellInv[i][j])
    for (int j = 0; j < 3; j++)
    {
        frac[j] = 0;
        for (int i = 0; i < 3; i++)
            frac[j] += cart[i] * cellInv[i][j];
    }
}

void fracToCart(const double frac[3], const double cell[3][3], double cart[3])
{
    // cart[j] = sum_i frac[i] * cell[i][j]
    for (int j = 0; j < 3; j++)
    {
        cart[j] = 0;
        for (int i = 0; i < 3; i++)
            cart[j] += frac[i] * cell[i][j];
    }
}

void invertCell(const double cell[3][3], double inv[3][3])
{
    M3 m;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m[i][j] = cell[i][j];
    M3 mi = inv3(m);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            inv[i][j] = mi[i][j];
}

double wrapFrac(double f)
{
    f -= std::floor(f);
    if (f >= 1.0 - 1e-8) f = 0.0;
    if (f < -1e-10) f += 1.0;
    return f;
}

// Compute lattice parameters (lengths and angles) from the cell matrix.
void cellToParameters(const double cell[3][3],
                      double& a, double& b, double& c,
                      double& alpha, double& beta, double& gamma)
{
    a = cellVecLen(cell[0]);
    b = cellVecLen(cell[1]);
    c = cellVecLen(cell[2]);

    double cosAlpha = (a > 1e-12 && b > 1e-12 && c > 1e-12)
        ? dot3(cell[1], cell[2]) / (b * c) : 0.0;
    double cosBeta  = (a > 1e-12 && c > 1e-12)
        ? dot3(cell[0], cell[2]) / (a * c) : 0.0;
    double cosGamma = (a > 1e-12 && b > 1e-12)
        ? dot3(cell[0], cell[1]) / (a * b) : 0.0;

    cosAlpha = std::max(-1.0, std::min(1.0, cosAlpha));
    cosBeta  = std::max(-1.0, std::min(1.0, cosBeta));
    cosGamma = std::max(-1.0, std::min(1.0, cosGamma));

    alpha = std::acos(cosAlpha) * 180.0 / kPi;
    beta  = std::acos(cosBeta)  * 180.0 / kPi;
    gamma = std::acos(cosGamma) * 180.0 / kPi;
}

// Build a canonical cell matrix from lattice parameters.
// Convention: a along x, b in the xy-plane, c general.
// (Same as pymatgen Lattice.from_parameters)
void parametersToCell(double a, double b, double c,
                      double alpha, double beta, double gamma,
                      double cell[3][3])
{
    double alphaR = alpha * kPi / 180.0;
    double betaR  = beta  * kPi / 180.0;
    double gammaR = gamma * kPi / 180.0;

    double cosAlpha = std::cos(alphaR);
    double cosBeta  = std::cos(betaR);
    double cosGamma = std::cos(gammaR);
    double sinGamma = std::sin(gammaR);

    std::memset(cell, 0, sizeof(double) * 9);

    cell[0][0] = a;

    cell[1][0] = b * cosGamma;
    cell[1][1] = b * sinGamma;

    cell[2][0] = c * cosBeta;
    double cy = (std::abs(sinGamma) > 1e-10)
        ? c * (cosAlpha - cosBeta * cosGamma) / sinGamma
        : 0.0;
    cell[2][1] = cy;
    double cz2 = c * c - cell[2][0] * cell[2][0] - cy * cy;
    cell[2][2] = (cz2 > 0) ? std::sqrt(cz2) : 0.0;
}

// Canonicalize a grain's cell: redefined to the canonical form
// (a along x, b in the xy-plane) while preserving fractional coordinates.
// This matches pymatgen's Lattice.from_parameters behaviour that aimsgb
// applies at the end of every make_supercell call.
void canonicalizeGrain(Grain& g)
{
    double oldCellInv[3][3];
    invertCell(g.cell, oldCellInv);

    double a, b, c, alpha, beta, gamma;
    cellToParameters(g.cell, a, b, c, alpha, beta, gamma);

    double newCell[3][3];
    parametersToCell(a, b, c, alpha, beta, gamma, newCell);

    // Re-express atoms: fractional coordinates stay the same,
    // Cartesian coordinates are recomputed in the canonical cell.
    for (size_t ai = 0; ai < g.atoms.size(); ai++)
    {
        double cart[3] = {g.atoms[ai].x, g.atoms[ai].y, g.atoms[ai].z};
        double frac[3];
        cartToFrac(cart, oldCellInv, frac);
        for (int d = 0; d < 3; d++) frac[d] = wrapFrac(frac[d]);
        double newCart[3];
        fracToCart(frac, newCell, newCart);
        g.atoms[ai].x = newCart[0];
        g.atoms[ai].y = newCart[1];
        g.atoms[ai].z = newCart[2];
    }

    std::memcpy(g.cell, newCell, sizeof(g.cell));
}

// Make supercell with 3x3 integer scaling matrix M.
// new_lattice = M @ old_lattice
// For each atom at old frac f, new frac = (f + t) @ M^{-1} for integer t
Grain makeSupercell(const Grain& grain, const int M[3][3])
{
    Grain result;

    // New cell = M * old_cell
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            result.cell[i][j] = 0;
            for (int k = 0; k < 3; k++)
                result.cell[i][j] += (double)M[i][k] * grain.cell[k][j];
        }

    // Inverse of M (double)
    M3 Md = fromInt3x3(M);
    M3 Minv = inv3(Md);

    // Inverse of old cell
    double cellInv[3][3];
    invertCell(grain.cell, cellInv);

    // Search range: compute safe bound from row sums of |M|
    int searchRange = 0;
    for (int i = 0; i < 3; i++)
    {
        int rowSum = 0;
        for (int j = 0; j < 3; j++)
            rowSum += std::abs(M[i][j]);
        if (rowSum > searchRange) searchRange = rowSum;
    }
    searchRange += 1;

    for (size_t ai = 0; ai < grain.atoms.size(); ai++)
    {
        const AtomSite& atom = grain.atoms[ai];
        double cart[3] = {atom.x, atom.y, atom.z};
        double oldFrac[3];
        cartToFrac(cart, cellInv, oldFrac);

        for (int ti = -searchRange; ti <= searchRange; ti++)
        for (int tj = -searchRange; tj <= searchRange; tj++)
        for (int tk = -searchRange; tk <= searchRange; tk++)
        {
            double shifted[3] = {oldFrac[0] + ti, oldFrac[1] + tj, oldFrac[2] + tk};
            // newFrac = shifted @ Minv
            double newFrac[3];
            for (int j2 = 0; j2 < 3; j2++)
            {
                newFrac[j2] = 0;
                for (int k2 = 0; k2 < 3; k2++)
                    newFrac[j2] += shifted[k2] * Minv[k2][j2];
            }

            // Check [0, 1) with tolerance
            bool inside = true;
            for (int c = 0; c < 3; c++)
            {
                if (newFrac[c] < -1e-8 || newFrac[c] >= 1.0 - 1e-8)
                    inside = false;
            }
            if (!inside) continue;

            for (int c = 0; c < 3; c++)
                newFrac[c] = wrapFrac(newFrac[c]);

            double newCart[3];
            fracToCart(newFrac, result.cell, newCart);

            AtomSite newAtom = atom;
            newAtom.x = newCart[0];
            newAtom.y = newCart[1];
            newAtom.z = newCart[2];
            result.atoms.push_back(newAtom);
        }
    }

    // Canonicalize: redefine the cell to canonical form (a along x, b in
    // xy-plane) to match aimsgb's Lattice.from_parameters behaviour.
    canonicalizeGrain(result);

    return result;
}

// Make supercell with diagonal scaling matrix
Grain makeSupercellDiag(const Grain& grain, int na, int nb, int nc)
{
    int M[3][3] = {{na,0,0},{0,nb,0},{0,0,nc}};
    return makeSupercell(grain, M);
}

bool isOrthogonal(const double cell[3][3], double tol = 0.01)
{
    double d01 = cell[0][0]*cell[1][0] + cell[0][1]*cell[1][1] + cell[0][2]*cell[1][2];
    double d02 = cell[0][0]*cell[2][0] + cell[0][1]*cell[2][1] + cell[0][2]*cell[2][2];
    double d12 = cell[1][0]*cell[2][0] + cell[1][1]*cell[2][1] + cell[1][2]*cell[2][2];
    double l0 = cellVecLen(cell[0]);
    double l1 = cellVecLen(cell[1]);
    double l2 = cellVecLen(cell[2]);
    if (l0 < 1e-10 || l1 < 1e-10 || l2 < 1e-10) return true;
    return (std::abs(d01)/(l0*l1) < tol &&
            std::abs(d02)/(l0*l2) < tol &&
            std::abs(d12)/(l1*l2) < tol);
}

// Force the grain's cell to be orthogonal in the given direction.
// Adapted from aimsgb Grain.set_orthogonal_grain
Grain setOrthogonalGrain(const Grain& g, int direction)
{
    Grain result;
    std::memcpy(result.cell, g.cell, sizeof(result.cell));

    // The other two vectors
    int i1 = (direction + 1) % 3;
    int i2 = (direction + 2) % 3;
    double* v1 = result.cell[i1];
    double* v2 = result.cell[i2];

    // new_dir = projection of cell[direction] onto cross(v1, v2)
    double cross[3] = {
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    };
    double cn = norm3(cross);
    if (cn < 1e-12) { result.atoms = g.atoms; return result; }
    for (int i = 0; i < 3; i++) cross[i] /= cn;

    double proj = dot3(g.cell[direction], cross);
    for (int i = 0; i < 3; i++)
        result.cell[direction][i] = proj * cross[i];

    // Re-express atoms in new fractional coords
    double oldInv[3][3], newInv[3][3];
    invertCell(g.cell, oldInv);
    invertCell(result.cell, newInv);

    result.atoms.reserve(g.atoms.size());
    for (size_t ai = 0; ai < g.atoms.size(); ai++)
    {
        double cart[3] = {g.atoms[ai].x, g.atoms[ai].y, g.atoms[ai].z};
        // Keep Cartesian positions, they'll be wrapped into new cell
        double frac[3];
        cartToFrac(cart, newInv, frac);
        for (int c = 0; c < 3; c++) frac[c] = wrapFrac(frac[c]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);
        AtomSite a = g.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }
    return result;
}

// Create grain B by rotating grain A by 180° around a Cartesian axis.
// Tries all three axes and picks the first one that creates a different structure.
Grain getBFromA(const Grain& grainA)
{
    double cellInv[3][3];
    invertCell(grainA.cell, cellInv);

    for (int tryAxis = 0; tryAxis < 3; tryAxis++)
    {
        Grain b;
        std::memcpy(b.cell, grainA.cell, sizeof(b.cell));

        bool different = false;
        b.atoms.reserve(grainA.atoms.size());
        for (size_t ai = 0; ai < grainA.atoms.size(); ai++)
        {
            double cart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
            // 180° rotation around axes: negate the other two components
            for (int c = 0; c < 3; c++)
                if (c != tryAxis)
                    cart[c] = -cart[c];

            // Wrap into cell
            double frac[3];
            cartToFrac(cart, cellInv, frac);
            for (int c = 0; c < 3; c++) frac[c] = wrapFrac(frac[c]);
            double newCart[3];
            fracToCart(frac, grainA.cell, newCart);

            // Check if position changed
            double origCart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
            double origFrac[3];
            cartToFrac(origCart, cellInv, origFrac);
            for (int c = 0; c < 3; c++) origFrac[c] = wrapFrac(origFrac[c]);

            double df = 0;
            for (int c = 0; c < 3; c++)
            {
                double dd = frac[c] - origFrac[c];
                if (dd > 0.5) dd -= 1.0;
                if (dd < -0.5) dd += 1.0;
                df += dd * dd;
            }
            if (df > 1e-6) different = true;

            AtomSite a = grainA.atoms[ai];
            a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
            b.atoms.push_back(a);
        }

        if (different)
            return b;
    }

    // Fallback: just return a copy rotated around axis 0
    Grain b;
    std::memcpy(b.cell, grainA.cell, sizeof(b.cell));
    b.atoms.reserve(grainA.atoms.size());
    double cellInv2[3][3];
    invertCell(grainA.cell, cellInv2);
    for (size_t ai = 0; ai < grainA.atoms.size(); ai++)
    {
        double cart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
        cart[1] = -cart[1]; cart[2] = -cart[2];
        double frac[3];
        cartToFrac(cart, cellInv2, frac);
        for (int c = 0; c < 3; c++) frac[c] = wrapFrac(frac[c]);
        double newCart[3];
        fracToCart(frac, grainA.cell, newCart);
        AtomSite a = grainA.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        b.atoms.push_back(a);
    }
    return b;
}

// Stack two grains along the given direction with vacuum and gap.
// Follows aimsgb's Grain.stack_grains: build a combined cell from lattice
// parameters, re-express grain A atoms as fracs in the new cell, shift
// grain B Cartesian coords by (lenA + gap) along the stacking direction,
// then express as fracs in the new cell.
Grain stackGrains(const Grain& grainA, const Grain& grainB, int dir, double vacuum, double gap)
{
    Grain result;

    // Compute lattice parameters for both grains.
    double aA, bA, cA, alphaA, betaA, gammaA;
    cellToParameters(grainA.cell, aA, bA, cA, alphaA, betaA, gammaA);

    double aB, bB, cB, alphaB, betaB, gammaB;
    cellToParameters(grainB.cell, aB, bB, cB, alphaB, betaB, gammaB);

    double abcA[3] = {aA, bA, cA};
    double abcB[3] = {aB, bB, cB};
    double angles[3] = {alphaB, betaB, gammaB};

    // Shift distance in Cartesian (matches aimsgb logic for direction != 1)
    double l;
    if (dir == 1)
    {
        double sinGamma = std::sin(gammaB * kPi / 180.0);
        l = (abcA[dir] + gap) * sinGamma;
    }
    else
    {
        l = abcA[dir] + gap;
    }

    // Combined cell parameters: extend the stacking direction
    double abcNew[3] = {abcA[0], abcA[1], abcA[2]};
    abcNew[dir] += abcB[dir] + 2.0 * gap + vacuum;

    // Build combined cell in canonical form
    parametersToCell(abcNew[0], abcNew[1], abcNew[2],
                     angles[0], angles[1], angles[2],
                     result.cell);

    double newCellInv[3][3];
    invertCell(result.cell, newCellInv);

    // Grain A atoms: convert Cartesian to frac in combined cell
    for (size_t ai = 0; ai < grainA.atoms.size(); ai++)
    {
        double cart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
        double frac[3];
        cartToFrac(cart, newCellInv, frac);
        for (int c = 0; c < 3; c++) frac[c] = wrapFrac(frac[c]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);

        AtomSite a = grainA.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }

    // Grain B atoms: shift in Cartesian by l_vector, then convert
    double lVec[3] = {0, 0, 0};
    lVec[dir] = l;

    for (size_t ai = 0; ai < grainB.atoms.size(); ai++)
    {
        double cart[3] = {grainB.atoms[ai].x + lVec[0],
                          grainB.atoms[ai].y + lVec[1],
                          grainB.atoms[ai].z + lVec[2]};
        double frac[3];
        cartToFrac(cart, newCellInv, frac);
        for (int c = 0; c < 3; c++) frac[c] = wrapFrac(frac[c]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);

        AtomSite a = grainB.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }

    return result;
}

// Remove overlapping atoms within a distance threshold.
int removeOverlaps(std::vector<AtomSite>& atoms, double distThreshold)
{
    if (distThreshold <= 0 || atoms.size() < 2) return 0;

    double distSq = distThreshold * distThreshold;
    std::vector<bool> keep(atoms.size(), true);
    int removed = 0;

    // Simple O(n^2) overlap removal - mark later atoms for removal
    for (size_t i = 0; i < atoms.size(); i++)
    {
        if (!keep[i]) continue;
        for (size_t j = i + 1; j < atoms.size(); j++)
        {
            if (!keep[j]) continue;
            double dx = atoms[i].x - atoms[j].x;
            double dy = atoms[i].y - atoms[j].y;
            double dz = atoms[i].z - atoms[j].z;
            if (dx*dx + dy*dy + dz*dz < distSq)
            {
                keep[j] = false;
                removed++;
            }
        }
    }

    if (removed > 0)
    {
        std::vector<AtomSite> kept;
        kept.reserve(atoms.size() - removed);
        for (size_t i = 0; i < atoms.size(); i++)
            if (keep[i]) kept.push_back(atoms[i]);
        atoms.swap(kept);
    }
    return removed;
}

// ── Build result ────────────────────────────────────────────────

struct BuildResult
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

void drawBuildResultSummary(const BuildResult& result)
{
    if (result.message.empty())
        return;

    if (result.success)
    {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.4f, 1.0f, 0.5f, 1.0f));
        ImGui::TextWrapped("OK: %s", result.message.c_str());
        ImGui::PopStyleColor();
    }
    else
    {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
        ImGui::TextWrapped("Error: %s", result.message.c_str());
        ImGui::PopStyleColor();
        return;
    }

    ImGui::Spacing();

    if (ImGui::BeginTable("##ResultSummary", 2, ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn("Value");

        auto row = [](const char* label, const char* fmt, ...) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label);
            ImGui::TableNextColumn();
            va_list args;
            va_start(args, fmt);
            ImGui::TextV(fmt, args);
            va_end(args);
        };

        row("Input atoms", "%d", result.inputAtoms);
        row("Output atoms", "%d", result.generatedAtoms);
        row("Removed overlaps", "%d", result.removedOverlap);
        row("Sigma", "%d", result.sigma);
        row("Misorientation angle", "%.4f deg", result.thetaDeg);
        row("Boundary type", "%s", result.boundaryType.empty() ? "unknown" : result.boundaryType.c_str());

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Separator();
        ImGui::TableNextColumn();
        ImGui::Separator();

        row("Axis [u v w]", "[%d %d %d]", result.axis[0], result.axis[1], result.axis[2]);
        row("GB plane (h k l)", "(%d %d %d)", result.gbPlane[0], result.gbPlane[1], result.gbPlane[2]);
        row("CSL matrix", "[%d %d %d]  [%d %d %d]  [%d %d %d]",
            result.csl[0][0], result.csl[0][1], result.csl[0][2],
            result.csl[1][0], result.csl[1][1], result.csl[1][2],
            result.csl[2][0], result.csl[2][1], result.csl[2][2]);
        row("Grain A / Grain B", "%d / %d unit cells", result.ucA, result.ucB);
        row("Vacuum padding", "%.3f A", result.vacuumPadding);
        row("Gap", "%.3f A", result.gap);
        row("Overlap distance", "%.3f A", result.overlapDist);

        ImGui::EndTable();
    }
}

} // namespace

// ── Dialog methods ──────────────────────────────────────────────

CSLGrainBoundaryDialog::CSLGrainBoundaryDialog() = default;

CSLGrainBoundaryDialog::~CSLGrainBoundaryDialog()
{
    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void CSLGrainBoundaryDialog::initRenderResources(Renderer& renderer)
{
    m_renderer        = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vao, m_previewCylinder->vao);
    m_previewShadow   = createShadowMap(1, 1);
    m_glReady         = true;
}

void CSLGrainBoundaryDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewW && h == m_previewH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    glGenTextures(1, &m_previewColorTex);
    glBindTexture(GL_TEXTURE_2D, m_previewColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_previewDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_previewDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_previewFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_previewColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_previewDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_previewW = w;
    m_previewH = h;
}

void CSLGrainBoundaryDialog::rebuildPreviewBuffers(const Structure& s,
                                                    const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    if (!m_glReady || s.atoms.empty())
        return;

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    StructureInstanceData data = buildStructureInstanceData(
        s, false, kIdent, radii, shininess);

    std::array<bool,119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    m_previewBufDirty = false;
}

void CSLGrainBoundaryDialog::renderPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensurePreviewFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw      = m_camYaw;
    cam.pitch    = m_camPitch;
    cam.distance = m_camDistance;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          *m_previewCylinder, m_previewBuffers.bondCount);

    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          m_previewShadow, *m_previewSphere,
                          m_previewBuffers.atomCount);

    m_renderer->drawBoxLines(frame.projection, frame.view,
                             m_previewBuffers.lineVAO,
                             m_previewBuffers.boxLines.size());

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

void CSLGrainBoundaryDialog::autoFitPreviewCamera()
{
    m_camYaw   = 45.0f;
    m_camPitch = 35.0f;

    if (m_previewBuffers.atomCount == 0) {
        m_camDistance = 10.0f;
        return;
    }

    float maxR = 0.0f;
    for (size_t i = 0; i < m_previewBuffers.atomPositions.size(); ++i) {
        float r = (i < m_previewBuffers.atomRadii.size())
                  ? m_previewBuffers.atomRadii[i] : 0.0f;
        float d = glm::length(m_previewBuffers.atomPositions[i]
                              - m_previewBuffers.orbitCenter) + r;
        maxR = std::max(maxR, d);
    }
    maxR = std::max(maxR, 1.0f);

    const float halfFov = glm::radians(22.5f);
    float dist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_camDistance = dist;
}

void CSLGrainBoundaryDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

void CSLGrainBoundaryDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("CSL Grain Boundary", NULL, false, enabled))
        m_openRequested = true;
}

void CSLGrainBoundaryDialog::drawDialog(Structure& structure,
                                        const std::vector<glm::vec3>& elementColors,
                                        const std::vector<float>& elementRadii,
                                        const std::vector<float>& elementShininess,
                                        const std::function<void(Structure&)>& updateBuffers)
{
    // ── Static dialog state ────────────────────────────────────
    static Structure inputStructure;
    static char statusMsg[256] = "(no structure loaded)";
    static int axis[3] = {0, 0, 1};
    static int sigmaMax = 200;
    static std::vector<SigmaCandidate> sigmaCandidates;
    static int sigmaSelection = 0;
    static int planeSelection = 0;
    static int lastAxisForSigma[3] = {0, 0, 0};
    static int lastSigmaMaxForSigma = 0;
    static int ucA = 1;
    static int ucB = 1;
    static float vacuumPadding = 0.0f;
    static float gapDist = 0.0f;
    static float overlapDist = 0.0f;
    static bool conventionalCell = false;
    static BuildResult lastResult;

    // ── Handle pending drop ────────────────────────────────────
    if (!m_pendingDropPath.empty())
    {
        Structure loaded;
        std::string err;
        if (loadStructureFromFile(m_pendingDropPath, loaded, err))
        {
            if (!loaded.hasUnitCell)
            {
                std::snprintf(statusMsg, sizeof(statusMsg), "Error: structure has no unit cell");
                std::cout << "[CSL] Dropped file has no unit cell: " << m_pendingDropPath << std::endl;
            }
            else
            {
                inputStructure = loaded;
                std::snprintf(statusMsg, sizeof(statusMsg), "Loaded: %d atoms", (int)inputStructure.atoms.size());
                std::cout << "[CSL] Loaded input structure: " << m_pendingDropPath
                          << " (" << inputStructure.atoms.size() << " atoms)" << std::endl;
                m_previewBufDirty = true;
                if (m_glReady) {
                    rebuildPreviewBuffers(inputStructure, elementRadii, elementShininess);
                    autoFitPreviewCamera();
                }
            }
        }
        else
        {
            std::snprintf(statusMsg, sizeof(statusMsg), "Error: %s", err.c_str());
            std::cout << "[CSL] Drop-load failed: " << m_pendingDropPath
                      << " (" << err << ")" << std::endl;
        }
        m_pendingDropPath.clear();
    }

    // ── Open popup ─────────────────────────────────────────────
    if (m_openRequested)
    {
        ImGui::OpenPopup("CSL Grain Boundary Builder");
        m_openRequested = false;
        m_isOpen = true;
    }

    ImGui::SetNextWindowSize(ImVec2(780.0f, 580.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(620.0f, 400.0f), ImVec2(FLT_MAX, FLT_MAX));
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("CSL Grain Boundary Builder", &dialogOpen, ImGuiWindowFlags_None))
    {
        // ── Input Structure ────────────────────────────────────
        ImGui::SeparatorText("Input Structure");
        {
            constexpr float kPrevH = 260.0f;

            ImGui::InvisibleButton("##cslDropZone", ImVec2(-1.0f, kPrevH));
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 dropMin = ImGui::GetItemRectMin();
            ImVec2 dropMax = ImGui::GetItemRectMax();
            const bool dropZoneHovered = ImGui::IsItemHovered();
            const bool dropZoneActive  = ImGui::IsItemActive();
            dl->AddRect(dropMin, dropMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

            if (inputStructure.hasUnitCell && !inputStructure.atoms.empty())
            {
                // ── 3D preview ──
                if (m_glReady) {
                    if (m_previewBufDirty)
                        rebuildPreviewBuffers(inputStructure, elementRadii, elementShininess);

                    const float pad = 5.0f;
                    const ImVec2 prevSize(dropMax.x - dropMin.x - 2.0f * pad,
                                          dropMax.y - dropMin.y - 2.0f * pad);
                    const int pw = std::max(1, (int)prevSize.x);
                    const int ph = std::max(1, (int)prevSize.y);

                    renderPreviewToFBO(pw, ph);

                    const ImVec2 prevMin(dropMin.x + pad, dropMin.y + pad);
                    const ImVec2 prevMax(prevMin.x + prevSize.x, prevMin.y + prevSize.y);
                    dl->AddImage((ImTextureID)(intptr_t)m_previewColorTex,
                                 prevMin, prevMax,
                                 ImVec2(0, 1), ImVec2(1, 0));

                    // Orbit with left-drag
                    if (dropZoneActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f)) {
                        ImVec2 delta = ImGui::GetIO().MouseDelta;
                        m_camYaw   -= delta.x * 0.5f;
                        m_camPitch += delta.y * 0.5f;
                    }

                    // Zoom with scroll
                    if (dropZoneHovered) {
                        float wheel = ImGui::GetIO().MouseWheel;
                        if (wheel != 0.0f) {
                            m_camDistance -= wheel * m_camDistance * 0.1f;
                            m_camDistance  = std::max(Camera::kMinDistance,
                                                      std::min(Camera::kMaxDistance, m_camDistance));
                        }
                    }
                }

                // Info line below preview
                const auto& cv = inputStructure.cellVectors;
                double la = std::sqrt(cv[0][0]*cv[0][0] + cv[0][1]*cv[0][1] + cv[0][2]*cv[0][2]);
                double lb = std::sqrt(cv[1][0]*cv[1][0] + cv[1][1]*cv[1][1] + cv[1][2]*cv[1][2]);
                double lc = std::sqrt(cv[2][0]*cv[2][0] + cv[2][1]*cv[2][1] + cv[2][2]*cv[2][2]);
                ImGui::TextDisabled("a=%.3f  b=%.3f  c=%.3f A    %d atoms    Left-drag=orbit  Scroll=zoom",
                                    la, lb, lc, (int)inputStructure.atoms.size());
            }
            else
            {
                // Show drop hint
                ImVec2 mid((dropMin.x + dropMax.x) * 0.5f, (dropMin.y + dropMax.y) * 0.5f);
                dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                            ImVec2(mid.x - 100.0f, mid.y - 10.0f),
                            ImGui::GetColorU32(ImGuiCol_TextDisabled),
                            "Drop a crystal structure file here");
            }
        }

        // ── Misorientation ─────────────────────────────────────
        ImGui::SeparatorText("Misorientation");
        {
            ImGui::SetNextItemWidth(220.0f);
            ImGui::InputInt3("Axis [u v w]", axis);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Crystallographic rotation axis in Miller indices.");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Max sigma", &sigmaMax);
            if (sigmaMax < 3) sigmaMax = 3;
            ImGui::SameLine();
            if (ImGui::Button("Find"))
            {
                sigmaCandidates = computeGBInfo(axis, sigmaMax);
                sigmaSelection = 0;
                planeSelection = 0;
            }

            // Clear sigma list when axis or max changes
            if (axis[0] != lastAxisForSigma[0] || axis[1] != lastAxisForSigma[1] ||
                axis[2] != lastAxisForSigma[2] || sigmaMax != lastSigmaMaxForSigma)
            {
                sigmaCandidates.clear();
                sigmaSelection = 0;
                planeSelection = 0;
                lastAxisForSigma[0] = axis[0];
                lastAxisForSigma[1] = axis[1];
                lastAxisForSigma[2] = axis[2];
                lastSigmaMaxForSigma = sigmaMax;
            }

            if (!sigmaCandidates.empty())
            {
                auto sigmaGetter = [](void* data, int idx) -> const char* {
                    static char label[96];
                    std::vector<SigmaCandidate>* vec = static_cast<std::vector<SigmaCandidate>*>(data);
                    if (idx < 0 || idx >= (int)vec->size()) return "";
                    const SigmaCandidate& c = (*vec)[idx];
                    std::snprintf(label, sizeof(label), "S%d  m=%d n=%d  %.2f deg",
                                  c.sigma, c.m, c.n, c.thetaDeg);
                    return label;
                };
                if (sigmaSelection >= (int)sigmaCandidates.size())
                    sigmaSelection = 0;
                ImGui::SetNextItemWidth(-1.0f);
                ImGui::Combo("##sigma", &sigmaSelection, sigmaGetter,
                             &sigmaCandidates, (int)sigmaCandidates.size());
            }
            else
            {
                ImGui::TextDisabled("Set axis and press Find.");
            }
        }

        // ── Boundary Geometry ──────────────────────────────────
        ImGui::SeparatorText("Boundary Geometry");
        {
            if (!sigmaCandidates.empty() && sigmaSelection >= 0 && sigmaSelection < (int)sigmaCandidates.size())
            {
                const SigmaCandidate& sel = sigmaCandidates[sigmaSelection];

                // Plane dropdown: 3 planes (columns of CSL)
                // Build labels with type classification
                char planeLabels[3][96];
                for (int i = 0; i < 3; i++)
                {
                    std::string type = classifyBoundaryType(axis, sel.plane[i].data());
                    std::snprintf(planeLabels[i], sizeof(planeLabels[i]), "(%d %d %d)  %s",
                                  sel.plane[i][0], sel.plane[i][1], sel.plane[i][2], type.c_str());
                }
                auto planeLabelGetter = [](void* data, int idx) -> const char* {
                    if (idx < 0 || idx >= 3) return "";
                    return static_cast<char(*)[96]>(data)[idx];
                };

                if (planeSelection >= 3) planeSelection = 0;
                ImGui::SetNextItemWidth(260.0f);
                ImGui::Combo("GB plane", &planeSelection, planeLabelGetter, planeLabels, 3);

                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("GB plane normals derived from CSL matrix columns.\n"
                                      "Twist = plane parallel to axis, Tilt = plane perpendicular to axis.");

                ImGui::SameLine();
                ImGui::AlignTextToFramePadding();
                std::string btype = classifyBoundaryType(axis, sel.plane[planeSelection].data());
                ImGui::TextDisabled("Type: %s", btype.c_str());

                // Show CSL matrix
                ImGui::TextDisabled("CSL: [%d %d %d] [%d %d %d] [%d %d %d]",
                                    sel.csl[0][0], sel.csl[0][1], sel.csl[0][2],
                                    sel.csl[1][0], sel.csl[1][1], sel.csl[1][2],
                                    sel.csl[2][0], sel.csl[2][1], sel.csl[2][2]);
            }
            else
            {
                ImGui::TextDisabled("Select a Sigma candidate first.");
            }

            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Grain A (uc)", &ucA);
            if (ucA < 1) ucA = 1;
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Grain B (uc)", &ucB);
            if (ucB < 1) ucB = 1;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Number of unit cells for each grain along the GB normal.");
        }

        // ── Options ────────────────────────────────────────────
        ImGui::SeparatorText("Options");
        {
            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Vacuum (A)", &vacuumPadding, 0.1f, 0.0f, 50.0f, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Vacuum space at the end of the bicrystal cell.");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Gap (A)", &gapDist, 0.05f, 0.0f, 10.0f, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Extra distance between the two grains at the interface.");

            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Overlap dist (A)", &overlapDist, 0.05f, 0.0f, 5.0f, "%.3f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Atoms closer than this distance will be removed (0 = off).");

            ImGui::Checkbox("Conventional cell", &conventionalCell);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("If unchecked (default), reduce to the primitive cell.\nIf checked, keep the conventional (unreduced) cell.");
        }

        // ── Generate ───────────────────────────────────────────
        ImGui::Spacing();
        if (ImGui::Button("Generate Grain Boundary", ImVec2(-1.0f, 0.0f)))
        {
            if (!inputStructure.hasUnitCell || inputStructure.atoms.empty())
            {
                lastResult.success = false;
                lastResult.message = "Please load a structure with a unit cell first.";
            }
            else if (sigmaCandidates.empty())
            {
                lastResult.success = false;
                lastResult.message = "Generate and select a Sigma candidate first.";
            }
            else
            {
                const SigmaCandidate& sel = sigmaCandidates[sigmaSelection];
                const int* plane = sel.plane[planeSelection].data();
                int direction = planeSelection; // column index = stacking direction

                // Build grains using CSL supercell (aimsgb algorithm)
                Grain inputGrain = structureToGrain(inputStructure);

                // Scaling matrix = CSL^T (transpose of CSL)
                int cslT[3][3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        cslT[i][j] = sel.csl[j][i];

                Grain grainA = makeSupercell(inputGrain, cslT);
                lastResult.inputAtoms = (int)inputStructure.atoms.size();

                // Force orthogonal if needed
                if (!isOrthogonal(grainA.cell))
                    grainA = setOrthogonalGrain(grainA, direction);

                // Scale grain A by ucB for creating grain B, then scale grain A by ucA
                int scaleB[3] = {1, 1, 1};
                scaleB[direction] = ucB;
                Grain tempA = makeSupercellDiag(grainA, scaleB[0], scaleB[1], scaleB[2]);

                Grain grainB = getBFromA(tempA);
                // Normalize grain B
                int ones[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
                grainB = makeSupercell(grainB, ones);

                int scaleA[3] = {1, 1, 1};
                scaleA[direction] = ucA;
                grainA = makeSupercellDiag(grainA, scaleA[0], scaleA[1], scaleA[2]);

                // Stack grains
                Grain gb = stackGrains(grainA, grainB, direction, vacuumPadding, gapDist);

                // Remove overlaps
                int removed = removeOverlaps(gb.atoms, overlapDist);

                // Convert to Structure
                structure = grainToStructure(gb);

                // Reduce to primitive cell unless "Conventional cell" is checked
                if (!conventionalCell)
                    reduceToPrimitive(structure);

                lastResult.success = true;
                lastResult.sigma = sel.sigma;
                lastResult.thetaDeg = sel.thetaDeg;
                lastResult.generatedAtoms = (int)structure.atoms.size();
                lastResult.removedOverlap = removed;
                lastResult.boundaryType = classifyBoundaryType(axis, plane);
                std::memcpy(lastResult.axis, axis, sizeof(lastResult.axis));
                lastResult.gbPlane[0] = plane[0];
                lastResult.gbPlane[1] = plane[1];
                lastResult.gbPlane[2] = plane[2];
                std::memcpy(lastResult.csl, sel.csl, sizeof(lastResult.csl));
                lastResult.ucA = ucA;
                lastResult.ucB = ucB;
                lastResult.vacuumPadding = vacuumPadding;
                lastResult.gap = gapDist;
                lastResult.overlapDist = overlapDist;

                lastResult.message = "CSL grain boundary generated.";

                updateBuffers(structure);

                std::cout << "[Operation] Built CSL grain boundary: "
                          << "Sigma=" << sel.sigma
                          << ", theta=" << sel.thetaDeg << " deg"
                          << ", plane=(" << plane[0] << " " << plane[1] << " " << plane[2] << ")"
                          << ", type=" << lastResult.boundaryType
                          << ", atoms=" << lastResult.generatedAtoms
                          << std::endl;
            }

            if (!lastResult.success)
            {
                std::cout << "[Operation] CSL grain boundary generation failed: "
                          << lastResult.message << std::endl;
            }
        }

        // ── Result ─────────────────────────────────────────────
        ImGui::Spacing();
        if (!lastResult.message.empty())
        {
            ImGui::SeparatorText("Result");
            ImGui::BeginChild("##ResultScroll", ImVec2(0.0f, 0.0f), ImGuiChildFlags_Borders);
            drawBuildResultSummary(lastResult);
            ImGui::EndChild();
        }

        ImGui::EndPopup();
    }

    if (!dialogOpen)
    {
        m_isOpen = false;
    }
}
