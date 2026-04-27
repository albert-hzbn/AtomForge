#include "algorithms/CSLComputation.h"

// CSL grain-boundary generation in AtomForge follows the AIMSGB workflow
// described in: Jianli Cheng, Jian Luo, and Kesong Yang, "Aimsgb: An
// Algorithm and Open-Source Python Library to Generate Periodic Grain
// Boundary Structures," Comput. Mater. Sci. 155, 92-103 (2018).

#include "ElementData.h"

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>

static const double kPi = 3.14159265358979323846;

// ── Helper math ─────────────────────────────────────────────────

static int gcd2(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0) { int t = a % b; a = b; b = t; }
    return a == 0 ? 1 : a;
}

static int gcdN(const int* v, int n)
{
    int g = 0;
    for (int i = 0; i < n; i++) g = gcd2(g, std::abs(v[i]));
    return g == 0 ? 1 : g;
}

static bool coPrime(int a, int b) { return gcd2(a, b) <= 1; }

static bool isInteger(double x, double tol = 1e-5)
{
    return std::abs(x - std::round(x)) < tol;
}

static bool isIntegerVec(const double* v, int n, double tol = 1e-5)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += (v[i] - std::round(v[i])) * (v[i] - std::round(v[i]));
    return std::sqrt(sum) < tol;
}

static int getSmallestMultiplier(const double* v, int n, int maxN = 10000)
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

// ── M3 ──────────────────────────────────────────────────────────

M3::M3() { std::memset(v, 0, sizeof(v)); }

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

M3 roundM3(const M3& m, int decimals)
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

void toIntegerVec3(const double in[3], int out[3])
{
    int mult = getSmallestMultiplier(in, 3);
    for (int i = 0; i < 3; i++)
        out[i] = (int)std::round(in[i] * mult);
    reduceIntVec(out, 3);
}

void reduceIntVec(int* v, int n)
{
    int g = gcdN(v, n);
    if (g > 1)
        for (int i = 0; i < n; i++) v[i] /= g;
}

// ── Unimodular matrices ─────────────────────────────────────────

static const int UNIMODULAR[5][3][3] = {
    {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
    {{ 1, 0, 1}, { 0, 1, 0}, { 0, 1, 1}},
    {{ 1, 0, 1}, { 0, 1, 0}, { 0, 1,-1}},
    {{ 1, 0, 1}, { 0, 1, 0}, {-1, 1, 0}},
    {{ 1, 0, 1}, { 1, 1, 0}, { 1, 1, 1}},
};

// ── Rotation and CSL functions ──────────────────────────────────

M3 getRotateMatrix(const int axis[3], double angleDeg)
{
    double ax[3] = {(double)axis[0], (double)axis[1], (double)axis[2]};
    double n = norm3(ax);
    if (n < 1e-12) return eye3();
    ax[0] /= n; ax[1] /= n; ax[2] /= n;

    double angle = angleDeg * kPi / 180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);

    M3 K;
    K[0][1] = -ax[2]; K[0][2] =  ax[1];
    K[1][0] =  ax[2]; K[1][2] = -ax[0];
    K[2][0] = -ax[1]; K[2][1] =  ax[0];

    M3 KK = mul3(K, K);
    M3 R = eye3();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] += K[i][j] * s + KK[i][j] * (1.0 - c);
    return R;
}

M3 oLatticeToCsl(M3 oLattice, double n)
{
    M3 csl = transpose3(oLattice);

    if (n < 0)
    {
        for (int j = 0; j < 3; j++) csl[0][j] *= -1.0;
        n *= -1.0;
    }

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
                    for (int ki = 1; ki < mB && !changed; ki++)
                    {
                        for (int sign = 0; sign < 2 && !changed; sign++)
                        {
                            int kk = (sign == 0) ? ki : -ki;
                            double handle[3];
                            for (int c2 = 0; c2 < 3; c2++)
                                handle[c2] = csl[a][c2] + kk * csl[b][c2];
                            if (getSmallestMultiplier(handle, 3) < m[a])
                            {
                                for (int c2 = 0; c2 < 3; c2++)
                                    csl[a][c2] += kk * csl[b][c2];
                                changed = true;
                            }
                        }
                    }
                }
            }
            if (!changed)
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        csl[i][j] *= m[i];
                break;
            }
        }
    }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            csl[i][j] = std::round(csl[i][j]);

    for (int iter = 0; iter < 200; iter++)
    {
        bool changed = false;
        for (int i = 0; i < 3 && !changed; i++)
        {
            for (int j = 0; j < 3 && !changed; j++)
            {
                if (i == j) continue;
                {
                    double sumAbs[3], origAbs[3];
                    double sumMax = 0, origMax = 0, sumSum = 0, origSum = 0;
                    for (int c2 = 0; c2 < 3; c2++)
                    {
                        sumAbs[c2] = std::abs(csl[i][c2] + csl[j][c2]);
                        origAbs[c2] = std::abs(csl[i][c2]);
                        if (sumAbs[c2] > sumMax) sumMax = sumAbs[c2];
                        if (origAbs[c2] > origMax) origMax = origAbs[c2];
                        sumSum += sumAbs[c2];
                        origSum += origAbs[c2];
                    }
                    bool go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    while (go)
                    {
                        for (int c2 = 0; c2 < 3; c2++)
                            csl[i][c2] += csl[j][c2];
                        changed = true;
                        sumMax = 0; origMax = 0; sumSum = 0; origSum = 0;
                        for (int c2 = 0; c2 < 3; c2++)
                        {
                            sumAbs[c2] = std::abs(csl[i][c2] + csl[j][c2]);
                            origAbs[c2] = std::abs(csl[i][c2]);
                            if (sumAbs[c2] > sumMax) sumMax = sumAbs[c2];
                            if (origAbs[c2] > origMax) origMax = origAbs[c2];
                            sumSum += sumAbs[c2];
                            origSum += origAbs[c2];
                        }
                        go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    }
                }
                if (!changed)
                {
                    double sumAbs[3], origAbs[3];
                    double sumMax = 0, origMax = 0, sumSum = 0, origSum = 0;
                    for (int c2 = 0; c2 < 3; c2++)
                    {
                        sumAbs[c2] = std::abs(csl[i][c2] - csl[j][c2]);
                        origAbs[c2] = std::abs(csl[i][c2]);
                        if (sumAbs[c2] > sumMax) sumMax = sumAbs[c2];
                        if (origAbs[c2] > origMax) origMax = origAbs[c2];
                        sumSum += sumAbs[c2];
                        origSum += origAbs[c2];
                    }
                    bool go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    while (go)
                    {
                        for (int c2 = 0; c2 < 3; c2++)
                            csl[i][c2] -= csl[j][c2];
                        changed = true;
                        sumMax = 0; origMax = 0; sumSum = 0; origSum = 0;
                        for (int c2 = 0; c2 < 3; c2++)
                        {
                            sumAbs[c2] = std::abs(csl[i][c2] - csl[j][c2]);
                            origAbs[c2] = std::abs(csl[i][c2]);
                            if (sumAbs[c2] > sumMax) sumMax = sumAbs[c2];
                            if (origAbs[c2] > origMax) origMax = origAbs[c2];
                            sumSum += sumAbs[c2];
                            origSum += origAbs[c2];
                        }
                        go = (sumMax < origMax) || (sumMax == origMax && sumSum < origSum);
                    }
                }
            }
        }
        if (!changed) break;
    }

    return transpose3(csl);
}

M3 reduceCsl(M3 csl)
{
    csl = transpose3(csl);
    for (int i = 0; i < 3; i++)
    {
        int iv[3];
        for (int j = 0; j < 3; j++) iv[j] = (int)std::round(csl[i][j]);
        reduceIntVec(iv, 3);
        for (int j = 0; j < 3; j++) csl[i][j] = (double)iv[j];
    }
    return transpose3(csl);
}

M3 orthogonalizeCsl(M3 csl, const int axis[3])
{
    csl = transpose3(csl);

    double axisD[3] = {(double)axis[0], (double)axis[1], (double)axis[2]};

    M3 cslT = transpose3(csl);
    double c[3];
    solve3x3(cslT, axisD, c);

    if (!isIntegerVec(c, 3))
    {
        int mult = getSmallestMultiplier(c, 3);
        for (int i = 0; i < 3; i++) c[i] *= mult;
    }
    for (int i = 0; i < 3; i++) c[i] = std::round(c[i]);

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

    if (ind != 2)
    {
        double tmp[3];
        std::memcpy(tmp, csl[2], sizeof(tmp));
        for (int j = 0; j < 3; j++) csl[2][j] = -csl[ind][j];
        std::memcpy(csl[ind], tmp, sizeof(tmp));
        double ct = c[2]; c[2] = -c[ind]; c[ind] = ct;
    }

    double newRow2[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            newRow2[j] += c[i] * csl[i][j];
    std::memcpy(csl[2], newRow2, sizeof(newRow2));

    if (c[2] < 0)
        for (int j = 0; j < 3; j++) csl[1][j] *= -1;

    double u1[3];
    double n1 = norm3(csl[2]);
    if (n1 < 1e-12) return transpose3(csl);
    for (int i = 0; i < 3; i++) u1[i] = csl[2][i] / n1;

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

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            csl[i][j] = std::round(csl[i][j]);

    return transpose3(csl);
}

M3 getCslMatrix(int sigma, const M3& rotateMatrix)
{
    M3 rInv = inv3(rotateMatrix);

    M3 t;
    for (int ui = 0; ui < 5; ui++)
    {
        M3 u = fromInt3x3(UNIMODULAR[ui]);
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

// ── Sigma / GB information ──────────────────────────────────────

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
    int axis[3] = {axisIn[0], axisIn[1], axisIn[2]};
    reduceIntVec(axis, 3);

    int axisNorm2 = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    if (axisNorm2 <= 0) return std::vector<SigmaCandidate>();

    float axisNorm = std::sqrt((float)axisNorm2);

    // Check if this is a principal axis: sorted abs values are {0,0,1}
    int sortedAbs[3] = {std::abs(axis[0]), std::abs(axis[1]), std::abs(axis[2])};
    std::sort(sortedAbs, sortedAbs + 3);
    bool isPrincipalAxis = (sortedAbs[0] == 0 && sortedAbs[1] == 0 && sortedAbs[2] == 1);

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

    // Helper: check if a CSL row is a positive standard basis vector
    auto isStdBasisRow = [](const int v[3]) -> bool {
        return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] == 1) &&
               (v[0] + v[1] + v[2] == 1);
    };

    // Helper: fill SigmaCandidate from CSL int matrix
    auto fillCandidate = [](SigmaCandidate& c, int sigma, int m, int n, float theta,
                            const int csl[3][3]) {
        c.sigma = sigma;
        c.m = m;
        c.n = n;
        c.thetaDeg = theta;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                c.csl[i][j] = csl[i][j];
        for (int col = 0; col < 3; col++)
        {
            c.plane[col][0] = csl[0][col];
            c.plane[col][1] = csl[1][col];
            c.plane[col][2] = csl[2][col];
        }
    };

    std::vector<SigmaCandidate> result;
    for (std::map<int, std::vector<TMN>>::iterator it = sigmaTMN.begin(); it != sigmaTMN.end(); ++it)
    {
        int sigma = it->first;
        std::vector<TMN>& entries = it->second;
        std::sort(entries.begin(), entries.end(), [](const TMN& a, const TMN& b) { return a.theta < b.theta; });

        double minTheta = entries[0].theta;

        M3 rotMatrix = getRotateMatrix(axis, minTheta);
        M3 cslRaw = getCslMatrix(sigma, rotMatrix);
        M3 cslOrtho = orthogonalizeCsl(cslRaw, axis);
        M3 cslReduced = reduceCsl(cslOrtho);

        int cslInt[3][3];
        roundToInt3x3(cslReduced, cslInt);

        // For principal axes, derive second CSL by modifying m,n in rows
        // (matches aimsgb ext_csl algorithm)
        if (isPrincipalAxis && (int)entries.size() > 1)
        {
            int ind = -1;
            int extCsl[3][3];

            for (int ei = 0; ei < (int)entries.size(); ei++)
            {
                int mVal = entries[ei].m;
                int nVal = entries[ei].n;

                int tryCsl[3][3];
                std::memcpy(tryCsl, cslInt, sizeof(tryCsl));

                for (int row = 0; row < 3; row++)
                {
                    if (isStdBasisRow(cslInt[row]))
                        continue;

                    int a = cslInt[row][0], b = cslInt[row][1];
                    int sa = (a >= 0) ? 1 : -1;
                    int sb = (b >= 0) ? 1 : -1;
                    if (std::abs(a) > std::abs(b))
                    {
                        tryCsl[row][0] = sa * mVal;
                        tryCsl[row][1] = sb * nVal;
                        tryCsl[row][2] = 0;
                    }
                    else
                    {
                        tryCsl[row][0] = sa * nVal;
                        tryCsl[row][1] = sb * mVal;
                        tryCsl[row][2] = 0;
                    }
                }

                if (std::memcmp(tryCsl, cslInt, sizeof(cslInt)) != 0)
                {
                    ind = ei;
                    std::memcpy(extCsl, tryCsl, sizeof(extCsl));
                    break;
                }
            }

            if (ind >= 0)
            {
                float theta2 = 90.0f - (float)minTheta;
                SigmaCandidate c1, c2;
                if (ind > 0)
                {
                    fillCandidate(c1, sigma, entries[0].m, entries[0].n, (float)minTheta, cslInt);
                    fillCandidate(c2, sigma, entries[ind].m, entries[ind].n, theta2, extCsl);
                }
                else
                {
                    fillCandidate(c1, sigma, entries[0].m, entries[0].n, theta2, cslInt);
                    fillCandidate(c2, sigma, entries[ind].m, entries[ind].n, (float)minTheta, extCsl);
                }
                result.push_back(c1);
                result.push_back(c2);
            }
            else
            {
                SigmaCandidate cand;
                fillCandidate(cand, sigma, entries[0].m, entries[0].n, (float)minTheta, cslInt);
                result.push_back(cand);
            }
        }
        else
        {
            SigmaCandidate cand;
            fillCandidate(cand, sigma, entries[0].m, entries[0].n, (float)minTheta, cslInt);
            result.push_back(cand);
        }
    }
    return result;
}

// ── Grain operations ────────────────────────────────────────────

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

bool reduceToPrimitive(Structure& s, double symprec)
{
#ifdef ATOMS_ENABLE_SPGLIB
    if (!s.hasUnitCell || s.atoms.empty()) return false;

    int natom = (int)s.atoms.size();

    double lattice[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            lattice[i][j] = s.cellVectors[i][j];

    M3 cellM;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cellM[i][j] = lattice[i][j];
    M3 cellInv = inv3(cellM);

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

    std::vector<int> types(natom);
    for (int ai = 0; ai < natom; ai++)
        types[ai] = s.atoms[ai].atomicNumber;

    int num_prim = spg_standardize_cell(lattice,
                                        reinterpret_cast<double(*)[3]>(positions.data()),
                                        types.data(),
                                        natom,
                                        1,
                                        0,
                                        symprec);
    if (num_prim <= 0 || num_prim > natom) return false;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            s.cellVectors[i][j] = lattice[i][j];

    struct ElemInfo { std::string sym; float r, g, b; };
    std::map<int, ElemInfo> zInfo;
    for (const auto& a : s.atoms)
        zInfo[a.atomicNumber] = {a.symbol, a.r, a.g, a.b};

    std::vector<AtomSite> newAtoms(num_prim);
    for (int ai = 0; ai < num_prim; ai++)
    {
        double frac[3] = {positions[ai*3], positions[ai*3+1], positions[ai*3+2]};
        double cart[3] = {0, 0, 0};
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
                cart[j] += frac[i] * lattice[i][j];
        newAtoms[ai].x = cart[0];
        newAtoms[ai].y = cart[1];
        newAtoms[ai].z = cart[2];
        newAtoms[ai].atomicNumber = types[ai];
        if (auto it = zInfo.find(types[ai]); it != zInfo.end())
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

// ── GB-aware primitive reduction ────────────────────────────────
// Finds in-plane translational symmetry and reduces the cell while
// keeping the stacking direction (GB normal) fixed, so the grain
// boundary stays at frac 0 / 0.5 along that axis.

bool reduceToPrimitiveGB(Structure& s, int stackDir, double tol)
{
    if (!s.hasUnitCell || s.atoms.empty() || s.atoms.size() < 2)
        return false;

    int natom = (int)s.atoms.size();

    double cell[3][3], cellInv[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cell[i][j] = s.cellVectors[i][j];
    invertCell(cell, cellInv);

    // Fractional coordinates
    std::vector<std::array<double,3>> fc(natom);
    for (int i = 0; i < natom; i++)
    {
        double cart[3] = {s.atoms[i].x, s.atoms[i].y, s.atoms[i].z};
        cartToFrac(cart, cellInv, fc[i].data());
        for (int d = 0; d < 3; d++)
            fc[i][d] = wrapFrac(fc[i][d]);
    }

    double tolSq = tol * tol;

    // ---------- helper: does translation t map ALL atoms? ----------
    auto mapsAll = [&](const double t[3]) -> bool {
        for (int i = 0; i < natom; i++)
        {
            double sf[3] = { wrapFrac(fc[i][0] + t[0]),
                             wrapFrac(fc[i][1] + t[1]),
                             wrapFrac(fc[i][2] + t[2]) };
            bool found = false;
            for (int j = 0; j < natom && !found; j++)
            {
                if (s.atoms[j].atomicNumber != s.atoms[i].atomicNumber) continue;
                double diff[3];
                for (int d = 0; d < 3; d++)
                {
                    diff[d] = sf[d] - fc[j][d];
                    diff[d] -= std::round(diff[d]);
                }
                double dCart[3] = {0,0,0};
                for (int d = 0; d < 3; d++)
                    for (int k = 0; k < 3; k++)
                        dCart[k] += diff[d] * cell[d][k];
                if (dCart[0]*dCart[0]+dCart[1]*dCart[1]+dCart[2]*dCart[2] < tolSq)
                    found = true;
            }
            if (!found) return false;
        }
        return true;
    };

    // ---------- find reference species (most common) ----------
    std::map<int,int> specCount;
    for (auto& a : s.atoms) specCount[a.atomicNumber]++;
    int refZ = s.atoms[0].atomicNumber;
    int maxC = 0;
    for (auto& [z, count] : specCount)
        if (count > maxC) { maxC = count; refZ = z; }
    int refIdx = 0;
    for (int i = 0; i < natom; i++)
        if (s.atoms[i].atomicNumber == refZ) { refIdx = i; break; }

    // ---------- collect valid in-plane translations ----------
    std::vector<std::array<double,3>> validT;
    for (int j = 0; j < natom; j++)
    {
        if (j == refIdx || s.atoms[j].atomicNumber != refZ) continue;
        double t[3];
        for (int d = 0; d < 3; d++)
        {
            t[d] = fc[j][d] - fc[refIdx][d];
            if (t[d] >  0.5) t[d] -= 1.0;
            if (t[d] < -0.5) t[d] += 1.0;
        }
        // Must be in-plane (no stacking component)
        if (std::abs(t[stackDir]) > 0.01) continue;
        t[stackDir] = 0.0;
        // Skip zero-length
        double ct[3] = {0,0,0};
        for (int d = 0; d < 3; d++)
            for (int k = 0; k < 3; k++)
                ct[k] += t[d] * cell[d][k];
        if (ct[0]*ct[0]+ct[1]*ct[1]+ct[2]*ct[2] < tolSq) continue;

        if (mapsAll(t))
            validT.push_back({t[0], t[1], t[2]});
    }

    if (validT.empty()) return false; // already primitive

    // ---------- build candidate Cartesian vectors ----------
    int d1 = (stackDir + 1) % 3;
    int d2 = (stackDir + 2) % 3;

    struct V3 { double v[3]; };
    std::vector<V3> cand;

    auto addFrac = [&](double f0, double f1, double f2) {
        double fv[3] = {f0, f1, f2};
        V3 c = {{0,0,0}};
        for (int d = 0; d < 3; d++)
            for (int k = 0; k < 3; k++)
                c.v[k] += fv[d] * cell[d][k];
        double len = c.v[0]*c.v[0]+c.v[1]*c.v[1]+c.v[2]*c.v[2];
        if (len > tolSq) cand.push_back(c);
    };

    // Original in-plane cell vectors
    double e1[3]={0,0,0}, e2[3]={0,0,0};
    e1[d1] = 1.0;
    e2[d2] = 1.0;
    addFrac(e1[0], e1[1], e1[2]);
    addFrac(e2[0], e2[1], e2[2]);
    addFrac(-e1[0], -e1[1], -e1[2]);
    addFrac(-e2[0], -e2[1], -e2[2]);

    // Valid translations, their negatives, and combinations with cell vectors
    for (auto& t : validT)
    {
        addFrac( t[0],  t[1],  t[2]);
        addFrac(-t[0], -t[1], -t[2]);
        addFrac(e1[0]+t[0], e1[1]+t[1], e1[2]+t[2]);
        addFrac(e1[0]-t[0], e1[1]-t[1], e1[2]-t[2]);
        addFrac(e2[0]+t[0], e2[1]+t[1], e2[2]+t[2]);
        addFrac(e2[0]-t[0], e2[1]-t[1], e2[2]-t[2]);
    }
    // Combinations among translations
    for (int i = 0; i < (int)validT.size(); i++)
        for (int j = i+1; j < (int)validT.size(); j++)
        {
            addFrac(validT[i][0]+validT[j][0],
                    validT[i][1]+validT[j][1],
                    validT[i][2]+validT[j][2]);
            addFrac(validT[i][0]-validT[j][0],
                    validT[i][1]-validT[j][1],
                    validT[i][2]-validT[j][2]);
        }

    // ---------- pick two shortest independent vectors ----------
    auto vlen = [](const V3& a) {
        return std::sqrt(a.v[0]*a.v[0]+a.v[1]*a.v[1]+a.v[2]*a.v[2]);
    };
    std::sort(cand.begin(), cand.end(),
              [&](const V3& a, const V3& b){ return vlen(a) < vlen(b); });

    V3 v1 = cand[0];
    V3 v2 = {{0,0,0}};
    bool gotV2 = false;
    for (int i = 1; i < (int)cand.size() && !gotV2; i++)
    {
        double cx = v1.v[1]*cand[i].v[2] - v1.v[2]*cand[i].v[1];
        double cy = v1.v[2]*cand[i].v[0] - v1.v[0]*cand[i].v[2];
        double cz = v1.v[0]*cand[i].v[1] - v1.v[1]*cand[i].v[0];
        double cl  = std::sqrt(cx*cx+cy*cy+cz*cz);
        if (cl > tol * vlen(v1))
        {
            v2 = cand[i];
            gotV2 = true;
        }
    }
    if (!gotV2) return false;

    // Ensure right-handed system relative to stacking direction
    {
        double cx = v1.v[1]*v2.v[2] - v1.v[2]*v2.v[1];
        double cy = v1.v[2]*v2.v[0] - v1.v[0]*v2.v[2];
        double cz = v1.v[0]*v2.v[1] - v1.v[1]*v2.v[0];
        double dotStack = cx*cell[stackDir][0]+cy*cell[stackDir][1]+cz*cell[stackDir][2];
        double origCx = cell[d1][1]*cell[d2][2] - cell[d1][2]*cell[d2][1];
        double origCy = cell[d1][2]*cell[d2][0] - cell[d1][0]*cell[d2][2];
        double origCz = cell[d1][0]*cell[d2][1] - cell[d1][1]*cell[d2][0];
        double origDot = origCx*cell[stackDir][0]+origCy*cell[stackDir][1]+origCz*cell[stackDir][2];
        if (dotStack * origDot < 0)
            for (int j = 0; j < 3; j++) v2.v[j] = -v2.v[j];
    }

    // ---------- build new cell ----------
    double newCell[3][3];
    for (int j = 0; j < 3; j++)
    {
        newCell[d1][j]       = v1.v[j];
        newCell[d2][j]       = v2.v[j];
        newCell[stackDir][j] = cell[stackDir][j];
    }

    // Check volume actually decreased
    M3 oldM, newM;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            oldM[i][j] = cell[i][j];
            newM[i][j] = newCell[i][j];
        }
    double oldVol = std::abs(det3(oldM));
    double newVol = std::abs(det3(newM));
    if (newVol >= oldVol - 1e-4) return false;

    // ---------- re-map atoms into smaller cell ----------
    double ncInv[3][3];
    invertCell(newCell, ncInv);

    std::vector<AtomSite> newAtoms;
    newAtoms.reserve(natom);
    for (int i = 0; i < natom; i++)
    {
        double cart[3] = {s.atoms[i].x, s.atoms[i].y, s.atoms[i].z};
        double f[3];
        cartToFrac(cart, ncInv, f);
        for (int d = 0; d < 3; d++) f[d] = wrapFrac(f[d]);
        double nc[3];
        fracToCart(f, newCell, nc);

        // Check for duplicate (PBC-aware)
        bool dup = false;
        for (auto& a : newAtoms)
        {
            double df[3];
            double ac[3] = {a.x, a.y, a.z};
            cartToFrac(ac, ncInv, df);
            for (int d = 0; d < 3; d++) df[d] = wrapFrac(df[d]);
            double dd[3];
            for (int d = 0; d < 3; d++)
            {
                dd[d] = f[d] - df[d];
                dd[d] -= std::round(dd[d]);
            }
            double dCart[3] = {0,0,0};
            for (int d = 0; d < 3; d++)
                for (int k = 0; k < 3; k++)
                    dCart[k] += dd[d] * newCell[d][k];
            if (dCart[0]*dCart[0]+dCart[1]*dCart[1]+dCart[2]*dCart[2] < tolSq)
            {
                dup = true;
                break;
            }
        }
        if (!dup)
        {
            AtomSite a = s.atoms[i];
            a.x = nc[0]; a.y = nc[1]; a.z = nc[2];
            newAtoms.push_back(a);
        }
    }

    // Update structure
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            s.cellVectors[i][j] = newCell[i][j];
    s.atoms = newAtoms;
    return true;
}

double cellVecLen(const double v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void cartToFrac(const double cart[3], const double cellInv[3][3], double frac[3])
{
    for (int j = 0; j < 3; j++)
    {
        frac[j] = 0;
        for (int i = 0; i < 3; i++)
            frac[j] += cart[i] * cellInv[i][j];
    }
}

void fracToCart(const double frac[3], const double cell[3][3], double cart[3])
{
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

void canonicalizeGrain(Grain& g)
{
    double oldCellInv[3][3];
    invertCell(g.cell, oldCellInv);

    double a, b, c, alpha, beta, gamma;
    cellToParameters(g.cell, a, b, c, alpha, beta, gamma);

    double newCell[3][3];
    parametersToCell(a, b, c, alpha, beta, gamma, newCell);

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

Grain makeSupercell(const Grain& grain, const int M[3][3])
{
    Grain result;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            result.cell[i][j] = 0;
            for (int k = 0; k < 3; k++)
                result.cell[i][j] += (double)M[i][k] * grain.cell[k][j];
        }

    M3 Md = fromInt3x3(M);
    M3 Minv = inv3(Md);

    double cellInv[3][3];
    invertCell(grain.cell, cellInv);

    int searchRange = 0;
    for (int i = 0; i < 3; i++)
    {
        int rowSum = 0;
        for (int j = 0; j < 3; j++)
            rowSum += std::abs(M[i][j]);
        if (rowSum > searchRange) searchRange = rowSum;
    }
    searchRange += 1;

    // Parallelize per-atom: each atom's search loop is independent.
    const int numAtoms = (int)grain.atoms.size();
    std::vector<std::vector<AtomSite>> atomCandidates((size_t)numAtoms);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int ai = 0; ai < numAtoms; ai++)
    {
        const AtomSite& atom = grain.atoms[(size_t)ai];
        double cart[3] = {atom.x, atom.y, atom.z};
        double oldFrac[3];
        cartToFrac(cart, cellInv, oldFrac);

        std::vector<AtomSite>& cands = atomCandidates[(size_t)ai];
        for (int ti = -searchRange; ti <= searchRange; ti++)
        for (int tj = -searchRange; tj <= searchRange; tj++)
        for (int tk = -searchRange; tk <= searchRange; tk++)
        {
            double shifted[3] = {oldFrac[0] + ti, oldFrac[1] + tj, oldFrac[2] + tk};
            double newFrac[3];
            for (int j2 = 0; j2 < 3; j2++)
            {
                newFrac[j2] = 0;
                for (int k2 = 0; k2 < 3; k2++)
                    newFrac[j2] += shifted[k2] * Minv[k2][j2];
            }

            bool inside = true;
            for (int c2 = 0; c2 < 3; c2++)
            {
                if (newFrac[c2] < -1e-8 || newFrac[c2] >= 1.0 - 1e-8)
                    inside = false;
            }
            if (!inside) continue;

            for (int c2 = 0; c2 < 3; c2++)
                newFrac[c2] = wrapFrac(newFrac[c2]);

            double newCart[3];
            fracToCart(newFrac, result.cell, newCart);

            AtomSite newAtom = atom;
            newAtom.x = newCart[0];
            newAtom.y = newCart[1];
            newAtom.z = newCart[2];
            cands.push_back(newAtom);
        }
    }
    for (auto& cands : atomCandidates)
        result.atoms.insert(result.atoms.end(), cands.begin(), cands.end());

    canonicalizeGrain(result);

    return result;
}

Grain makeSupercellDiag(const Grain& grain, int na, int nb, int nc)
{
    int M[3][3] = {{na,0,0},{0,nb,0},{0,0,nc}};
    return makeSupercell(grain, M);
}

bool isOrthogonal(const double cell[3][3], double tol)
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

Grain setOrthogonalGrain(const Grain& g, int direction)
{
    Grain result;
    std::memcpy(result.cell, g.cell, sizeof(result.cell));

    int i1 = (direction + 1) % 3;
    int i2 = (direction + 2) % 3;
    double* v1 = result.cell[i1];
    double* v2 = result.cell[i2];

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

    double oldInv[3][3], newInv[3][3];
    invertCell(g.cell, oldInv);
    invertCell(result.cell, newInv);

    result.atoms.reserve(g.atoms.size());
    for (size_t ai = 0; ai < g.atoms.size(); ai++)
    {
        double cart[3] = {g.atoms[ai].x, g.atoms[ai].y, g.atoms[ai].z};
        double frac[3];
        cartToFrac(cart, newInv, frac);
        for (int c2 = 0; c2 < 3; c2++) frac[c2] = wrapFrac(frac[c2]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);
        AtomSite a = g.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }
    return result;
}

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
            for (int c2 = 0; c2 < 3; c2++)
                if (c2 != tryAxis)
                    cart[c2] = -cart[c2];

            double frac[3];
            cartToFrac(cart, cellInv, frac);
            for (int c2 = 0; c2 < 3; c2++) frac[c2] = wrapFrac(frac[c2]);
            double newCart[3];
            fracToCart(frac, grainA.cell, newCart);

            double origCart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
            double origFrac[3];
            cartToFrac(origCart, cellInv, origFrac);
            for (int c2 = 0; c2 < 3; c2++) origFrac[c2] = wrapFrac(origFrac[c2]);

            double df = 0;
            for (int c2 = 0; c2 < 3; c2++)
            {
                double dd = frac[c2] - origFrac[c2];
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
        for (int c2 = 0; c2 < 3; c2++) frac[c2] = wrapFrac(frac[c2]);
        double newCart[3];
        fracToCart(frac, grainA.cell, newCart);
        AtomSite a = grainA.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        b.atoms.push_back(a);
    }
    return b;
}

Grain stackGrains(const Grain& grainA, const Grain& grainB, int dir, double vacuum, double gap)
{
    Grain result;

    double aA, bA, cA, alphaA, betaA, gammaA;
    cellToParameters(grainA.cell, aA, bA, cA, alphaA, betaA, gammaA);

    double aB, bB, cB, alphaB, betaB, gammaB;
    cellToParameters(grainB.cell, aB, bB, cB, alphaB, betaB, gammaB);

    double abcA[3] = {aA, bA, cA};
    double abcB[3] = {aB, bB, cB};
    double angles[3] = {alphaB, betaB, gammaB};

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

    double abcNew[3] = {abcA[0], abcA[1], abcA[2]};
    abcNew[dir] += abcB[dir] + 2.0 * gap + vacuum;

    parametersToCell(abcNew[0], abcNew[1], abcNew[2],
                     angles[0], angles[1], angles[2],
                     result.cell);

    double newCellInv[3][3];
    invertCell(result.cell, newCellInv);

    // Place atoms directly — no centering shift.
    // The natural stacking already puts grainA at frac [0, ~0.5)
    // and grainB at [~0.5, ~1.0), matching the aimsgb convention
    // where the GB plane sits at frac 0.5 and the periodic boundary at 0.0.
    double lVec[3] = {0, 0, 0};
    lVec[dir] = l;

    for (size_t ai = 0; ai < grainA.atoms.size(); ai++)
    {
        double cart[3] = {grainA.atoms[ai].x, grainA.atoms[ai].y, grainA.atoms[ai].z};
        double frac[3];
        cartToFrac(cart, newCellInv, frac);
        for (int c2 = 0; c2 < 3; c2++) frac[c2] = wrapFrac(frac[c2]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);

        AtomSite a = grainA.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }

    for (size_t ai = 0; ai < grainB.atoms.size(); ai++)
    {
        double cart[3] = {grainB.atoms[ai].x + lVec[0],
                          grainB.atoms[ai].y + lVec[1],
                          grainB.atoms[ai].z + lVec[2]};
        double frac[3];
        cartToFrac(cart, newCellInv, frac);
        for (int c2 = 0; c2 < 3; c2++) frac[c2] = wrapFrac(frac[c2]);
        double newCart[3];
        fracToCart(frac, result.cell, newCart);

        AtomSite a = grainB.atoms[ai];
        a.x = newCart[0]; a.y = newCart[1]; a.z = newCart[2];
        result.atoms.push_back(a);
    }

    return result;
}

int removeOverlaps(std::vector<AtomSite>& atoms, double distThreshold)
{
    if (distThreshold <= 0 || atoms.size() < 2) return 0;

    double distSq = distThreshold * distThreshold;
    std::vector<bool> keep(atoms.size(), true);
    int removed = 0;

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
