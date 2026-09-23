#include "algorithms/ElasticConstants.h"

#include <cmath>

namespace
{
int voigtIndex(int i, int j)
{
    // Standard Voigt pair index: xx=0, yy=1, zz=2, yz=3, xz=4, xy=5.
    if (i == j) return i;
    if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
    if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
    return 5; // (0,1) or (1,0)
}
}

ElasticTensor makeVoigtElasticTensor(const double voigt[6][6])
{
    ElasticTensor t;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
        const int a = voigtIndex(i, j);
        const int b = voigtIndex(k, l);
        // Voigt matrices are conventionally given as upper-triangular by the
        // caller; symmetrize so either triangle works.
        t.c[i][j][k][l] = (voigt[a][b] != 0.0) ? voigt[a][b] : voigt[b][a];
    }
    return t;
}

ElasticTensor makeCubicElasticTensor(double c11, double c12, double c44)
{
    double v[6][6] = {};
    v[0][0] = v[1][1] = v[2][2] = c11;
    v[0][1] = v[1][0] = v[0][2] = v[2][0] = v[1][2] = v[2][1] = c12;
    v[3][3] = v[4][4] = v[5][5] = c44;
    return makeVoigtElasticTensor(v);
}

ElasticTensor makeHexagonalElasticTensor(double c11, double c12, double c13, double c33, double c44)
{
    double v[6][6] = {};
    v[0][0] = v[1][1] = c11;
    v[2][2] = c33;
    v[0][1] = v[1][0] = c12;
    v[0][2] = v[2][0] = v[1][2] = v[2][1] = c13;
    v[3][3] = v[4][4] = c44;
    v[5][5] = 0.5 * (c11 - c12); // C66, hexagonal symmetry constraint
    return makeVoigtElasticTensor(v);
}

ElasticTensor rotateElasticTensor(const ElasticTensor& c, const glm::dmat3& axes)
{
    // R[p][i] = component p (crystal frame) of new axis i, i.e. axes[i] are
    // the columns; C'_ijkl = R_pi R_qj R_rk R_sl C_pqrs.
    double R[3][3];
    for (int p = 0; p < 3; ++p)
        for (int i = 0; i < 3; ++i)
            R[p][i] = axes[i][p];

    ElasticTensor out;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
        double sum = 0.0;
        for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
        for (int r = 0; r < 3; ++r)
        for (int s = 0; s < 3; ++s)
            sum += R[p][i] * R[q][j] * R[r][k] * R[s][l] * c.c[p][q][r][s];
        out.c[i][j][k][l] = sum;
    }
    return out;
}

bool isNearlyIsotropic(const ElasticTensor& c, double relativeTolerance)
{
    // Fit the best isotropic tensor C_ijkl = lambda*d_ij*d_kl + mu*(d_ik*d_jl+d_il*d_jk)
    // via the standard invariant projections, then compare residual norm.
    double trace = 0.0; // sum_i C_iiii-ish bulk-type contraction
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            trace += c.c[i][i][j][j];
    const double lambdaPlus2muOver3 = trace / 9.0; // (3*lambda+2*mu)/3 averaged form, see below

    double shearSum = 0.0;
    int shearCount = 0;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        if (i == j) continue;
        shearSum += c.c[i][j][i][j];
        shearCount++;
    }
    const double mu = shearSum / shearCount;
    const double lambda = lambdaPlus2muOver3 - (2.0 / 3.0) * mu;

    double residual = 0.0;
    double norm = 0.0;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
        const double dij = (i == j) ? 1.0 : 0.0;
        const double dkl = (k == l) ? 1.0 : 0.0;
        const double dik = (i == k) ? 1.0 : 0.0;
        const double djl = (j == l) ? 1.0 : 0.0;
        const double dil = (i == l) ? 1.0 : 0.0;
        const double djk = (j == k) ? 1.0 : 0.0;
        const double iso = lambda * dij * dkl + mu * (dik * djl + dil * djk);
        const double diff = c.c[i][j][k][l] - iso;
        residual += diff * diff;
        norm += c.c[i][j][k][l] * c.c[i][j][k][l];
    }

    if (norm <= 1e-12)
        return true;
    return std::sqrt(residual / norm) <= relativeTolerance;
}
