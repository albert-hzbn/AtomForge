#include "algorithms/InterfaceBuilder.h"

#include "graphics/StructureInstanceBuilder.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_set>

// -- 2x2 matrix utilities ----------------------------------------------------

int mat2DetInt(const int m[2][2])
{
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

void applyMat(const int n[2][2], const double basis[2][2], double out[2][2])
{
    out[0][0] = n[0][0]*basis[0][0] + n[0][1]*basis[1][0];
    out[0][1] = n[0][0]*basis[0][1] + n[0][1]*basis[1][1];
    out[1][0] = n[1][0]*basis[0][0] + n[1][1]*basis[1][0];
    out[1][1] = n[1][0]*basis[0][1] + n[1][1]*basis[1][1];
}

Mat2 rotation2D(double theta)
{
    double c = std::cos(theta), s = std::sin(theta);
    Mat2 r;
    r.m[0][0] = c; r.m[0][1] = -s;
    r.m[1][0] = s; r.m[1][1] = c;
    return r;
}

double angleOf(double x, double y)
{
    return std::atan2(y, x);
}

double wrapDegPm180(double deg)
{
    return std::fmod(deg + 540.0, 360.0) - 180.0;
}

double vecLen(double x, double y)
{
    return std::sqrt(x*x + y*y);
}

LatticeKey equivalentLatticeKey(const double v[2][2])
{
    double l1 = vecLen(v[0][0], v[0][1]);
    double l2 = vecLen(v[1][0], v[1][1]);
    double area = std::abs(v[0][0]*v[1][1] - v[0][1]*v[1][0]);
    if (l1 > l2) std::swap(l1, l2);
    LatticeKey k;
    k.l1 = std::round(l1 * 1e6) / 1e6;
    k.l2 = std::round(l2 * 1e6) / 1e6;
    k.area = std::round(area * 1e6) / 1e6;
    return k;
}

// -- Supercell generation and strain matching ---------------------------------

std::vector<SupercellEntry> generateUniqueSupercells(
    const double basis[2][2], int nmax, int maxCells)
{
    std::vector<SupercellEntry> result;
    std::unordered_set<LatticeKey, LatticeKeyHash> seen;

    for (int n11 = -nmax; n11 <= nmax; ++n11)
    for (int n12 = -nmax; n12 <= nmax; ++n12)
    for (int n21 = -nmax; n21 <= nmax; ++n21)
    for (int n22 = -nmax; n22 <= nmax; ++n22)
    {
        int m[2][2] = {{n11, n12}, {n21, n22}};
        int det = mat2DetInt(m);
        if (det <= 0 || det > maxCells)
            continue;
        double v[2][2];
        applyMat(m, basis, v);
        double detV = v[0][0]*v[1][1] - v[0][1]*v[1][0];
        if (std::abs(detV) < 1e-10)
            continue;
        LatticeKey key = equivalentLatticeKey(v);
        if (seen.count(key))
            continue;
        seen.insert(key);
        SupercellEntry e;
        std::memcpy(e.mat, m, sizeof(m));
        std::memcpy(e.vecs, v, sizeof(v));
        e.det = det;
        result.push_back(e);
    }
    return result;
}

bool strainComponents(const double v[2][2], const double u[2][2],
                      double& exx, double& eyy, double& exy)
{
    double alpha = angleOf(v[0][0], v[0][1]);
    Mat2 rLocal = rotation2D(-alpha);

    double vt[2][2], ut[2][2];
    for (int i = 0; i < 2; ++i)
    {
        vt[i][0] = v[i][0]*rLocal.m[0][0] + v[i][1]*rLocal.m[1][0];
        vt[i][1] = v[i][0]*rLocal.m[0][1] + v[i][1]*rLocal.m[1][1];
        ut[i][0] = u[i][0]*rLocal.m[0][0] + u[i][1]*rLocal.m[1][0];
        ut[i][1] = u[i][0]*rLocal.m[0][1] + u[i][1]*rLocal.m[1][1];
    }

    if (std::abs(ut[0][0]) < 1e-12 || std::abs(ut[1][1]) < 1e-12)
        return false;

    exx = vt[0][0] / ut[0][0] - 1.0;
    eyy = vt[1][1] / ut[1][1] - 1.0;
    exy = 0.5 * (vt[1][0] - (vt[0][0] / ut[0][0]) * ut[1][0]) / ut[1][1];
    return true;
}

double meanAbsStrain(double exx, double eyy, double exy)
{
    return (std::abs(exx) + std::abs(eyy) + std::abs(exy)) / 3.0;
}

double cubicElasticDensity(double exx, double eyy, double exy, const float c[6][6])
{
    double c11 = c[0][0];
    double c12 = c[0][1];
    double c44 = (c[3][3] + c[4][4] + c[5][5]) / 3.0;
    return exx*exx*c11 + exx*eyy*c12 + 0.5*exy*exy*c44;
}

// -- Structure manipulation ---------------------------------------------------

void get2DBasis(const Structure& s, double basis[2][2])
{
    basis[0][0] = s.cellVectors[0][0];
    basis[0][1] = s.cellVectors[0][1];
    basis[1][0] = s.cellVectors[1][0];
    basis[1][1] = s.cellVectors[1][1];
}

void buildMat3From2x2(const int m2[2][2], int m3[3][3])
{
    m3[0][0] = m2[0][0]; m3[0][1] = m2[0][1]; m3[0][2] = 0;
    m3[1][0] = m2[1][0]; m3[1][1] = m2[1][1]; m3[1][2] = 0;
    m3[2][0] = 0;        m3[2][1] = 0;        m3[2][2] = 1;
}

Structure makeSupercell2D(const Structure& base, const int mat2[2][2])
{
    int mat3[3][3];
    buildMat3From2x2(mat2, mat3);
    return buildSupercell(base, mat3);
}

Structure repeatLayersZ(const Structure& base, int layers)
{
    if (layers <= 1) return base;

    Structure result = base;
    double cz[3] = {base.cellVectors[2][0], base.cellVectors[2][1], base.cellVectors[2][2]};

    std::vector<AtomSite> origAtoms = base.atoms;

    // Each layer is independent: build in parallel, then merge in order.
    std::vector<std::vector<AtomSite>> layerAtoms((size_t)(layers - 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int l = 1; l < layers; ++l)
    {
        std::vector<AtomSite>& la = layerAtoms[(size_t)(l - 1)];
        la.resize(origAtoms.size());
        for (size_t k = 0; k < origAtoms.size(); ++k)
        {
            la[k] = origAtoms[k];
            la[k].x += l * cz[0];
            la[k].y += l * cz[1];
            la[k].z += l * cz[2];
        }
    }

    result.atoms.reserve((size_t)layers * origAtoms.size());
    for (auto& la : layerAtoms)
        result.atoms.insert(result.atoms.end(), la.begin(), la.end());

    result.cellVectors[2][0] *= layers;
    result.cellVectors[2][1] *= layers;
    result.cellVectors[2][2] *= layers;

    return result;
}

Structure applyTransform2D(const Structure& b, const double F[2][2])
{
    Structure result = b;

    for (auto& atom : result.atoms)
    {
        double ox = atom.x, oy = atom.y;
        atom.x = F[0][0]*ox + F[0][1]*oy;
        atom.y = F[1][0]*ox + F[1][1]*oy;
    }

    if (result.hasUnitCell)
    {
        for (int i = 0; i < 2; ++i)
        {
            double ox = result.cellVectors[i][0];
            double oy = result.cellVectors[i][1];
            result.cellVectors[i][0] = F[0][0]*ox + F[0][1]*oy;
            result.cellVectors[i][1] = F[1][0]*ox + F[1][1]*oy;
        }
    }

    return result;
}

Structure assembleInterface(const Structure& aSuper, const Structure& bStrained,
                            double zGap, double vacuum)
{
    double aTop = -1e30, bBottom = 1e30;
    for (const auto& a : aSuper.atoms)
        aTop = std::max(aTop, a.z);
    for (const auto& a : bStrained.atoms)
        bBottom = std::min(bBottom, a.z);

    double bShift = aTop - bBottom + zGap;

    Structure result;
    result.hasUnitCell = aSuper.hasUnitCell;
    if (result.hasUnitCell)
    {
        result.cellVectors[0] = aSuper.cellVectors[0];
        result.cellVectors[1] = aSuper.cellVectors[1];
    }

    for (const auto& atom : aSuper.atoms)
        result.atoms.push_back(atom);

    for (auto atom : bStrained.atoms)
    {
        atom.z += bShift;
        result.atoms.push_back(atom);
    }

    double zMin = 1e30, zMax = -1e30;
    for (const auto& a : result.atoms)
    {
        zMin = std::min(zMin, a.z);
        zMax = std::max(zMax, a.z);
    }
    double height = (zMax - zMin) + vacuum;
    result.cellVectors[2] = {{0.0, 0.0, height}};

    double shift = -zMin + vacuum * 0.5;
    for (auto& a : result.atoms)
        a.z += shift;

    result.cellOffset = {{0.0, 0.0, 0.0}};
    return result;
}

Structure repeatInterfaceXY(const Structure& iface, int rx, int ry)
{
    if (rx <= 1 && ry <= 1) return iface;

    Structure result;
    result.hasUnitCell = iface.hasUnitCell;
    result.cellVectors = iface.cellVectors;
    result.cellOffset = iface.cellOffset;

    double ax = iface.cellVectors[0][0], ay = iface.cellVectors[0][1];
    double bx = iface.cellVectors[1][0], by = iface.cellVectors[1][1];

    // Each tile is independent: build in parallel, merge in order.
    const int totalTiles = rx * ry;
    std::vector<std::vector<AtomSite>> tileAtoms((size_t)totalTiles);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int t = 0; t < totalTiles; ++t)
    {
        const int ix = t / ry;
        const int iy = t % ry;
        const double dx = ix * ax + iy * bx;
        const double dy = ix * ay + iy * by;
        std::vector<AtomSite>& tile = tileAtoms[(size_t)t];
        tile.resize(iface.atoms.size());
        for (size_t k = 0; k < iface.atoms.size(); ++k)
        {
            tile[k] = iface.atoms[k];
            tile[k].x += dx;
            tile[k].y += dy;
        }
    }

    result.atoms.reserve((size_t)totalTiles * iface.atoms.size());
    for (auto& tile : tileAtoms)
        result.atoms.insert(result.atoms.end(), tile.begin(), tile.end());

    result.cellVectors[0][0] *= rx;
    result.cellVectors[0][1] *= rx;
    result.cellVectors[1][0] *= ry;
    result.cellVectors[1][1] *= ry;

    return result;
}

// -- Orientation relationship -------------------------------------------------

float orientationAngleFromPlaneDir(const Structure& s,
                                   const float hkl[3],
                                   const float uvw[3])
{
    if (!s.hasUnitCell) return 0.0f;

    glm::dvec3 a(s.cellVectors[0][0], s.cellVectors[0][1], s.cellVectors[0][2]);
    glm::dvec3 b(s.cellVectors[1][0], s.cellVectors[1][1], s.cellVectors[1][2]);
    glm::dvec3 c(s.cellVectors[2][0], s.cellVectors[2][1], s.cellVectors[2][2]);

    double vol = glm::dot(a, glm::cross(b, c));
    if (std::abs(vol) < 1e-30) return 0.0f;

    glm::dvec3 ra = glm::cross(b, c) / vol;
    glm::dvec3 rb = glm::cross(c, a) / vol;
    glm::dvec3 rc = glm::cross(a, b) / vol;

    glm::dvec3 normal = (double)hkl[0] * ra + (double)hkl[1] * rb + (double)hkl[2] * rc;
    glm::dvec3 d = (double)uvw[0] * a + (double)uvw[1] * b + (double)uvw[2] * c;

    double nLen = glm::length(normal);
    double dLen = glm::length(d);
    if (nLen < 1e-12 || dLen < 1e-12) return 0.0f;

    glm::dvec3 dInPlane = d - (glm::dot(d, normal) / glm::dot(normal, normal)) * normal;
    if (glm::length(dInPlane) < 1e-12) return 0.0f;

    if (std::hypot(dInPlane.x, dInPlane.y) >= 1e-12)
        return (float)(std::atan2(dInPlane.y, dInPlane.x) * 180.0 / M_PI);

    return 0.0f;
}
