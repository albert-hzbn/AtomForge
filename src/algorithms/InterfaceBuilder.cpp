#include "algorithms/InterfaceBuilder.h"

#include <vector>

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
