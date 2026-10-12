#pragma once
// Structure helpers for the interface builder (algorithms/InterfaceMatching.h).

#include "model/Structure.h"

// The structure stacked `layers` times along its third cell vector.
Structure repeatLayersZ(const Structure& base, int layers);

// Atom positions and the in-plane cell vectors mapped by the 2x2 matrix F
// (rows: x' = F[0][0] x + F[0][1] y, y' = F[1][0] x + F[1][1] y).
Structure applyTransform2D(const Structure& b, const double F[2][2]);

// The interface repeated rx times along its first and ry times along its
// second cell vector.
Structure repeatInterfaceXY(const Structure& iface, int rx, int ry);
