#pragma once

#include "model/Structure.h"

#include <vector>
#include <glm/glm.hpp>

// ---------------------------------------------------------------------------
// CellSlabPlane
// One half-open slab = two parallel planes at [d1, d2] (Å) along the (h,k,l)
// unit normal measured from the structure reference point.
// When usePeriodic=true the bounds are integer multiples of d_hkl:
//   d1 = startPlane * d_hkl
//   d2 = (startPlane + nPeriods) * d_hkl
// ---------------------------------------------------------------------------
struct CellSlabPlane
{
    int   h = 1, k = 0, l = 0;
    float d1          = -5.0f;   // lower bound (Å from ref) – manual mode
    float d2          =  5.0f;   // upper bound (Å from ref) – manual mode
    bool  usePeriodic = true;    // snap to crystal planes (default on)
    int   startPlane  = 0;       // lower plane index (periodic mode)
    int   nPeriods    = 1;       // number of planes to include
    float color[3]    = {0.4f, 0.7f, 1.0f};
};

// ---------------------------------------------------------------------------
// Pure-compute free functions (no ImGui, no OpenGL, no global state).
// `source` is the un-tiled unit-cell structure; it supplies cell vectors for
// reciprocal-lattice calculations.
// ---------------------------------------------------------------------------

// Reference point for a slab's distance measurement.
// Periodic slabs:  phase-aligned origin — shifted from cellOffset so that
//                  d = n * d_hkl coincides with the n-th actual equivalent
//                  atom plane (A→A repeat, not A→B).  Works for any basis
//                  including substituted/doped structures.
// Manual slabs:    supercell centroid (unchanged, user-relative behaviour).
glm::vec3 cscSlabRef(const CellSlabPlane& slab, const Structure& source,
                     const glm::vec3& centroid);

// Tolerance (Å) added to slab bounds so atoms that sit exactly on a boundary
// plane are kept despite floating-point rounding.
inline constexpr float kCscSlabTol = 0.02f;

// Unit normal in Cartesian (Å^-1 space) for the (h,k,l) family.
// Uses the reciprocal lattice B = (A^-T) when a unit cell is present,
// otherwise falls back to the normalised (h,k,l) Cartesian vector.
glm::vec3 cscNormal(const CellSlabPlane& slab, const Structure& source);

// Interplanar spacing d_hkl (Å) from the reciprocal lattice (lattice planes only).
// Returns 1.0 if no unit cell.
float cscDhkl(const CellSlabPlane& slab, const Structure& source);

// Structural period (Å) = N·d_hkl where N is the smallest integer such that
// N·d_hkl·n̂ is a real-space lattice vector.  Both slab faces at multiples of
// this distance pass through the same crystallographic plane type.
//   FCC (111) → 3·d_hkl     BCC (110) → 2·d_hkl     SC (any) → d_hkl
float cscStructuralDhkl(const CellSlabPlane& slab, const Structure& source);

// Actual atom-layer spacing (Å): smallest gap between distinct atom planes
// within one structural period.  For FCC(111): equals d_hkl (= structural_d/3).
// For UI display only.
float cscAtomLayerSpacing(const CellSlabPlane& slab, const Structure& source);

// Effective lower bound (Å from reference).
// Periodic: startPlane * cscStructuralDhkl.  Manual: slab.d1.
float cscD1(const CellSlabPlane& slab, const Structure& source);

// Effective upper bound (Å from reference).
// Periodic: (startPlane+nPeriods)*cscStructuralDhkl.  Manual: slab.d2.
float cscD2(const CellSlabPlane& slab, const Structure& source);

// Centroid of all atoms in s (arithmetic mean of Cartesian coordinates).
glm::vec3 cscCentroid(const Structure& s);

// Build an (nx × ny × nz) supercell with symmetric tiling around the
// source centroid.  Tiling ranges: ia ∈ [-(nx/2), nx-(nx/2)), etc.
// The returned structure has unit-cell vectors scaled by nx/ny/nz.
Structure cscBuildSupercell(const Structure& source, int nx, int ny, int nz);

// Keep only atoms in `sc` that satisfy all enabled slab constraints.
// The reference origin is the centroid of `sc` (the supercell).
Structure cscApplySlabs(const Structure& sc,
                        const std::vector<CellSlabPlane>& slabs,
                        const Structure& source);

// True when the enabled slab normals span ℝ³ (i.e. the intersection is
// a closed, bounded convex region).
bool cscIsBounded(const std::vector<CellSlabPlane>& slabs,
                  const Structure& source);

// Sutherland–Hodgman single-plane clip.
// Returns vertices on the side where dot(p, n) >= plane_d.
std::vector<glm::vec3> cscClipHalfspace(
    const std::vector<glm::vec3>& poly,
    const glm::vec3& n, float plane_d);

// Remove spatially coincident atoms (distance < tol Å), keeping the first
// occurrence of each unique position.  Eliminates duplicate boundary atoms
// produced by periodic tiling (corner/edge/face equivalents).
Structure cscDeduplicateAtoms(const Structure& s, float tol = 0.01f);
