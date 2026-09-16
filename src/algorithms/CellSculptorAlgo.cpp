#include "algorithms/CellSculptorAlgo.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include <glm/gtc/matrix_inverse.hpp>

// ---------------------------------------------------------------------------
// Internal helper: reciprocal-lattice matrix B = (A^-T).
// Columns of A are the lattice vectors a, b, c.
// B * hkl gives the (hkl) plane normal in Cartesian space (un-normalised).
static glm::mat3 recipMatrix(const Structure& source)
{
    const auto& cv = source.cellVectors;
    const glm::mat3 A(
        (float)cv[0][0], (float)cv[0][1], (float)cv[0][2],
        (float)cv[1][0], (float)cv[1][1], (float)cv[1][2],
        (float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
    return glm::transpose(glm::inverse(A));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Find the smallest N ≥ 1 such that translating the entire atom arrangement
// by N·d_hkl·n̂ maps every source atom onto another source atom (mod the
// unit cell).  This is the *structural period multiplier*; nPeriods then
// counts in units of N·d_hkl so that both slab faces always cut through
// crystallographically equivalent planes.
//
// This is atom-arrangement aware, not just lattice-vector aware.  The key
// difference:
//   FCC (111)  → conventional cell check: N=3  (atom check: N=3) ✓
//   FCC (110)  → conventional cell check: N=2  (atom check: N=1) ✓
//                because (a/2,a/2,0) maps FCC atoms to FCC atoms even
//                though it is not an integer combination of (a,0,0)…
//   BCC (110)  → N=1 (already a body-center translation)
//   SC  (any)  → N=1
static int structuralPeriodMultiplier(const CellSlabPlane& slab,
                                      const Structure& source)
{
    if (!source.hasUnitCell || source.atoms.empty()) return 1;

    const glm::vec3 hkl((float)slab.h, (float)slab.k, (float)slab.l);
    const glm::mat3 B    = recipMatrix(source);
    const glm::vec3 Ghkl = B * hkl;
    const float     Gsq  = glm::dot(Ghkl, Ghkl);
    if (Gsq < 1e-18f) return 1;

    // Single-layer step vector in Cartesian: d_hkl · n̂ = G_hkl / |G_hkl|²
    const glm::vec3 step = Ghkl / Gsq;

    // Build A (columns = cell vectors) and its inverse for fractional coords.
    const auto& cv = source.cellVectors;
    const glm::mat3 A(
        (float)cv[0][0], (float)cv[0][1], (float)cv[0][2],
        (float)cv[1][0], (float)cv[1][1], (float)cv[1][2],
        (float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
    const glm::mat3 Ainv = glm::inverse(A);

    // Fractional positions of all source atoms, reduced to [0,1), paired with
    // their element symbol so that species-breaking translations are rejected.
    struct FracAtom { glm::vec3 f; std::string sym; };
    std::vector<FracAtom> fracAtoms;
    fracAtoms.reserve(source.atoms.size());
    for (const auto& a : source.atoms)
    {
        glm::vec3 f = Ainv * glm::vec3((float)a.x, (float)a.y, (float)a.z);
        for (int i = 0; i < 3; ++i)
            f[i] -= std::floor(f[i]);
        fracAtoms.push_back({f, a.symbol});
    }

    // Fractional representation of one step.
    const glm::vec3 fracStep = Ainv * step;
    const float tol = 0.01f;

    for (int N = 1; N <= 24; ++N)
    {
        const glm::vec3 fT = fracStep * (float)N;
        bool isSymmetry = true;

        for (const auto& fa : fracAtoms)
        {
            // Translate and wrap to [0,1).
            glm::vec3 ft;
            for (int i = 0; i < 3; ++i)
            {
                ft[i] = fa.f[i] + fT[i];
                ft[i] -= std::floor(ft[i]);
            }

            // Check that ft maps onto an atom of the SAME species.
            bool found = false;
            for (const auto& g : fracAtoms)
            {
                if (g.sym != fa.sym) continue; // species must match
                glm::vec3 d = ft - g.f;
                for (int i = 0; i < 3; ++i)
                    d[i] -= std::round(d[i]);
                if (std::abs(d.x) < tol && std::abs(d.y) < tol && std::abs(d.z) < tol)
                { found = true; break; }
            }
            if (!found) { isSymmetry = false; break; }
        }

        if (isSymmetry) return N;
    }
    return 1; // fallback: treat any d_hkl as the period
}

// ---------------------------------------------------------------------------
// Collect all unique atom-plane projections (mod d_period) by tiling ±1 in
// every lattice direction.  Using d_struct = N·d_hkl as the modulus exposes
// all crystallographically distinct layer positions within one structural
// period (e.g. three positions for FCC(111): A, B, C at 0, d_hkl, 2·d_hkl).
// Returns a sorted, deduplicated list in [0, d_period).
static std::vector<float> allLayerProjections(const CellSlabPlane& slab,
                                              const Structure& source,
                                              float d_period)
{
    if (d_period < 1e-9f || source.atoms.empty()) return {};

    const glm::vec3 n   = cscNormal(slab, source);
    const float     tol = d_period * 0.005f;
    const auto&     cv  = source.cellVectors;

    std::vector<float> projs;
    projs.reserve(source.atoms.size() * 27);

    for (int ia = -1; ia <= 1; ++ia)
    for (int ib = -1; ib <= 1; ++ib)
    for (int ic = -1; ic <= 1; ++ic)
    {
        const glm::vec3 t(
            (float)(cv[0][0]*ia + cv[1][0]*ib + cv[2][0]*ic),
            (float)(cv[0][1]*ia + cv[1][1]*ib + cv[2][1]*ic),
            (float)(cv[0][2]*ia + cv[1][2]*ib + cv[2][2]*ic));
        for (const auto& a : source.atoms)
        {
            const glm::vec3 p((float)a.x + t.x, (float)a.y + t.y, (float)a.z + t.z);
            float proj = glm::dot(p, n);
            proj = std::fmod(proj, d_period);
            if (proj < 0.0f) proj += d_period;
            projs.push_back(proj);
        }
    }

    std::sort(projs.begin(), projs.end());
    projs.erase(std::unique(projs.begin(), projs.end(),
        [tol](float a, float b){ return b - a < tol; }), projs.end());
    while (!projs.empty() && projs.back() >= d_period - tol)
        projs.pop_back();
    return projs;
}

// Phase of the slab relative to cellOffset, measured mod d_struct.
// Shifting ref back by this amount aligns startPlane=0 with the first
// actual atom plane in the structural repeat.
static float slabPhase(const CellSlabPlane& slab, const Structure& source,
                       float d_struct)
{
    const std::vector<float> projs = allLayerProjections(slab, source, d_struct);
    return projs.empty() ? 0.0f : projs.front();
}

glm::vec3 cscSlabRef(const CellSlabPlane& slab, const Structure& source,
                     const glm::vec3& centroid)
{
    // Manual slabs: centroid-relative (user-intuitive).
    if (!slab.usePeriodic || !source.hasUnitCell)
        return centroid;

    // Periodic slabs: shift ref so that integer multiples of the structural
    // period (N·d_hkl) from `ref` land on atom planes of the same type.
    //   startPlane=0 → first atom plane at or just above the cell origin
    //   startPlane=1 → next crystallographically equivalent plane (N·d_hkl further)
    const glm::vec3 base((float)source.cellOffset[0],
                         (float)source.cellOffset[1],
                         (float)source.cellOffset[2]);
    if (source.atoms.empty()) return base;

    const glm::vec3 n       = cscNormal(slab, source);
    const float     d_struct = (float)structuralPeriodMultiplier(slab, source)
                               * cscDhkl(slab, source);
    const float     phase   = slabPhase(slab, source, d_struct);
    return base - n * phase;
}

// ---------------------------------------------------------------------------
glm::vec3 cscNormal(const CellSlabPlane& slab, const Structure& source)
{
    const glm::vec3 hkl((float)slab.h, (float)slab.k, (float)slab.l);
    const glm::vec3 n = source.hasUnitCell ? recipMatrix(source) * hkl : hkl;
    const float len = glm::length(n);
    return len > 1e-9f ? n / len : glm::vec3(0, 0, 1);
}

// ---------------------------------------------------------------------------
float cscDhkl(const CellSlabPlane& slab, const Structure& source)
{
    if (!source.hasUnitCell) return 1.0f;
    const glm::vec3 hkl((float)slab.h, (float)slab.k, (float)slab.l);
    const float len = glm::length(recipMatrix(source) * hkl);
    return len > 1e-9f ? 1.0f / len : 1.0f;
}

// Actual atom-layer spacing (Å): smallest gap between distinct atom planes
// within one structural period.  E.g. for FCC(111): d_hkl (= d_struct/3).
float cscAtomLayerSpacing(const CellSlabPlane& slab, const Structure& source)
{
    const float d_struct = (float)structuralPeriodMultiplier(slab, source)
                           * cscDhkl(slab, source);
    if (d_struct < 1e-9f) return d_struct;

    const std::vector<float> projs = allLayerProjections(slab, source, d_struct);
    if (projs.size() < 2) return d_struct;

    const float tol = d_struct * 0.005f;
    float minGap = d_struct - projs.back() + projs.front(); // wrap-around gap
    for (size_t i = 1; i < projs.size(); ++i)
        minGap = std::min(minGap, projs[i] - projs[i-1]);
    return (minGap > tol) ? minGap : d_struct;
}

// ---------------------------------------------------------------------------
float cscStructuralDhkl(const CellSlabPlane& slab, const Structure& source)
{
    return (float)structuralPeriodMultiplier(slab, source) * cscDhkl(slab, source);
}

// ---------------------------------------------------------------------------
float cscD1(const CellSlabPlane& slab, const Structure& source)
{
    if (!slab.usePeriodic) return slab.d1;
    // nPeriods counts in units of structural_d = N·d_hkl so that each period
    // spans a full A→A repeat (same crystallographic plane type at both ends).
    return (float)slab.startPlane * cscStructuralDhkl(slab, source);
}

float cscD2(const CellSlabPlane& slab, const Structure& source)
{
    if (!slab.usePeriodic) return slab.d2;
    return (float)(slab.startPlane + slab.nPeriods) * cscStructuralDhkl(slab, source);
}

// ---------------------------------------------------------------------------
glm::vec3 cscCentroid(const Structure& s)
{
    if (s.atoms.empty()) return glm::vec3(0.0f);
    glm::vec3 sum(0.0f);
    for (const auto& a : s.atoms)
        sum += glm::vec3((float)a.x, (float)a.y, (float)a.z);
    return sum / (float)s.atoms.size();
}

// ---------------------------------------------------------------------------
Structure cscBuildSupercell(const Structure& source, int nx, int ny, int nz)
{
    if (source.atoms.empty()) return {};
    if (nx == 1 && ny == 1 && nz == 1) return source;

    Structure sc;
    sc.hasUnitCell = source.hasUnitCell;

    const auto& cv = source.cellVectors;
    for (int i = 0; i < 3; ++i)
    {
        sc.cellVectors[0][i] = cv[0][i] * nx;
        sc.cellVectors[1][i] = cv[1][i] * ny;
        sc.cellVectors[2][i] = cv[2][i] * nz;
    }
    // Symmetric tiling: ia ∈ [-(nx/2), nx-(nx/2))
    const int ia_min = -(nx / 2), ia_max = nx - (nx / 2);
    const int ib_min = -(ny / 2), ib_max = ny - (ny / 2);
    const int ic_min = -(nz / 2), ic_max = nz - (nz / 2);

    // Shift the cell offset to align the box with the actual atom positions.
    sc.cellOffset[0] = source.cellOffset[0]
                      + cv[0][0]*ia_min + cv[1][0]*ib_min + cv[2][0]*ic_min;
    sc.cellOffset[1] = source.cellOffset[1]
                      + cv[0][1]*ia_min + cv[1][1]*ib_min + cv[2][1]*ic_min;
    sc.cellOffset[2] = source.cellOffset[2]
                      + cv[0][2]*ia_min + cv[1][2]*ib_min + cv[2][2]*ic_min;

    // Pre-compute A matrix and its inverse for fractional-coordinate wrapping.
    // Wrapping every source atom to fractional [0,1) before tiling ensures that
    // atoms sitting at face/edge/corner positions (or slightly outside the cell
    // due to editor moves or floating-point drift) end up in the correct tiled
    // cell rather than doubling up or disappearing.
    const glm::mat3 A(
        (float)cv[0][0], (float)cv[0][1], (float)cv[0][2],
        (float)cv[1][0], (float)cv[1][1], (float)cv[1][2],
        (float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
    const glm::mat3 Ainv   = glm::inverse(A);
    const glm::vec3 origin((float)source.cellOffset[0],
                            (float)source.cellOffset[1],
                            (float)source.cellOffset[2]);

    sc.atoms.reserve(source.atoms.size() * (size_t)(nx * ny * nz));
    for (int ia = ia_min; ia < ia_max; ++ia)
    for (int ib = ib_min; ib < ib_max; ++ib)
    for (int ic = ic_min; ic < ic_max; ++ic)
    {
        const double tx = cv[0][0]*ia + cv[1][0]*ib + cv[2][0]*ic;
        const double ty = cv[0][1]*ia + cv[1][1]*ib + cv[2][1]*ic;
        const double tz = cv[0][2]*ia + cv[1][2]*ib + cv[2][2]*ic;
        for (std::size_t index=0;index<source.atoms.size();++index)
        {
            AtomSite a = source.atoms[index];

            // Wrap the source atom to fractional [0,1) relative to cellOffset
            // so that atoms outside the canonical cell (e.g. face atoms at
            // exactly 1.0, or atoms moved by the editor) are placed correctly.
            if (source.hasUnitCell)
            {
                const glm::vec3 cart((float)a.x, (float)a.y, (float)a.z);
                glm::vec3 frac = Ainv * (cart - origin);
                for (int i = 0; i < 3; ++i)
                    frac[i] -= std::floor(frac[i]);
                const glm::vec3 wrapped = A * frac + origin;
                a.x = (double)wrapped.x;
                a.y = (double)wrapped.y;
                a.z = (double)wrapped.z;
            }

            a.x += tx; a.y += ty; a.z += tz;
            sc.atoms.push_back(a);
            if (source.grainColors.size()==source.atoms.size()) sc.grainColors.push_back(source.grainColors[index]);
            if (source.grainRegionIds.size()==source.atoms.size()) sc.grainRegionIds.push_back(source.grainRegionIds[index]);
        }
    }
    return sc;
}

// ---------------------------------------------------------------------------
Structure cscApplySlabs(const Structure& sc,
                        const std::vector<CellSlabPlane>& slabs,
                        const Structure& source)
{
    if (sc.atoms.empty() || slabs.empty()) return sc;

    const glm::vec3 centroid = cscCentroid(sc);

    // Pre-compute per-slab data for the enabled slabs.
    // `periodic` flag controls the upper-bound check:
    //   periodic  → half-open [d1-tol, d2-tol): the atom plane sitting exactly
    //               at d2 belongs to the NEXT period and must be excluded.
    //               (For FCC (111) this is the B layer when the slab spans A→B.)
    //   manual    → closed   [d1-tol, d2+tol]: user may want to include the
    //               upper boundary.
    struct SlabData { glm::vec3 n, ref; float d1, d2; bool periodic; };
    std::vector<SlabData> enabled;
    for (const auto& slab : slabs)
        enabled.push_back({ cscNormal(slab, source),
                             cscSlabRef(slab, source, centroid),
                             cscD1(slab, source),
                             cscD2(slab, source),
                             slab.usePeriodic && source.hasUnitCell });

    Structure result = sc;
    result.atoms.clear();
    result.grainColors.clear(); result.grainRegionIds.clear();
    result.dislocationLoopPoints.clear(); result.dislocationDetectionDone=false;
    result.atoms.reserve(sc.atoms.size());

    for (std::size_t index=0;index<sc.atoms.size();++index)
    {
        const auto& atom=sc.atoms[index];
        const glm::vec3 p((float)atom.x, (float)atom.y, (float)atom.z);
        bool keep = true;
        for (const auto& sd : enabled)
        {
            const float proj = glm::dot(p - sd.ref, sd.n);
            if (proj < sd.d1 - kCscSlabTol) { keep = false; break; }
            // Periodic: half-open upper bound — atoms at exactly d2 (the next
            // period's first layer) are excluded.  The kCscSlabTol margin
            // absorbs float-to-double conversion errors (~1e-4 Å) without
            // accidentally including the boundary layer (min real layer gap
            // > 0.5 Å >> kCscSlabTol for any physical crystal).
            if (sd.periodic  ? (proj >= sd.d2 - kCscSlabTol)
                             : (proj >  sd.d2 + kCscSlabTol))
            { keep = false; break; }
        }
        if (keep) {
            result.atoms.push_back(atom);
            if (sc.grainColors.size()==sc.atoms.size()) result.grainColors.push_back(sc.grainColors[index]);
            if (sc.grainRegionIds.size()==sc.atoms.size()) result.grainRegionIds.push_back(sc.grainRegionIds[index]);
        }
    }

    // Build cell vectors from the cutting planes when exactly 3 slabs are
    // enabled.  The cell is the parallelepiped whose three edge directions
    // are the slab normals and whose corner is the intersection of the three
    // d1 halfspace planes (solved as a 3×3 linear system).
    if (enabled.size() == 3)
    {
        const glm::vec3& n0 = enabled[0].n;
        const glm::vec3& n1 = enabled[1].n;
        const glm::vec3& n2 = enabled[2].n;

        // Coefficient matrix A where each row is one slab normal.
        // glm::mat3(n0,n1,n2) puts them as columns → transpose gives rows.
        const glm::mat3 A = glm::transpose(glm::mat3(n0, n1, n2));
        const float det = glm::determinant(A);

        if (std::abs(det) > 1e-6f)
        {
            // rhs[i] = d1_i + dot(ref_i, n_i)  (plane equation offset)
            const glm::vec3 rhs(
                enabled[0].d1 + glm::dot(enabled[0].ref, n0),
                enabled[1].d1 + glm::dot(enabled[1].ref, n1),
                enabled[2].d1 + glm::dot(enabled[2].ref, n2));

            const glm::vec3 corner = glm::inverse(A) * rhs;

            result.hasUnitCell    = true;
            result.cellOffset[0]  = (double)corner.x;
            result.cellOffset[1]  = (double)corner.y;
            result.cellOffset[2]  = (double)corner.z;

            // Each cell vector = slab normal × slab width (Å).
            for (int j = 0; j < 3; ++j)
            {
                const glm::vec3 cv = glm::inverse(A)[j] * (enabled[j].d2 - enabled[j].d1);
                result.cellVectors[j][0] = (double)cv.x;
                result.cellVectors[j][1] = (double)cv.y;
                result.cellVectors[j][2] = (double)cv.z;
            }
        }
        else
        {
            result.hasUnitCell = false;
        }
    }
    else
    {
        result.hasUnitCell = false;
    }

    // Remove duplicate atoms after final cell assignment.
    // For 3-slab bounded cuts this also removes periodic-equivalent boundary
    // pairs on opposite faces (fractional 0 ≡ 1) by canonical wrapping.
    result = cscDeduplicateAtoms(result);

    return result;
}

// ---------------------------------------------------------------------------
Structure cscDeduplicateAtoms(const Structure& s, float tol)
{
    if (s.atoms.empty()) return s;

    const float tol2 = tol * tol;
    const float inv  = 1.0f / tol;

    // If a valid unit cell exists, canonicalize each position by wrapping to
    // fractional [0,1). This lets us collapse periodic-equivalent boundary
    // duplicates (e.g. frac 0 and frac 1 on opposite faces).
    bool usePeriodicCanonical = false;
    glm::mat3 cell(1.0f), invCell(1.0f);
    glm::vec3 origin(0.0f);
    if (s.hasUnitCell)
    {
        cell = glm::mat3(
            glm::vec3((float)s.cellVectors[0][0], (float)s.cellVectors[0][1], (float)s.cellVectors[0][2]),
            glm::vec3((float)s.cellVectors[1][0], (float)s.cellVectors[1][1], (float)s.cellVectors[1][2]),
            glm::vec3((float)s.cellVectors[2][0], (float)s.cellVectors[2][1], (float)s.cellVectors[2][2]));
        const float det = glm::determinant(cell);
        if (std::abs(det) > 1e-8f)
        {
            invCell = glm::inverse(cell);
            origin = glm::vec3((float)s.cellOffset[0],
                               (float)s.cellOffset[1],
                               (float)s.cellOffset[2]);
            usePeriodicCanonical = true;
        }
    }

    // Spatial grid: each cell covers (tol × tol × tol) Å³.
    // For each candidate atom we only check the 3³ = 27 neighbour cells.
    using Key = std::tuple<int,int,int>;
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<int>{}(std::get<0>(k));
            h ^= std::hash<int>{}(std::get<1>(k)) + 0x9e3779b9u + (h<<6) + (h>>2);
            h ^= std::hash<int>{}(std::get<2>(k)) + 0x9e3779b9u + (h<<6) + (h>>2);
            return h;
        }
    };
    // grid cell → indices into result.atoms
    std::unordered_map<Key, std::vector<int>, KeyHash> grid;
    grid.reserve(s.atoms.size());

    Structure result = s;
    result.atoms.clear();
    result.grainColors.clear(); result.grainRegionIds.clear();
    result.dislocationLoopPoints.clear(); result.dislocationDetectionDone=false;
    result.atoms.reserve(s.atoms.size());

    std::vector<glm::vec3> canonicalPositions;
    canonicalPositions.reserve(s.atoms.size());

    for (std::size_t index=0;index<s.atoms.size();++index)
    {
        const auto& a=s.atoms[index];
        glm::vec3 canonical((float)a.x, (float)a.y, (float)a.z);
        if (usePeriodicCanonical)
        {
            glm::vec3 frac = invCell * (canonical - origin);
            for (int i = 0; i < 3; ++i)
                frac[i] -= std::floor(frac[i]);
            canonical = origin + cell * frac;
        }

        const float ax = canonical.x;
        const float ay = canonical.y;
        const float az = canonical.z;
        const int   cx = (int)std::floor(ax * inv);
        const int   cy = (int)std::floor(ay * inv);
        const int   cz = (int)std::floor(az * inv);

        bool dup = false;
        for (int dx = -1; dx <= 1 && !dup; ++dx)
        for (int dy = -1; dy <= 1 && !dup; ++dy)
        for (int dz = -1; dz <= 1 && !dup; ++dz)
        {
            auto it = grid.find({cx+dx, cy+dy, cz+dz});
            if (it == grid.end()) continue;
            for (int idx : it->second)
            {
                const glm::vec3& b = canonicalPositions[(size_t)idx];
                const float ex = ax - b.x;
                const float ey = ay - b.y;
                const float ez = az - b.z;
                if (ex*ex + ey*ey + ez*ez < tol2) { dup = true; break; }
            }
        }
        if (!dup)
        {
            grid[{cx, cy, cz}].push_back((int)result.atoms.size());
            result.atoms.push_back(a);
            if (s.grainColors.size()==s.atoms.size()) result.grainColors.push_back(s.grainColors[index]);
            if (s.grainRegionIds.size()==s.atoms.size()) result.grainRegionIds.push_back(s.grainRegionIds[index]);
            canonicalPositions.push_back(canonical);
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
bool cscIsBounded(const std::vector<CellSlabPlane>& slabs,
                  const Structure& source)
{
    std::vector<glm::vec3> normals;
    normals.reserve(slabs.size());
    for (const auto& s : slabs)
        normals.push_back(cscNormal(s, source));

    if (normals.size() < 3) return false;

    glm::vec3 n0(0.0f), n1(0.0f);
    int found = 0;
    for (const auto& n : normals)
    {
        if (found == 0)
        {
            if (glm::length(n) < 1e-6f) continue;
            n0 = n; ++found;
        }
        else if (found == 1)
        {
            if (glm::length(glm::cross(n0, n)) > 1e-4f) { n1 = n; ++found; }
        }
        else
        {
            if (std::abs(glm::dot(glm::cross(n0, n1), n)) > 1e-4f)
                return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
std::vector<glm::vec3> cscClipHalfspace(
    const std::vector<glm::vec3>& poly,
    const glm::vec3& n, float plane_d)
{
    if (poly.empty()) return {};
    std::vector<glm::vec3> result;
    result.reserve(poly.size() + 1);
    for (size_t i = 0; i < poly.size(); ++i)
    {
        const glm::vec3& a = poly[i];
        const glm::vec3& b = poly[(i + 1) % poly.size()];
        const float da = glm::dot(a, n) - plane_d;
        const float db = glm::dot(b, n) - plane_d;
        if (da >= 0.0f) result.push_back(a);
        if ((da > 0.0f && db < 0.0f) || (da < 0.0f && db > 0.0f))
            result.push_back(a + (da / (da - db)) * (b - a));
    }
    return result;
}
