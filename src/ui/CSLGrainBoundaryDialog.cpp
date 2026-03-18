#include "ui/CSLGrainBoundaryDialog.h"

#include "ElementData.h"
#include "imgui.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace
{
struct WorkingAtom
{
    AtomSite atom;
    glm::vec3 pos;
    int grain = 0; // 1 or 2
    bool keep = true;
};

struct BuildResult
{
    bool success = false;
    std::string message;
    int sigma = 0;
    float thetaDeg = 0.0f;
    int inputAtoms = 0;
    int generatedAtoms = 0;
    int removedOverlap = 0;

    // Parameters used for the most recent successful generation.
    int basisIndex = 2;
    float latticeParameter = 0.0f;
    int atomicNumber = 0;
    int axis[3] = {0, 0, 1};
    int m = 1;
    int n = 1;
    int gbPlane[3] = {1, 0, 0};
    int dims[3] = {1, 1, 1};
    int reps[3] = {1, 1, 1};
    float overlapFraction = 0.0f;
    int removeFromGrain = 2;
    bool rigidTranslationEnabled = false;
    float rigidTx = 0.0f;
    float rigidTy = 0.0f;
    float vacuumPadding = 0.0f;
};

struct SigmaCandidate
{
    int sigma = 0;
    int m = 1;
    int n = 1;
    float thetaDeg = 0.0f;
};

enum class CubicBasisType
{
    SC = 0,
    BCC = 1,
    FCC = 2,
    Diamond = 3,
};

const char* cubicBasisName(CubicBasisType basis)
{
    switch (basis)
    {
        case CubicBasisType::SC: return "sc";
        case CubicBasisType::BCC: return "bcc";
        case CubicBasisType::FCC: return "fcc";
        case CubicBasisType::Diamond: return "diamond";
        default: return "unknown";
    }
}

int gcd3(int a, int b, int c)
{
    auto gcd2 = [](int p, int q) {
        p = std::abs(p);
        q = std::abs(q);
        while (q != 0)
        {
            int t = p % q;
            p = q;
            q = t;
        }
        return (p == 0) ? 1 : p;
    };

    int x = std::abs(a);
    int y = std::abs(b);
    int z = std::abs(c);
    int g = gcd2(gcd2(x, y), z);
    return (g == 0) ? 1 : g;
}

bool areCoprime(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    if (a == 0 && b == 0)
        return false;
    if (a == 0)
        return b == 1;
    if (b == 0)
        return a == 1;

    while (b != 0)
    {
        int t = a % b;
        a = b;
        b = t;
    }
    return a == 1;
}

void normalizeMillerTriplet(int& a, int& b, int& c)
{
    if (a == 0 && b == 0 && c == 0)
        return;

    int g = gcd3(a, b, c);
    a /= g;
    b /= g;
    c /= g;

    // (h k l) and (-h -k -l) represent the same plane.
    // Canonicalize by forcing the first non-zero component to be positive.
    if ((a < 0) ||
        (a == 0 && b < 0) ||
        (a == 0 && b == 0 && c < 0))
    {
        a = -a;
        b = -b;
        c = -c;
    }
}

glm::vec3 chooseStableInPlaneVector(const glm::vec3& normal)
{
    // Pick the world axis least aligned with normal for numerical stability.
    const glm::vec3 refs[3] = {
        glm::vec3(1.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 1.0f)
    };

    int best = 0;
    float bestDot = std::abs(glm::dot(normal, refs[0]));
    for (int i = 1; i < 3; ++i)
    {
        float d = std::abs(glm::dot(normal, refs[i]));
        if (d < bestDot)
        {
            bestDot = d;
            best = i;
        }
    }
    glm::vec3 candidate = glm::cross(normal, refs[best]);
    float l = glm::length(candidate);
    if (l < 1e-8f)
        return glm::vec3(1.0f, 0.0f, 0.0f);
    return candidate / l;
}

std::vector<SigmaCandidate> buildSigmaCandidates(int axisU, int axisV, int axisW, int sigmaMax)
{
    std::vector<SigmaCandidate> out;
    if (sigmaMax < 3 || (axisU == 0 && axisV == 0 && axisW == 0))
        return out;

    int g = gcd3(axisU, axisV, axisW);
    int ru = axisU / g;
    int rv = axisV / g;
    int rw = axisW / g;
    int axisNorm2 = ru * ru + rv * rv + rw * rw;
    if (axisNorm2 <= 0)
        return out;

    float axisNorm = std::sqrt((float)axisNorm2);
    std::map<int, SigmaCandidate> bestBySigma;

    // AIMSGB/Grimmer style search domain.
    int mMax = (int)std::ceil(std::sqrt(4.0f * (float)sigmaMax));
    for (int m = 0; m <= mMax; ++m)
    {
        for (int n = 1; n <= mMax; ++n)
        {
            if (!areCoprime(m, n))
                continue;

            int sigmaRaw = m * m + axisNorm2 * n * n;
            if (sigmaRaw <= 1)
                continue;

            int sigma = sigmaRaw;
            while (sigma % 2 == 0)
                sigma /= 2;

            if (sigma <= 1 || sigma > sigmaMax)
                continue;

            float theta = 0.0f;
            if (m == 0)
                theta = 3.14159265358979323846f;
            else
                theta = 2.0f * std::atan((float)n * axisNorm / (float)m);
            float thetaDeg = theta * 180.0f / 3.14159265358979323846f;

            SigmaCandidate c;
            c.sigma = sigma;
            c.m = m;
            c.n = n;
            c.thetaDeg = thetaDeg;

            std::map<int, SigmaCandidate>::iterator it = bestBySigma.find(sigma);
            if (it == bestBySigma.end() ||
                c.thetaDeg < it->second.thetaDeg ||
                (std::abs(c.thetaDeg - it->second.thetaDeg) < 1e-6f && (c.m + c.n) < (it->second.m + it->second.n)))
                bestBySigma[sigma] = c;
        }
    }

    out.reserve(bestBySigma.size());
    for (std::map<int, SigmaCandidate>::const_iterator it = bestBySigma.begin(); it != bestBySigma.end(); ++it)
        out.push_back(it->second);
    return out;
}

std::vector<glm::vec3> basisFractionalPositions(CubicBasisType basis)
{
    if (basis == CubicBasisType::SC)
    {
        return std::vector<glm::vec3>(1, glm::vec3(0.0f, 0.0f, 0.0f));
    }

    if (basis == CubicBasisType::BCC)
    {
        std::vector<glm::vec3> out;
        out.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
        out.push_back(glm::vec3(0.5f, 0.5f, 0.5f));
        return out;
    }

    if (basis == CubicBasisType::FCC)
    {
        std::vector<glm::vec3> out;
        out.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
        out.push_back(glm::vec3(0.0f, 0.5f, 0.5f));
        out.push_back(glm::vec3(0.5f, 0.0f, 0.5f));
        out.push_back(glm::vec3(0.5f, 0.5f, 0.0f));
        return out;
    }

    // Diamond = FCC + basis offset (1/4,1/4,1/4)
    std::vector<glm::vec3> out;
    out.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    out.push_back(glm::vec3(0.0f, 0.5f, 0.5f));
    out.push_back(glm::vec3(0.5f, 0.0f, 0.5f));
    out.push_back(glm::vec3(0.5f, 0.5f, 0.0f));
    out.push_back(glm::vec3(0.25f, 0.25f, 0.25f));
    out.push_back(glm::vec3(0.25f, 0.75f, 0.75f));
    out.push_back(glm::vec3(0.75f, 0.25f, 0.75f));
    out.push_back(glm::vec3(0.75f, 0.75f, 0.25f));
    return out;
}

float safeLength(const glm::vec3& v)
{
    return std::sqrt(glm::dot(v, v));
}

glm::vec3 normalizeOr(const glm::vec3& v, const glm::vec3& fallback)
{
    float l = safeLength(v);
    if (l < 1e-8f)
        return fallback;
    return v / l;
}

glm::vec3 rotateAroundAxis(const glm::vec3& p,
                           const glm::vec3& origin,
                           const glm::vec3& axisUnit,
                           float angleRad)
{
    glm::vec3 r = p - origin;
    float c = std::cos(angleRad);
    float s = std::sin(angleRad);
    glm::vec3 rotated = r * c + glm::cross(axisUnit, r) * s + axisUnit * glm::dot(axisUnit, r) * (1.0f - c);
    return origin + rotated;
}

BuildResult buildCslBoundary(Structure& structure,
                             CubicBasisType basis,
                             float latticeParameter,
                             int atomicNumber,
                             int axisU,
                             int axisV,
                             int axisW,
                             int selectedSigma,
                             float selectedThetaDeg,
                             int m,
                             int n,
                             int h,
                             int k,
                             int l,
                             int dimNormal,
                             int dimInplane1,
                             int dimInplane2,
                             float overlapFraction,
                             int removeFromGrain,
                             bool rigidTranslationEnabled,
                             float rigidTx,
                             float rigidTy,
                             float vacuumPadding,
                             int repU,
                             int repV,
                             int repW)
{
    BuildResult result;
    result.inputAtoms = 0;

    if (latticeParameter <= 0.0f)
    {
        result.message = "Lattice parameter must be positive.";
        return result;
    }
    if (atomicNumber <= 0 || atomicNumber > 118)
    {
        result.message = "Atomic number must be in range [1, 118].";
        return result;
    }
    if (m <= 0 || n <= 0)
    {
        result.message = "m and n must be positive integers.";
        return result;
    }
    if (selectedSigma <= 1)
    {
        result.message = "Selected Sigma must be greater than 1.";
        return result;
    }
    if (dimNormal <= 0 || dimInplane1 <= 0 || dimInplane2 <= 0)
    {
        result.message = "All dimension values must be positive.";
        return result;
    }
    if (repV <= 0 || repW <= 0)
    {
        result.message = "In-plane replication counts must be positive integers.";
        return result;
    }

    // Replication along GB-normal creates repeated GB images.
    // Keep one period only in the normal direction by design.
    (void)repU;
    repU = 1;
    if (axisU == 0 && axisV == 0 && axisW == 0)
    {
        result.message = "Rotation axis [u v w] cannot be [0 0 0].";
        return result;
    }
    if (h == 0 && k == 0 && l == 0)
    {
        result.message = "GB plane (h k l) cannot be (0 0 0).";
        return result;
    }

    int hNorm = h;
    int kNorm = k;
    int lNorm = l;
    normalizeMillerTriplet(hNorm, kNorm, lNorm);

    std::vector<glm::vec3> basisFrac = basisFractionalPositions(basis);
    if (basisFrac.empty())
    {
        result.message = "Invalid cubic basis selection.";
        return result;
    }

    std::vector<AtomSite> baseAtoms;
    baseAtoms.reserve(basisFrac.size());
    float cr = 0.5f, cg = 0.5f, cb = 0.5f;
    getDefaultElementColor(atomicNumber, cr, cg, cb);
    for (int i = 0; i < (int)basisFrac.size(); ++i)
    {
        const glm::vec3& f = basisFrac[i];
        AtomSite atom;
        atom.symbol = elementSymbol(atomicNumber);
        atom.atomicNumber = atomicNumber;
        atom.x = f.x * latticeParameter;
        atom.y = f.y * latticeParameter;
        atom.z = f.z * latticeParameter;
        atom.r = cr;
        atom.g = cg;
        atom.b = cb;
        baseAtoms.push_back(atom);
    }
    result.inputAtoms = (int)baseAtoms.size();

    glm::vec3 axis = normalizeOr(glm::vec3((float)axisU, (float)axisV, (float)axisW), glm::vec3(0.0f, 0.0f, 1.0f));
    glm::vec3 gbNormal = normalizeOr(glm::vec3((float)hNorm, (float)kNorm, (float)lNorm), glm::vec3(1.0f, 0.0f, 0.0f));

    // Build an in-plane basis from normal.
    glm::vec3 t1 = glm::cross(gbNormal, axis);
    if (safeLength(t1) < 1e-8f)
        t1 = chooseStableInPlaneVector(gbNormal);
    else
        t1 = normalizeOr(t1, chooseStableInPlaneVector(gbNormal));
    glm::vec3 t2 = normalizeOr(glm::cross(gbNormal, t1), glm::vec3(0.0f, 1.0f, 0.0f));

    int sigma = selectedSigma;
    float theta = selectedThetaDeg * 3.14159265358979323846f / 180.0f;
    float halfTheta = 0.5f * theta;

    result.sigma = sigma;
    result.thetaDeg = selectedThetaDeg;
    // Target bicrystal box in local GB coordinates:
    // u along normal (thickness), v/w in-plane.
    float halfU = (float)dimNormal * latticeParameter;
    float halfV = 0.5f * (float)dimInplane1 * latticeParameter;
    float halfW = 0.5f * (float)dimInplane2 * latticeParameter;

    float searchRadius = std::sqrt(halfU * halfU + halfV * halfV + halfW * halfW) + 2.0f * latticeParameter;
    int searchCells = std::max(1, (int)std::ceil(searchRadius / latticeParameter));

    std::vector<WorkingAtom> generated;
    generated.reserve((size_t)baseAtoms.size() * (size_t)(2 * searchCells + 1) * (size_t)(2 * searchCells + 1) * (size_t)(2 * searchCells + 1));

    glm::vec3 rigidShift = rigidTranslationEnabled ? (rigidTx * t1 + rigidTy * t2) : glm::vec3(0.0f);

    for (int i = -searchCells; i <= searchCells; ++i)
    {
        for (int j = -searchCells; j <= searchCells; ++j)
        {
            for (int kCell = -searchCells; kCell <= searchCells; ++kCell)
            {
                glm::vec3 latticeShift((float)i * latticeParameter,
                                       (float)j * latticeParameter,
                                       (float)kCell * latticeParameter);

                for (int ib = 0; ib < (int)baseAtoms.size(); ++ib)
                {
                    glm::vec3 p0 = glm::vec3((float)baseAtoms[ib].x, (float)baseAtoms[ib].y, (float)baseAtoms[ib].z) + latticeShift;

                    // Grain 1 (negative side, -theta/2)
                    glm::vec3 p1 = rotateAroundAxis(p0, glm::vec3(0.0f), axis, -halfTheta);
                    float u1 = glm::dot(p1, gbNormal);
                    float v1 = glm::dot(p1, t1);
                    float w1 = glm::dot(p1, t2);
                    if (u1 <= 0.0f && u1 >= -halfU && std::abs(v1) <= halfV && std::abs(w1) <= halfW)
                    {
                        WorkingAtom wa;
                        wa.atom = baseAtoms[ib];
                        wa.pos = p1;
                        wa.grain = 1;
                        generated.push_back(wa);
                    }

                    // Grain 2 (positive side, +theta/2) + optional rigid translation
                    glm::vec3 p2 = rotateAroundAxis(p0, glm::vec3(0.0f), axis, +halfTheta) + rigidShift;
                    float u2 = glm::dot(p2, gbNormal);
                    float v2 = glm::dot(p2, t1);
                    float w2 = glm::dot(p2, t2);
                    if (u2 >= 0.0f && u2 <= halfU && std::abs(v2) <= halfV && std::abs(w2) <= halfW)
                    {
                        WorkingAtom wa;
                        wa.atom = baseAtoms[ib];
                        wa.pos = p2;
                        wa.grain = 2;
                        generated.push_back(wa);
                    }
                }
            }
        }
    }

    if (generated.empty())
    {
        result.message = "No atoms generated. Try larger dimensions or different axis/plane.";
        return result;
    }

    float overlapDist = std::max(0.0f, overlapFraction) * latticeParameter;

    int removedOverlap = 0;
    if (overlapDist > 0.0f)
    {
        float overlapSq = overlapDist * overlapDist;
        for (int i = 0; i < (int)generated.size(); ++i)
        {
            if (!generated[i].keep || generated[i].grain == removeFromGrain)
                continue;

            for (int j = 0; j < (int)generated.size(); ++j)
            {
                if (!generated[j].keep || generated[j].grain != removeFromGrain)
                    continue;

                glm::vec3 d = generated[i].pos - generated[j].pos;
                if (glm::dot(d, d) < overlapSq)
                {
                    generated[j].keep = false;
                    ++removedOverlap;
                }
            }
        }
    }

    std::vector<AtomSite> outAtoms;
    outAtoms.reserve(generated.size());

    for (int i = 0; i < (int)generated.size(); ++i)
    {
        if (!generated[i].keep)
            continue;

        generated[i].atom.x = generated[i].pos.x;
        generated[i].atom.y = generated[i].pos.y;
        generated[i].atom.z = generated[i].pos.z;
        outAtoms.push_back(generated[i].atom);
    }

    if (outAtoms.empty())
    {
        result.message = "All atoms were removed by overlap settings.";
        return result;
    }

    float pad = std::max(0.0f, vacuumPadding);

    float lengthU = 2.0f * halfU + pad;
    float lengthV = 2.0f * halfV + pad;
    float lengthW = 2.0f * halfW + pad;

    glm::vec3 vecU = gbNormal * lengthU;
    glm::vec3 vecV = t1 * lengthV;
    glm::vec3 vecW = t2 * lengthW;
    glm::vec3 origin = -0.5f * vecU - 0.5f * vecV - 0.5f * vecW;

    // Robustly place every atom inside the final periodic box by mapping to
    // fractional coordinates and wrapping to [0,1).
    glm::mat3 cell(vecU, vecV, vecW);
    float det = glm::determinant(cell);
    if (std::abs(det) < 1e-10f)
    {
        result.message = "Generated cell is singular. Try different axis/plane or larger dimensions.";
        return result;
    }
    glm::mat3 invCell = glm::inverse(cell);

    const double eps = 1e-8;
    for (int i = 0; i < (int)outAtoms.size(); ++i)
    {
        glm::vec3 p((float)outAtoms[i].x, (float)outAtoms[i].y, (float)outAtoms[i].z);
        glm::vec3 frac = invCell * (p - origin);

        frac.x -= std::floor(frac.x);
        frac.y -= std::floor(frac.y);
        frac.z -= std::floor(frac.z);

        // Avoid edge-case rendering ambiguity for atoms that land exactly on 1.0.
        if (frac.x >= 1.0f - (float)eps) frac.x = 0.0f;
        if (frac.y >= 1.0f - (float)eps) frac.y = 0.0f;
        if (frac.z >= 1.0f - (float)eps) frac.z = 0.0f;

        glm::vec3 wrapped = origin + cell * frac;
        outAtoms[i].x = wrapped.x;
        outAtoms[i].y = wrapped.y;
        outAtoms[i].z = wrapped.z;
    }

    // Optional supercell replication of the generated GB box.
    std::vector<AtomSite> replicated;
    replicated.reserve((size_t)outAtoms.size() * (size_t)repU * (size_t)repV * (size_t)repW);
    for (int iu = 0; iu < repU; ++iu)
    {
        for (int iv = 0; iv < repV; ++iv)
        {
            for (int iw = 0; iw < repW; ++iw)
            {
                glm::vec3 shift = (float)iu * vecU + (float)iv * vecV + (float)iw * vecW;
                for (int i = 0; i < (int)outAtoms.size(); ++i)
                {
                    AtomSite a = outAtoms[i];
                    a.x += shift.x;
                    a.y += shift.y;
                    a.z += shift.z;
                    replicated.push_back(a);
                }
            }
        }
    }

    glm::vec3 vecUFinal = vecU * (float)repU;
    glm::vec3 vecVFinal = vecV * (float)repV;
    glm::vec3 vecWFinal = vecW * (float)repW;

    structure.atoms.swap(replicated);
    structure.hasUnitCell = true;
    structure.cellOffset = { origin.x, origin.y, origin.z };
    structure.cellVectors = {{
        {{ vecUFinal.x, vecUFinal.y, vecUFinal.z }},
        {{ vecVFinal.x, vecVFinal.y, vecVFinal.z }},
        {{ vecWFinal.x, vecWFinal.y, vecWFinal.z }}
    }};

    result.generatedAtoms = (int)structure.atoms.size();
    result.removedOverlap = removedOverlap;
    result.success = true;

    result.basisIndex = (int)basis;
    result.latticeParameter = latticeParameter;
    result.atomicNumber = atomicNumber;
    result.axis[0] = axisU; result.axis[1] = axisV; result.axis[2] = axisW;
    result.m = m; result.n = n;
    result.gbPlane[0] = hNorm; result.gbPlane[1] = kNorm; result.gbPlane[2] = lNorm;
    result.dims[0] = dimNormal; result.dims[1] = dimInplane1; result.dims[2] = dimInplane2;
    result.reps[0] = repU; result.reps[1] = repV; result.reps[2] = repW;
    result.overlapFraction = overlapFraction;
    result.removeFromGrain = removeFromGrain;
    result.rigidTranslationEnabled = rigidTranslationEnabled;
    result.rigidTx = rigidTx;
    result.rigidTy = rigidTy;
    result.vacuumPadding = vacuumPadding;

    result.message = std::string("CSL grain boundary generated from ") + cubicBasisName(basis) + " cubic lattice.";
    return result;
}

} // namespace

void CSLGrainBoundaryDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("CSL Grain Boundary...", NULL, false, enabled))
        m_openRequested = true;
}

void CSLGrainBoundaryDialog::drawDialog(Structure& structure,
                                        const std::function<void(Structure&)>& updateBuffers)
{
    static int axis[3] = {0, 0, 1};
    static int basisIndex = (int)CubicBasisType::FCC;
    static float latticeParameter = 3.60f;
    static int atomicNumber = 29;
    static int sigmaMax = 200;
    static std::vector<SigmaCandidate> sigmaCandidates;
    static int sigmaSelection = 0;
    static int lastAxisForSigma[3] = {0, 0, 1};
    static int lastSigmaMaxForSigma = 200;
    static int gbPlane[3] = {1, 0, 0};
    static int dims[3] = {2, 2, 2}; // normal, inplane1, inplane2
    static int replicationsInPlane[2] = {1, 1};
    static float overlapFraction = 0.20f;
    static int removeFromGrain = 2; // 1 or 2
    static bool rigidTrans = false;
    static float rigidT[2] = {0.0f, 0.0f};
    static float vacuumPadding = 0.0f;
    static BuildResult lastResult;

    if (m_openRequested)
    {
        ImGui::OpenPopup("CSL Grain Boundary Builder");
        m_openRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(860.0f, 760.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("CSL Grain Boundary Builder", &dialogOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::TextWrapped("Inspired by CSL grain-boundary workflows for cubic systems. "
                           "This in-app builder creates a bicrystal from an ideal cubic lattice (sc/bcc/fcc/diamond) and applies "
                           "misorientation, optional rigid translation, and overlap removal.");

        ImGui::Separator();
        ImGui::Text("Cubic Lattice Source");
        const char* basisItems[] = { "sc", "bcc", "fcc", "diamond" };
        ImGui::Combo("Basis", &basisIndex, basisItems, 4);
        ImGui::DragFloat("Lattice parameter (A)", &latticeParameter, 0.01f, 0.1f, 20.0f, "%.4f");
        ImGui::InputInt("Atomic number (Z)", &atomicNumber);
        if (atomicNumber < 1) atomicNumber = 1;
        if (atomicNumber > 118) atomicNumber = 118;
        ImGui::Text("Element: %s (%s)", elementName(atomicNumber), elementSymbol(atomicNumber));

        ImGui::Separator();
        ImGui::Text("Misorientation Parameters");
        ImGui::InputInt3("Axis [u v w]", axis);
        ImGui::InputInt("Sigma search max", &sigmaMax);
        if (sigmaMax < 3) sigmaMax = 3;

        if (axis[0] != lastAxisForSigma[0] || axis[1] != lastAxisForSigma[1] || axis[2] != lastAxisForSigma[2] ||
            sigmaMax != lastSigmaMaxForSigma)
        {
            sigmaCandidates.clear();
            sigmaSelection = 0;
            lastAxisForSigma[0] = axis[0];
            lastAxisForSigma[1] = axis[1];
            lastAxisForSigma[2] = axis[2];
            lastSigmaMaxForSigma = sigmaMax;
        }

        ImGui::SameLine();
        if (ImGui::Button("Generate Sigma List"))
        {
            sigmaCandidates = buildSigmaCandidates(axis[0], axis[1], axis[2], sigmaMax);
            sigmaSelection = 0;
        }

        if (!sigmaCandidates.empty())
        {
            auto sigmaGetter = [](void* data, int idx) -> const char* {
                static char label[96];
                std::vector<SigmaCandidate>* vec = static_cast<std::vector<SigmaCandidate>*>(data);
                if (idx < 0 || idx >= (int)vec->size()) return "";
                const SigmaCandidate& c = (*vec)[idx];
                std::snprintf(label, sizeof(label), "Sigma %d  |  m=%d n=%d  |  theta=%.4f deg", c.sigma, c.m, c.n, c.thetaDeg);
                return label;
            };
            if (sigmaSelection >= (int)sigmaCandidates.size())
                sigmaSelection = 0;
            ImGui::Combo("Sigma candidate", &sigmaSelection, sigmaGetter, &sigmaCandidates, (int)sigmaCandidates.size());

            const SigmaCandidate& selected = sigmaCandidates[sigmaSelection];
            ImGui::Text("Selected Sigma: %d", selected.sigma);
            ImGui::Text("Selected m,n: %d, %d", selected.m, selected.n);
            ImGui::Text("Misorientation angle: %.4f deg", selected.thetaDeg);
        }
        else
        {
            ImGui::TextDisabled("No Sigma candidates generated yet.");
            ImGui::TextDisabled("Use 'Generate Sigma List' after setting axis and search max.");
        }

        ImGui::Separator();
        ImGui::Text("Boundary / Cell Parameters");
        ImGui::InputInt3("GB Plane (h k l)", gbPlane);
        ImGui::InputInt3("Dimensions [l1 l2 l3]", dims);
        ImGui::TextDisabled("l1=normal replication per grain, l2/l3=in-plane replication.");
        ImGui::Text("Replication along normal is fixed to 1 (to avoid multiple GB images).");
        ImGui::InputInt2("In-plane replications [r2 r3]", replicationsInPlane);
        if (replicationsInPlane[0] < 1) replicationsInPlane[0] = 1;
        if (replicationsInPlane[1] < 1) replicationsInPlane[1] = 1;
        ImGui::TextDisabled("Replicate only parallel to GB plane.");

        ImGui::DragFloat("Overlap distance fraction", &overlapFraction, 0.01f, 0.0f, 2.0f, "%.3f");
        ImGui::TextDisabled("Atoms closer than fraction*lattice_scale are removed from selected grain.");

        int removeChoice = (removeFromGrain == 1) ? 0 : 1;
        const char* removeItems[] = { "Remove from grain 1", "Remove from grain 2" };
        if (ImGui::Combo("Overlap removal", &removeChoice, removeItems, 2))
            removeFromGrain = (removeChoice == 0) ? 1 : 2;

        ImGui::Checkbox("Enable rigid translation of grain 2", &rigidTrans);
        if (rigidTrans)
            ImGui::InputFloat2("Rigid translation (in-plane)", rigidT);

        ImGui::DragFloat("Vacuum/box padding (A)", &vacuumPadding, 0.05f, 0.0f, 20.0f, "%.2f");

        ImGui::Separator();
        if (ImGui::Button("Generate Grain Boundary (Replace Structure)", ImVec2(-1.0f, 0.0f)))
        {
            if (sigmaCandidates.empty())
            {
                lastResult.success = false;
                lastResult.message = "Generate and select a Sigma candidate first.";
            }
            else
            {
            const SigmaCandidate& selected = sigmaCandidates[sigmaSelection];
            CubicBasisType basis = CubicBasisType::FCC;
            if (basisIndex == 0) basis = CubicBasisType::SC;
            else if (basisIndex == 1) basis = CubicBasisType::BCC;
            else if (basisIndex == 2) basis = CubicBasisType::FCC;
            else if (basisIndex == 3) basis = CubicBasisType::Diamond;

            lastResult = buildCslBoundary(structure,
                                          basis,
                                          latticeParameter,
                                          atomicNumber,
                                          axis[0], axis[1], axis[2],
                                          selected.sigma,
                                          selected.thetaDeg,
                                          selected.m, selected.n,
                                          gbPlane[0], gbPlane[1], gbPlane[2],
                                          dims[0], dims[1], dims[2],
                                          overlapFraction,
                                          removeFromGrain,
                                          rigidTrans,
                                          rigidT[0], rigidT[1],
                                          vacuumPadding,
                                          1, replicationsInPlane[0], replicationsInPlane[1]);
            if (lastResult.success)
                updateBuffers(structure);
            }
        }

        ImGui::Spacing();
        if (!lastResult.message.empty())
        {
            ImGui::TextWrapped("Status: %s", lastResult.message.c_str());
            if (lastResult.success)
            {
                ImGui::Text("Input atoms: %d", lastResult.inputAtoms);
                ImGui::Text("Output atoms: %d", lastResult.generatedAtoms);
                ImGui::Text("Removed overlaps: %d", lastResult.removedOverlap);
                ImGui::Text("Sigma: %d", lastResult.sigma);
                ImGui::Text("Misorientation angle: %.4f deg", lastResult.thetaDeg);
                ImGui::Separator();
                ImGui::Text("GB Parameters Used");
                ImGui::Text("Basis: %s", cubicBasisName((CubicBasisType)lastResult.basisIndex));
                ImGui::Text("Lattice parameter: %.4f A", lastResult.latticeParameter);
                ImGui::Text("Element: Z=%d (%s)", lastResult.atomicNumber, elementSymbol(lastResult.atomicNumber));
                ImGui::Text("Axis [u v w]: [%d %d %d]", lastResult.axis[0], lastResult.axis[1], lastResult.axis[2]);
                ImGui::Text("Sigma / m / n: %d / %d / %d", lastResult.sigma, lastResult.m, lastResult.n);
                ImGui::Text("GB plane (h k l): (%d %d %d)", lastResult.gbPlane[0], lastResult.gbPlane[1], lastResult.gbPlane[2]);
                ImGui::Text("Dimensions [l1 l2 l3]: [%d %d %d]", lastResult.dims[0], lastResult.dims[1], lastResult.dims[2]);
                ImGui::Text("Replications [r1 r2 r3]: [%d %d %d]", lastResult.reps[0], lastResult.reps[1], lastResult.reps[2]);
                ImGui::Text("Overlap fraction: %.3f", lastResult.overlapFraction);
                ImGui::Text("Overlap removal target grain: %d", lastResult.removeFromGrain);
                ImGui::Text("Rigid translation enabled: %s", lastResult.rigidTranslationEnabled ? "yes" : "no");
                ImGui::Text("Rigid translation (in-plane): [%.4f %.4f]", lastResult.rigidTx, lastResult.rigidTy);
                ImGui::Text("Vacuum/box padding: %.3f A", lastResult.vacuumPadding);
            }
        }

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();
}
