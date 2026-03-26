#include "ui/BulkCrystalBuilderDialog.h"

#include "ElementData.h"
#include "imgui.h"
#include "ui/PeriodicTableDialog.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

namespace
{
enum class CrystalSystem
{
    Triclinic = 0,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
};

struct BulkBuildResult
{
    bool success = false;
    std::string message;
    int generatedAtoms = 0;
    int spaceGroupNumber = 0;
    std::string spaceGroupSymbol;
};

struct SpaceGroupRange
{
    CrystalSystem system;
    const char* label;
    int first;
    int last;
};

struct LatticeParameters
{
    double a = 4.0;
    double b = 4.0;
    double c = 4.0;
    double alpha = 90.0;
    double beta = 90.0;
    double gamma = 90.0;
};

const SpaceGroupRange kSpaceGroupRanges[] = {
    { CrystalSystem::Triclinic,    "Triclinic",    1,   2 },
    { CrystalSystem::Monoclinic,   "Monoclinic",   3,  15 },
    { CrystalSystem::Orthorhombic, "Orthorhombic", 16, 74 },
    { CrystalSystem::Tetragonal,   "Tetragonal",   75, 142 },
    { CrystalSystem::Trigonal,     "Trigonal",     143, 167 },
    { CrystalSystem::Hexagonal,    "Hexagonal",    168, 194 },
    { CrystalSystem::Cubic,        "Cubic",        195, 230 },
};

bool hasElementColor(int atomicNumber, const std::vector<glm::vec3>& elementColors)
{
    return atomicNumber >= 0 && atomicNumber < (int)elementColors.size();
}

void applyElementToAtom(AtomSite& atom,
                        int atomicNumber,
                        const std::vector<glm::vec3>& elementColors)
{
    atom.atomicNumber = atomicNumber;
    atom.symbol = elementSymbol(atomicNumber);

    if (hasElementColor(atomicNumber, elementColors))
    {
        atom.r = elementColors[atomicNumber].r;
        atom.g = elementColors[atomicNumber].g;
        atom.b = elementColors[atomicNumber].b;
    }
    else
    {
        getDefaultElementColor(atomicNumber, atom.r, atom.g, atom.b);
    }
}

const char* crystalSystemLabel(CrystalSystem system)
{
    for (size_t i = 0; i < sizeof(kSpaceGroupRanges) / sizeof(kSpaceGroupRanges[0]); ++i)
    {
        if (kSpaceGroupRanges[i].system == system)
            return kSpaceGroupRanges[i].label;
    }
    return "Unknown";
}

const SpaceGroupRange& currentRange(CrystalSystem system)
{
    for (size_t i = 0; i < sizeof(kSpaceGroupRanges) / sizeof(kSpaceGroupRanges[0]); ++i)
    {
        if (kSpaceGroupRanges[i].system == system)
            return kSpaceGroupRanges[i];
    }
    return kSpaceGroupRanges[0];
}

glm::mat3 buildLatticeFromParameters(const LatticeParameters& params)
{
    const double degToRad = 3.14159265358979323846 / 180.0;
    const double alpha = params.alpha * degToRad;
    const double beta = params.beta * degToRad;
    const double gamma = params.gamma * degToRad;

    const double cosAlpha = std::cos(alpha);
    const double cosBeta = std::cos(beta);
    const double cosGamma = std::cos(gamma);
    const double sinGamma = std::sin(gamma);

    glm::vec3 a((float)params.a, 0.0f, 0.0f);
    glm::vec3 b((float)(params.b * cosGamma),
                (float)(params.b * sinGamma),
                0.0f);

    double cx = params.c * cosBeta;
    double cy = 0.0;
    if (std::abs(sinGamma) > 1e-10)
        cy = params.c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
    double cz2 = params.c * params.c - cx * cx - cy * cy;
    if (cz2 < 0.0 && cz2 > -1e-8)
        cz2 = 0.0;
    glm::vec3 c((float)cx, (float)cy, (float)((cz2 > 0.0) ? std::sqrt(cz2) : 0.0));

    return glm::mat3(a, b, c);
}

bool validateParameters(CrystalSystem system, const LatticeParameters& params, std::string& error)
{
    auto near = [](double lhs, double rhs) {
        return std::abs(lhs - rhs) < 1e-6;
    };

    if (params.a <= 0.0 || params.b <= 0.0 || params.c <= 0.0)
    {
        error = "Lattice lengths must be positive.";
        return false;
    }

    if (params.alpha <= 0.0 || params.alpha >= 180.0 ||
        params.beta <= 0.0 || params.beta >= 180.0 ||
        params.gamma <= 0.0 || params.gamma >= 180.0)
    {
        error = "Angles must be in the range (0, 180).";
        return false;
    }

    switch (system)
    {
        case CrystalSystem::Monoclinic:
            if (!near(params.alpha, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Monoclinic builder uses alpha=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Orthorhombic:
            if (!near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Orthorhombic builder uses alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Tetragonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Tetragonal builder uses a=b and alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Trigonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 120.0))
            {
                error = "Trigonal builder uses the hexagonal setting: a=b, alpha=beta=90 deg, gamma=120 deg.";
                return false;
            }
            break;
        case CrystalSystem::Hexagonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 120.0))
            {
                error = "Hexagonal builder uses a=b, alpha=beta=90 deg, gamma=120 deg.";
                return false;
            }
            break;
        case CrystalSystem::Cubic:
            if (!near(params.a, params.b) || !near(params.b, params.c) ||
                !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Cubic builder uses a=b=c and alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Triclinic:
        default:
            break;
    }

    glm::mat3 lattice = buildLatticeFromParameters(params);
    if (std::abs(glm::determinant(lattice)) < 1e-10f)
    {
        error = "Lattice parameters produce a zero-volume cell.";
        return false;
    }

    return true;
}

void applySystemConstraints(CrystalSystem system, LatticeParameters& params)
{
    switch (system)
    {
        case CrystalSystem::Monoclinic:
            params.alpha = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Orthorhombic:
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Tetragonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Trigonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 120.0;
            break;
        case CrystalSystem::Hexagonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 120.0;
            break;
        case CrystalSystem::Cubic:
            params.b = params.a;
            params.c = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Triclinic:
        default:
            break;
    }
}

void drawLatticeParameterInputs(CrystalSystem system, LatticeParameters& latticeParams)
{
    const float fieldWidth = 100.0f;
    const float spacing = 8.0f;
    const float hintGap = 20.0f;

    auto inputField = [&](const char* label, double* val) {
        ImGui::PushItemWidth(fieldWidth);
        ImGui::InputDouble(label, val, 0.0, 0.0, "%.4f");
        ImGui::PopItemWidth();
    };

    switch (system)
    {
        case CrystalSystem::Triclinic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            inputField("alpha", &latticeParams.alpha);
            ImGui::SameLine(0.0f, spacing); inputField("beta", &latticeParams.beta);
            ImGui::SameLine(0.0f, spacing); inputField("gamma", &latticeParams.gamma);
            break;
        case CrystalSystem::Monoclinic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            inputField("beta", &latticeParams.beta);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("alpha=gamma=90");
            latticeParams.alpha = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Orthorhombic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("b", &latticeParams.b);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("alpha=beta=gamma=90");
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Tetragonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=a, angles=90");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
        case CrystalSystem::Trigonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("hex: b=a, gamma=120");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Hexagonal:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, spacing); inputField("c", &latticeParams.c);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=a, gamma=120");
            latticeParams.b = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 120.0;
            break;
        case CrystalSystem::Cubic:
            inputField("a", &latticeParams.a);
            ImGui::SameLine(0.0f, hintGap); ImGui::AlignTextToFramePadding(); ImGui::TextDisabled("b=c=a, angles=90");
            latticeParams.b = latticeParams.a;
            latticeParams.c = latticeParams.a;
            latticeParams.alpha = 90.0;
            latticeParams.beta = 90.0;
            latticeParams.gamma = 90.0;
            break;
    }
}

std::vector<int>& hallBySpaceGroup()
{
    static std::vector<int> hallByNumber(231, 0);
    static bool initialized = false;
    if (!initialized)
    {
#ifdef ATOMS_ENABLE_SPGLIB
        for (int hall = 1; hall <= 530; ++hall)
        {
            SpglibSpacegroupType type = spg_get_spacegroup_type(hall);
            if (type.number >= 1 && type.number <= 230 && hallByNumber[type.number] == 0)
                hallByNumber[type.number] = hall;
        }
#endif
        initialized = true;
    }
    return hallByNumber;
}

std::string spaceGroupLabel(int spaceGroupNumber)
{
#ifdef ATOMS_ENABLE_SPGLIB
    const std::vector<int>& hallMap = hallBySpaceGroup();
    if (spaceGroupNumber >= 1 && spaceGroupNumber <= 230 && hallMap[spaceGroupNumber] > 0)
    {
        SpglibSpacegroupType type = spg_get_spacegroup_type(hallMap[spaceGroupNumber]);
        std::ostringstream out;
        out << type.number << " " << type.international_short;
        return out.str();
    }
#endif
    std::ostringstream out;
    out << spaceGroupNumber;
    return out.str();
}

glm::vec3 wrapFractional(const glm::vec3& frac)
{
    glm::vec3 wrapped = frac;
    wrapped.x -= std::floor(wrapped.x);
    wrapped.y -= std::floor(wrapped.y);
    wrapped.z -= std::floor(wrapped.z);
    if (wrapped.x >= 1.0f - 1e-6f) wrapped.x = 0.0f;
    if (wrapped.y >= 1.0f - 1e-6f) wrapped.y = 0.0f;
    if (wrapped.z >= 1.0f - 1e-6f) wrapped.z = 0.0f;
    return wrapped;
}

bool sameFractionalPosition(const glm::vec3& lhs, const glm::vec3& rhs, float tol)
{
    glm::vec3 delta = lhs - rhs;
    delta.x -= std::round(delta.x);
    delta.y -= std::round(delta.y);
    delta.z -= std::round(delta.z);
    return std::abs(delta.x) <= tol && std::abs(delta.y) <= tol && std::abs(delta.z) <= tol;
}

BulkBuildResult buildBulkCrystal(Structure& structure,
                                 CrystalSystem system,
                                 int spaceGroupNumber,
                                 const LatticeParameters& params,
                                 const std::vector<AtomSite>& asymmetricAtoms,
                                 const std::vector<glm::vec3>& elementColors)
{
    BulkBuildResult result;
    result.spaceGroupNumber = spaceGroupNumber;
    result.spaceGroupSymbol = spaceGroupLabel(spaceGroupNumber);

#ifndef ATOMS_ENABLE_SPGLIB
    (void)system;
    (void)params;
    (void)asymmetricAtoms;
    (void)elementColors;
    result.message = "Bulk crystal builder requires spglib at build time.";
    return result;
#else
    std::string validationError;
    if (!validateParameters(system, params, validationError))
    {
        result.message = validationError;
        return result;
    }
    if (asymmetricAtoms.empty())
    {
        result.message = "Add at least one asymmetric-unit atom before building.";
        return result;
    }

    const std::vector<int>& hallMap = hallBySpaceGroup();
    if (spaceGroupNumber < 1 || spaceGroupNumber > 230 || hallMap[spaceGroupNumber] == 0)
    {
        result.message = "Selected space group is not available in the symmetry database.";
        return result;
    }

    const int hallNumber = hallMap[spaceGroupNumber];
    int rotations[192][3][3];
    double translations[192][3];
    int operationCount = spg_get_symmetry_from_database(rotations, translations, hallNumber);
    if (operationCount <= 0)
    {
        result.message = "Could not load symmetry operations for the selected space group.";
        return result;
    }

    glm::mat3 lattice = buildLatticeFromParameters(params);
    glm::vec3 a = lattice[0];
    glm::vec3 b = lattice[1];
    glm::vec3 c = lattice[2];

    std::vector<AtomSite> generatedAtoms;
    generatedAtoms.reserve(asymmetricAtoms.size() * (size_t)operationCount);
    // Keep fractional coords in a separate list for deduplication.
    // generatedAtoms.x/y/z are overwritten with Cartesian before push_back,
    // so they cannot be used for the fractional duplicate check.
    std::vector<glm::vec3> generatedFrac;
    generatedFrac.reserve(asymmetricAtoms.size() * (size_t)operationCount);

    for (size_t atomIndex = 0; atomIndex < asymmetricAtoms.size(); ++atomIndex)
    {
        const AtomSite& basisAtom = asymmetricAtoms[atomIndex];
        glm::vec3 frac((float)basisAtom.x, (float)basisAtom.y, (float)basisAtom.z);
        frac = wrapFractional(frac);

        for (int op = 0; op < operationCount; ++op)
        {
            glm::vec3 transformed(0.0f);
            for (int row = 0; row < 3; ++row)
            {
                transformed[row] = (float)translations[op][row];
                for (int col = 0; col < 3; ++col)
                    transformed[row] += (float)rotations[op][row][col] * frac[col];
            }
            transformed = wrapFractional(transformed);

            bool duplicate = false;
            for (size_t existingIndex = 0; existingIndex < generatedAtoms.size(); ++existingIndex)
            {
                if (generatedAtoms[existingIndex].atomicNumber != basisAtom.atomicNumber)
                    continue;
                if (sameFractionalPosition(transformed, generatedFrac[existingIndex], 1e-5f))
                {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate)
                continue;

            glm::vec3 cart = a * transformed.x + b * transformed.y + c * transformed.z;

            AtomSite generated;
            generated.atomicNumber = basisAtom.atomicNumber;
            generated.symbol = elementSymbol(generated.atomicNumber);
            generated.x = transformed.x;
            generated.y = transformed.y;
            generated.z = transformed.z;
            if (generated.atomicNumber >= 0 && generated.atomicNumber < (int)elementColors.size())
            {
                generated.r = elementColors[generated.atomicNumber].r;
                generated.g = elementColors[generated.atomicNumber].g;
                generated.b = elementColors[generated.atomicNumber].b;
            }
            else
            {
                getDefaultElementColor(generated.atomicNumber, generated.r, generated.g, generated.b);
            }

            generated.x = cart.x;
            generated.y = cart.y;
            generated.z = cart.z;
            generatedFrac.push_back(transformed);
            generatedAtoms.push_back(generated);
        }
    }

    if (generatedAtoms.empty())
    {
        result.message = "Symmetry expansion produced no atoms.";
        return result;
    }

    structure.atoms.swap(generatedAtoms);
    structure.hasUnitCell = true;
    structure.cellOffset = {{0.0, 0.0, 0.0}};
    structure.cellVectors = {{
        {{ a.x, a.y, a.z }},
        {{ b.x, b.y, b.z }},
        {{ c.x, c.y, c.z }}
    }};

    result.generatedAtoms = (int)structure.atoms.size();
    result.success = true;
    result.message = "Bulk crystal generated from asymmetric-unit atoms and selected space-group symmetry.";
    return result;
#endif
}

void addDefaultAsymmetricAtom(std::vector<AtomSite>& atoms,
                             const std::vector<glm::vec3>& elementColors)
{
    AtomSite atom;
    applyElementToAtom(atom, 1, elementColors);
    atom.x = 0.0;
    atom.y = 0.0;
    atom.z = 0.0;
    atoms.push_back(atom);
}
} // namespace

BulkCrystalBuilderDialog::BulkCrystalBuilderDialog() = default;

void BulkCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Bulk Crystal", NULL, false, enabled))
        m_openRequested = true;
}

void BulkCrystalBuilderDialog::drawDialog(Structure& structure,
                                          const std::vector<glm::vec3>& elementColors,
                                          const std::function<void(Structure&)>& updateBuffers)
{
    static int crystalSystemIndex = (int)CrystalSystem::Cubic;
    static int selectedSpaceGroup = 225;
    static LatticeParameters latticeParams;
    static std::vector<AtomSite> asymmetricAtoms;
    static BulkBuildResult lastResult;
    static int selectedEditElement = 6;
    static int targetAtomIndex = -1;
    static bool showElementPicker = false;
    static bool reopenBulkDialogAfterPicker = false;
    static bool scrollRowsToBottom = false;
    static int lastCrystalSystemIndex = crystalSystemIndex;

    if (asymmetricAtoms.empty())
        addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);

    if (m_openRequested)
    {
        ImGui::OpenPopup("Build Bulk Crystal");
        m_openRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(720.0f, 540.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("Build Bulk Crystal", &dialogOpen, ImGuiWindowFlags_None))
    {
#ifndef ATOMS_ENABLE_SPGLIB
        ImGui::TextWrapped("spglib is not available in this build, so symmetry expansion cannot be generated.");
#else
        const char* systemLabels[] = {
            "Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic"
        };

        const float halfWidth = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.5f;

        ImGui::PushItemWidth(halfWidth);
        if (ImGui::Combo("##system", &crystalSystemIndex, systemLabels, 7))
        {
            const SpaceGroupRange& range = currentRange((CrystalSystem)crystalSystemIndex);
            selectedSpaceGroup = range.first;
        }
        ImGui::PopItemWidth();
        if (crystalSystemIndex != lastCrystalSystemIndex)
        {
            applySystemConstraints((CrystalSystem)crystalSystemIndex, latticeParams);
            lastCrystalSystemIndex = crystalSystemIndex;
        }

        const SpaceGroupRange& range = currentRange((CrystalSystem)crystalSystemIndex);
        if (selectedSpaceGroup < range.first || selectedSpaceGroup > range.last)
            selectedSpaceGroup = range.first;
        int sgIndex = selectedSpaceGroup - range.first;
        auto sgGetter = [](void* data, int idx) -> const char* {
            static char label[64];
            const SpaceGroupRange* sgRange = static_cast<const SpaceGroupRange*>(data);
            int sgNumber = sgRange->first + idx;
            std::snprintf(label, sizeof(label), "%s", spaceGroupLabel(sgNumber).c_str());
            return label;
        };
        ImGui::SameLine();
        ImGui::PushItemWidth(halfWidth);
        if (ImGui::Combo("##spacegroup", &sgIndex, sgGetter, (void*)&range, range.last - range.first + 1))
            selectedSpaceGroup = range.first + sgIndex;
        ImGui::PopItemWidth();

        ImGui::Separator();
        drawLatticeParameterInputs((CrystalSystem)crystalSystemIndex, latticeParams);

        ImGui::Separator();
        ImGui::Text("Asymmetric Unit  ");
        ImGui::SameLine();
        if (ImGui::Button("+ Add"))
        {
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            scrollRowsToBottom = true;
            const AtomSite& atom = asymmetricAtoms.back();
            std::cout << "[Operation] Added asymmetric atom (bulk builder): "
                      << atom.symbol << "(" << atom.atomicNumber << ")"
                      << " at [" << atom.x << ", " << atom.y << ", " << atom.z << "]"
                      << std::endl;
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear"))
        {
            asymmetricAtoms.clear();
            addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            std::cout << "[Operation] Cleared asymmetric atoms (bulk builder); reset to default H(1)" << std::endl;
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(fractional coords)");

        int pendingDelete = -1;
        float listHeight = ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 2 - ImGui::GetStyle().ItemSpacing.y * 2;
        if (listHeight < 60.0f) listHeight = 60.0f;
        if (ImGui::BeginChild("##bulk-atom-rows", ImVec2(-1.0f, listHeight), true))
        {
            for (int i = 0; i < (int)asymmetricAtoms.size(); ++i)
            {
                AtomSite& atom = asymmetricAtoms[i];
                ImGui::PushID(i);

                if (ImGui::Button("Delete"))
                    pendingDelete = i;
                ImGui::SameLine();
                ImGui::Text("#%d", i);
                ImGui::SameLine();

                float editPos[3] = { (float)atom.x, (float)atom.y, (float)atom.z };

                ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 50.0f);
                if (ImGui::DragFloat3("##bulk-direct", editPos, 0.005f, -10.0f, 10.0f, "%.4f"))
                {
                    atom.x = editPos[0];
                    atom.y = editPos[1];
                    atom.z = editPos[2];
                }
                ImGui::PopItemWidth();

                ImGui::SameLine();
                std::string elemButton = std::string(atom.symbol) + "##bulk-element";
                if (ImGui::Button(elemButton.c_str()))
                {
                    targetAtomIndex = i;
                    selectedEditElement = atom.atomicNumber;
                    showElementPicker = true;
                    reopenBulkDialogAfterPicker = true;
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Click to substitute element");

                ImGui::PopID();
            }

            if (scrollRowsToBottom)
            {
                ImGui::SetScrollY(ImGui::GetScrollMaxY());
                scrollRowsToBottom = false;
            }
            ImGui::EndChild();
        }

        if (pendingDelete >= 0 && pendingDelete < (int)asymmetricAtoms.size())
        {
            const AtomSite deletedAtom = asymmetricAtoms[pendingDelete];
            asymmetricAtoms.erase(asymmetricAtoms.begin() + pendingDelete);
            if (asymmetricAtoms.empty())
                addDefaultAsymmetricAtom(asymmetricAtoms, elementColors);
            std::cout << "[Operation] Deleted asymmetric atom (bulk builder): "
                      << deletedAtom.symbol << "(" << deletedAtom.atomicNumber << ")"
                      << " row=" << pendingDelete << std::endl;
        }

        ImGui::Separator();
            if (ImGui::Button("Build"))
        {
            applySystemConstraints((CrystalSystem)crystalSystemIndex, latticeParams);
            lastResult = buildBulkCrystal(structure,
                                          (CrystalSystem)crystalSystemIndex,
                                          selectedSpaceGroup,
                                          latticeParams,
                                          asymmetricAtoms,
                                          elementColors);
            if (lastResult.success)
            {
                updateBuffers(structure);
                std::cout << "[Operation] Built bulk crystal: system="
                          << crystalSystemLabel((CrystalSystem)crystalSystemIndex)
                          << ", SG=" << lastResult.spaceGroupNumber
                          << " (" << lastResult.spaceGroupSymbol << ")"
                          << ", asym_atoms=" << asymmetricAtoms.size()
                          << ", generated_atoms=" << lastResult.generatedAtoms
                          << std::endl;
            }
            else
            {
                std::cout << "[Operation] Bulk crystal build failed: " << lastResult.message << std::endl;
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Close"))
            dialogOpen = false;

        if (!lastResult.message.empty())
        {
            if (lastResult.success)
                ImGui::Text("OK: %d atoms, SG %d (%s)", lastResult.generatedAtoms, lastResult.spaceGroupNumber, lastResult.spaceGroupSymbol.c_str());
            else
                ImGui::TextWrapped("Error: %s", lastResult.message.c_str());
        }
#endif

        ImGui::EndPopup();
    }

    if (!dialogOpen)
        ImGui::CloseCurrentPopup();

    if (showElementPicker)
    {
        ImGui::OpenPopup("Bulk Crystal Element Picker");
        showElementPicker = false;
    }

    ImGui::SetNextWindowSize(ImVec2(940.0f, 560.0f), ImGuiCond_Appearing);
    bool pickerOpen = true;
    if (ImGui::BeginPopupModal("Bulk Crystal Element Picker", &pickerOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select replacement element from the periodic table.");
        ImGui::Separator();
        drawPeriodicTableInlineSelector(selectedEditElement);
        ImGui::Separator();
        if (ImGui::Button("Apply Element##bulk-builder"))
        {
            if (targetAtomIndex >= 0 && targetAtomIndex < (int)asymmetricAtoms.size())
            {
                AtomSite& atom = asymmetricAtoms[targetAtomIndex];
                const int oldAtomicNumber = atom.atomicNumber;
                const std::string oldSymbol = atom.symbol;
                applyElementToAtom(atom, selectedEditElement, elementColors);
                std::cout << "[Operation] Substituted asymmetric atom (bulk builder): "
                          << "row=" << targetAtomIndex << " "
                          << oldSymbol << "(" << oldAtomicNumber << ") -> "
                          << atom.symbol << "(" << atom.atomicNumber << ")"
                          << std::endl;
            }
            targetAtomIndex = -1;
            if (reopenBulkDialogAfterPicker)
                m_openRequested = true;
            reopenBulkDialogAfterPicker = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }
    if (!pickerOpen)
    {
        targetAtomIndex = -1;
        if (reopenBulkDialogAfterPicker)
            m_openRequested = true;
        reopenBulkDialogAfterPicker = false;
        ImGui::CloseCurrentPopup();
    }
}