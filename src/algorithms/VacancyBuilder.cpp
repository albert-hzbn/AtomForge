#include "algorithms/VacancyBuilder.h"
#include "science/ScienceCore.h"

#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>
#include <random>
#include <vector>

namespace atomforge
{

namespace
{
double distance(const AtomSite& a, const AtomSite& b)
{
    return glm::length(glm::dvec3(a.x - b.x, a.y - b.y, a.z - b.z));
}
}

VacancyResult buildVacancies(const Structure& source, const VacancyParams& params)
{
    VacancyResult result;
    result.structure = source;

    std::vector<int> eligible;
    eligible.reserve(source.atoms.size());
    for (int i = 0; i < (int)source.atoms.size(); ++i)
        if (params.element.empty() || source.atoms[(size_t)i].symbol == params.element)
            eligible.push_back(i);
    result.eligibleCount = (int)eligible.size();

    if (eligible.empty())
    {
        result.message = params.element.empty()
            ? "Structure has no atoms."
            : "No atoms of element '" + params.element + "' found.";
        return result;
    }

    int requested = params.targetCount > 0
        ? params.targetCount
        : (int)std::lround(params.targetPercentage / 100.0 * (double)eligible.size());
    result.requestedCount = requested;
    if (requested < 0 || requested > (int)eligible.size())
    {
        result.message = "Requested " + std::to_string(requested) + " vacancies but only "
            + std::to_string(eligible.size()) + " atoms are eligible.";
        return result;
    }

    if (requested == 0)
    {
        result.success = true;
        result.message = "Requested zero vacancies; structure unchanged.";
        return result;
    }

    std::mt19937 rng(params.seed);
    std::shuffle(eligible.begin(), eligible.end(), rng);

    science::MinimumImage mic;
    if (source.hasUnitCell)
    {
        science::Mat3 cell{};
        for (int r = 0; r < 3; ++r) cell[r] = source.cellVectors[r];
        mic = science::MinimumImage(cell, {true, true, true});
    }
    std::vector<int> removedIndices;
    removedIndices.reserve((size_t)requested);
    for (int idx : eligible)
    {
        if ((int)removedIndices.size() >= requested)
            break;
        if (params.minSeparation > 0.0)
        {
            bool tooClose = false;
            for (int chosen : removedIndices)
            {
                const auto& p = source.atoms[(size_t)idx];
                const auto& q = source.atoms[(size_t)chosen];
                // Periodic images count: separation uses the minimum image.
                const double separation = source.hasUnitCell
                    ? science::norm(mic({p.x - q.x, p.y - q.y, p.z - q.z}))
                    : distance(p, q);
                if (separation < params.minSeparation)
                {
                    tooClose = true;
                    break;
                }
            }
            if (tooClose)
                continue;
        }
        removedIndices.push_back(idx);
    }

    std::vector<bool> remove(source.atoms.size(), false);
    for (int idx : removedIndices) remove[(size_t)idx] = true;

    Structure out = source;
    out.atoms.clear();
    out.grainColors.clear();
    out.grainRegionIds.clear();
    const bool keepGrainColors = source.grainColors.size() == source.atoms.size();
    const bool keepGrainRegions = source.grainRegionIds.size() == source.atoms.size();
    for (size_t i = 0; i < source.atoms.size(); ++i)
    {
        if (remove[i]) continue;
        out.atoms.push_back(source.atoms[i]);
        if (keepGrainColors) out.grainColors.push_back(source.grainColors[i]);
        if (keepGrainRegions) out.grainRegionIds.push_back(source.grainRegionIds[i]);
    }

    result.structure = out;
    result.removedCount = (int)removedIndices.size();
    result.success = true;
    result.message = "Removed " + std::to_string(result.removedCount) + " of "
        + std::to_string(requested) + " requested vacancies ("
        + std::to_string(eligible.size()) + " eligible sites)."
        + (result.removedCount < requested
               ? " Minimum separation prevented reaching the full request."
               : "");
    return result;
}

}
