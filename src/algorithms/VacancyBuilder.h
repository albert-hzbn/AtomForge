#pragma once

#include "model/Structure.h"

#include <string>

// Point-defect tool: randomly removes atoms to create vacancies at a target
// concentration, the inverse operation of interstitial insertion
// (InterstitialVoidAnalysis.h) and substitutional doping
// (SubstitutionalSolidSolutionBuilder.h). Unlike those two, no BABEL/literature
// program is being matched here; this is a straightforward geometric utility.
namespace atomforge
{

struct VacancyParams
{
    // Chemical symbol restricting which atoms are eligible for removal.
    // Empty means every atom in the structure is eligible.
    std::string element;

    // Fraction of eligible atoms to remove, in percent (0-100). Ignored when
    // targetCount > 0.
    double targetPercentage = 0.0;

    // Exact number of atoms to remove. Takes precedence over
    // targetPercentage when positive.
    int targetCount = 0;

    // Minimum Cartesian separation (Angstrom) enforced between removed
    // sites, so vacancies do not cluster into a single large void unless
    // requested. 0 disables the check.
    double minSeparation = 0.0;

    unsigned seed = 1;
};

struct VacancyResult
{
    bool success = false;
    std::string message;
    Structure structure;

    int eligibleCount = 0;
    int requestedCount = 0;
    int removedCount = 0;
};

VacancyResult buildVacancies(const Structure& source, const VacancyParams& params);

}
