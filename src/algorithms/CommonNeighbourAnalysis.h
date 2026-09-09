#pragma once

#include "model/Structure.h"
#include <map>
#include <string>
#include <vector>

namespace atomforge::analysis
{
struct Signature
{
    int common = 0;
    int bonds = 0;
    int chain = 0;

    bool operator<(const Signature& other) const
    {
        if (common != other.common) return common < other.common;
        if (bonds != other.bonds) return bonds < other.bonds;
        return chain < other.chain;
    }

    bool operator==(const Signature& other) const
    {
        return common == other.common && bonds == other.bonds && chain == other.chain;
    }
};

struct AtomRow
{
    int index = -1;
    int atomicNumber = 0;
    std::string symbol;
    int coordination = 0;
    Signature dominantSignature;
    int dominantSignatureCount = 0;
    std::string environment;
};

struct CnaResult
{
    bool valid = false;
    std::string message;

    int atomCount = 0;
    int pairCount = 0;
    bool pbcUsed = false;

    std::map<Signature, int> signatureCounts;
    std::map<std::string, int> environmentCounts;
    std::vector<AtomRow> atomRows;
};

struct CnaParams
{
    float cutoffScale = 1.18f;
    bool usePbc = true;
};

[[nodiscard]] CnaResult computeCna(const Structure& structure, const CnaParams& params = {});
} // namespace atomforge::analysis
