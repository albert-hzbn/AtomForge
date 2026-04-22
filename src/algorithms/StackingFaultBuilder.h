#pragma once

#include "io/StructureLoader.h"

#include <string>
#include <vector>

enum class StackingFaultFamily
{
    Unknown = 0,
    Fcc,
    Hcp,
    Bcc,
};

enum class StackingFaultPlane
{
    Auto = 0,
    Fcc111,
    HcpBasal,
    HcpPrismatic,
    HcpPyramidal,
    Bcc110,
    Bcc112,
};

enum class StackingFaultCellMode
{
    SmallestUnitCell = 0,
    OrthogonalCell,
};

struct StackingFaultSequenceItem
{
    Structure   structure;
    std::string label;
    float       displacementFactor = 0.0f;
};

struct StackingFaultDetectionResult
{
    bool                success = false;
    std::string         message;
    StackingFaultFamily family = StackingFaultFamily::Unknown;
    int                 fccCount = 0;
    int                 hcpCount = 0;
    int                 bccCount = 0;
    int                 recognizedCount = 0;
};

struct StackingFaultParams
{
    StackingFaultPlane    plane = StackingFaultPlane::Auto;
    StackingFaultCellMode cellMode = StackingFaultCellMode::SmallestUnitCell;
    int                   layerCount = 9;
    float                 interval = 0.1f;
    float                 maxDisplacementFactor = 2.0f;
    bool                  usePbcForDetection = true;
};

struct StackingFaultResult
{
    bool                             success = false;
    std::string                      message;
    StackingFaultFamily              family = StackingFaultFamily::Unknown;
    std::string                      detectedPhase;
    std::vector<StackingFaultSequenceItem> sequence;
    int                              detectedLayerCount = 0;
    int                              shiftedAtomCount = 0;
    float                            partialDisplacement = 0.0f;
};

const char* stackingFaultFamilyName(StackingFaultFamily family);

StackingFaultDetectionResult detectStackingFaultFamily(const Structure& structure,
                                                       bool usePbcForDetection = true);

StackingFaultResult buildStackingFaultSequence(const Structure& base,
                                               const StackingFaultParams& params);