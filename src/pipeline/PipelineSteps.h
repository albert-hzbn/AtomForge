#pragma once
// The families of built-in modifier types, gathered by the registry in
// Modifiers.cpp. Adding a family: declare its list here, define it in its own
// file and append it in registry().

#include "pipeline/Pipeline.h"

#include <vector>

namespace atomforge::pipeline
{
// Per-atom analyses, statistics, binning, colouring and derived selections.
std::vector<ModifierType> analysisStepTypes();
// Runs a user's Python script on the structure.
ModifierType pythonScriptStepType();
}
