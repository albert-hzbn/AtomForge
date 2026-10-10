#pragma once

#include "pipeline/Pipeline.h"

#include <functional>
#include <string>

// Widgets for a modifier's parameters, built from its description (kinds int,
// float, bool, vector, ivec3, choice and text). Shared by the pipeline panel and
// the one-step operation dialogs; `changed(name, value)` receives each edit.
void drawModifierParameters(const atomforge::pipeline::ModifierType& type, const atomforge::pipeline::Json& parameters,
                            const std::function<void(const std::string&, const atomforge::pipeline::Json&)>& changed);
