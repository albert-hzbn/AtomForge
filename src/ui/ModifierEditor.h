#pragma once

#include "pipeline/Pipeline.h"

#include <functional>
#include <string>
#include <vector>

// Widgets for a modifier's parameters, built from its description (kinds int,
// float, bool, vector, ivec3, choice and text). Shared by the pipeline panel and
// the one-step operation dialogs; `changed(name, value)` receives each edit.
// Command-line options (kind "options", the CLI builders) get one labelled
// field per flag of the usage in the modifier's help, plus a field for others.
void drawModifierParameters(const atomforge::pipeline::ModifierType& type, const atomforge::pipeline::Json& parameters,
                            const std::function<void(const std::string&, const atomforge::pipeline::Json&)>& changed);

// The modifier's help for display: without the command-line usage lines when
// its options are shown as fields.
std::string modifierDescription(const atomforge::pipeline::ModifierType& type);

namespace modifierFlags
{
struct Flag
{
    std::string name;   // without "--"
    std::string value;  // placeholder or default ("" for a switch)
    bool optional = true;
    bool repeatable = false;
};
struct Usage
{
    std::vector<Flag> flags;  // --input and --output excluded
    std::string description;
};
Usage parse(const std::string& help);
// Options string <-> values per flag (switches: "1" when given; repeatable
// flags: entries joined by "; "), unknown tokens in `other`.
std::vector<std::string> values(const Usage& usage, const std::string& options, std::string& other);
std::string compose(const Usage& usage, const std::vector<std::string>& values, const std::string& other);
}
