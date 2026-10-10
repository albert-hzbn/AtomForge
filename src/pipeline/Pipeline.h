#pragma once
// Non-destructive modification pipeline (in the spirit of OVITO): an input
// structure flows through an ordered list of modifiers, like a shell pipe
//     replicate 2 2 1 | select-expression "fz > 0.5" | delete-selected
// The input is never changed; editing, disabling, moving or deleting a
// modifier re-evaluates the output from the input. Modifiers act on the whole
// structure or on the current selection, which selection modifiers set.
//
// Adding a modifier: describe it (id, parameters with defaults) and give its
// apply function in modifierTypes() (Modifiers.cpp). The desktop panel, the
// text syntax and project files use that description; no other code changes.
#include "model/Structure.h"
#include "science/Json.h"

#include <functional>
#include <string>
#include <vector>

namespace atomforge::pipeline
{
using science::Json;

// What flows between modifiers.
struct PipelineData
{
    Structure structure;
    std::vector<char> selected;         // per atom, 1 = selected (same size as atoms)
    std::vector<std::string> notes;     // messages from the last modifier
    std::size_t selectedCount() const;
    void clearSelection() { selected.assign(structure.atoms.size(), 0); }
};

// Parameter kinds: int, float, bool, vector (3 numbers), ivec3 (3 integers),
// string, element (one symbol), elements (a list; takes the rest of its step in
// the text syntax), expression (likewise), options (command-line flags; likewise),
// choice (one of `options`, '|'-separated).
struct ModifierParameter
{
    const char* name;
    const char* label;
    const char* kind;
    const char* value;    // default, in the text syntax ("2 2 1", "true", "x > 0")
    const char* options;  // choice values, else ""
};

struct ModifierType
{
    const char* id;
    const char* title;
    const char* category;  // Selection, Modification, Cell or Analysis
    const char* help;
    std::vector<ModifierParameter> parameters;
    std::function<void(const Json& parameters, PipelineData& data)> apply;
};

const std::vector<ModifierType>& modifierTypes();
const ModifierType* findModifierType(const std::string& id);
// Adds (or replaces, by id) a modifier type. The desktop and command line
// register the Build and Edit operations this way; call it at start-up, before
// pipelines are evaluated.
void registerModifierType(ModifierType type);

// One step of a pipeline: its type, parameter values (JSON, keyed by name),
// whether it is applied and an optional display name.
struct Modifier
{
    std::string type;
    Json parameters = Json::object();
    bool enabled = true;
    std::string label;
};

// A modifier of `type` with default parameters.
Modifier makeModifier(const std::string& type);
// Parameter value or its default; throws for unknown names.
Json parameter(const Json& parameters, const ModifierType& type, const std::string& name);

// Outcome of one modifier during evaluation.
struct StageResult
{
    std::size_t atoms = 0;
    std::size_t selected = 0;
    double milliseconds = 0;
    std::string error;
    std::vector<std::string> notes;
    bool skipped = false;  // disabled, or after an error
};

// Outputs of each stage kept between evaluations: editing modifier i only
// re-runs modifiers i and later.
struct PipelineCache
{
    std::vector<PipelineData> outputs;  // after each modifier
    std::vector<StageResult> stages;
    std::size_t valid = 0;              // leading outputs that are still correct
    void invalidateFrom(std::size_t index) { valid = std::min(valid, index); }
    void clear() { outputs.clear(); stages.clear(); valid = 0; }
};

class Pipeline
{
public:
    std::vector<Modifier> modifiers;

    // Runs the enabled modifiers on `input`. A failing modifier stops the
    // pipeline: its stage records the error and the output is the data before it.
    PipelineData evaluate(const Structure& input, std::vector<StageResult>* stages = nullptr) const;
    PipelineData evaluate(const Structure& input, PipelineCache& cache) const;

    Json toJson() const;
    static Pipeline fromJson(const Json& json);
    // The shell-pipe syntax: modifiers separated by '|', each "id positional... key=value...";
    // vector parameters take three numbers, an expression takes the rest of its step.
    static Pipeline parse(const std::string& text);
    std::string toText() const;
};
}
