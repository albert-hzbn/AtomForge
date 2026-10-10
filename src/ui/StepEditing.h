#pragma once
// A Build or Edit dialog opened from a pipeline step edits the step, not the
// structure: it is filled from the step's parameters (and, for steps that take
// a structure, from the structure entering the step), and its main button
// ("Update step") hands the new parameters back to the pipeline.
#include "model/Structure.h"
#include "pipeline/Options.h"
#include "pipeline/Pipeline.h"

#include <functional>
#include <string>

struct StepEdit
{
    std::string step;                                                  // modifier type, e.g. "build-gb"
    atomforge::pipeline::Json parameters = atomforge::pipeline::Json::object();
    Structure input;                                                   // the structure entering the step
    std::function<void(const atomforge::pipeline::Json&)> commit;      // writes the step's new parameters

    bool active() const { return static_cast<bool>(commit); }
    // Build steps keep their command-line flags in "options".
    std::string options() const
    {
        return parameters.contains("options") && parameters.at("options").isString() ? parameters.at("options").string() : "";
    }
    atomforge::pipeline::options::Reader reader() const { return atomforge::pipeline::options::Reader(options()); }
    // Hands back new options (the step's other parameters are kept) or parameters.
    void commitOptions(const std::string& text)
    {
        atomforge::pipeline::Json updated = parameters;
        updated["options"] = text;
        commitParameters(updated);
    }
    void commitParameters(const atomforge::pipeline::Json& updated)
    {
        parameters = updated;
        if (commit) commit(updated);
    }
    // The dialog closed: it builds normally again.
    void finish() { commit = {}; }
};
