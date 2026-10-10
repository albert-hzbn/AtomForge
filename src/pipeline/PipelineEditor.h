#pragma once
// Editing state behind the desktop pipeline panel, without any GUI: the
// captured input structure, the modifier list and the per-stage cache. Each
// edit invalidates exactly the stages it affects.
#include "pipeline/Pipeline.h"

#include <cstdint>
#include <optional>

namespace atomforge::pipeline
{
class PipelineEditor
{
public:
    bool active() const { return m_input.has_value(); }
    // Starts (or restarts) the pipeline on a copy of `structure`; the modifiers are kept.
    void setInput(const Structure& structure);
    const Structure& input() const { return *m_input; }
    // Ends the pipeline; returns its last output (the structure to keep).
    Structure bake();
    void reset();

    const Pipeline& pipeline() const { return m_pipeline; }
    std::size_t size() const { return m_pipeline.modifiers.size(); }
    const Modifier& modifier(std::size_t index) const { return m_pipeline.modifiers.at(index); }

    // Inserts a modifier of `type` with default parameters at `index` (end when
    // past the end); returns its position.
    std::size_t add(const std::string& type, std::size_t index = static_cast<std::size_t>(-1));
    void remove(std::size_t index);
    // Moves the modifier at `from` so that it ends up at position `to`.
    void move(std::size_t from, std::size_t to);
    std::size_t duplicate(std::size_t index);
    void setEnabled(std::size_t index, bool enabled);
    void setLabel(std::size_t index, const std::string& label);
    void setParameter(std::size_t index, const std::string& name, const Json& value);
    // Replaces the whole modifier list (text syntax or loaded JSON).
    void setPipeline(Pipeline pipeline);

    // Re-runs the stages that changed since the last call; returns the output.
    const PipelineData& evaluate();
    const std::vector<StageResult>& stages() const { return m_cache.stages; }
    bool dirty() const { return m_dirty; }
    // Output of stage `index` (the data after that modifier), once evaluated.
    const PipelineData& stageOutput(std::size_t index) const { return m_cache.outputs.at(index); }

    // Whether `structure` is still the last output (it differs after edits made
    // outside the pipeline, which the next evaluation would overwrite).
    bool isCurrentOutput(const Structure& structure) const;
    static std::uint64_t fingerprint(const Structure& structure);

    // Pipeline and options for project files (the input is stored separately).
    Json toJson() const;
    void fromJson(const Json& json);

private:
    void changed(std::size_t from);
    std::optional<Structure> m_input;
    Pipeline m_pipeline;
    PipelineCache m_cache;
    PipelineData m_output;
    std::uint64_t m_outputPrint = 0;
    bool m_dirty = true;
};
}
