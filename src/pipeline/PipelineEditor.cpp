#include "pipeline/PipelineEditor.h"

#include <cstring>
#include <stdexcept>

namespace atomforge::pipeline
{
void PipelineEditor::setInput(const Structure& structure)
{
    m_input = structure;
    changed(0);
}

Structure PipelineEditor::bake()
{
    Structure out = active() ? evaluate().structure : Structure{};
    reset();
    return out;
}

void PipelineEditor::reset()
{
    m_input.reset();
    m_pipeline.modifiers.clear();
    m_cache.clear();
    m_output = {};
    m_outputPrint = 0;
    m_dirty = true;
}

void PipelineEditor::changed(std::size_t from)
{
    m_cache.invalidateFrom(from);
    m_dirty = true;
}

std::size_t PipelineEditor::add(const std::string& type, std::size_t index)
{
    index = std::min(index, size());
    m_pipeline.modifiers.insert(m_pipeline.modifiers.begin() + static_cast<std::ptrdiff_t>(index), makeModifier(type));
    changed(index);
    return index;
}

void PipelineEditor::remove(std::size_t index)
{
    if (index >= size()) throw std::out_of_range("No modifier " + std::to_string(index));
    m_pipeline.modifiers.erase(m_pipeline.modifiers.begin() + static_cast<std::ptrdiff_t>(index));
    changed(index);
}

void PipelineEditor::move(std::size_t from, std::size_t to)
{
    if (from >= size() || to >= size()) throw std::out_of_range("Modifier position out of range");
    if (from == to) return;
    Modifier moving = std::move(m_pipeline.modifiers[from]);
    m_pipeline.modifiers.erase(m_pipeline.modifiers.begin() + static_cast<std::ptrdiff_t>(from));
    m_pipeline.modifiers.insert(m_pipeline.modifiers.begin() + static_cast<std::ptrdiff_t>(to), std::move(moving));
    changed(std::min(from, to));
}

std::size_t PipelineEditor::duplicate(std::size_t index)
{
    if (index >= size()) throw std::out_of_range("No modifier " + std::to_string(index));
    m_pipeline.modifiers.insert(m_pipeline.modifiers.begin() + static_cast<std::ptrdiff_t>(index) + 1, m_pipeline.modifiers[index]);
    changed(index + 1);
    return index + 1;
}

void PipelineEditor::setEnabled(std::size_t index, bool enabled)
{
    Modifier& m = m_pipeline.modifiers.at(index);
    if (m.enabled == enabled) return;
    m.enabled = enabled;
    changed(index);
}

void PipelineEditor::setLabel(std::size_t index, const std::string& label)
{
    m_pipeline.modifiers.at(index).label = label;  // display only
}

void PipelineEditor::setParameter(std::size_t index, const std::string& name, const Json& value)
{
    Modifier& m = m_pipeline.modifiers.at(index);
    const ModifierType* type = findModifierType(m.type);
    bool known = false;
    for (const auto& p : type->parameters) known = known || name == p.name;
    if (!known) throw std::runtime_error(m.type + " has no parameter " + name);
    m.parameters[name] = value;
    changed(index);
}

void PipelineEditor::setPipeline(Pipeline pipeline)
{
    m_pipeline = std::move(pipeline);
    m_cache.clear();
    changed(0);
}

const PipelineData& PipelineEditor::evaluate()
{
    if (!active()) throw std::runtime_error("The pipeline has no input structure");
    if (m_dirty) {
        m_output = m_pipeline.evaluate(*m_input, m_cache);
        m_outputPrint = fingerprint(m_output.structure);
        m_dirty = false;
    }
    return m_output;
}

bool PipelineEditor::isCurrentOutput(const Structure& structure) const
{
    return !m_dirty && fingerprint(structure) == m_outputPrint;
}

std::uint64_t PipelineEditor::fingerprint(const Structure& structure)
{
    // FNV-1a over species, positions and cell.
    std::uint64_t hash = 1469598103934665603ull;
    const auto mix = [&](const void* data, std::size_t bytes) {
        const auto* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < bytes; ++i) { hash ^= p[i]; hash *= 1099511628211ull; }
    };
    const std::size_t count = structure.atoms.size();
    mix(&count, sizeof(count));
    for (const auto& a : structure.atoms) {
        mix(a.symbol.data(), a.symbol.size());
        mix(&a.x, sizeof(double)); mix(&a.y, sizeof(double)); mix(&a.z, sizeof(double));
    }
    mix(&structure.hasUnitCell, sizeof(bool));
    for (const auto& row : structure.cellVectors) mix(row.data(), sizeof(double) * 3);
    return hash;
}

Json PipelineEditor::toJson() const
{
    Json json = m_pipeline.toJson();
    json["active"] = active();
    return json;
}

void PipelineEditor::fromJson(const Json& json)
{
    setPipeline(Pipeline::fromJson(json));
}
}
