#pragma once

#include "io/StructureLoader.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>

struct EditorSnapshot
{
    Structure structure;
    std::vector<float> elementRadii;
    std::vector<glm::vec3> elementColors;
    std::vector<float> elementShininess;
};

inline bool operator==(const EditorSnapshot& lhs, const EditorSnapshot& rhs)
{
    if (lhs.structure.hasUnitCell != rhs.structure.hasUnitCell)
        return false;
    if (lhs.structure.cellOffset != rhs.structure.cellOffset)
        return false;
    if (lhs.structure.cellVectors != rhs.structure.cellVectors)
        return false;
    if (lhs.structure.atoms.size() != rhs.structure.atoms.size())
        return false;

    for (size_t index = 0; index < lhs.structure.atoms.size(); ++index)
    {
        const AtomSite& a = lhs.structure.atoms[index];
        const AtomSite& b = rhs.structure.atoms[index];
        if (a.symbol != b.symbol ||
            a.atomicNumber != b.atomicNumber ||
            a.x != b.x || a.y != b.y || a.z != b.z ||
            a.r != b.r || a.g != b.g || a.b != b.b)
            return false;
    }

    if (lhs.elementRadii != rhs.elementRadii)
        return false;
    if (lhs.elementColors.size() != rhs.elementColors.size())
        return false;

    for (size_t index = 0; index < lhs.elementColors.size(); ++index)
    {
        const glm::vec3& a = lhs.elementColors[index];
        const glm::vec3& b = rhs.elementColors[index];
        if (a.x != b.x || a.y != b.y || a.z != b.z)
            return false;
    }

    if (lhs.elementShininess != rhs.elementShininess)
        return false;

    return true;
}

class UndoRedoManager
{
public:
    void reset(const EditorSnapshot& initialSnapshot)
    {
        m_history.clear();
        m_history.resize(1);
        m_history[0] = initialSnapshot;
        m_index = 0;
    }

    void commit(const EditorSnapshot& snapshot)
    {
        if (m_history.empty())
        {
            reset(snapshot);
            return;
        }

        if (m_index >= m_history.size())
        {
            reset(snapshot);
            return;
        }

        if (m_index + 1 < m_history.size())
            m_history.erase(m_history.begin() + (std::ptrdiff_t)m_index + 1, m_history.end());

        m_history.push_back(snapshot);
        if (m_history.size() > kMaxEntries)
        {
            m_history.erase(m_history.begin());
        }
        else
        {
            ++m_index;
        }

        if (m_index >= m_history.size())
            m_index = m_history.size() - 1;
    }

    bool canUndo() const { return !m_history.empty() && m_index > 0; }
    bool canRedo() const { return !m_history.empty() && m_index + 1 < m_history.size(); }

    const EditorSnapshot& undo()
    {
        if (canUndo())
            --m_index;
        return m_history[m_index];
    }

    const EditorSnapshot& redo()
    {
        if (canRedo())
            ++m_index;
        return m_history[m_index];
    }

private:
    static const size_t kMaxEntries = 128;

    std::vector<EditorSnapshot> m_history;
    size_t m_index = 0;
};