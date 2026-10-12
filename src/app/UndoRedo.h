#pragma once

#include <cstdint>

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
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

// What changed from `before` to `after`, as a short phrase for the Undo and
// Redo commands ("delete 3 atoms").
inline std::string describeChange(const EditorSnapshot& before, const EditorSnapshot& after)
{
    const std::size_t a = before.structure.atoms.size(), b = after.structure.atoms.size();
    const auto atoms = [](std::size_t n) { return std::to_string(n) + (n == 1 ? " atom" : " atoms"); };
    if (a == 0 && b > 0) return "load the structure";
    if (a != b) {
        // Adding or deleting keeps the other atoms, in order; anything else
        // (another frame, a builder's result) replaces the structure.
        const auto& small = a < b ? before.structure.atoms : after.structure.atoms;
        const auto& large = a < b ? after.structure.atoms : before.structure.atoms;
        std::size_t matched = 0;
        for (std::size_t i = 0; i < large.size() && matched < small.size(); ++i) {
            const AtomSite& p = small[matched];
            const AtomSite& q = large[i];
            if (p.symbol == q.symbol && p.x == q.x && p.y == q.y && p.z == q.z) ++matched;
        }
        if (matched < small.size()) return "replace the structure (" + atoms(b) + ")";
        return b > a ? "add " + atoms(b - a) : "delete " + atoms(a - b);
    }
    std::size_t moved = 0, changed = 0;
    for (std::size_t i = 0; i < a; ++i) {
        const AtomSite& p = before.structure.atoms[i];
        const AtomSite& q = after.structure.atoms[i];
        if (p.x != q.x || p.y != q.y || p.z != q.z) ++moved;
        if (p.symbol != q.symbol || p.atomicNumber != q.atomicNumber) ++changed;
    }
    const bool cell = before.structure.hasUnitCell != after.structure.hasUnitCell ||
                      before.structure.cellVectors != after.structure.cellVectors;
    if (cell) return "change the cell";
    if (changed) return "change the element of " + atoms(changed);
    if (moved) return "move " + atoms(moved);
    return "change display settings";
}

class UndoRedoManager
{
public:
    void reset(const EditorSnapshot& initialSnapshot)
    {
        m_history.clear();
        m_history.resize(1);
        m_history[0] = initialSnapshot;
        m_ids.assign(1, ++m_nextId);
        m_descriptions.assign(1, std::string());
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

        // Nothing changed (e.g. a view update): no step to undo, and the
        // document stays as saved.
        if (m_history[m_index] == snapshot)
            return;

        if (m_index + 1 < m_history.size())
        {
            m_history.erase(m_history.begin() + (std::ptrdiff_t)m_index + 1, m_history.end());
            m_ids.erase(m_ids.begin() + (std::ptrdiff_t)m_index + 1, m_ids.end());
            m_descriptions.erase(m_descriptions.begin() + (std::ptrdiff_t)m_index + 1, m_descriptions.end());
        }

        // Display-only changes (colours, radii) are undoable but keep the id:
        // they are not saved in structure files, so the document is unchanged.
        const std::string description = describeChange(m_history.back(), snapshot);
        const bool displayOnly = description == "change display settings";
        m_descriptions.push_back(description);
        m_history.push_back(snapshot);
        m_ids.push_back(displayOnly ? m_ids.back() : ++m_nextId);
        if (m_history.size() > kMaxEntries)
        {
            m_history.erase(m_history.begin());
            m_ids.erase(m_ids.begin());
            m_descriptions.erase(m_descriptions.begin());
        }
        else
        {
            ++m_index;
        }

        if (m_index >= m_history.size())
            m_index = m_history.size() - 1;
    }

    bool canUndo() const { return !m_history.empty() && m_index > 0; }
    // Identifies the current state: unchanged by undo/redo back to it, so a
    // document is unmodified while this equals the id recorded when saved.
    // The history: entries 0..size()-1, the current one at index(); entry i
    // was reached by description(i) ("" for the first).
    std::size_t size() const { return m_history.size(); }
    std::size_t index() const { return m_index; }
    const std::string& description(std::size_t i) const { static const std::string none; return i < m_descriptions.size() ? m_descriptions[i] : none; }
    const EditorSnapshot& jumpTo(std::size_t i)
    {
        if (!m_history.empty()) m_index = std::min(i, m_history.size() - 1);
        return m_history[m_index];
    }

    // What Undo and Redo would change ("" when unavailable).
    std::string undoDescription() const { return canUndo() && m_index < m_descriptions.size() ? m_descriptions[m_index] : std::string(); }
    std::string redoDescription() const { return canRedo() && m_index + 1 < m_descriptions.size() ? m_descriptions[m_index + 1] : std::string(); }
    std::uint64_t currentId() const { return m_ids.empty() ? 0 : m_ids[std::min(m_index, m_ids.size() - 1)]; }
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
    std::vector<std::uint64_t> m_ids;  // one per history entry
    std::vector<std::string> m_descriptions;  // the change into each entry
    std::uint64_t m_nextId = 0;
    size_t m_index = 0;
};