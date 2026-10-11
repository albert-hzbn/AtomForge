// Undo history: no-op commits, change descriptions, and the id that tells
// whether a document differs from its saved state.
#include "app/UndoRedo.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void check(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

EditorSnapshot snapshot(int atoms)
{
    EditorSnapshot s;
    for (int i = 0; i < atoms; ++i) s.structure.atoms.push_back({"Cu", 29, double(i), 0.0, 0.0});
    return s;
}
}

int main()
{
    try {
        UndoRedoManager history;
        EditorSnapshot state = snapshot(4);
        history.reset(state);
        const auto saved = history.currentId();

        history.commit(state);
        check(!history.canUndo() && history.currentId() == saved, "An unchanged commit must not add an undo step");

        EditorSnapshot coloured = state;
        coloured.structure.atoms[0].r = 0.2f;
        history.commit(coloured);
        check(history.canUndo() && history.undoDescription() == "change display settings", "Colour change description");
        check(history.currentId() == saved, "A display-only change must leave the document unmodified");

        EditorSnapshot moved = coloured;
        moved.structure.atoms[1].x += 0.5;
        moved.structure.atoms[2].y += 0.5;
        history.commit(moved);
        check(history.undoDescription() == "move 2 atoms", "Move description: " + history.undoDescription());
        check(history.currentId() != saved, "Moving atoms modifies the document");

        EditorSnapshot substituted = moved;
        substituted.structure.atoms[3].symbol = "Ni";
        substituted.structure.atoms[3].atomicNumber = 28;
        history.commit(substituted);
        check(history.undoDescription() == "change the element of 1 atom", "Element description: " + history.undoDescription());

        EditorSnapshot fewer = substituted;
        fewer.structure.atoms.resize(1);
        history.commit(fewer);
        check(history.undoDescription() == "delete 3 atoms", "Delete description: " + history.undoDescription());

        history.undo();
        check(history.redoDescription() == "delete 3 atoms" && history.undoDescription() == "change the element of 1 atom",
              "Descriptions follow undo and redo");
        history.undo();
        history.undo();
        history.undo();
        check(history.currentId() == saved, "Undoing back to the saved state leaves the document unmodified");
        check(!history.canUndo(), "Back at the start of the history");

        std::cout << "Undo history checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
