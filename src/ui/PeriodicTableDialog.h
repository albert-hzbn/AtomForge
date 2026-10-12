#pragma once

#include <string>
#include <vector>

struct ElementSelection {
    std::string symbol;
    int atomicNumber;
};

// Call once (per substitution request) to open the modal picker.
void openPeriodicTable();

// Draw the periodic-table picker modal.  Call every ImGui frame.
// Returns true exactly once, when the user clicks "Apply".
// On true: outSelections contains the single selected element.
// Click an element to replace selection, then click "Apply" to confirm.
bool drawPeriodicTable(std::vector<ElementSelection>& outSelections);

// Draw an inline periodic table selector (non-modal), useful for edit panels.
// Returns true when the selected element changes via click (or typing).
// `activated`, when given, is set on a double-click. `typing` lets typed
// symbols ("cu") pick an element; use it only where no text field is active.
bool drawPeriodicTableInlineSelector(int& selectedAtomicNumber, bool* activated = nullptr, bool typing = false);
