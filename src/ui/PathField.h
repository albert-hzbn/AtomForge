#pragma once

#include <cstddef>

// The folder field of the file dialogs. A path too long for the field shows
// its end (".../Documents/AtomForge") while not being edited, since the
// end names the folder; clicking it edits the whole path. Returns true when
// Enter confirms an edit, as InputText with EnterReturnsTrue does.
bool pathField(const char* id, char* buffer, std::size_t size);
