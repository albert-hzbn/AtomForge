#pragma once

#include "StructureLoader.h"
#include "ui/TransformAtomsDialog.h"
#include "ui/EditMenuDialogs.h"

#include <array>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

struct FileBrowser
{
    FileBrowser();

    // Initialize browser state from a starting path.
    void initFromPath(const std::string& initialPath);

    // Draw the unified menu bar (File / Edit / About) and all related popups.
    // updateBuffers is called whenever a new structure is loaded or edited.
    void draw(Structure& structure,
              EditMenuDialogs& editMenuDialogs,
              const std::function<void(Structure&)>& updateBuffers);

    bool isTransformMatrixEnabled() const { return transformDialog.isEnabled(); }
    const int (&getTransformMatrix() const)[3][3] { return transformDialog.getMatrix(); }
    void clearTransformMatrix() { transformDialog.clearTransform(); }
    bool isShowElementEnabled() const { return showElementLabels; }

    // Programmatically trigger file open dialog (for keyboard shortcuts)
    void openFileDialog() { openStructurePopup = true; }
    // Trigger save-as dialog (e.g. from Ctrl+S shortcut)
    void saveFileDialog() { saveStructurePopup = true; }

private:
    void pushHistory(const std::string& dir);
    static std::string toLower(const std::string& s);
    bool isAllowedFile(const std::string& name) const;

    bool showAbout;
    bool showEditColors;
    bool showElementLabels;
    bool openStructurePopup;
    bool saveStructurePopup;

    std::string openDir;
    std::vector<std::string> dirHistory;
    int historyIndex;

    // Save-as dialog state
    std::string saveDir;
    std::vector<std::string> saveDirHistory;
    int saveHistoryIndex;
    int selectedSaveFormat;
    char saveFilename[1024];
    char saveStatusMsg[256];

    int selectedAtomicNumber;

    // Map of atomic number to user-adjusted color
    std::unordered_map<int, std::array<float, 3>> elementColorOverrides;

    std::vector<std::string> driveRoots;
    std::vector<std::string> allowedExtensions;

    char openFilename[1024];

    TransformAtomsDialog transformDialog;
};
