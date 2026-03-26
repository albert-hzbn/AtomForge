#pragma once

#include "StructureLoader.h"
#include "ui/BulkCrystalBuilderDialog.h"
#include "ui/CSLGrainBoundaryDialog.h"
#include "ui/NanoCrystalBuilderDialog.h"
#include "ui/InterfaceBuilderDialog.h"
#include "ui/CommonNeighbourAnalysis.h"
#include "ui/RadialDistributionAnalysis.h"
#include "ui/TransformAtomsDialog.h"
#include "ui/EditMenuDialogs.h"

#include <array>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

enum class ViewMode
{
    Isometric,
    Orthographic,
};

enum class ImageExportFormat
{
    Png,
    Jpg,
    Svg,
};

struct ImageExportRequest
{
    std::string outputPath;
    ImageExportFormat format = ImageExportFormat::Png;
    bool includeBackground = true;
};

struct LatticePlane
{
    int h = 1;
    int k = 0;
    int l = 0;
    float offset = 1.0f;
    float opacity = 0.20f;
    std::array<float, 3> color = {0.95f, 0.62f, 0.20f};
    bool visible = true;
};

struct FileBrowser
{
    FileBrowser();

    // Initialize browser state from a starting path.
    void initFromPath(const std::string& initialPath);

    // Draw the unified menu bar (File / Edit / About) and all related popups.
    // updateBuffers is called whenever a new structure is loaded or edited.
    void draw(Structure& structure,
              EditMenuDialogs& editMenuDialogs,
              const std::function<void(Structure&)>& updateBuffers,
              bool canUndo,
              bool canRedo);

    bool isTransformMatrixEnabled() const { return transformDialog.isEnabled(); }
    const int (&getTransformMatrix() const)[3][3] { return transformDialog.getMatrix(); }
    void clearTransformMatrix() { transformDialog.clearTransform(); }
    bool isShowElementEnabled() const { return showElementLabels; }
    bool isShowBondsEnabled() const { return showBonds; }
    bool isBondElementFilterEnabled() const { return bondElementFilterEnabled; }
    const std::array<bool, 119>& getBondElementFilterMask() const { return bondElementFilterMask; }
    bool isOrthographicViewEnabled() const { return viewMode == ViewMode::Orthographic; }
    bool isBoxSelectModeEnabled() const { return boxSelectMode; }
    bool isShowLatticePlanesEnabled() const { return showLatticePlanes; }
    const std::vector<LatticePlane>& getLatticePlanes() const { return latticePlanes; }
    void clearLatticePlanes() { latticePlanes.clear(); }
    bool consumeMeasureDistanceRequest()
    {
        bool requested = requestMeasureDistance;
        requestMeasureDistance = false;
        return requested;
    }
    bool consumeMeasureAngleRequest()
    {
        bool requested = requestMeasureAngle;
        requestMeasureAngle = false;
        return requested;
    }
    bool consumeAtomInfoRequest()
    {
        bool requested = requestAtomInfo;
        requestAtomInfo = false;
        return requested;
    }
    bool consumeViewAxisXRequest()
    {
        bool requested = requestViewAxisX;
        requestViewAxisX = false;
        return requested;
    }
    bool consumeViewAxisYRequest()
    {
        bool requested = requestViewAxisY;
        requestViewAxisY = false;
        return requested;
    }
    bool consumeViewAxisZRequest()
    {
        bool requested = requestViewAxisZ;
        requestViewAxisZ = false;
        return requested;
    }
    bool consumeResetDefaultViewRequest()
    {
        bool requested = requestResetDefaultView;
        requestResetDefaultView = false;
        return requested;
    }
    bool consumeStructureInfoRequest()
    {
        bool requested = requestStructureInfo;
        requestStructureInfo = false;
        return requested;
    }
    bool consumeUndoRequest()
    {
        bool requested = requestUndo;
        requestUndo = false;
        return requested;
    }
    bool consumeRedoRequest()
    {
        bool requested = requestRedo;
        requestRedo = false;
        return requested;
    }
    bool consumeCloseStructureRequest()
    {
        bool requested = requestCloseStructure;
        requestCloseStructure = false;
        return requested;
    }
    bool consumeImageExportRequest(ImageExportRequest& request)
    {
        if (!requestImageExport)
            return false;

        request = pendingImageExport;
        requestImageExport = false;
        return true;
    }

    // Programmatically trigger file open dialog (for keyboard shortcuts)
    void openFileDialog() { openStructurePopup = true; }
    // Trigger save-as dialog (e.g. from Ctrl+S shortcut)
    void saveFileDialog() { saveStructurePopup = true; }
    // Trigger image-export dialog (e.g. from Ctrl+Shift+S shortcut)
    void exportImageDialog() { exportImagePopup = true; }
    // Trigger unload of current structure (e.g. from Ctrl+W shortcut)
    void closeStructure() { requestCloseStructure = true; }

    // Apply persistent user element color overrides to a loaded structure.
    void applyElementColorOverrides(Structure& structure) const;

    // Initialise the nanocrystal dialog's own GL preview resources.
    void initNanoCrystalRenderResources(Renderer& renderer);
    bool isNanoCrystalDialogOpen() const;
    void feedDropToNanoCrystalDialog(const std::string& path);

    // CSL grain boundary dialog GL resources and drop routing.
    void initCSLGrainBoundaryRenderResources(Renderer& renderer);
    bool isCSLGrainBoundaryDialogOpen() const;
    void feedDropToCSLGrainBoundaryDialog(const std::string& path);

    // Interface builder dialog GL resources and drop routing.
    void initInterfaceBuilderRenderResources(Renderer& renderer);
    bool isInterfaceBuilderDialogOpen() const;
    void feedDropToInterfaceBuilderDialog(const std::string& path);

    // Show a modal error popup for file-load failures.
    void showLoadError(const std::string& message);

private:
    void updateBondElementFilterMask();
    void pushHistory(const std::string& dir);
    static std::string toLower(const std::string& s);
    bool isAllowedFile(const std::string& name) const;

    bool showAbout;
    bool showManual;
    bool showEditColors;
    bool showElementLabels;
    bool showBonds;
    bool showLatticePlanes;
    bool showLatticePlanesDialog;
    bool bondElementFilterEnabled;
    ViewMode viewMode;
    bool boxSelectMode;
    bool requestMeasureDistance;
    bool requestMeasureAngle;
    bool requestAtomInfo;
    bool requestViewAxisX;
    bool requestViewAxisY;
    bool requestViewAxisZ;
    bool requestResetDefaultView;
    bool requestStructureInfo;
    bool requestUndo;
    bool requestRedo;
    bool requestCloseStructure;
    bool requestImageExport;
    bool openStructurePopup;
    bool saveStructurePopup;
    bool exportImagePopup;
    bool loadErrorPopupRequested;

    std::string openDir;
    std::vector<std::string> dirHistory;
    int historyIndex;

    // Save-as dialog state
    std::string saveDir;
    std::vector<std::string> saveDirHistory;
    int saveHistoryIndex;
    int selectedSaveFormat;

    // Image-export dialog state
    std::string exportDir;
    std::vector<std::string> exportDirHistory;
    int exportHistoryIndex;
    int selectedExportFormat;
    bool exportIncludeBackground;
    ImageExportRequest pendingImageExport;

    char openStatusMsg[256];
    char loadErrorMsg[512];
    char saveFilename[1024];
    char saveStatusMsg[256];
    char exportFilename[1024];
    char exportStatusMsg[256];

    int selectedAtomicNumber;

    // Map of atomic number to user-adjusted color
    std::unordered_map<int, std::array<float, 3>> elementColorOverrides;

    std::vector<std::string> driveRoots;
    std::vector<std::string> allowedExtensions;

    char openFilename[1024];
    char bondElementFilterInput[256];
    std::array<bool, 119> bondElementFilterMask;
    int latticePlaneInputH;
    int latticePlaneInputK;
    int latticePlaneInputL;
    float latticePlaneInputOffset;
    float latticePlaneInputOpacity;
    float latticePlaneInputColor[3];
    std::vector<LatticePlane> latticePlanes;

    BulkCrystalBuilderDialog bulkCrystalDialog;
    CSLGrainBoundaryDialog cslDialog;
    NanoCrystalBuilderDialog nanoCrystalDialog;
    InterfaceBuilderDialog interfaceBuilderDialog;
    CommonNeighbourAnalysisDialog cnaDialog;
    RadialDistributionAnalysisDialog rdfDialog;
    TransformAtomsDialog transformDialog;
};
