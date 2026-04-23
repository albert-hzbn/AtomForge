#pragma once

#include "StructureLoader.h"
#include "ui/BulkCrystalBuilderDialog.h"
#include "ui/CSLGrainBoundaryDialog.h"
#include "ui/NanoCrystalBuilderDialog.h"
#include "ui/CustomStructureDialog.h"
#include "ui/MergeStructuresDialog.h"
#include "ui/InterfaceBuilderDialog.h"
#include "ui/PolyCrystalBuilderDialog.h"
#include "ui/StackingFaultBuilderDialog.h"
#include "ui/SubstitutionalSolidSolutionDialog.h"
#include "ui/CommonNeighbourAnalysis.h"
#include "ui/RadialDistributionAnalysis.h"
#include "ui/PolyhedralOverlay.h"
#include "ui/TransformAtomsDialog.h"
#include "ui/EditMenuDialogs.h"

#include <array>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

struct Renderer;

struct ToastNotification
{
    std::string message;
    double expiresAt = 0.0;
    bool isError = false;
};

enum class ViewMode
{
    Isometric,
    Orthographic,
};

enum class AtomColorMode
{
    ElementType = 0,
    CrystalOrientation,
    GrainBoundary,
};

enum class AtomDisplayMode
{
    Balls = 0,
    BallAndStick,
    SpaceFilling,
    Polyhedral,
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
    int resolutionScale = 1;
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

struct MillerDirection
{
    int u = 1;
    int v = 0;
    int w = 0;
    float length = 5.0f;  // display length in Angstroms
    std::array<float, 3> color = {0.25f, 0.75f, 1.0f};
    bool visible = true;
};

struct FileBrowser
{
    FileBrowser();

    // Initialize browser state from a starting path.
    void initFromPath(const std::string& initialPath);

    // Draw the unified menu bar (File / Edit / About) and all related popups.
    // updateBuffers is called whenever a new structure is loaded or edited.
    // openInNewTab is called by builder dialogs to place results in a new tab.
    void draw(Structure& structure,
              EditMenuDialogs& editMenuDialogs,
              const std::function<void(Structure&)>& updateBuffers,
              const std::function<void(Structure)>& openInNewTab,
              bool canUndo,
              bool canRedo);

    bool isTransformMatrixEnabled() const { return transformDialog.isEnabled(); }
    const int (&getTransformMatrix() const)[3][3] { return transformDialog.getMatrix(); }
    void clearTransformMatrix() { transformDialog.clearTransform(); }
    bool isShowElementEnabled() const { return showElementLabels; }
    bool isShowBondsEnabled() const
    {
        if (atomDisplayMode == AtomDisplayMode::BallAndStick || atomDisplayMode == AtomDisplayMode::Polyhedral)
            return true;
        return showBonds;
    }
    bool isShowAtomsEnabled() const { return showAtoms; }
    bool isShowBoundingBoxEnabled() const { return showBoundingBox; }
    bool isLightThemeEnabled() const { return useLightTheme; }
    bool isBondElementFilterEnabled() const { return bondElementFilterEnabled; }
    const std::array<bool, 119>& getBondElementFilterMask() const { return bondElementFilterMask; }
    bool isOrthographicViewEnabled() const { return viewMode == ViewMode::Orthographic; }
    bool isBoxSelectModeEnabled() const { return boxSelectMode; }
    bool isLassoSelectModeEnabled() const { return lassoSelectMode; }
    AtomColorMode getAtomColorMode() const { return atomColorMode; }
    bool atomColorModeChanged();
    void setAtomColorMode(AtomColorMode mode) { atomColorMode = mode; atomColorModeJustChanged = true; }
    AtomDisplayMode getAtomDisplayMode() const { return atomDisplayMode; }
    float getAtomRadiusScale() const
    {
        switch (atomDisplayMode)
        {
            case AtomDisplayMode::BallAndStick: return 1.0f;
            case AtomDisplayMode::SpaceFilling: return 1.8f;
            case AtomDisplayMode::Polyhedral:   return 0.7f;
            default:                            return 1.0f;
        }
    }
    bool isShowLatticePlanesEnabled() const { return showLatticePlanes; }
    const std::vector<LatticePlane>& getLatticePlanes() const { return latticePlanes; }
    void clearLatticePlanes() { latticePlanes.clear(); }
    bool isShowMillerDirectionsEnabled() const { return showMillerDirections; }
    const std::vector<MillerDirection>& getMillerDirections() const { return millerDirections; }
    void clearMillerDirections() { millerDirections.clear(); }
    bool isShowVoronoiEnabled() const { return showVoronoi; }
    bool isShowPolyhedralViewerEnabled() const { return showPolyhedralViewer || atomDisplayMode == AtomDisplayMode::Polyhedral; }
    const PolyhedralOverlaySettings& getPolyhedralOverlaySettings() const { return polyhedralSettings; }
    // Returns the path of the most recently loaded structure (empty if no new load since last call).
    std::string consumeLastLoadedPath()
    {
        std::string p = std::move(lastLoadedPath);
        lastLoadedPath.clear();
        return p;
    }
    // Returns a file path the user selected for opening (not yet loaded).
    std::string consumePendingOpenPath()
    {
        std::string p = std::move(pendingOpenPath);
        pendingOpenPath.clear();
        return p;
    }
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
    bool consumeViewLatticeARequest()
    {
        bool requested = requestViewLatticeA;
        requestViewLatticeA = false;
        return requested;
    }
    bool consumeViewLatticeBRequest()
    {
        bool requested = requestViewLatticeB;
        requestViewLatticeB = false;
        return requested;
    }
    bool consumeViewLatticeCRequest()
    {
        bool requested = requestViewLatticeC;
        requestViewLatticeC = false;
        return requested;
    }
    float getRotateCrystalAngle() const { return rotateCrystalAngle; }
    bool consumeRotateCrystalXRequest()
    {
        bool requested = requestRotateCrystalX;
        requestRotateCrystalX = false;
        return requested;
    }
    bool consumeRotateCrystalYRequest()
    {
        bool requested = requestRotateCrystalY;
        requestRotateCrystalY = false;
        return requested;
    }
    bool consumeRotateCrystalZRequest()
    {
        bool requested = requestRotateCrystalZ;
        requestRotateCrystalZ = false;
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

    // Custom structure dialog GL resources and drop routing.
    void initCustomStructureRenderResources(Renderer& renderer);
    bool isCustomStructureDialogOpen() const;
    void feedDropToCustomStructureDialog(const std::string& path);

    // Merge structures dialog resources and drop routing.
    void initMergeStructuresRenderResources(Renderer& renderer);
    bool isMergeStructuresDialogOpen() const;
    void feedDropToMergeStructuresDialog(const std::string& path);

    // CSL grain boundary dialog GL resources and drop routing.
    void initCSLGrainBoundaryRenderResources(Renderer& renderer);
    bool isCSLGrainBoundaryDialogOpen() const;
    void feedDropToCSLGrainBoundaryDialog(const std::string& path);

    // Interface builder dialog GL resources and drop routing.
    void initInterfaceBuilderRenderResources(Renderer& renderer);
    bool isInterfaceBuilderDialogOpen() const;
    void feedDropToInterfaceBuilderDialog(const std::string& path);

    // Polycrystal builder dialog GL resources and drop routing.
    void initPolyCrystalRenderResources(Renderer& renderer);
    bool isPolyCrystalDialogOpen() const;
    void feedDropToPolyCrystalDialog(const std::string& path);

    // Stacking faults builder dialog GL resources and drop routing.
    void initStackingFaultRenderResources(Renderer& renderer);
    bool isStackingFaultDialogOpen() const;
    void feedDropToStackingFaultDialog(const std::string& path);

    // Substitutional solid solution dialog GL resources and drop routing.
    void initSubstitutionalSolidSolutionRenderResources(Renderer& renderer);
    bool isSubstitutionalSolidSolutionDialogOpen() const;
    void feedDropToSubstitutionalSolidSolutionDialog(const std::string& path);

    // Returns true if any builder/analysis/file dialog is currently open.
    bool isAnyDialogOpen() const;

    // Show a modal error popup for file-load failures.
    void showLoadError(const std::string& message);
    // Show a modal informational popup for successful load status.
    void showLoadInfo(const std::string& message);
    // Show a transient notification toast.
    void showNotification(const std::string& message, bool isError = false);

private:
    void updateBondElementFilterMask();
    void updatePolyhedralElementFilterMask(const char* input, std::array<bool, 119>& mask);
    void updatePolyhedralCenterAtomIndexFilter(const char* input);
    void pushHistory(const std::string& dir);
    static std::string toLower(const std::string& s);
    bool isAllowedFile(const std::string& name) const;
    void drawNotifications();

    bool useLightTheme;
    bool showAbout;
    bool showManual;
    bool showEditColors;
    bool showElementLabels;
    bool showBonds;
    bool showAtoms;
    bool showBoundingBox;
    bool showLatticePlanes;
    bool showLatticePlanesDialog;
    bool showMillerDirections;
    bool showMillerDirectionsDialog;
    bool showVoronoi;
    bool showPolyhedralViewer;
    bool showPolyhedralSettingsDialog;
    bool bondElementFilterEnabled;
    std::string lastLoadedPath;
    ViewMode viewMode;
    AtomColorMode atomColorMode;
    bool atomColorModeJustChanged;
    AtomDisplayMode atomDisplayMode;
    bool boxSelectMode;
    bool lassoSelectMode;
    bool requestMeasureDistance;
    bool requestMeasureAngle;
    bool requestAtomInfo;
    bool requestViewAxisX;
    bool requestViewAxisY;
    bool requestViewAxisZ;
    bool requestViewLatticeA;
    float rotateCrystalAngle;
    bool requestRotateCrystalX;
    bool requestRotateCrystalY;
    bool requestRotateCrystalZ;
    bool requestViewLatticeB;
    bool requestViewLatticeC;
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
    std::string pendingOpenPath;

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
    int exportResolutionScale;
    ImageExportRequest pendingImageExport;

    char openStatusMsg[256];
    char loadPopupTitle[64];
    char loadErrorMsg[512];
    char saveFilename[1024];
    char saveStatusMsg[256];
    char exportFilename[1024];
    char exportStatusMsg[256];
    std::vector<ToastNotification> notifications;

    int selectedAtomicNumber;

    // Map of atomic number to user-adjusted color
    std::unordered_map<int, std::array<float, 3>> elementColorOverrides;

    std::vector<std::string> driveRoots;
    std::vector<std::string> allowedExtensions;

    char openFilename[1024];
    char bondElementFilterInput[256];
    char polyhedralCenterAtomIndexInput[256];
    char polyhedralCenterFilterInput[256];
    char polyhedralLigandFilterInput[256];
    std::array<bool, 119> bondElementFilterMask;
    PolyhedralOverlaySettings polyhedralSettings;
    int latticePlaneInputH;
    int latticePlaneInputK;
    int latticePlaneInputL;
    float latticePlaneInputOffset;
    float latticePlaneInputOpacity;
    float latticePlaneInputColor[3];
    std::vector<LatticePlane> latticePlanes;

    int millerDirInputU;
    int millerDirInputV;
    int millerDirInputW;
    float millerDirInputLength;
    float millerDirInputColor[3];
    std::vector<MillerDirection> millerDirections;

    BulkCrystalBuilderDialog bulkCrystalDialog;
    CSLGrainBoundaryDialog cslDialog;
    NanoCrystalBuilderDialog nanoCrystalDialog;
    CustomStructureDialog customStructureDialog;
    MergeStructuresDialog mergeStructuresDialog;
    InterfaceBuilderDialog interfaceBuilderDialog;
    PolyCrystalBuilderDialog polyCrystalDialog;
    StackingFaultBuilderDialog stackingFaultDialog;
    SubstitutionalSolidSolutionDialog substitutionalSolidSolutionDialog;
    CommonNeighbourAnalysisDialog cnaDialog;
    RadialDistributionAnalysisDialog rdfDialog;
    TransformAtomsDialog transformDialog;
};
