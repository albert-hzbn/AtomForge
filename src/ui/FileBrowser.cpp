#include "ui/ResponsiveLayout.h"
#include "ui/UiScript.h"
#include "FileBrowser.h"
#include "io/StructureLoader.h"
#include "ElementData.h"
#include "ImGuiSetup.h"
#include "util/PathUtils.h"
#include "util/ApplicationPaths.h"
#include "app/StructureFileService.h"
#include "ui/PeriodicTableDialog.h"
#include "ui/DialogLayout.h"

#include "imgui.h"
#include <GLFW/glfw3.h>
#include <openbabel3/openbabel/elements.h>

#include <algorithm>
#include <filesystem>
#include <ctime>
#include <chrono>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#ifndef ATOMFORGE_VERSION
#define ATOMFORGE_VERSION "0.2.0"
#endif

// ---------------------------------------------------------------------------
// Save-format table (used by the Save As dialog)
// ---------------------------------------------------------------------------
struct SaveFormat { const char* label; const char* ext; const char* fmt; };
static const SaveFormat kSaveFormats[] = {
    { "XYZ (.xyz)",               ".xyz",  "xyz"   },
    { "CIF (.cif)",               ".cif",  "cif"   },
    { "VASP POSCAR (.vasp)",      ".vasp", "vasp"  },
    { "PDB (.pdb)",               ".pdb",  "pdb"   },
    { "SDF (.sdf)",               ".sdf",  "sdf"   },
    { "Mol2 (.mol2)",             ".mol2", "mol2"  },
    { "Quantum ESPRESSO (.pwi)",  ".pwi",  "pwscf" },
    { "Gaussian Input (.gjf)",    ".gjf",  "gjf"   },
};
static constexpr int kNumSaveFormats = (int)(sizeof(kSaveFormats) / sizeof(kSaveFormats[0]));
static constexpr int kDefaultSaveFormatIndex = 1; // CIF (.cif)

struct ImageExportFormatOption { const char* label; const char* ext; ImageExportFormat fmt; };
static const ImageExportFormatOption kImageExportFormats[] = {
    { "PNG (.png)", ".png", ImageExportFormat::Png },
    { "JPEG (.jpg)", ".jpg", ImageExportFormat::Jpg },
    { "SVG (.svg)", ".svg", ImageExportFormat::Svg },
};
static constexpr int kNumImageExportFormats = (int)(sizeof(kImageExportFormats) / sizeof(kImageExportFormats[0]));
static constexpr double kNotificationLifetimeSeconds = 3.5;
static constexpr std::size_t kMaxRecentFiles = 5;

namespace
{
std::vector<std::string>& sharedRecentFiles()
{
    static std::vector<std::string> files;
    return files;
}

bool& sharedRecentFilesLoaded()
{
    static bool loaded = false;
    return loaded;
}

std::string recentFilesStorePath()
{
    return joinPath(detectHomePath(), ".atomforge_recent_files.txt");
}

void loadSharedRecentFilesIfNeeded()
{
    if (sharedRecentFilesLoaded())
        return;

    sharedRecentFilesLoaded() = true;

    std::ifstream in(recentFilesStorePath());
    if (!in)
        return;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        const std::string normalized = normalizePathSeparators(line);
        auto& files = sharedRecentFiles();
        if (std::find(files.begin(), files.end(), normalized) != files.end())
            continue;

        files.push_back(normalized);
        if (files.size() >= kMaxRecentFiles)
            break;
    }
}

void saveSharedRecentFiles()
{
    std::ofstream out(recentFilesStorePath(), std::ios::trunc);
    if (!out)
        return;

    for (const std::string& path : sharedRecentFiles())
        out << path << '\n';
}
}

// "CIF structure, 42.1 KB, modified 2026-09-09 14:02" ("" if not a file).
std::string describeFile(const std::string& path)
{
    std::error_code error;
    const auto file = std::filesystem::u8path(path);
    if (!std::filesystem::is_regular_file(file, error)) return {};
    std::string extension = file.extension().u8string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const std::string name = file.filename().u8string();
    std::string kind = extension.size() > 1 ? extension.substr(1) : name;
    std::transform(kind.begin(), kind.end(), kind.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (extension == ".afcloud") kind = "Large dataset (opens in the large-data viewer)";
    else if (name.rfind("POSCAR", 0) == 0 || name.rfind("CONTCAR", 0) == 0 || extension == ".vasp") kind = "VASP structure";
    else kind += " structure";
    const double size = static_cast<double>(std::filesystem::file_size(file, error));
    char sizeText[32];
    if (size < 1024.0) std::snprintf(sizeText, sizeof(sizeText), "%.0f bytes", size);
    else if (size < 1024.0 * 1024.0) std::snprintf(sizeText, sizeof(sizeText), "%.1f KB", size / 1024.0);
    else if (size < 1024.0 * 1024.0 * 1024.0) std::snprintf(sizeText, sizeof(sizeText), "%.1f MB", size / (1024.0 * 1024.0));
    else std::snprintf(sizeText, sizeof(sizeText), "%.2f GB", size / (1024.0 * 1024.0 * 1024.0));
    std::string text = kind + ", " + sizeText;
    const auto written = std::filesystem::last_write_time(file, error);
    if (!error) {
        const auto system = std::chrono::time_point_cast<std::chrono::system_clock::duration>(
            written - std::filesystem::file_time_type::clock::now() + std::chrono::system_clock::now());
        const std::time_t time = std::chrono::system_clock::to_time_t(system);
        if (const std::tm* local = std::localtime(&time)) {
            char date[32];
            std::strftime(date, sizeof(date), "%Y-%m-%d %H:%M", local);
            text += std::string(", modified ") + date;
        }
    }
    return text;
}

// Persistent sidebar width shared across all three file-browser dialogs.
static float s_sidebarW = 130.0f;
 
FileBrowser::FileBrowser()
        : showAbout(false),
            showManual(false),
            showEditColors(false),
            showElementLabels(false),
            showBonds(false),
            showAtoms(true),
            showBoundingBox(true),
            showDislocationLines(false),
            showLatticePlanes(false),
            showLatticePlanesDialog(false),
            showMillerDirections(false),
            showMillerDirectionsDialog(false),
            showVoronoi(false),
            showPolyhedralViewer(false),
            showPolyhedralSettingsDialog(false),
            bondElementFilterEnabled(false),
            viewMode(ViewMode::Orthographic),
            atomColorMode(AtomColorMode::ElementType),
            atomColorModeJustChanged(false),
            atomDisplayMode(AtomDisplayMode::Balls),
            boxSelectMode(false),
            lassoSelectMode(false),
            requestMeasureDistance(false),
            requestMeasureAngle(false),
            requestAtomInfo(false),
            requestViewAxisX(false),
            requestViewAxisY(false),
            requestViewAxisZ(false),
            requestViewLatticeA(false),
            rotateCrystalAngle(5.0f),
            requestRotateCrystalX(false),
            requestRotateCrystalY(false),
            requestRotateCrystalZ(false),
            requestViewLatticeB(false),
            requestViewLatticeC(false),
            requestResetDefaultView(false),
            requestStructureInfo(false),
            requestUndo(false),
            requestRedo(false),
            requestCloseStructure(false),
            requestSaveFile(false),
            requestImageExport(false),
            openStructurePopup(false),
            saveStructurePopup(false),
            exportImagePopup(false),
            loadErrorPopupRequested(false),
            openDir("."),
      historyIndex(-1),
      saveDir("."),
      saveHistoryIndex(-1),
    selectedSaveFormat(kDefaultSaveFormatIndex),
        exportDir("."),
        exportHistoryIndex(-1),
        selectedExportFormat(0),
        exportIncludeBackground(true),
        exportIncludeGizmo(false),
        exportResolutionScale(1),
            selectedAtomicNumber(1),
            latticePlaneInputH(1),
            latticePlaneInputK(0),
            latticePlaneInputL(0),
            latticePlaneInputOffset(1.0f),
            latticePlaneInputOpacity(0.20f),
            millerDirInputU(1),
            millerDirInputV(0),
            millerDirInputW(0),
            millerDirInputLength(5.0f)
{
    allowedExtensions = {".cif", ".mol", ".pdb", ".xyz", ".sdf", ".vasp", ".mol2", ".pwi", ".gjf"};

    const std::string homePath = detectHomePath();

#ifdef _WIN32
    const char* systemDrive = std::getenv("SystemDrive");
    appendUniquePath(driveRoots, systemDrive ? std::string(systemDrive) + "/" : "C:/");
    appendUniquePath(driveRoots, homePath);

    DWORD driveMask = GetLogicalDrives();
    if (driveMask != 0)
    {
        for (int i = 0; i < 26; ++i)
        {
            if ((driveMask & (1u << i)) == 0)
                continue;

            std::string drive;
            drive.push_back((char)('A' + i));
            drive += ":/";
            appendUniquePath(driveRoots, drive);
        }
    }
#else
    appendUniquePath(driveRoots, "/");
    appendUniquePath(driveRoots, homePath);
    appendUniquePath(driveRoots, "/mnt");
    appendUniquePath(driveRoots, "/media");
#endif

    if (driveRoots.empty())
        appendUniquePath(driveRoots, ".");

    openFilename[0] = '\0';
    openStatusMsg[0] = '\0';
    std::snprintf(loadPopupTitle, sizeof(loadPopupTitle), "Load Error");
    loadErrorMsg[0] = '\0';
    saveFilename[0] = '\0';
    saveStatusMsg[0] = '\0';
    std::snprintf(exportFilename, sizeof(exportFilename), "%s", "structure.png");
    exportStatusMsg[0] = '\0';
    std::snprintf(bondElementFilterInput, sizeof(bondElementFilterInput), "%s", "O,F");
    polyhedralCenterAtomIndexInput[0] = '\0';
    polyhedralCenterFilterInput[0] = '\0';
    polyhedralLigandFilterInput[0] = '\0';
    latticePlaneInputColor[0] = 0.95f;
    latticePlaneInputColor[1] = 0.62f;
    latticePlaneInputColor[2] = 0.20f;
    millerDirInputColor[0] = 0.25f;
    millerDirInputColor[1] = 0.75f;
    millerDirInputColor[2] = 1.00f;
    bondElementFilterMask.fill(false);
    polyhedralSettings.centerElementMask.fill(false);
    polyhedralSettings.ligandElementMask.fill(false);
    updateBondElementFilterMask();
    loadSharedRecentFilesIfNeeded();
}

void FileBrowser::rememberRecentFile(const std::string& path)
{
    if (path.empty())
        return;

    loadSharedRecentFilesIfNeeded();

    const std::string normalized = normalizePathSeparators(path);
    auto& files = sharedRecentFiles();
    files.erase(
        std::remove(files.begin(), files.end(), normalized),
        files.end());
    files.insert(files.begin(), normalized);

    if (files.size() > kMaxRecentFiles)
        files.resize(kMaxRecentFiles);

    saveSharedRecentFiles();
}

void FileBrowser::updateBondElementFilterMask()
{
    bondElementFilterMask.fill(false);

    std::string token;
    const size_t len = std::strlen(bondElementFilterInput);
    for (size_t i = 0; i <= len; ++i)
    {
        const char ch = (i < len) ? bondElementFilterInput[i] : ',';
        if (ch == ',' || std::isspace((unsigned char)ch))
        {
            if (!token.empty())
            {
                const int z = atomicNumberFromSymbol(token);
                if (z >= 1 && z <= 118)
                    bondElementFilterMask[(size_t)z] = true;
                token.clear();
            }
            continue;
        }

        token.push_back(ch);
    }
}

void FileBrowser::updatePolyhedralElementFilterMask(const char* input, std::array<bool, 119>& mask)
{
    mask.fill(false);
    if (input == nullptr)
        return;

    std::string token;
    const size_t len = std::strlen(input);
    for (size_t i = 0; i <= len; ++i)
    {
        const char ch = (i < len) ? input[i] : ',';
        if (ch == ',' || std::isspace((unsigned char)ch))
        {
            if (!token.empty())
            {
                const int z = atomicNumberFromSymbol(token);
                if (z >= 1 && z <= 118)
                    mask[(size_t)z] = true;
                token.clear();
            }
            continue;
        }

        token.push_back(ch);
    }
}

void FileBrowser::updatePolyhedralCenterAtomIndexFilter(const char* input)
{
    polyhedralSettings.centerAtomIndices.clear();
    if (input == nullptr)
        return;

    std::string token;
    const size_t len = std::strlen(input);

    auto appendIndex = [&](int oneBased) {
        if (oneBased <= 0)
            return;
        const int zeroBased = oneBased - 1;
        if (std::find(polyhedralSettings.centerAtomIndices.begin(),
                      polyhedralSettings.centerAtomIndices.end(),
                      zeroBased) == polyhedralSettings.centerAtomIndices.end())
        {
            polyhedralSettings.centerAtomIndices.push_back(zeroBased);
        }
    };

    auto parseToken = [&](const std::string& value) {
        if (value.empty())
            return;

        const std::size_t dashPos = value.find('-');
        if (dashPos == std::string::npos)
        {
            appendIndex(std::atoi(value.c_str()));
            return;
        }

        const std::string left = value.substr(0, dashPos);
        const std::string right = value.substr(dashPos + 1);
        if (left.empty() || right.empty())
            return;

        int a = std::atoi(left.c_str());
        int b = std::atoi(right.c_str());
        if (a <= 0 || b <= 0)
            return;
        if (a > b)
            std::swap(a, b);

        const int count = std::min(b - a + 1, 20000);
        for (int i = 0; i < count; ++i)
            appendIndex(a + i);
    };

    for (size_t i = 0; i <= len; ++i)
    {
        const char ch = (i < len) ? input[i] : ',';
        if (ch == ',')
        {
            parseToken(token);
            token.clear();
            continue;
        }

        if (!std::isspace((unsigned char)ch))
            token.push_back(ch);
    }
}

void FileBrowser::initFromPath(const std::string& initialPath)
{
    saveFilename[0] = '\0';

    openDir = ".";
    auto pos = initialPath.find_last_of("/\\");
    if (pos != std::string::npos)
        openDir = normalizePathSeparators(initialPath.substr(0, pos));

    if (isDriveRootPath(openDir) && openDir.size() == 2)
        openDir += "/";

    if (openDir.empty())
        openDir = ".";

    pushHistory(openDir);

    std::string initialOpenFilename = initialPath;
    if (pos != std::string::npos)
        initialOpenFilename = initialPath.substr(pos + 1);

    std::snprintf(openFilename, sizeof(openFilename), "%s", initialOpenFilename.c_str());
    currentStructurePath = initialPath;
    if (!initialPath.empty())
    {
        rememberRecentFile(initialPath);
        lastLoadedPath = initialPath;
    }
}

void FileBrowser::triggerSaveAsDialog()
{
    saveStructurePopup = true;
    saveDir = openDir.empty() ? "." : openDir;
    saveDirHistory = dirHistory;
    saveHistoryIndex = historyIndex;
    if (saveDirHistory.empty())
    {
        saveDirHistory.push_back(saveDir);
        saveHistoryIndex = 0;
    }
    saveStatusMsg[0] = '\0';
}

void FileBrowser::performQuickSave(const Structure& structure)
{
    if (structure.atoms.empty())
    {
        showNotification("Error: no atoms to save.", true);
        return;
    }

    if (currentStructurePath.empty())
    {
        triggerSaveAsDialog();
        return;
    }

    std::string fullPath = normalizePathSeparators(currentStructurePath);
    const std::string lowerPath = toLower(fullPath);
    int formatIndex = -1;
    for (int i = 0; i < kNumSaveFormats; ++i)
    {
        if (hasExtension(lowerPath, toLower(kSaveFormats[i].ext)))
        {
            formatIndex = i;
            break;
        }
    }

    if (formatIndex < 0)
    {
        formatIndex = kDefaultSaveFormatIndex;
        fullPath = replaceFileExtension(fullPath, kSaveFormats[formatIndex].ext, "structure");
    }

    std::size_t savedAtomCount = 0;
    const bool ok = saveStructureWithOptionalSupercell(
        structure,
        isTransformMatrixEnabled(),
        getTransformMatrix(),
        fullPath,
        kSaveFormats[formatIndex].fmt,
        savedAtomCount);

    if (ok)
    {
        currentStructurePath = fullPath;
        lastLoadedPath = fullPath;

        const std::size_t slash = fullPath.find_last_of("/\\");
        const std::string shownName = (slash == std::string::npos)
            ? fullPath
            : fullPath.substr(slash + 1);

        std::ostringstream msg;
        msg << "Structure saved: " << shownName << " (" << savedAtomCount << " atoms)";
        showNotification(msg.str(), false);
        std::cout << "[Operation] Saved structure: " << fullPath
                  << " (format=" << kSaveFormats[formatIndex].fmt
                  << ", atoms=" << savedAtomCount << ")" << std::endl;
    }
    else
    {
        showNotification("Error: failed to save (format may not support this structure).", true);
        std::cout << "[Operation] Save failed: " << fullPath
                  << " (format=" << kSaveFormats[formatIndex].fmt << ")" << std::endl;
    }
}

void FileBrowser::drawMainMenuBar(Structure& structure,
                                  EditMenuDialogs& editMenuDialogs,
                                  const std::function<void(Structure&)>& updateBuffers,
                                  bool canUndo,
                                  bool canRedo)
{
    if (!ImGui::BeginMainMenuBar())
        return;

    forceTourMenu("File");
    if (ImGui::BeginMenu("File"))
    {
        if (ImGui::MenuItem("Open",  "Ctrl+O"))
            openStructurePopup = true;

        loadSharedRecentFilesIfNeeded();
        const auto& recent = sharedRecentFiles();
        forceTourSubmenu("Open Recent");
        if (ImGui::BeginMenu("Open Recent", !recent.empty()))
        {
            for (const std::string& recentPath : recent)
            {
                // The name, with its folder (shortened from the left) muted beside it.
                const std::size_t slash = recentPath.find_last_of("/\\");
                const std::string name = slash == std::string::npos ? recentPath : recentPath.substr(slash + 1);
                std::string folder = slash == std::string::npos ? std::string() : recentPath.substr(0, slash + 1);
                if (folder.size() > 48) folder = "\xe2\x80\xa6" + folder.substr(folder.size() - 47);
                std::error_code missing;
                const bool exists = std::filesystem::is_regular_file(std::filesystem::u8path(recentPath), missing);
                ImGui::PushID(recentPath.c_str());
                const bool chosen = ImGui::MenuItem(name.c_str(), folder.c_str(), false, exists);
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                    ImGui::SetTooltip("%s%s", recentPath.c_str(), exists ? "" : "\n(the file is no longer there)");
                ImGui::PopID();
                if (chosen)
                {
                    pendingOpenPath = recentPath;
                    openStatusMsg[0] = '\0';
                }
            }
            ImGui::EndMenu();
        }

        if (ImGui::MenuItem("Close", "Ctrl+W", false, !structure.atoms.empty()))
        {
            requestCloseStructure = true;
            std::cout << "[Operation] Close structure requested" << std::endl;
        }

        if (ImGui::MenuItem("Save", "Ctrl+S", false, !structure.atoms.empty()))
            requestSaveFile = true;

        if (ImGui::MenuItem("Save As", "Ctrl+Shift+S", false, !structure.atoms.empty()))
            triggerSaveAsDialog();

        if (ImGui::MenuItem("Export Image", "Ctrl+Alt+S", false, !structure.atoms.empty()))
        {
            exportImagePopup = true;
            exportDir = openDir;
            exportDirHistory = dirHistory;
            exportHistoryIndex = historyIndex;

            if (exportFilename[0] == '\0')
                std::snprintf(exportFilename, sizeof(exportFilename), "%s", "structure.png");

            exportStatusMsg[0] = '\0';
        }

        ImGui::Separator();

        if (ImGui::MenuItem("Open Project")) { projectAction=1; projectPicker.open("Open AtomForge project",false,"project.afproject"); }
        if (ImGui::MenuItem("Save Project As")) { projectAction=2; projectPicker.open("Save AtomForge project",true,"project.afproject"); }
        ImGui::MenuItem("Autosave project every minute",nullptr,&autosaveEnabled);
        if (ImGui::MenuItem("Recover Last Autosave")) projectAction=3;
        ImGui::Separator();
        if (ImGui::MenuItem("Quit"))
            requestQuit = true;

        ImGui::EndMenu();
    }

    forceTourMenu("Edit");
    if (ImGui::BeginMenu("Edit"))
    {
        const std::string undoLabel = (undoDescription.empty() ? std::string("Undo") : "Undo " + undoDescription) + "###undo";
        const std::string redoLabel = (redoDescription.empty() ? std::string("Redo") : "Redo " + redoDescription) + "###redo";
        if (ImGui::MenuItem(undoLabel.c_str(), "Ctrl+Z", false, canUndo))
            requestUndo = true;
        if (ImGui::MenuItem(redoLabel.c_str(), "Ctrl+Y", false, canRedo))
            requestRedo = true;
        // Every step of the history, newest first; choosing one undoes or
        // redoes up to it.
        forceTourSubmenu("History");
        if (ImGui::BeginMenu("History", historySteps.size() > 1)) {
            for (int step = static_cast<int>(historySteps.size()) - 1; step >= 0; --step) {
                const std::string& what = historySteps[static_cast<std::size_t>(step)];
                std::string label = step == 0 ? std::string("The opened structure") : what;
                if (!label.empty()) label[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(label[0])));
                label += "##history" + std::to_string(step);
                if (ImGui::MenuItem(label.c_str(), nullptr, step == historyCurrent)) historyJump = step;
            }
            ImGui::EndMenu();
        }
        ImGui::Separator();
        if (ImGui::MenuItem("Box Select Mode", nullptr, &boxSelectMode) && boxSelectMode)
            lassoSelectMode = false;
        if (ImGui::MenuItem("Lasso Select Mode", nullptr, &lassoSelectMode) && lassoSelectMode)
            boxSelectMode = false;
        ImGui::Separator();
        editMenuDialogs.drawMenuItems();
        dislocationDialog.drawMenuItem(true);
        interstitialAtomsDialog.drawMenuItem(!structure.atoms.empty());
        ImGui::Separator();
        transformDialog.drawMenuItem(structure.hasUnitCell);
        ImGui::Separator();
        mergeStructuresDialog.drawMenuItem(true);
        ImGui::Separator();
        cellSculptorDialog.drawMenuItem(!structure.atoms.empty() || true);
        ImGui::Separator();
        // Operations that are also pipeline steps (one-step dialogs).
        if (ImGui::BeginMenu("Structure Operations")) {
            operationDialog.drawMenuItems("Edit", "Structure Operations", !structure.atoms.empty());
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Select Atoms")) {
            operationDialog.drawMenuItems("Edit", "Select Atoms", !structure.atoms.empty());
            ImGui::EndMenu();
        }
        pipelineDialog.drawMenuItem(!structure.atoms.empty());

        ImGui::EndMenu();
    }

    forceTourMenu("Build");
    if (ImGui::BeginMenu("Build"))
    {
        bulkCrystalDialog.drawMenuItem(true);
#if ATOMFORGE_ENABLE_SSS_BUILDER
        substitutionalSolidSolutionDialog.drawMenuItem(true);
#endif
#if ATOMFORGE_ENABLE_SFE_BUILDER
        stackingFaultDialog.drawMenuItem(true);
#endif
        cslDialog.drawMenuItem(true);
        interfaceBuilderDialog.drawMenuItem(true);
        nanoCrystalDialog.drawMenuItem(true);
        customStructureDialog.drawMenuItem(true);
        polyCrystalDialog.drawMenuItem(true);
        amorphousBuilderDialog.drawMenuItem(true);
        ImGui::Separator();
        operationDialog.drawMenuItems("Build", "", !structure.atoms.empty());
        ImGui::EndMenu();
    }

    forceTourMenu("View");
    if (ImGui::BeginMenu("View"))
    {
        const bool hasAtoms = !structure.atoms.empty();
        const auto setColourMode = [&](AtomColorMode mode) {
            if (atomColorMode == mode) return;
            atomColorMode = mode;
            atomColorModeJustChanged = true;
            updateBuffers(structure);
        };
        const auto setDisplayMode = [&](AtomDisplayMode mode) {
            if (atomDisplayMode == mode) return;
            atomDisplayMode = mode;
            updateBuffers(structure);
        };
        forceTourSubmenu("Colour Atoms By");
        if (ImGui::BeginMenu("Colour Atoms By"))
        {
            if (ImGui::MenuItem("Element", nullptr, atomColorMode == AtomColorMode::ElementType)) setColourMode(AtomColorMode::ElementType);
            if (ImGui::MenuItem("Crystal Orientation", nullptr, atomColorMode == AtomColorMode::CrystalOrientation)) setColourMode(AtomColorMode::CrystalOrientation);
            if (ImGui::MenuItem("Grain Boundaries", nullptr, atomColorMode == AtomColorMode::GrainBoundary)) setColourMode(AtomColorMode::GrainBoundary);
            if (ImGui::MenuItem("Atom Property", nullptr, atomColorMode == AtomColorMode::AtomProperty,
                                !structure.atomProperty.empty() && structure.atomProperty.size() == structure.atoms.size()))
                setColourMode(AtomColorMode::AtomProperty);
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled) && structure.atomProperty.empty())
                ImGui::SetTooltip("Available after a per-atom analysis (Analysis menu or a compute-property step).");
            ImGui::EndMenu();
        }
        forceTourSubmenu("Atom Style");
        if (ImGui::BeginMenu("Atom Style"))
        {
            if (ImGui::MenuItem("Balls", nullptr, atomDisplayMode == AtomDisplayMode::Balls)) setDisplayMode(AtomDisplayMode::Balls);
            if (ImGui::MenuItem("Ball and Stick", nullptr, atomDisplayMode == AtomDisplayMode::BallAndStick)) setDisplayMode(AtomDisplayMode::BallAndStick);
            if (ImGui::MenuItem("Space Filling", nullptr, atomDisplayMode == AtomDisplayMode::SpaceFilling)) setDisplayMode(AtomDisplayMode::SpaceFilling);
            if (!hasAtoms && atomDisplayMode == AtomDisplayMode::Polyhedral) atomDisplayMode = AtomDisplayMode::BallAndStick;
            if (ImGui::MenuItem("Polyhedral", nullptr, atomDisplayMode == AtomDisplayMode::Polyhedral, hasAtoms)) setDisplayMode(AtomDisplayMode::Polyhedral);
            if (!hasAtoms && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("Open a structure to use the polyhedral style.");
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Show"))
        {
            ImGui::MenuItem("Atoms", nullptr, &showAtoms, hasAtoms);
            ImGui::MenuItem("Element Labels", nullptr, &showElementLabels);
            ImGui::MenuItem("Element Legend", nullptr, &showElementLegend);
            ImGui::MenuItem("Bounding Box", nullptr, &showBoundingBox, structure.hasUnitCell);
            ImGui::Separator();
            ImGui::MenuItem("Lattice Planes", nullptr, &showLatticePlanes, structure.hasUnitCell);
            ImGui::MenuItem("Miller Directions", nullptr, &showMillerDirections, structure.hasUnitCell);
            ImGui::MenuItem("Voronoi Volumes", nullptr, &showVoronoi, hasAtoms);
            ImGui::MenuItem("Dislocation Lines", nullptr, &showDislocationLines, hasAtoms);
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Projection"))
        {
            if (ImGui::MenuItem("Perspective", nullptr, viewMode == ViewMode::Isometric) && viewMode != ViewMode::Isometric)
            {
                viewMode = ViewMode::Isometric;
                requestResetDefaultView = true;
            }
            if (ImGui::MenuItem("Orthographic", nullptr, viewMode == ViewMode::Orthographic) && viewMode != ViewMode::Orthographic)
            {
                viewMode = ViewMode::Orthographic;
                requestResetDefaultView = true;
            }
            ImGui::EndMenu();
        }
        if (ImGui::MenuItem("Reset View", "R"))
            requestResetDefaultView = true;
        ImGui::Separator();
        if (ImGui::MenuItem("Lattice Planes Settings", nullptr, false, structure.hasUnitCell))
            showLatticePlanesDialog = true;
        if (ImGui::MenuItem("Miller Directions Settings", nullptr, false, structure.hasUnitCell))
            showMillerDirectionsDialog = true;
        ImGui::MenuItem("Polyhedral Viewer", nullptr, &showPolyhedralViewer, hasAtoms);
        ImGui::Separator();
        if (ImGui::MenuItem("Structure Info", nullptr, false, hasAtoms))
            requestStructureInfo = true;
        if (ImGui::MenuItem("Atom Info", nullptr, false, hasAtoms))
            requestAtomInfo = true;
        if (ImGui::MenuItem("Measure Distance", nullptr, false, hasAtoms))
            requestMeasureDistance = true;
        if (ImGui::MenuItem("Measure Angle", nullptr, false, hasAtoms))
            requestMeasureAngle = true;
        ImGui::Separator();
        forceTourSubmenu("Theme");
        if (ImGui::BeginMenu("Theme"))
        {
            if (ImGui::MenuItem("Light", nullptr, themeMode() == ThemeMode::Light)) setThemeMode(ThemeMode::Light);
            if (ImGui::MenuItem("Dark", nullptr, themeMode() == ThemeMode::Dark)) setThemeMode(ThemeMode::Dark);
            if (ImGui::MenuItem("Match the system", nullptr, themeMode() == ThemeMode::System)) setThemeMode(ThemeMode::System);
            ImGui::EndMenu();
        }
        ImGui::EndMenu();
    }

    forceTourMenu("Analysis");
    if (ImGui::BeginMenu("Analysis"))
    {
        if (ImGui::BeginMenu("Structure and defects")) {
            cnaDialog.drawMenuItem(!structure.atoms.empty());
            rdfDialog.drawMenuItem(!structure.atoms.empty());
            angularDistributionDialog.drawMenuItem(!structure.atoms.empty());
            drawShortRangeOrderMenuItem(!structure.atoms.empty(), shortRangeOrderDialog);
            dislocationAnalysisDialog.drawMenuItem();
            ImGui::Separator();
            scientificToolsDialog.drawMenuItems("Structure and defects");
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Trajectories and transport")) {
            trajectoryDialog.drawMenuItem();
            ImGui::Separator();
            scientificToolsDialog.drawMenuItems("Trajectories and transport");
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Electronic structure")) {
            electronicDialog.drawMenuItem();
            wannierDialog.drawMenuItem();
            lobsterDialog.drawMenuItem();
            ImGui::Separator();
            scientificToolsDialog.drawMenuItems("Electronic structure");
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Reciprocal space")) {
            scientificToolsDialog.drawMenuItems("Reciprocal space");
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Thermal and mechanical")) {
            scientificToolsDialog.drawMenuItems("Thermal and mechanical");
            ImGui::EndMenu();
        }
        ImGui::Separator();
        if (ImGui::BeginMenu("Simulation")) {
            scientificToolsDialog.drawMenuItems("Simulation");
            ImGui::EndMenu();
        }
        ImGui::EndMenu();
    }

    forceTourMenu("Settings");
    if (ImGui::BeginMenu("Settings"))
    {
        forceTourSubmenu("Interface size");
        if (ImGui::BeginMenu("Interface size"))
        {
            if (ImGui::MenuItem("Automatic (display scaling)", nullptr, interfaceScale() == 1.0f)) setInterfaceScale(1.0f);
            if (ImGui::MenuItem("Compact (85%)", nullptr, interfaceScale() == .85f)) setInterfaceScale(.85f);
            if (ImGui::MenuItem("Larger (125%)", nullptr, interfaceScale() == 1.25f)) setInterfaceScale(1.25f);
            if (ImGui::MenuItem("Large (150%)", nullptr, interfaceScale() == 1.5f)) setInterfaceScale(1.5f);
            if (ImGui::MenuItem("Extra large (200%)", nullptr, interfaceScale() == 2.0f)) setInterfaceScale(2.0f);
            ImGui::EndMenu();
        }
        ImGui::Separator();
        editMenuDialogs.drawSettingsMenuItems();
        if (ImGui::MenuItem("Polyhedral Settings"))
            showPolyhedralSettingsDialog = true;
        ImGui::Separator();
        if (ImGui::BeginMenu("Theme"))
        {
            if (ImGui::MenuItem("Light", nullptr, themeMode() == ThemeMode::Light)) setThemeMode(ThemeMode::Light);
            if (ImGui::MenuItem("Dark", nullptr, themeMode() == ThemeMode::Dark)) setThemeMode(ThemeMode::Dark);
            if (ImGui::MenuItem("Match the system", nullptr, themeMode() == ThemeMode::System)) setThemeMode(ThemeMode::System);
            ImGui::EndMenu();
        }
        ImGui::EndMenu();
    }

    forceTourMenu("Help");
    if (ImGui::BeginMenu("Help"))
    {
        if (ImGui::MenuItem("Manual"))
            showManual = true;

        if (ImGui::MenuItem("Find a Command", "Ctrl+K"))
            showCommandPalette = true;

        if (ImGui::MenuItem("Keyboard Shortcuts", "F1"))
            showShortcuts = true;

        if (ImGui::MenuItem("About"))
            showAbout = true;

        ImGui::EndMenu();
    }

    ImGui::EndMainMenuBar();
}

void FileBrowser::drawOpenStructureDialog()
{
    if (openStructurePopup)
    {
        ImGui::OpenPopup("Open Structure");
        openStructurePopup = false;
        openStatusMsg[0] = '\0';
    }

    responsive::windowSize(ImVec2(720.0f, 510.0f), ImGuiCond_Appearing);
    bool openStructureOpen = true;
    if (responsive::beginModal("Open Structure", &openStructureOpen, 0))
    {
        // --- Compact navigation bar with editable path ---
        static char s_openPathBuf[2048] = {};
        static std::string s_prevOpenDir;
        if (s_prevOpenDir != openDir)
        {
            std::snprintf(s_openPathBuf, sizeof(s_openPathBuf), "%s", openDir.c_str());
            s_prevOpenDir = openDir;
        }

        const float navBtnW = 32.0f;
        if (responsive::button("\xe2\x86\x90##openBack", ImVec2(navBtnW, 0.0f)) && historyIndex > 0)
        {
            historyIndex--;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x92##openFwd", ImVec2(navBtnW, 0.0f)) && historyIndex + 1 < (int)dirHistory.size())
        {
            historyIndex++;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x91##openUp", ImVec2(navBtnW, 0.0f)))
        {
            openDir = parentPath(openDir);
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(responsive::dp(0.0f), 8.0f);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
        if (ImGui::InputText("##openPathBar", s_openPathBuf, sizeof(s_openPathBuf),
                             ImGuiInputTextFlags_EnterReturnsTrue))
        {
            std::string newDir = normalizePathSeparators(std::string(s_openPathBuf));
            if (isDriveRootPath(newDir) && newDir.size() == 2)
                newDir += "/";
            openDir = newDir;
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }

        ImGui::Separator();

        // --- Sidebar + file list ---
        // The lists take the height the fields below leave.
        const float listH = std::max(responsive::dp(200.0f),
                                     ImGui::GetContentRegionAvail().y - 2.0f * ImGui::GetFrameHeightWithSpacing() - ImGui::GetTextLineHeightWithSpacing() - responsive::dp(30.0f));

        if (responsive::beginChild("##opensidebar", ImVec2(s_sidebarW, listH), true))
        {
            ImGui::TextDisabled("Locations");
            ImGui::Separator();
            if (ImGui::Selectable("Home##openHome"))
            {
                openDir = detectHomePath();
                pushHistory(openDir);
                openStatusMsg[0] = '\0';
            }
            ImGui::Spacing();
            ImGui::TextDisabled("Drives");
            ImGui::Separator();
            for (const auto& root : driveRoots)
            {
                if (ImGui::Selectable(root.c_str()))
                {
                    openDir = root;
                    pushHistory(openDir);
                    openStatusMsg[0] = '\0';
                }
            }
            ImGui::EndChild();
        }

        // Draggable splitter
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);
        ImGui::InvisibleButton("##openSplitter", ImVec2(6.0f, listH));
        if (ImGui::IsItemHovered() || ImGui::IsItemActive())
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        if (ImGui::IsItemActive())
        {
            s_sidebarW += ImGui::GetIO().MouseDelta.x;
            s_sidebarW = std::max(80.0f, std::min(s_sidebarW, 300.0f));
        }
        {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImVec2 r0 = ImGui::GetItemRectMin();
            const ImVec2 r1 = ImGui::GetItemRectMax();
            const float  cx = (r0.x + r1.x) * 0.5f;
            const ImU32  lc = ImGui::IsItemHovered() || ImGui::IsItemActive()
                              ? IM_COL32(150, 190, 230, 200)
                              : IM_COL32(120, 120, 120, 100);
            dl->AddLine(ImVec2(cx, r0.y + 4.0f), ImVec2(cx, r1.y - 4.0f), lc, 1.5f);
        }
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);

        bool openFileDoubleClicked = false;
        if (responsive::beginChild("##filebrowser", ImVec2(0.0f, listH), true))
        {
            std::vector<DirectoryEntry> entries;
            bool listed = loadDirectoryEntries(
                openDir,
                true,
                [&](const std::string& name) { return isAllowedFile(name); },
                entries);

            if (!listed)
            {
                ImGui::TextDisabled("Unable to open folder");
            }
            else
            {
                drawDirectoryEntries(entries, openFilename, 0,
                    [&](const std::string& name) {
                        openDir = joinPath(openDir, name);
                        pushHistory(openDir);
                        openStatusMsg[0] = '\0';
                    }, &openFileDoubleClicked);
            }

            ImGui::EndChild();
        }

        {
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted("File name");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::InputTextWithHint("##openFilename", "CIF, POSCAR/.vasp, XYZ, PDB, SDF, MOL2, PWI, GJF, LAMMPS dump or .afcloud",
                                         openFilename, sizeof(openFilename), ImGuiInputTextFlags_EnterReturnsTrue))
                openFileDoubleClicked = true;
        }
        // The selected file: format, size and date.
        if (openFilename[0] != '\0') {
            const std::string details = describeFile(joinPath(openDir, openFilename));
            if (!details.empty()) ImGui::TextDisabled("%s", details.c_str());
        }
        if (openStatusMsg[0] != '\0')
            dialogLayout::note(dialogLayout::Note::Error, openStatusMsg);

        bool doLoad = openFileDoubleClicked;
        switch (dialogLayout::footer("Open", "Cancel", openFilename[0] != '\0'))
        {
            case 1: doLoad = true; break;
            case 2: ImGui::CloseCurrentPopup(); break;
            default: break;
        }

        if (doLoad)
        {
            std::string fullPath = joinPath(openDir, openFilename);
            if (!isSupportedStructureFile(openFilename) && !isLargeDataFile(openFilename))
            {
                std::snprintf(openStatusMsg, sizeof(openStatusMsg),
                              "Unsupported file format.");
            }
            else
            {
                // Signal the main loop to load this path into a new tab.
                pendingOpenPath = fullPath;
                openStatusMsg[0] = '\0';
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    if (!openStructureOpen)
        ImGui::CloseCurrentPopup();
}

void FileBrowser::drawSaveAsDialog(Structure& structure)
{
    if (saveStructurePopup)
    {
        if (saveDir.empty())
            saveDir = openDir.empty() ? "." : openDir;

        if (saveDirHistory.empty())
        {
            saveDirHistory = dirHistory;
            saveHistoryIndex = historyIndex;
            if (saveDirHistory.empty())
            {
                saveDirHistory.push_back(saveDir);
                saveHistoryIndex = 0;
            }
        }

        if (saveFilename[0] == '\0')
        {
            selectedSaveFormat = kDefaultSaveFormatIndex;
            std::string fallbackBase = "structure";
            if (!currentStructurePath.empty())
            {
                const std::string path = normalizePathSeparators(currentStructurePath);
                const std::size_t slash = path.find_last_of('/');
                std::string fileName = (slash == std::string::npos) ? path : path.substr(slash + 1);
                const std::size_t dot = fileName.find_last_of('.');
                if (dot != std::string::npos)
                    fileName = fileName.substr(0, dot);
                if (!fileName.empty())
                    fallbackBase = fileName;
            }
            const std::string initialName = replaceFileExtension(
                fallbackBase,
                kSaveFormats[selectedSaveFormat].ext,
                "structure");
            std::snprintf(saveFilename, sizeof(saveFilename), "%s", initialName.c_str());
        }

        ImGui::OpenPopup("Save As");
        saveStructurePopup = false;
    }

    responsive::windowSize(ImVec2(720.0f, 560.0f), ImGuiCond_Appearing);
    bool saveAsOpen = true;
    if (responsive::beginModal("Save As", &saveAsOpen, 0))
    {
        // Helper: navigate to a new directory and record it in history.
        auto pushSaveDir = [&](const std::string& dir) {
            saveDir = dir;
            pushDirectoryHistory(saveDirHistory, saveHistoryIndex, saveDir);
        };

        // --- Compact navigation bar with editable path ---
        static char s_savePathBuf[2048] = {};
        static std::string s_prevSaveDir;
        if (s_prevSaveDir != saveDir)
        {
            std::snprintf(s_savePathBuf, sizeof(s_savePathBuf), "%s", saveDir.c_str());
            s_prevSaveDir = saveDir;
        }

        const float navBtnW = 32.0f;
        if (responsive::button("\xe2\x86\x90##saveBack", ImVec2(navBtnW, 0.0f)) && saveHistoryIndex > 0)
        {
            saveHistoryIndex--;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x92##saveFwd", ImVec2(navBtnW, 0.0f)) && saveHistoryIndex + 1 < (int)saveDirHistory.size())
        {
            saveHistoryIndex++;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x91##saveUp", ImVec2(navBtnW, 0.0f)))
            pushSaveDir(parentPath(saveDir));

        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(responsive::dp(0.0f), 8.0f);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
        if (ImGui::InputText("##savePathBar", s_savePathBuf, sizeof(s_savePathBuf),
                             ImGuiInputTextFlags_EnterReturnsTrue))
        {
            std::string newDir = normalizePathSeparators(std::string(s_savePathBuf));
            if (isDriveRootPath(newDir) && newDir.size() == 2)
                newDir += "/";
            pushSaveDir(newDir);
        }

        ImGui::Separator();

        // --- Sidebar + file list ---
        const float listH = 200.0f;

        if (responsive::beginChild("##savesidebar", ImVec2(s_sidebarW, listH), true))
        {
            ImGui::TextDisabled("Locations");
            ImGui::Separator();
            if (ImGui::Selectable("Home##saveHome"))
                pushSaveDir(detectHomePath());
            ImGui::Spacing();
            ImGui::TextDisabled("Drives");
            ImGui::Separator();
            for (const auto& root : driveRoots)
            {
                if (ImGui::Selectable(root.c_str()))
                    pushSaveDir(root);
            }
            ImGui::EndChild();
        }

        // Draggable splitter
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);
        ImGui::InvisibleButton("##saveSplitter", ImVec2(6.0f, listH));
        if (ImGui::IsItemHovered() || ImGui::IsItemActive())
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        if (ImGui::IsItemActive())
        {
            s_sidebarW += ImGui::GetIO().MouseDelta.x;
            s_sidebarW = std::max(80.0f, std::min(s_sidebarW, 300.0f));
        }
        {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImVec2 r0 = ImGui::GetItemRectMin();
            const ImVec2 r1 = ImGui::GetItemRectMax();
            const float  cx = (r0.x + r1.x) * 0.5f;
            const ImU32  lc = ImGui::IsItemHovered() || ImGui::IsItemActive()
                              ? IM_COL32(150, 190, 230, 200)
                              : IM_COL32(120, 120, 120, 100);
            dl->AddLine(ImVec2(cx, r0.y + 4.0f), ImVec2(cx, r1.y - 4.0f), lc, 1.5f);
        }
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);

        if (responsive::beginChild("##savefilebrowser", ImVec2(0.0f, listH), true))
        {
            const std::string saveExtFilter = toLower(kSaveFormats[selectedSaveFormat].ext);

            std::vector<DirectoryEntry> entries;
            bool listed = loadDirectoryEntries(
                saveDir,
                true,
                [&](const std::string& name) {
                    const std::string lowerName = toLower(name);
                    return hasExtension(lowerName, saveExtFilter);
                },
                entries);

            if (!listed)
            {
                ImGui::TextDisabled("Unable to open folder");
            }
            else
            {
                drawDirectoryEntries(entries, saveFilename, 10000,
                    [&](const std::string& name) {
                        pushSaveDir(joinPath(saveDir, name));
                    });
            }
            ImGui::EndChild();
        }

        dialogLayout::fieldLabel("File name", ImGui::GetFontSize() * 5.5f);
        ImGui::SetNextItemWidth(-FLT_MIN);
        bool saveEnterPressed = ImGui::InputText("##saveFilename", saveFilename, sizeof(saveFilename),
                                                 ImGuiInputTextFlags_EnterReturnsTrue);

        // Format selector
        dialogLayout::fieldLabel("Format", ImGui::GetFontSize() * 5.5f);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Combo("##saveFormat", &selectedSaveFormat,
                         [](void* d, int i) -> const char* {
                             return static_cast<const SaveFormat*>(d)[i].label;
                         }, (void*)kSaveFormats, kNumSaveFormats))
        {
            updateFilenameWithExtension(
                saveFilename,
                sizeof(saveFilename),
                kSaveFormats[selectedSaveFormat].ext,
                "structure");
            saveStatusMsg[0] = '\0';
        }

        if (saveStatusMsg[0] != '\0')
            dialogLayout::note(dialogLayout::Note::Error, saveStatusMsg);

        bool doSave = saveEnterPressed;
        switch (dialogLayout::footer("Save", "Cancel", saveFilename[0] != '\0'))
        {
            case 1: doSave = true; break;
            case 2: ImGui::CloseCurrentPopup(); break;
            default: break;
        }

        if (doSave)
        {
            if (structure.atoms.empty())
            {
                std::snprintf(saveStatusMsg, sizeof(saveStatusMsg), "Error: no atoms to save.");
            }
            else if (saveFilename[0] == '\0')
            {
                std::snprintf(saveStatusMsg, sizeof(saveStatusMsg), "Error: please enter a filename.");
            }
            else
            {
                std::string fullPath = joinPath(saveDir, saveFilename);
                std::size_t savedAtomCount = 0;
                bool ok = saveStructureWithOptionalSupercell(
                    structure,
                    isTransformMatrixEnabled(),
                    getTransformMatrix(),
                    fullPath,
                    kSaveFormats[selectedSaveFormat].fmt,
                    savedAtomCount);
                if (ok)
                {
                    currentStructurePath = fullPath;
                    lastLoadedPath = fullPath;

                    std::ostringstream msg;
                    msg << "Structure saved: " << saveFilename
                        << " (" << savedAtomCount << " atoms)";
                    showNotification(msg.str(), false);
                    std::cout << "[Operation] Saved structure: " << fullPath
                              << " (format=" << kSaveFormats[selectedSaveFormat].fmt
                              << ", atoms=" << savedAtomCount << ")" << std::endl;
                    ImGui::CloseCurrentPopup();
                }
                else
                {
                    std::snprintf(saveStatusMsg, sizeof(saveStatusMsg),
                                  "Error: failed to save (format may not support this structure).");
                    showNotification(saveStatusMsg, true);
                    std::cout << "[Operation] Save failed: " << fullPath
                              << " (format=" << kSaveFormats[selectedSaveFormat].fmt << ")" << std::endl;
                }
            }
        }
        ImGui::EndPopup();
    }
    if (!saveAsOpen)
        ImGui::CloseCurrentPopup();
}

void FileBrowser::drawExportImageDialog(Structure& structure)
{
    if (exportImagePopup)
    {
        if (exportDir.empty())
            exportDir = openDir.empty() ? "." : openDir;

        if (exportDirHistory.empty())
        {
            exportDirHistory = dirHistory;
            exportHistoryIndex = historyIndex;
            if (exportDirHistory.empty())
            {
                exportDirHistory.push_back(exportDir);
                exportHistoryIndex = 0;
            }
        }

        ImGui::OpenPopup("Export Image");
        exportImagePopup = false;
    }

    responsive::windowSize(ImVec2(760.0f, 680.0f), ImGuiCond_Appearing);
    bool exportImageOpen = true;
    if (responsive::beginModal("Export Image", &exportImageOpen, 0))
    {
        auto pushExportDir = [&](const std::string& dir) {
            exportDir = dir;
            pushDirectoryHistory(exportDirHistory, exportHistoryIndex, exportDir);
        };

        // --- Compact navigation bar with editable path ---
        static char s_exportPathBuf[2048] = {};
        static std::string s_prevExportDir;
        if (s_prevExportDir != exportDir)
        {
            std::snprintf(s_exportPathBuf, sizeof(s_exportPathBuf), "%s", exportDir.c_str());
            s_prevExportDir = exportDir;
        }

        const float navBtnW = 32.0f;
        if (responsive::button("\xe2\x86\x90##exportBack", ImVec2(navBtnW, 0.0f)) && exportHistoryIndex > 0)
        {
            exportHistoryIndex--;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x92##exportFwd", ImVec2(navBtnW, 0.0f)) && exportHistoryIndex + 1 < (int)exportDirHistory.size())
        {
            exportHistoryIndex++;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(responsive::dp(0.0f), 2.0f);
        if (responsive::button("\xe2\x86\x91##exportUp", ImVec2(navBtnW, 0.0f)))
            pushExportDir(parentPath(exportDir));
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(responsive::dp(0.0f), 8.0f);
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
        if (ImGui::InputText("##exportPathBar", s_exportPathBuf, sizeof(s_exportPathBuf),
                             ImGuiInputTextFlags_EnterReturnsTrue))
        {
            std::string newDir = normalizePathSeparators(std::string(s_exportPathBuf));
            if (isDriveRootPath(newDir) && newDir.size() == 2)
                newDir += "/";
            pushExportDir(newDir);
        }

        ImGui::Separator();

        // --- Sidebar + file list (taking the height the fields below leave) ---
        const float listH = std::max(responsive::dp(200.0f),
                                     ImGui::GetContentRegionAvail().y - 6.0f * ImGui::GetFrameHeightWithSpacing() - responsive::dp(28.0f));

        if (responsive::beginChild("##exportsidebar", ImVec2(s_sidebarW, listH), true))
        {
            ImGui::TextDisabled("Locations");
            ImGui::Separator();
            if (ImGui::Selectable("Home##exportHome"))
                pushExportDir(detectHomePath());
            ImGui::Spacing();
            ImGui::TextDisabled("Drives");
            ImGui::Separator();
            for (const auto& root : driveRoots)
            {
                if (ImGui::Selectable(root.c_str()))
                    pushExportDir(root);
            }
            ImGui::EndChild();
        }

        // Draggable splitter
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);
        ImGui::InvisibleButton("##exportSplitter", ImVec2(6.0f, listH));
        if (ImGui::IsItemHovered() || ImGui::IsItemActive())
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
        if (ImGui::IsItemActive())
        {
            s_sidebarW += ImGui::GetIO().MouseDelta.x;
            s_sidebarW = std::max(80.0f, std::min(s_sidebarW, 300.0f));
        }
        {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImVec2 r0 = ImGui::GetItemRectMin();
            const ImVec2 r1 = ImGui::GetItemRectMax();
            const float  cx = (r0.x + r1.x) * 0.5f;
            const ImU32  lc = ImGui::IsItemHovered() || ImGui::IsItemActive()
                              ? IM_COL32(150, 190, 230, 200)
                              : IM_COL32(120, 120, 120, 100);
            dl->AddLine(ImVec2(cx, r0.y + 4.0f), ImVec2(cx, r1.y - 4.0f), lc, 1.5f);
        }
        ImGui::SameLine(responsive::dp(0.0f), 0.0f);

        if (responsive::beginChild("##exportfilebrowser", ImVec2(0.0f, listH), true))
        {
            const std::string exportExtFilter = toLower(kImageExportFormats[selectedExportFormat].ext);

            std::vector<DirectoryEntry> entries;
            bool listed = loadDirectoryEntries(
                exportDir,
                true,
                [&](const std::string& name) {
                    const std::string lowerName = toLower(name);
                    return hasExtension(lowerName, exportExtFilter);
                },
                entries);

            if (!listed)
            {
                ImGui::TextDisabled("Unable to open folder");
            }
            else
            {
                drawDirectoryEntries(entries, exportFilename, 20000,
                    [&](const std::string& name) {
                        pushExportDir(joinPath(exportDir, name));
                    });
            }
            ImGui::EndChild();
        }

        const float labelColumn = ImGui::GetFontSize() * 7.5f;
        dialogLayout::fieldLabel("File name", labelColumn);
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputText("##exportFilename", exportFilename, sizeof(exportFilename));

        dialogLayout::fieldLabel("Format", labelColumn);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::Combo("##exportFormat", &selectedExportFormat,
                         [](void* d, int i) -> const char* {
                             return static_cast<const ImageExportFormatOption*>(d)[i].label;
                         }, (void*)kImageExportFormats, kNumImageExportFormats))
        {
            updateFilenameWithExtension(
                exportFilename,
                sizeof(exportFilename),
                kImageExportFormats[selectedExportFormat].ext,
                "structure");
            exportStatusMsg[0] = '\0';
        }

        ImGui::Checkbox("Include background", &exportIncludeBackground);
        if (kImageExportFormats[selectedExportFormat].fmt == ImageExportFormat::Jpg && !exportIncludeBackground)
            ImGui::TextDisabled("JPEG has no transparency; transparent areas will be blended on white.");
        ImGui::Checkbox("Include orientation gizmo", &exportIncludeGizmo);

        dialogLayout::fieldLabel("Resolution", labelColumn);
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::SliderInt("##exportScale", &exportResolutionScale, 1, 8, "%dx");
        if (exportResolutionScale > 1)
            ImGui::TextDisabled("Output will be %dx the current window size.", exportResolutionScale);

        if (exportStatusMsg[0] != '\0')
            dialogLayout::note(dialogLayout::Note::Error, exportStatusMsg);

        const int exportPressed = dialogLayout::footer("Export", "Cancel");
        if (exportPressed == 2)
            ImGui::CloseCurrentPopup();
        if (exportPressed == 1)
        {
            if (structure.atoms.empty())
            {
                std::snprintf(exportStatusMsg, sizeof(exportStatusMsg), "Error: no atoms to export.");
            }
            else if (exportFilename[0] == '\0')
            {
                std::snprintf(exportStatusMsg, sizeof(exportStatusMsg), "Error: please enter a filename.");
            }
            else
            {
                std::string finalName(exportFilename);
                std::string selectedExt = kImageExportFormats[selectedExportFormat].ext;
                std::string selectedExtLower = toLower(selectedExt);
                std::string currentExt;
                std::size_t dot = finalName.find_last_of('.');
                if (dot != std::string::npos)
                    currentExt = toLower(finalName.substr(dot));

                if (currentExt != selectedExtLower)
                {
                    finalName = replaceFileExtension(finalName, selectedExt, "structure");
                    std::snprintf(exportFilename, sizeof(exportFilename), "%s", finalName.c_str());
                }

                pendingImageExport.outputPath = joinPath(exportDir, finalName);
                pendingImageExport.format = kImageExportFormats[selectedExportFormat].fmt;
                pendingImageExport.includeBackground = exportIncludeBackground;
                pendingImageExport.includeGizmo = exportIncludeGizmo;
                pendingImageExport.resolutionScale = exportResolutionScale;
                requestImageExport = true;

                std::cout << "[Operation] Image export requested: " << pendingImageExport.outputPath
                          << " (format=" << selectedExt << ", background="
                          << (exportIncludeBackground ? "on" : "off")
                          << ", scale=" << exportResolutionScale << "x)" << std::endl;

                showNotification(std::string("Image export started: ") + finalName, false);
                exportStatusMsg[0] = '\0';
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    if (!exportImageOpen)
        ImGui::CloseCurrentPopup();
}

void FileBrowser::draw(Structure& structure,
                       EditMenuDialogs& editMenuDialogs,
                       const std::function<void(Structure&)>& updateBuffers,
                       const std::function<void(Structure)>& openInNewTab,
                       bool canUndo,
                       bool canRedo)
{
    if (requestSaveFile)
    {
        requestSaveFile = false;
        performQuickSave(structure);
    }

    if (pendingColourMode >= 0) {
        if (atomColorMode != static_cast<AtomColorMode>(pendingColourMode)) {
            atomColorMode = static_cast<AtomColorMode>(pendingColourMode);
            atomColorModeJustChanged = true;
            updateBuffers(structure);
        }
        pendingColourMode = -1;
    }
    drawMainMenuBar(structure, editMenuDialogs, updateBuffers, canUndo, canRedo);

    // The main window around the view: welcome screen, toolbar and status bar.
    drawWelcome(structure);
    drawToolbar(structure, canUndo, canRedo);
    drawStatusBar(structure);

    if (StepEdit request; pipelineDialog.consumeDialogRequest(request)) openStepDialog(std::move(request), editMenuDialogs);
    editMenuDialogs.drawPopups(structure, updateBuffers);

    transformDialog.drawDialog([&]() { updateBuffers(structure); });

    // Callback for builder dialogs: captures the built result and opens it in a new tab,
    // restoring the current tab's structure to its pre-build state.
    Structure preBuilderState = structure;
    auto updateFromBuilderToNewTab = [&](Structure& s) {
        Structure result = std::move(s);
        s = std::move(preBuilderState);
        // Say where the result went (a new tab, or this empty one).
        showNotification("Built " + std::to_string(result.atoms.size()) + (result.atoms.size() == 1 ? " atom" : " atoms") +
                         (s.atoms.empty() ? "." : " in a new tab."));
        openInNewTab(std::move(result));
    };

    bulkCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors, updateFromBuilderToNewTab);
    cslDialog.drawDialog(structure, editMenuDialogs.elementColors,
                         editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                         updateFromBuilderToNewTab);
    nanoCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateFromBuilderToNewTab);
    customStructureDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                      updateFromBuilderToNewTab);
    mergeStructuresDialog.drawDialog(structure, updateFromBuilderToNewTab);
    interfaceBuilderDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                      updateFromBuilderToNewTab);
    polyCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateFromBuilderToNewTab);
    dislocationDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii,
                                 editMenuDialogs.elementShininess,
                                 [&](Structure& s) { updateBuffers(s); });
#if ATOMFORGE_ENABLE_SFE_BUILDER
    stackingFaultDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                   editMenuDialogs.elementRadii,
                                   editMenuDialogs.elementShininess,
                                   updateFromBuilderToNewTab);
#endif
#if ATOMFORGE_ENABLE_SSS_BUILDER
    substitutionalSolidSolutionDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                                  editMenuDialogs.elementRadii,
                                                  editMenuDialogs.elementShininess,
                                                  updateFromBuilderToNewTab);
#endif
    amorphousBuilderDialog.drawDialog(editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii,
                                      updateFromBuilderToNewTab);
    interstitialAtomsDialog.drawDialog(structure,
                                       editMenuDialogs.elementColors,
                                       editMenuDialogs.elementRadii,
                                       editMenuDialogs.elementShininess,
                                       [&](Structure& s) { updateBuffers(s); });
    cnaDialog.drawDialog(structure);
    rdfDialog.drawDialog(structure);
    electronicDialog.drawDialog();
    wannierDialog.drawDialog();
    lobsterDialog.drawDialog();
    dislocationAnalysisDialog.drawDialog();
    const auto showAtomProperty = [&] {
        propertyDisplay.autoRange = true;
        if (atomColorMode != AtomColorMode::AtomProperty) atomColorModeJustChanged = true;
        atomColorMode = AtomColorMode::AtomProperty;
    };
    trajectoryDialog.draw(structure, updateBuffers, showAtomProperty);
    if (trajectoryDialog.consumeLoaded()) requestResetDefaultView = true;  // frame the trajectory
    pipelineDialog.draw(structure, updateBuffers, showAtomProperty);
    OperationDialog::Callbacks operations;
    operations.update = updateBuffers;
    operations.select = [&](const std::vector<int>& atoms) { requestedSelection = atoms; selectionRequested = true; };
    operations.addToPipeline = [&](const std::vector<atomforge::pipeline::Modifier>& steps) { pipelineDialog.addSteps(steps, structure); };
    operations.showAtomProperty = showAtomProperty;
    operationDialog.draw(structure, viewSelectionAtoms, operations);
    scientificToolsDialog.setTrajectorySource(trajectoryDialog.loadedPath());
    scientificToolsDialog.draw(structure, updateFromBuilderToNewTab,
        [&](const std::string& name, const std::vector<double>& values) {
            structure.atomProperty = values;
            structure.atomPropertyName = name;
            propertyDisplay.autoRange = true;
            atomColorMode = AtomColorMode::AtomProperty;
            atomColorModeJustChanged = true;
            updateBuffers(structure);
        });
    drawAtomPropertySettings(structure, updateBuffers);
    drawShortRangeOrderDialog(shortRangeOrderDialog, structure);
    angularDistributionDialog.drawDialog(structure);
    cellSculptorDialog.drawDialog(structure, updateBuffers);
    runStepDialogCheck(editMenuDialogs);
    runUiTour(structure, editMenuDialogs, updateBuffers);
    runTutorialCapture(structure, editMenuDialogs, updateBuffers);

    if (loadErrorPopupRequested)
    {
        ImGui::OpenPopup(loadPopupTitle);
        loadErrorPopupRequested = false;
    }

    responsive::windowSize(ImVec2(720.0f, 0.0f), ImGuiCond_Appearing);
    bool loadErrorOpen = true;
    if (responsive::beginModal(loadPopupTitle, &loadErrorOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::TextUnformatted(loadErrorMsg);
        ImGui::Spacing();
        if (responsive::button("OK", responsive::size(120,0)))
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (!loadErrorOpen)
        ImGui::CloseCurrentPopup();

    drawOpenStructureDialog();
    drawSaveAsDialog(structure);
    drawExportImageDialog(structure);

    if (showLatticePlanesDialog)
    {
        ImGui::OpenPopup("Lattice Planes");
        showLatticePlanesDialog = false;
    }

    ImGui::SetNextWindowSizeConstraints(ImVec2(responsive::dp(780.0f), 0.0f), ImVec2(responsive::dp(780.0f), FLT_MAX));
    bool latticePlanesOpen = true;
    if (responsive::beginModal("Lattice Planes", &latticePlanesOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("Current structure has no unit cell. Lattice planes are unavailable.");
            if (responsive::button("Close", responsive::size(120,0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
        else
        {
            ImGui::Checkbox("Show lattice planes", &showLatticePlanes);
            ImGui::TextDisabled("Live edit enabled: HKL and offset updates are shown immediately.");
            ImGui::Separator();

            ImGui::Text("New plane:");
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("H##plane", &latticePlaneInputH, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("K##plane", &latticePlaneInputK, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("L##plane", &latticePlaneInputL, 0.2f, -50, 50);

            ImGui::SetNextItemWidth(responsive::dp(210.0f));
            ImGui::DragFloat("Offset n##plane", &latticePlaneInputOffset, 0.01f, -50.0f, 50.0f, "%.3f");
            ImGui::SetNextItemWidth(responsive::dp(210.0f));
            ImGui::SliderFloat("Opacity##plane", &latticePlaneInputOpacity, 0.0f, 1.0f, "%.2f");
            ImGui::SameLine();
            ImGui::ColorEdit3("Color##plane", latticePlaneInputColor, ImGuiColorEditFlags_NoInputs);

            const bool invalidMiller =
                (latticePlaneInputH == 0 && latticePlaneInputK == 0 && latticePlaneInputL == 0);
            if (invalidMiller)
                ImGui::TextDisabled("(h, k, l) cannot all be zero.");

            if ((responsive::button("Add Plane", responsive::size(140.0f,0.0f)) || uiScript::take("lattice.add")) && !invalidMiller)
            {
                LatticePlane plane;
                plane.h = latticePlaneInputH;
                plane.k = latticePlaneInputK;
                plane.l = latticePlaneInputL;
                plane.offset = latticePlaneInputOffset;
                plane.opacity = latticePlaneInputOpacity;
                plane.color = {
                    latticePlaneInputColor[0],
                    latticePlaneInputColor[1],
                    latticePlaneInputColor[2]
                };
                plane.visible = true;
                latticePlanes.push_back(plane);
                showLatticePlanes = true;
            }

            ImGui::Separator();
            ImGui::Text("Saved planes:");

            if (latticePlanes.empty())
            {
                ImGui::TextDisabled("No lattice planes added.");
            }
            else if (responsive::beginChild("##lattice-plane-list", responsive::size(0.0f,220.0f), true))
            {
                int deleteIndex = -1;
                for (size_t i = 0; i < latticePlanes.size(); ++i)
                {
                    ImGui::PushID((int)i + 30000);
                    LatticePlane& plane = latticePlanes[i];

                    ImGui::Checkbox("##visible", &plane.visible);
                    ImGui::SameLine();
                    ImGui::Text("Plane %d", (int)i + 1);

                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("H##plane-h", &plane.h, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("K##plane-k", &plane.k, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("L##plane-l", &plane.l, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(120.0f));
                    ImGui::DragFloat("n##plane-offset", &plane.offset, 0.01f, -50.0f, 50.0f, "%.3f");
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(95.0f));
                    ImGui::SliderFloat("a##plane-opacity", &plane.opacity, 0.0f, 1.0f, "%.2f");

                    if (plane.h == 0 && plane.k == 0 && plane.l == 0)
                        plane.h = 1;

                    float color[3] = {plane.color[0], plane.color[1], plane.color[2]};
                    if (ImGui::ColorEdit3("##plane-color", color, ImGuiColorEditFlags_NoInputs))
                    {
                        plane.color[0] = color[0];
                        plane.color[1] = color[1];
                        plane.color[2] = color[2];
                    }

                    ImGui::SameLine();
                    if (responsive::button("Delete"))
                        deleteIndex = (int)i;

                    ImGui::Separator();

                    ImGui::PopID();
                }

                if (deleteIndex >= 0)
                    latticePlanes.erase(latticePlanes.begin() + deleteIndex);

                ImGui::EndChild();
            }

            ImGui::Spacing();
            if (responsive::button("Close", responsive::size(120,0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }
    if (!latticePlanesOpen)
        ImGui::CloseCurrentPopup();

    // ---- Miller Directions dialog ----------------------------------------
    if (showMillerDirectionsDialog)
    {
        ImGui::OpenPopup("Miller Directions");
        showMillerDirectionsDialog = false;
    }

    ImGui::SetNextWindowSizeConstraints(ImVec2(responsive::dp(680.0f), 0.0f), ImVec2(responsive::dp(680.0f), FLT_MAX));
    bool millerDirOpen = true;
    if (responsive::beginModal("Miller Directions", &millerDirOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("Current structure has no unit cell. Miller directions are unavailable.");
            if (responsive::button("Close", responsive::size(120,0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
        else
        {
            ImGui::Checkbox("Show Miller directions", &showMillerDirections);
            ImGui::TextDisabled("Arrows are drawn from the cell origin along the [uvw] direction.");
            ImGui::Separator();

            ImGui::Text("New direction:");
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("U##dir", &millerDirInputU, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("V##dir", &millerDirInputV, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(responsive::dp(85.0f));
            ImGui::DragInt("W##dir", &millerDirInputW, 0.2f, -50, 50);

            ImGui::SetNextItemWidth(responsive::dp(210.0f));
            ImGui::DragFloat("Length (A)##dir", &millerDirInputLength, 0.1f, 0.5f, 100.0f, "%.1f");
            ImGui::SameLine();
            ImGui::ColorEdit3("Color##dir", millerDirInputColor, ImGuiColorEditFlags_NoInputs);

            const bool invalidDir =
                (millerDirInputU == 0 && millerDirInputV == 0 && millerDirInputW == 0);
            if (invalidDir)
                ImGui::TextDisabled("(u, v, w) cannot all be zero.");

            if ((responsive::button("Add Direction", responsive::size(140.0f,0.0f)) || uiScript::take("miller.add")) && !invalidDir)
            {
                MillerDirection md;
                md.u = millerDirInputU;
                md.v = millerDirInputV;
                md.w = millerDirInputW;
                md.length = millerDirInputLength;
                md.color = { millerDirInputColor[0], millerDirInputColor[1], millerDirInputColor[2] };
                md.visible = true;
                millerDirections.push_back(md);
                showMillerDirections = true;
            }

            ImGui::Separator();
            ImGui::Text("Saved directions:");

            if (millerDirections.empty())
            {
                ImGui::TextDisabled("No Miller directions added.");
            }
            else if (responsive::beginChild("##miller-dir-list", responsive::size(0.0f,200.0f), true))
            {
                int deleteIndex = -1;
                for (size_t i = 0; i < millerDirections.size(); ++i)
                {
                    ImGui::PushID((int)i + 40000);
                    MillerDirection& md = millerDirections[i];

                    ImGui::Checkbox("##vis", &md.visible);
                    ImGui::SameLine();
                    ImGui::Text("Dir %d", (int)i + 1);

                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("U##du", &md.u, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("V##dv", &md.v, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(80.0f));
                    ImGui::DragInt("W##dw", &md.w, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(responsive::dp(100.0f));
                    ImGui::DragFloat("L##dl", &md.length, 0.1f, 0.5f, 100.0f, "%.1f");

                    if (md.u == 0 && md.v == 0 && md.w == 0)
                        md.u = 1;

                    float col[3] = { md.color[0], md.color[1], md.color[2] };
                    ImGui::SameLine();
                    if (ImGui::ColorEdit3("##dc", col, ImGuiColorEditFlags_NoInputs))
                    {
                        md.color[0] = col[0];
                        md.color[1] = col[1];
                        md.color[2] = col[2];
                    }
                    ImGui::SameLine();
                    if (responsive::button("Delete##ddel"))
                        deleteIndex = (int)i;

                    ImGui::Separator();
                    ImGui::PopID();
                }

                if (deleteIndex >= 0)
                    millerDirections.erase(millerDirections.begin() + deleteIndex);

                ImGui::EndChild();
            }

            ImGui::Spacing();
            if (responsive::button("Close##dclose", responsive::size(120,0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }
    if (!millerDirOpen)
        ImGui::CloseCurrentPopup();
    // ---- end Miller Directions dialog ------------------------------------

    if (showPolyhedralSettingsDialog)
    {
        ImGui::OpenPopup("Polyhedral Settings");
        showPolyhedralSettingsDialog = false;
    }

    ImGui::SetNextWindowSizeConstraints(ImVec2(responsive::dp(620.0f), 0.0f), ImVec2(responsive::dp(620.0f), FLT_MAX));
    bool polyhedralSettingsOpen = true;
    if (responsive::beginModal("Polyhedral Settings", &polyhedralSettingsOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::SetNextItemWidth(responsive::dp(240.0f));
        ImGui::SliderInt("Max displayed centers", &polyhedralSettings.maxDisplayedCenters, 1, 2000);
        ImGui::SetNextItemWidth(responsive::dp(240.0f));
        ImGui::SliderInt("Nearest neighbors", &polyhedralSettings.maxNeighborCandidates, 4, 64);

        ImGui::Separator();

        if (ImGui::Checkbox("Specify center atom indices", &polyhedralSettings.centerAtomIndexFilterEnabled)
            && polyhedralSettings.centerAtomIndexFilterEnabled)
        {
            updatePolyhedralCenterAtomIndexFilter(polyhedralCenterAtomIndexInput);
        }
        if (polyhedralSettings.centerAtomIndexFilterEnabled)
        {
            ImGui::SetNextItemWidth(responsive::dp(420.0f));
            if (ImGui::InputText("Center atom IDs##polyIndexFilter",
                                 polyhedralCenterAtomIndexInput,
                                 sizeof(polyhedralCenterAtomIndexInput)))
            {
                updatePolyhedralCenterAtomIndexFilter(polyhedralCenterAtomIndexInput);
            }
            ImGui::TextDisabled("1-based IDs, comma-separated, ranges allowed (e.g. 1,4,10-20)");
            ImGui::TextDisabled("Parsed centers: %d", (int)polyhedralSettings.centerAtomIndices.size());
        }

        ImGui::Separator();

        ImGui::Checkbox("Show polyhedral edges", &polyhedralSettings.showEdges);
        ImGui::SetNextItemWidth(responsive::dp(220.0f));
        ImGui::SliderFloat("Face opacity", &polyhedralSettings.faceOpacity, 0.00f, 1.00f, "%.2f");
        ImGui::SetNextItemWidth(responsive::dp(220.0f));
        ImGui::SliderFloat("Edge opacity", &polyhedralSettings.edgeOpacity, 0.00f, 1.00f, "%.2f");

        ImGui::Separator();

        if (ImGui::Checkbox("Filter center elements", &polyhedralSettings.centerElementFilterEnabled)
            && polyhedralSettings.centerElementFilterEnabled)
        {
            updatePolyhedralElementFilterMask(polyhedralCenterFilterInput,
                                              polyhedralSettings.centerElementMask);
        }
        if (polyhedralSettings.centerElementFilterEnabled)
        {
            ImGui::SetNextItemWidth(responsive::dp(340.0f));
            if (ImGui::InputText("Centers##polyFilter", polyhedralCenterFilterInput,
                                 sizeof(polyhedralCenterFilterInput)))
            {
                updatePolyhedralElementFilterMask(polyhedralCenterFilterInput,
                                                  polyhedralSettings.centerElementMask);
            }
            ImGui::TextDisabled("Comma-separated symbols, e.g. Ti,Zr,Fe");
        }

        ImGui::Spacing();

        if (ImGui::Checkbox("Filter ligand elements", &polyhedralSettings.ligandElementFilterEnabled)
            && polyhedralSettings.ligandElementFilterEnabled)
        {
            updatePolyhedralElementFilterMask(polyhedralLigandFilterInput,
                                              polyhedralSettings.ligandElementMask);
        }
        if (polyhedralSettings.ligandElementFilterEnabled)
        {
            ImGui::SetNextItemWidth(responsive::dp(340.0f));
            if (ImGui::InputText("Ligands##polyFilter", polyhedralLigandFilterInput,
                                 sizeof(polyhedralLigandFilterInput)))
            {
                updatePolyhedralElementFilterMask(polyhedralLigandFilterInput,
                                                  polyhedralSettings.ligandElementMask);
            }
            ImGui::TextDisabled("Comma-separated symbols, e.g. O,N,S,F");
        }

        ImGui::Separator();
        if (responsive::button("Reset Defaults", responsive::size(150.0f,0.0f)))
        {
            polyhedralSettings = PolyhedralOverlaySettings{};
            polyhedralCenterAtomIndexInput[0] = '\0';
            polyhedralCenterFilterInput[0] = '\0';
            polyhedralLigandFilterInput[0] = '\0';
        }
        ImGui::SameLine(responsive::dp(0.0f), 8.0f);
        if (responsive::button("Close##polySettings", responsive::size(120.0f,0.0f)))
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
    if (!polyhedralSettingsOpen)
        ImGui::CloseCurrentPopup();

    if (showManual)
    {
        showManual = false;
        openPdfManual();
    }

    auto wrappedBullet = [](const char* text) {
        ImGui::Bullet();
        ImGui::SameLine();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextWrapped("%s", text);
        ImGui::PopTextWrapPos();
    };

    if (ImGui::IsKeyPressed(ImGuiKey_F1, false) && !ImGui::GetIO().WantTextInput)
        showShortcuts = true;
    drawShortcuts();
    drawCommandPalette(structure, editMenuDialogs, updateBuffers);

    if (showAbout)
    {
        ImGui::OpenPopup("About");
        showAbout = false;
    }

    responsive::windowSize(ImVec2(860.0f, 620.0f), ImGuiCond_Appearing);
    bool aboutOpen = true;
    if (responsive::beginModal("About", &aboutOpen, ImGuiWindowFlags_NoResize))
    {
        ImDrawList* drawList = ImGui::GetWindowDrawList();
        const ImVec2 iconPos = ImGui::GetCursorScreenPos();
        const float iconSize = 52.0f;
        const ImU32 boxCol  = IM_COL32(46, 94, 119, 166);
        const ImU32 nodeCol = IM_COL32(63, 110, 134, 255);
        const ImU32 coreCol = IM_COL32(18, 63, 87, 255);

        const float p0 = iconPos.x + 6.0f;
        const float p1 = iconPos.x + iconSize * 0.5f;
        const float p2 = iconPos.x + iconSize - 6.0f;
        const float q0 = iconPos.y + 6.0f;
        const float q1 = iconPos.y + iconSize * 0.5f;
        const float q2 = iconPos.y + iconSize - 6.0f;

        // Box outline (4 sides only).
        drawList->AddLine(ImVec2(p0, q0), ImVec2(p2, q0), boxCol, 3.0f);
        drawList->AddLine(ImVec2(p0, q2), ImVec2(p2, q2), boxCol, 3.0f);
        drawList->AddLine(ImVec2(p0, q0), ImVec2(p0, q2), boxCol, 3.0f);
        drawList->AddLine(ImVec2(p2, q0), ImVec2(p2, q2), boxCol, 3.0f);

        // Corner atoms and central atom (no plus cross).
        drawList->AddCircleFilled(ImVec2(p0, q0), 8.5f, nodeCol);
        drawList->AddCircleFilled(ImVec2(p2, q0), 8.5f, nodeCol);
        drawList->AddCircleFilled(ImVec2(p0, q2), 8.5f, nodeCol);
        drawList->AddCircleFilled(ImVec2(p2, q2), 8.5f, nodeCol);
        drawList->AddCircleFilled(ImVec2(p1, q1), 10.2f, coreCol);

        ImGui::Dummy(ImVec2(iconSize + 8.0f, iconSize));
        ImGui::SameLine();

        ImGui::BeginGroup();
        ImGui::Text("AtomForge");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Version " ATOMFORGE_VERSION " - Atomistic structure builder, viewer, and editor with periodic-cell and finite-geometry workflows");
        ImGui::PopTextWrapPos();
        ImGui::EndGroup();
        ImGui::Separator();

        if (responsive::beginChild("##about-scroll", ImVec2(0.0f, -ImGui::GetFrameHeightWithSpacing() - 6.0f), false))
        {
            ImGui::Text("Creator");
            wrappedBullet("Albert Linda");

            ImGui::Spacing();
            ImGui::Text("Core Libraries");
            wrappedBullet("ImGui: immediate-mode GUI framework (by Omar Cornut)");
            wrappedBullet("OpenGL 3.3+: GPU rendering and visualization");
            wrappedBullet("GLFW 3: window and input management");
            wrappedBullet("GLEW: OpenGL extension loader");
            wrappedBullet("GLM: mathematics library for transformations");

            ImGui::Spacing();
            ImGui::Text("Chemistry Libraries");
            wrappedBullet("Open Babel: molecular file I/O (.cif, .mol, .pdb, .xyz, .sdf, .mol2, etc.)");
            wrappedBullet("spglib: crystallographic space group and symmetry operations");

            ImGui::Spacing();
            ImGui::Text("Data & References");
            wrappedBullet("Periodic table element data (atomic numbers, masses, default colors)");
            wrappedBullet("Cordero et al. (Dalton Trans. 2008): Van der Waals radii");
            wrappedBullet("Standard crystallographic conventions for cell vectors and positions");

            ImGui::Spacing();
            ImGui::Text("Algorithm References");
            wrappedBullet("CSL Grain Boundary Builder: Jianli Cheng, Jian Luo, and Kesong Yang, Aimsgb: An Algorithm and Open-Source Python Library to Generate Periodic Grain Boundary Structures, Comput. Mater. Sci. 155, 92-103 (2018).");
            wrappedBullet("Nanocrystal Builder (Wulff construction): Georgios D. Barmparis, Zbigniew Lodziana, Nuria Lopez, and Ioannis N. Remediakis, Nanoparticle shapes by using Wulff constructions and first-principles calculations, Beilstein J. Nanotechnol. 6, 361-368 (2015).");
            wrappedBullet("Dislocation Builder and Post-Processing: native reimplementations of the physical models and algorithms used by the BABEL dislocation toolkit (Douin, Clouet et al.). Isotropic and anisotropic (Stroh sextic formalism) displacement fields: A. N. Stroh, Philos. Mag. 3, 625 (1958); J. Math. Phys. 41, 77 (1962); D. M. Barnett and J. Lothe, Phys. Norv. 7, 13 (1973); J. P. Hirth and J. Lothe, Theory of Dislocations, 2nd ed. (1982), ch. 13. Nye tensor: J. F. Nye, Acta Metall. 1, 153 (1953); Hartley and Mishin, Acta Mater. 53, 1313 (2005), with the index-order correction documented in BABEL's own source (attributed there to David Olmsted). Vitek differential-displacement map: V. Vitek, Philos. Mag. 18, 773 (1968).");

            ImGui::Spacing();
            ImGui::Text("Features");
            wrappedBullet("Interactive 3D visualization with real-time rendering");
            wrappedBullet("Full undo/redo support for all editing operations");
            wrappedBullet("Multi-format structure loading and exporting");
            wrappedBullet("Crystallographic analysis: space groups, symmetry, cell metrics");
            wrappedBullet("Bulk crystal, CSL grain boundary, nanocrystal, polycrystal, substitutional solid solution, and custom mesh-filled structure builders");
            wrappedBullet("Wulff-construction nanocrystal builder with symmetry-expanded facet preview for periodic reference crystals, including non-cubic systems");
            wrappedBullet("Interactive Merge Structures workflow with per-structure 3D gizmo transforms");
            wrappedBullet("Customizable atom colors, sizes, and materials");
            wrappedBullet("Editable Wulff facet colors with dedicated 3D preview and high-contrast face boundaries");
            wrappedBullet("Crystal-orientation IPF coloring with in-view legend and saved sidecar metadata");
            wrappedBullet("Distance and angle measurement tools");
            wrappedBullet("Light and dark UI themes with matching overlay rendering");
            wrappedBullet("High-resolution image export with adjustable scale and HiDPI-aware interface scaling");
            wrappedBullet("Supercell transformation with periodic boundary visualization");
            wrappedBullet("Headless CLI builder workflows via --build modes (bulk, gb, poly, nano, amorphous, sss, custom)");

            ImGui::EndChild();
        }

        if (responsive::button("Manual##about", responsive::size(120.0f,0.0f)))
            openPdfManual();
        ImGui::SameLine();
        if (responsive::button("Close##about", responsive::size(120.0f,0.0f)))
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
    if (!aboutOpen)
        ImGui::CloseCurrentPopup();

    if (showEditColors)
    {
        ImGui::OpenPopup("Edit Element Colors");
        showEditColors = false;

        if (selectedAtomicNumber < 1 || selectedAtomicNumber > 118)
            selectedAtomicNumber = 1;
    }

    responsive::windowSize(ImVec2(950, 600), ImGuiCond_FirstUseEver);
    bool editElementColorsOpen = true;
    if (responsive::beginModal("Edit Element Colors", &editElementColorsOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select an element to edit its color.");
        ImGui::Separator();

        drawPeriodicTableInlineSelector(selectedAtomicNumber);

        if (selectedAtomicNumber >= 1 && selectedAtomicNumber <= 118)
        {
            const char* selectedElementSymbol = elementSymbol(selectedAtomicNumber);

            // Determine the current color (override > structure atoms > default CPK)
            float color[3] = {0.0f, 0.0f, 0.0f};
            bool hasAtoms = false;

            auto overrideIt = elementColorOverrides.find(selectedAtomicNumber);
            if (overrideIt != elementColorOverrides.end())
            {
                color[0] = overrideIt->second[0];
                color[1] = overrideIt->second[1];
                color[2] = overrideIt->second[2];
            }
            else
            {
                auto it = std::find_if(structure.atoms.begin(), structure.atoms.end(),
                                       [&](const AtomSite& a) { return a.symbol == selectedElementSymbol; });
                if (it != structure.atoms.end())
                {
                    color[0] = it->r;
                    color[1] = it->g;
                    color[2] = it->b;
                    hasAtoms = true;
                }
                else
                {
                    getDefaultElementColor(selectedAtomicNumber, color[0], color[1], color[2]);
                }
            }

            ImGui::Text("Element: %s", selectedElementSymbol);
            if (!hasAtoms)
                ImGui::TextDisabled("(Not present in current structure)");

            if (ImGui::ColorEdit3("Color", color))
            {
                elementColorOverrides[selectedAtomicNumber] = { color[0], color[1], color[2] };

                for (auto& a : structure.atoms)
                {
                    if (a.symbol == selectedElementSymbol)
                    {
                        a.r = color[0];
                        a.g = color[1];
                        a.b = color[2];
                    }
                }

                updateBuffers(structure);
            }
        }
        else
        {
            ImGui::TextDisabled("No element selected.");
        }

        ImGui::EndPopup();
    }
    if (!editElementColorsOpen)
        ImGui::CloseCurrentPopup();

    drawNotifications();
}

void FileBrowser::drawAtomPropertySettings(Structure& structure, const std::function<void(Structure&)>& updateBuffers)
{
    if (atomColorMode != AtomColorMode::AtomProperty) return;
    if (structure.atomProperty.size() != structure.atoms.size() || structure.atoms.empty()) {
        // The property no longer matches the atoms (edited or replaced structure).
        atomColorMode = AtomColorMode::ElementType;
        atomColorModeJustChanged = true;
        updateBuffers(structure);
        return;
    }
    responsive::windowSize(ImVec2(360, 250), ImGuiCond_FirstUseEver);
    // At the top left of the view, clear of the dialog that produced the values.
    ImGui::SetNextWindowPos(ImVec2(responsive::dp(16), ImGui::GetFrameHeight() * 4.5f), ImGuiCond_FirstUseEver);
    if (!responsive::begin("Atom property colouring", nullptr, ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    ImGui::TextWrapped("%s", structure.atomPropertyName.empty() ? "Per-atom property" : structure.atomPropertyName.c_str());
    const auto range = atomforge::science::finiteRange(structure.atomProperty);
    std::size_t valid = 0;
    for (double v : structure.atomProperty) valid += std::isfinite(v);
    ImGui::TextDisabled("%zu valid of %zu atoms; data range %.4g to %.4g", valid, structure.atomProperty.size(), range.low, range.high);
    bool changed = ImGui::Checkbox("Automatic colour range", &propertyDisplay.autoRange);
    if (propertyDisplay.autoRange) propertyDisplay.range = range;
    ImGui::BeginDisabled(propertyDisplay.autoRange);
    changed |= ImGui::InputDouble("Minimum", &propertyDisplay.range.low, 0, 0, "%.6g");
    changed |= ImGui::InputDouble("Maximum", &propertyDisplay.range.high, 0, 0, "%.6g");
    ImGui::EndDisabled();
    changed |= ImGui::Checkbox("Hide atoms outside the range", &propertyDisplay.hideOutside);
    changed |= ImGui::Checkbox("Hide atoms without a value", &propertyDisplay.hideInvalid);
    if (responsive::button("Show element colours")) {
        atomColorMode = AtomColorMode::ElementType;
        atomColorModeJustChanged = true;
        changed = true;
    }
    ImGui::End();
    if (changed) updateBuffers(structure);
}

void FileBrowser::applyElementColorOverrides(Structure& structure) const
{
    for (std::size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const int atomicNumber = structure.atoms[i].atomicNumber;
        std::unordered_map<int, std::array<float, 3>>::const_iterator it = elementColorOverrides.find(atomicNumber);
        if (it == elementColorOverrides.end())
            continue;

        structure.atoms[i].r = it->second[0];
        structure.atoms[i].g = it->second[1];
        structure.atoms[i].b = it->second[2];
    }
}

void FileBrowser::initNanoCrystalRenderResources(Renderer& renderer)
{
    nanoCrystalDialog.initRenderResources(renderer);
}

bool FileBrowser::isNanoCrystalDialogOpen() const
{
    return nanoCrystalDialog.isOpen();
}

void FileBrowser::feedDropToNanoCrystalDialog(const std::string& path)
{
    nanoCrystalDialog.feedDroppedFile(path);
}

bool FileBrowser::isCustomStructureDialogOpen() const
{
    return customStructureDialog.isOpen();
}

void FileBrowser::feedDropToCustomStructureDialog(const std::string& path)
{
    customStructureDialog.feedDroppedFile(path);
}

bool FileBrowser::isMergeStructuresDialogOpen() const
{
    return mergeStructuresDialog.isOpen();
}

void FileBrowser::feedDropToMergeStructuresDialog(const std::string& path)
{
    mergeStructuresDialog.feedDroppedFile(path);
}

void FileBrowser::initMergeStructuresRenderResources(Renderer& renderer)
{
    mergeStructuresDialog.initRenderResources(renderer);
}

void FileBrowser::initCustomStructureRenderResources(Renderer& renderer)
{
    customStructureDialog.initRenderResources(renderer);
}

void FileBrowser::initInterfaceBuilderRenderResources(Renderer& renderer)
{
    interfaceBuilderDialog.initRenderResources(renderer);
}

void FileBrowser::initCSLGrainBoundaryRenderResources(Renderer& renderer)
{
    cslDialog.initRenderResources(renderer);
}

bool FileBrowser::isCSLGrainBoundaryDialogOpen() const
{
    return cslDialog.isOpen();
}

void FileBrowser::feedDropToCSLGrainBoundaryDialog(const std::string& path)
{
    cslDialog.feedDroppedFile(path);
}

bool FileBrowser::isInterfaceBuilderDialogOpen() const
{
    return interfaceBuilderDialog.isOpen();
}

void FileBrowser::feedDropToInterfaceBuilderDialog(const std::string& path)
{
    interfaceBuilderDialog.feedDroppedFile(path);
}

void FileBrowser::initPolyCrystalRenderResources(Renderer& renderer)
{
    polyCrystalDialog.initRenderResources(renderer);
}

bool FileBrowser::isPolyCrystalDialogOpen() const
{
    return polyCrystalDialog.isOpen();
}

void FileBrowser::feedDropToPolyCrystalDialog(const std::string& path)
{
    polyCrystalDialog.feedDroppedFile(path);
}

void FileBrowser::initStackingFaultRenderResources(Renderer& renderer)
{
    stackingFaultDialog.initRenderResources(renderer);
}

void FileBrowser::initDislocationRenderResources(Renderer& renderer)
{
    dislocationDialog.initRenderResources(renderer);
}

bool FileBrowser::isStackingFaultDialogOpen() const
{
    return stackingFaultDialog.isOpen();
}

bool FileBrowser::isDislocationDialogOpen() const
{
    return dislocationDialog.isOpen();
}

void FileBrowser::feedDropToStackingFaultDialog(const std::string& path)
{
    stackingFaultDialog.feedDroppedFile(path);
}

void FileBrowser::feedDropToDislocationDialog(const std::string& path)
{
    dislocationDialog.feedDroppedFile(path);
}

void FileBrowser::initSubstitutionalSolidSolutionRenderResources(Renderer& renderer)
{
    substitutionalSolidSolutionDialog.initRenderResources(renderer);
}

bool FileBrowser::isSubstitutionalSolidSolutionDialogOpen() const
{
    return substitutionalSolidSolutionDialog.isOpen();
}

bool FileBrowser::isElectronicDialogOpen() const
{
    return electronicDialog.isOpen();
}

void FileBrowser::feedDropToElectronicDialog(const std::string& path)
{
    electronicDialog.feedDroppedFile(path);
}

bool FileBrowser::isAnyDialogOpen() const
{
    return openStructurePopup
        || saveStructurePopup
        || exportImagePopup
        || loadErrorPopupRequested
        || isNanoCrystalDialogOpen()
        || isCustomStructureDialogOpen()
        || isMergeStructuresDialogOpen()
        || isCSLGrainBoundaryDialogOpen()
        || isInterfaceBuilderDialogOpen()
        || isPolyCrystalDialogOpen()
        || isDislocationDialogOpen()
    #if ATOMFORGE_ENABLE_SFE_BUILDER
        || isStackingFaultDialogOpen()
    #endif
    #if ATOMFORGE_ENABLE_SSS_BUILDER
        || isSubstitutionalSolidSolutionDialogOpen()
    #endif
        || isAmorphousBuilderDialogOpen()
        || isInterstitialAtomsDialogOpen()
        || isCellSculptorDialogOpen()
        || isElectronicDialogOpen()
        || showAbout
        || showManual
        || showEditColors
        || showLatticePlanesDialog
        || showMillerDirectionsDialog
        || showPolyhedralSettingsDialog;
}

void FileBrowser::feedDropToSubstitutionalSolidSolutionDialog(const std::string& path)
{
    substitutionalSolidSolutionDialog.feedDroppedFile(path);
}

bool FileBrowser::isAmorphousBuilderDialogOpen() const
{
    return amorphousBuilderDialog.isOpen();
}

bool FileBrowser::isInterstitialAtomsDialogOpen() const
{
    return interstitialAtomsDialog.isOpen();
}

void FileBrowser::feedDropToInterstitialAtomsDialog(const std::string& path)
{
    interstitialAtomsDialog.feedDroppedFile(path);
}

void FileBrowser::drawInterstitialVoidOverlay(ImDrawList* drawList,
                                              const glm::mat4& projection,
                                              const glm::mat4& view,
                                              int framebufferWidth,
                                              int framebufferHeight,
                                              const Structure& structure) const
{
    interstitialAtomsDialog.drawVoidOverlay(drawList,
                                            projection,
                                            view,
                                            framebufferWidth,
                                            framebufferHeight,
                                            structure);
}

void FileBrowser::initCellSculptorRenderResources(Renderer& renderer)
{
    cellSculptorDialog.initRenderResources(renderer);
}

void FileBrowser::initInterstitialAtomsRenderResources(Renderer& renderer)
{
    interstitialAtomsDialog.initRenderResources(renderer);
}

bool FileBrowser::isCellSculptorDialogOpen() const
{
    return cellSculptorDialog.isOpen();
}

void FileBrowser::feedDropToCellSculptorDialog(const std::string& path)
{
    cellSculptorDialog.feedDroppedFile(path);
}

void FileBrowser::showLoadError(const std::string& message)
{
    std::snprintf(loadPopupTitle, sizeof(loadPopupTitle), "%s", "Load Error");
    std::snprintf(
        loadErrorMsg,
        sizeof(loadErrorMsg),
        "%s",
        message.empty() ? "Failed to load file." : message.c_str());
    loadErrorPopupRequested = true;
    showNotification(loadErrorMsg, true);
}

void FileBrowser::showLoadInfo(const std::string& message)
{
    showNotification(message.empty() ? "Load completed." : message, false);
}

void FileBrowser::showNotification(const std::string& message, bool isError)
{
    if (message.empty())
        return;

    ToastNotification notification;
    notification.message = message;
    notification.expiresAt = ImGui::GetTime() + kNotificationLifetimeSeconds;
    notification.isError = isError;
    notifications.push_back(notification);
}

void FileBrowser::drawNotifications()
{
    const double now = ImGui::GetTime();
    notifications.erase(
        std::remove_if(notifications.begin(), notifications.end(),
            [&](const ToastNotification& notification) {
                return notification.expiresAt <= now;
            }),
        notifications.end());

    if (notifications.empty())
        return;

    // Cards in the theme's popup colours, stacked upwards from above the
    // status bar, with a status-coloured accent bar; each fades in and out.
    ImDrawList* drawList = ImGui::GetForegroundDrawList();
    const ImGuiIO& io = ImGui::GetIO();
    const float margin = responsive::dp(16.0f);
    const float padX = responsive::dp(14.0f), padY = responsive::dp(10.0f), accentW = responsive::dp(4.0f);
    const float gap = responsive::dp(8.0f);
    const float wrapWidth = responsive::dp(380.0f);
    const float rounding = responsive::dp(8.0f);
    float bottom = io.DisplaySize.y - ImGui::GetFrameHeight() - margin;
    for (auto it = notifications.rbegin(); it != notifications.rend(); ++it)
    {
        const ToastNotification& notification = *it;
        const double age = kNotificationLifetimeSeconds - (notification.expiresAt - now);
        const float alpha = static_cast<float>(std::clamp(std::min(age / 0.15, (notification.expiresAt - now) / 0.4), 0.0, 1.0));
        const ImVec2 textSize = ImGui::CalcTextSize(notification.message.c_str(), nullptr, false, wrapWidth);
        const float boxW = textSize.x + accentW + padX * 2.0f;
        const float boxH = textSize.y + padY * 2.0f;
        const float x = io.DisplaySize.x - margin - boxW;
        const float y = bottom - boxH;
        const auto colour = [alpha](ImVec4 c, float a = 1.0f) { c.w *= a * alpha; return ImGui::ColorConvertFloat4ToU32(c); };
        const ImVec4 accent = notification.isError ? themeStatusError() : themeStatusGood();
        drawList->AddRectFilled(ImVec2(x + 2, y + 3), ImVec2(x + boxW + 2, y + boxH + 3), colour(ImVec4(0, 0, 0, 1), 0.12f), rounding);
        drawList->AddRectFilled(ImVec2(x, y), ImVec2(x + boxW, y + boxH), colour(ImGui::GetStyleColorVec4(ImGuiCol_PopupBg)), rounding);
        drawList->AddRectFilled(ImVec2(x, y), ImVec2(x + accentW, y + boxH), colour(accent), rounding, ImDrawFlags_RoundCornersLeft);
        drawList->AddRect(ImVec2(x, y), ImVec2(x + boxW, y + boxH), colour(ImGui::GetStyleColorVec4(ImGuiCol_Border)), rounding);
        drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize(), ImVec2(x + accentW + padX, y + padY),
                          colour(ImGui::GetStyleColorVec4(ImGuiCol_Text)), notification.message.c_str(), nullptr, wrapWidth);
        bottom = y - gap;
    }
}

bool FileBrowser::atomColorModeChanged()
{
    bool changed = atomColorModeJustChanged;
    atomColorModeJustChanged = false;
    return changed;
}

void FileBrowser::pushHistory(const std::string& dir)
{
    pushDirectoryHistory(dirHistory, historyIndex, dir);
}

std::string FileBrowser::toLower(const std::string& s)
{
    return toLowerStr(s);
}

bool FileBrowser::isAllowedFile(const std::string& name) const
{
    auto dot = name.find_last_of('.');
    if (dot == std::string::npos)
        return false;

    std::string ext = toLower(name.substr(dot));
    for (const auto& allowed : allowedExtensions)
    {
        if (ext == allowed)
            return true;
    }
    return false;
}

void FileBrowser::openPdfManual()
{
    const auto path = findManualPdf(applicationDirectory());
    if (path.empty())
    {
        showNotification("Manual PDF not found. Restore AtomForge-manual.pdf from the application package.", true);
        return;
    }
    const auto open = ImGui::GetPlatformIO().Platform_OpenInShellFn;
    if (!open || !open(ImGui::GetCurrentContext(), path.u8string().c_str()))
        showNotification("Unable to open the manual. Install or select a default PDF viewer.", true);
}

atomforge::Workspace FileBrowser::workspace(const Structure& structure) const
{
    auto saved=electronicDialog.snapshot(); saved.structure=structure;
    saved.settings["desktop.useLightTheme"]=isLightThemeEnabled();
    saved.settings["desktop.showElementLabels"]=showElementLabels;
    saved.settings["desktop.showBonds"]=showBonds;
    saved.settings["desktop.showAtoms"]=showAtoms;
    saved.settings["desktop.showBoundingBox"]=showBoundingBox;
    saved.settings["desktop.showDislocationLines"]=showDislocationLines;
    saved.settings["desktop.showLatticePlanes"]=showLatticePlanes;
    saved.settings["desktop.showMillerDirections"]=showMillerDirections;
    saved.settings["desktop.showVoronoi"]=showVoronoi;
    saved.settings["desktop.showPolyhedralViewer"]=showPolyhedralViewer;
    saved.settings["desktop.bondElementFilterEnabled"]=bondElementFilterEnabled;
    saved.settings["desktop.viewMode"]=static_cast<int>(viewMode);
    saved.settings["desktop.atomColorMode"]=static_cast<int>(atomColorMode);
    saved.settings["desktop.atomDisplayMode"]=static_cast<int>(atomDisplayMode);
    saved.settings["desktop.property.autoRange"]=propertyDisplay.autoRange;
    saved.settings["desktop.property.low"]=propertyDisplay.range.low;
    saved.settings["desktop.property.high"]=propertyDisplay.range.high;
    saved.settings["desktop.property.hideOutside"]=propertyDisplay.hideOutside;
    saved.settings["desktop.property.hideInvalid"]=propertyDisplay.hideInvalid;
    saved.science=scientificToolsDialog.snapshot();
    saved.pipeline=pipelineDialog.snapshot();
    if (pipelineDialog.isActive()) {
        saved.hasPipelineInput=true;
        saved.pipelineInput=pipelineDialog.editor().input();
    }
    return saved;
}

void FileBrowser::restoreWorkspace(const atomforge::Workspace& saved, Structure& structure)
{
    electronicDialog.restore(saved); structure=saved.structure;
    const auto read=[&](const char* key,double fallback) {
        const auto found=saved.settings.find(key); return found==saved.settings.end() ? fallback : found->second;
    };
    showElementLabels=read("desktop.showElementLabels",showElementLabels)!=0;
    showBonds=read("desktop.showBonds",showBonds)!=0;
    showAtoms=read("desktop.showAtoms",showAtoms)!=0;
    showBoundingBox=read("desktop.showBoundingBox",showBoundingBox)!=0;
    showDislocationLines=read("desktop.showDislocationLines",showDislocationLines)!=0;
    showLatticePlanes=read("desktop.showLatticePlanes",showLatticePlanes)!=0;
    showMillerDirections=read("desktop.showMillerDirections",showMillerDirections)!=0;
    showVoronoi=read("desktop.showVoronoi",showVoronoi)!=0;
    showPolyhedralViewer=read("desktop.showPolyhedralViewer",showPolyhedralViewer)!=0;
    bondElementFilterEnabled=read("desktop.bondElementFilterEnabled",bondElementFilterEnabled)!=0;
    viewMode=static_cast<ViewMode>(static_cast<int>(std::clamp(read("desktop.viewMode",0),0.0,1.0)));
    atomColorMode=static_cast<AtomColorMode>(static_cast<int>(std::clamp(read("desktop.atomColorMode",0),0.0,3.0)));
    if (atomColorMode==AtomColorMode::AtomProperty && structure.atomProperty.size()!=structure.atoms.size())
        atomColorMode=AtomColorMode::ElementType;
    propertyDisplay.autoRange=read("desktop.property.autoRange",propertyDisplay.autoRange)!=0;
    propertyDisplay.range.low=read("desktop.property.low",propertyDisplay.range.low);
    propertyDisplay.range.high=read("desktop.property.high",propertyDisplay.range.high);
    propertyDisplay.hideOutside=read("desktop.property.hideOutside",propertyDisplay.hideOutside)!=0;
    propertyDisplay.hideInvalid=read("desktop.property.hideInvalid",propertyDisplay.hideInvalid)!=0;
    scientificToolsDialog.restore(saved.science);
    pipelineDialog.restore(saved.pipeline, saved.hasPipelineInput ? &saved.pipelineInput : nullptr);
    atomDisplayMode=static_cast<AtomDisplayMode>(static_cast<int>(std::clamp(read("desktop.atomDisplayMode",0),0.0,3.0)));
}

ProjectRequest FileBrowser::drawProjectPicker()
{
    if (projectAction==3) { projectAction=0; return {3,{}}; }
    if (auto path=projectPicker.draw()) { const int action=projectAction; projectAction=0; return {action,*path}; }
    return {};
}

void FileBrowser::openStepDialog(StepEdit request, EditMenuDialogs& editMenuDialogs)
{
    // The same dialog as the menu entry of the step (ui/MenuParity), opened on the step.
    const std::string step = request.step;
    if (step == "build-bulk") bulkCrystalDialog.editStep(std::move(request));
    else if (step == "build-gb") cslDialog.editStep(std::move(request));
    else if (step == "build-nano") nanoCrystalDialog.editStep(std::move(request));
    else if (step == "build-custom") customStructureDialog.editStep(std::move(request));
    else if (step == "build-interface") interfaceBuilderDialog.editStep(std::move(request));
    else if (step == "build-poly") polyCrystalDialog.editStep(std::move(request));
    else if (step == "build-amorphous") amorphousBuilderDialog.editStep(std::move(request));
    else if (step == "build-dislocation") dislocationDialog.editStep(std::move(request));
#if ATOMFORGE_ENABLE_SFE_BUILDER
    else if (step == "build-stacking-fault") stackingFaultDialog.editStep(std::move(request));
#endif
#if ATOMFORGE_ENABLE_SSS_BUILDER
    else if (step == "build-sss") substitutionalSolidSolutionDialog.editStep(std::move(request));
#endif
    else if (step == "add-atom" || step == "set-cell") editMenuDialogs.editStep(std::move(request));
    else if (step == "insert-interstitials") interstitialAtomsDialog.editStep(std::move(request));
    else if (step == "supercell") transformDialog.editStep(std::move(request));
    else if (step == "merge") mergeStructuresDialog.editStep(std::move(request));
    else if (step == "build-sculpt") cellSculptorDialog.editStep(std::move(request));
    else if (step == "relax" || step == "nvt-dynamics" || step == "npt-dynamics" || step == "standardize-cell")
        scientificToolsDialog.editStep(std::move(request));
    // Steps without a dialog of their own use the operation window of their menu entry.
    else operationDialog.editStep(std::move(request));
}

std::vector<std::string> FileBrowser::recentFiles()
{
    loadSharedRecentFilesIfNeeded();
    return sharedRecentFiles();
}

bool FileBrowser::isLightThemeEnabled() const { return isLightTheme(); }
