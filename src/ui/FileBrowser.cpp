#include "FileBrowser.h"
#include "ElementData.h"
#include "util/PathUtils.h"
#include "app/StructureFileService.h"
#include "ui/PeriodicTableDialog.h"

#include "imgui.h"
#include <GLFW/glfw3.h>
#include <openbabel3/openbabel/elements.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
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

struct ImageExportFormatOption { const char* label; const char* ext; ImageExportFormat fmt; };
static const ImageExportFormatOption kImageExportFormats[] = {
    { "PNG (.png)", ".png", ImageExportFormat::Png },
    { "JPEG (.jpg)", ".jpg", ImageExportFormat::Jpg },
    { "SVG (.svg)", ".svg", ImageExportFormat::Svg },
};
static constexpr int kNumImageExportFormats = (int)(sizeof(kImageExportFormats) / sizeof(kImageExportFormats[0]));
static constexpr double kNotificationLifetimeSeconds = 3.5;

FileBrowser::FileBrowser()
        : showAbout(false),
            showManual(false),
            showEditColors(false),
            showElementLabels(false),
            showBonds(false),
            showLatticePlanes(false),
            showLatticePlanesDialog(false),
            bondElementFilterEnabled(false),
            viewMode(ViewMode::Orthographic),
            atomColorMode(AtomColorMode::ElementType),
            atomColorModeJustChanged(false),
            boxSelectMode(false),
            requestMeasureDistance(false),
            requestMeasureAngle(false),
            requestAtomInfo(false),
            requestViewAxisX(false),
            requestViewAxisY(false),
            requestViewAxisZ(false),
            requestViewLatticeA(false),
            requestViewLatticeB(false),
            requestViewLatticeC(false),
            requestResetDefaultView(false),
            requestStructureInfo(false),
            requestUndo(false),
            requestRedo(false),
            requestCloseStructure(false),
            requestImageExport(false),
            openStructurePopup(false),
            saveStructurePopup(false),
            exportImagePopup(false),
            loadErrorPopupRequested(false),
            openDir("."),
      historyIndex(-1),
      saveDir("."),
      saveHistoryIndex(-1),
      selectedSaveFormat(0),
        exportDir("."),
        exportHistoryIndex(-1),
        selectedExportFormat(0),
        exportIncludeBackground(true),
            selectedAtomicNumber(1),
            latticePlaneInputH(1),
            latticePlaneInputK(0),
            latticePlaneInputL(0),
            latticePlaneInputOffset(1.0f),
            latticePlaneInputOpacity(0.20f)
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
    latticePlaneInputColor[0] = 0.95f;
    latticePlaneInputColor[1] = 0.62f;
    latticePlaneInputColor[2] = 0.20f;
    bondElementFilterMask.fill(false);
    updateBondElementFilterMask();
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

void FileBrowser::initFromPath(const std::string& initialPath)
{
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
}

void FileBrowser::draw(Structure& structure,
                       EditMenuDialogs& editMenuDialogs,
                       const std::function<void(Structure&)>& updateBuffers,
                       bool canUndo,
                       bool canRedo)
{
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("Open",  "Ctrl+O"))
                openStructurePopup = true;

            if (ImGui::MenuItem("Close", "Ctrl+W", false, !structure.atoms.empty()))
            {
                requestCloseStructure = true;
                std::cout << "[Operation] Close structure requested" << std::endl;
            }

            if (ImGui::MenuItem("Save As", "Ctrl+S"))
            {
                saveStructurePopup = true;
                saveDir = openDir;
                saveDirHistory = dirHistory;
                saveHistoryIndex = historyIndex;
                saveStatusMsg[0] = '\0';
            }

            if (ImGui::MenuItem("Export Image", "Ctrl+Shift+S", false, !structure.atoms.empty()))
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

            if (ImGui::MenuItem("Quit"))
                glfwSetWindowShouldClose(glfwGetCurrentContext(), true);

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Edit"))
        {
            if (ImGui::MenuItem("Undo", "Ctrl+Z", false, canUndo))
                requestUndo = true;
            if (ImGui::MenuItem("Redo", "Ctrl+Y", false, canRedo))
                requestRedo = true;
            ImGui::Separator();
            ImGui::MenuItem("Box Select Mode", nullptr, &boxSelectMode);
            ImGui::Separator();
            editMenuDialogs.drawMenuItems();
            ImGui::Separator();
            transformDialog.drawMenuItem(structure.hasUnitCell);

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Build"))
        {
            bulkCrystalDialog.drawMenuItem(true);
            cslDialog.drawMenuItem(true);
            nanoCrystalDialog.drawMenuItem(true);
            interfaceBuilderDialog.drawMenuItem(true);
            polyCrystalDialog.drawMenuItem(true);
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("View"))
        {
            ImGui::MenuItem("Show Element", nullptr, &showElementLabels);
            bool bondFilterChanged = false;
            ImGui::MenuItem("Show Bonds", nullptr, &showBonds);
            if (ImGui::MenuItem("Filter Bonds By Elements", nullptr, &bondElementFilterEnabled))
                bondFilterChanged = true;
            if (bondElementFilterEnabled)
            {
                ImGui::SetNextItemWidth(240.0f);
                if (ImGui::InputText("Bond Elements##filter", bondElementFilterInput, sizeof(bondElementFilterInput)))
                {
                    updateBondElementFilterMask();
                    bondFilterChanged = true;
                }
                ImGui::TextDisabled("Comma-separated symbols (e.g. O,F,Cl)");
            }
            if (bondFilterChanged)
                updateBuffers(structure);
            ImGui::Separator();

            if (ImGui::BeginMenu("Color Structure By"))
            {
                if (ImGui::MenuItem("Element Type", nullptr,
                                    atomColorMode == AtomColorMode::ElementType))
                {
                    if (atomColorMode != AtomColorMode::ElementType)
                    {
                        atomColorMode = AtomColorMode::ElementType;
                        atomColorModeJustChanged = true;
                        updateBuffers(structure);
                    }
                }
                if (ImGui::MenuItem("Crystal Orientation", nullptr,
                                    atomColorMode == AtomColorMode::CrystalOrientation))
                {
                    if (atomColorMode != AtomColorMode::CrystalOrientation)
                    {
                        atomColorMode = AtomColorMode::CrystalOrientation;
                        atomColorModeJustChanged = true;
                        updateBuffers(structure);
                    }
                }
                if (ImGui::MenuItem("Grain Boundary", nullptr,
                                    atomColorMode == AtomColorMode::GrainBoundary))
                {
                    if (atomColorMode != AtomColorMode::GrainBoundary)
                    {
                        atomColorMode = AtomColorMode::GrainBoundary;
                        atomColorModeJustChanged = true;
                        updateBuffers(structure);
                    }
                }
                ImGui::EndMenu();
            }

            ImGui::Separator();

            ImGui::MenuItem("Show Lattice Planes", nullptr, &showLatticePlanes, structure.hasUnitCell);
            if (ImGui::MenuItem("Lattice Planes", nullptr, false, structure.hasUnitCell))
                showLatticePlanesDialog = true;

            ImGui::Separator();

            const bool isIsometricView = (viewMode == ViewMode::Isometric);
            if (ImGui::MenuItem("Isometric View", nullptr, isIsometricView) && !isIsometricView)
            {
                viewMode = ViewMode::Isometric;
                requestResetDefaultView = true;
            }

            const bool isOrthographicView = (viewMode == ViewMode::Orthographic);
            if (ImGui::MenuItem("Orthographic View", nullptr, isOrthographicView) && !isOrthographicView)
            {
                viewMode = ViewMode::Orthographic;
                requestResetDefaultView = true;
            }

            ImGui::Separator();
            if (ImGui::MenuItem("Structure Info"))
                requestStructureInfo = true;
            if (ImGui::MenuItem("Measure Distance"))
                requestMeasureDistance = true;
            if (ImGui::MenuItem("Measure Angle"))
                requestMeasureAngle = true;
            if (ImGui::MenuItem("Atom Info"))
                requestAtomInfo = true;
            if (ImGui::MenuItem("Reset Default View"))
                requestResetDefaultView = true;
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Analysis"))
        {
            cnaDialog.drawMenuItem(!structure.atoms.empty());
            rdfDialog.drawMenuItem(!structure.atoms.empty());
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Help"))
        {
            if (ImGui::MenuItem("Manual"))
                showManual = true;

            if (ImGui::MenuItem("About"))
                showAbout = true;

            ImGui::EndMenu();
        }

        ImGui::EndMainMenuBar();
    }

    // Draw toolbar below menu bar with axis view and measurement options
    {
        ImGui::SetNextWindowPos(ImVec2(0.0f, ImGui::GetFrameHeight()), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, 0.0f), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("##ViewToolbar", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(12.0f, 4.0f));
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(8.0f, 4.0f));
            
            // View Axis section
            ImGui::AlignTextToFramePadding();
            ImGui::Text("View Axis:");
            ImGui::SameLine(0.0f, 8.0f);
            if (ImGui::Button("X##axis", ImVec2(40.0f, 0.0f)))
                requestViewAxisX = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("Y##axis", ImVec2(40.0f, 0.0f)))
                requestViewAxisY = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("Z##axis", ImVec2(40.0f, 0.0f)))
                requestViewAxisZ = true;

            ImGui::SameLine(0.0f, 16.0f);

            const bool hasInputCell = structure.hasUnitCell && !structure.atoms.empty();
            if (!hasInputCell) ImGui::BeginDisabled();
            if (ImGui::Button("a##latview", ImVec2(34.0f, 0.0f)))
                requestViewLatticeA = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("b##latview", ImVec2(34.0f, 0.0f)))
                requestViewLatticeB = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("c##latview", ImVec2(34.0f, 0.0f)))
                requestViewLatticeC = true;
            if (!hasInputCell) ImGui::EndDisabled();

            ImGui::SameLine(0.0f, 16.0f);
            
            // Measure section
            ImGui::AlignTextToFramePadding();
            ImGui::Text("Measure:");
            ImGui::SameLine(0.0f, 8.0f);
            if (ImGui::Button("Distance##measure", ImVec2(80.0f, 0.0f)))
                requestMeasureDistance = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("Angle##measure", ImVec2(80.0f, 0.0f)))
                requestMeasureAngle = true;
            
            ImGui::PopStyleVar(2);
            ImGui::End();
        }
    }

    editMenuDialogs.drawPopups(structure, updateBuffers);

    transformDialog.drawDialog([&]() { updateBuffers(structure); });
    bulkCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors, updateBuffers);
    cslDialog.drawDialog(structure, editMenuDialogs.elementColors,
                         editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                         updateBuffers);
    nanoCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateBuffers);
    interfaceBuilderDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                      updateBuffers);
    polyCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateBuffers);
    cnaDialog.drawDialog(structure);
    rdfDialog.drawDialog(structure);

    if (loadErrorPopupRequested)
    {
        ImGui::OpenPopup(loadPopupTitle);
        loadErrorPopupRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(720.0f, 0.0f), ImGuiCond_Appearing);
    bool loadErrorOpen = true;
    if (ImGui::BeginPopupModal(loadPopupTitle, &loadErrorOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::TextUnformatted(loadErrorMsg);
        ImGui::Spacing();
        if (ImGui::Button("OK", ImVec2(120, 0)))
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (!loadErrorOpen)
        ImGui::CloseCurrentPopup();

    if (openStructurePopup)
    {
        ImGui::OpenPopup("Open Structure");
        openStructurePopup = false;
        openStatusMsg[0] = '\0';
    }

    bool openStructureOpen = true;
    if (ImGui::BeginPopupModal("Open Structure", &openStructureOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Current folder: %s", openDir.c_str());
        ImGui::SameLine();
        if (ImGui::Button(".."))
        {
            openDir = parentPath(openDir);
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }

        ImGui::SameLine();
        if (ImGui::Button("Back") && historyIndex > 0)
        {
            historyIndex--;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }
        ImGui::SameLine();
        if (ImGui::Button("Forward") && historyIndex + 1 < (int)dirHistory.size())
        {
            historyIndex++;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }

        ImGui::SameLine();
        if (ImGui::Button("Root"))
        {
            openDir = !driveRoots.empty() ? driveRoots.front() : "/";
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }
        ImGui::SameLine();
        if (ImGui::Button("Home"))
        {
            openDir = detectHomePath();
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }

        ImGui::Separator();

        if (ImGui::TreeNodeEx("Supported Formats", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::BulletText("XYZ (.xyz) - XYZ structure format");
            ImGui::BulletText("CIF (.cif) - Crystallographic Information File");
            ImGui::BulletText("PDB (.pdb) - Protein Data Bank format");
            ImGui::BulletText("SDF (.sdf) - Structure Data File");
            ImGui::BulletText("MOL (.mol) - MDL MOL format");
            ImGui::BulletText("VASP (.vasp) - VASP POSCAR format");
            ImGui::BulletText("Mol2 (.mol2) - Sybyl Mol2 format");
            ImGui::BulletText("Quantum ESPRESSO (.pwi) - PWscf input");
            ImGui::BulletText("Gaussian (.gjf) - Gaussian input");
            ImGui::TreePop();
        }

        if (ImGui::BeginChild("##filebrowser", ImVec2(500, 300), true))
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
                    });
            }

            ImGui::EndChild();
        }

        ImGui::InputText("Filename", openFilename, sizeof(openFilename));
        if (openStatusMsg[0] != '\0')
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", openStatusMsg);

        if (ImGui::Button("Load"))
        {
            std::string fullPath = joinPath(openDir, openFilename);
            Structure newStructure;
            std::string loadError;
            if (loadStructureFromPath(fullPath, newStructure, loadError))
            {
                structure = std::move(newStructure);
                latticePlanes.clear();
                showLatticePlanes = false;
                applyElementColorOverrides(structure);
                showLoadInfo(std::string("Structure loaded. ") + structure.ipfLoadStatus);

                updateBuffers(structure);
                requestResetDefaultView = true;
                openStatusMsg[0] = '\0';
                std::cout << "[Operation] Loaded structure: " << fullPath
                          << " (atoms=" << structure.atoms.size() << ")" << std::endl;
                ImGui::CloseCurrentPopup();
            }
            else
            {
                std::snprintf(
                    openStatusMsg,
                    sizeof(openStatusMsg),
                    "%s",
                    loadError.empty() ? "Failed to load file." : loadError.c_str());
                std::cout << "[Operation] Load failed: " << fullPath
                          << " (" << (loadError.empty() ? "Failed to load file." : loadError) << ")" << std::endl;
            }
        }
        ImGui::EndPopup();
    }
    if (!openStructureOpen)
        ImGui::CloseCurrentPopup();

    // ---- Save As dialog ------------------------------------------------
    if (saveStructurePopup)
    {
        ImGui::OpenPopup("Save As");
        saveStructurePopup = false;
    }

    ImGui::SetNextWindowSize(ImVec2(560, 480), ImGuiCond_FirstUseEver);
    bool saveAsOpen = true;
    if (ImGui::BeginPopupModal("Save As", &saveAsOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        // Helper: navigate to a new directory and record it in history.
        auto pushSaveDir = [&](const std::string& dir) {
            saveDir = dir;
            pushDirectoryHistory(saveDirHistory, saveHistoryIndex, saveDir);
        };

        ImGui::Text("Current folder: %s", saveDir.c_str());
        ImGui::SameLine();
        if (ImGui::Button("..##save"))
        {
            pushSaveDir(parentPath(saveDir));
        }
        ImGui::SameLine();
        if (ImGui::Button("Back##save") && saveHistoryIndex > 0)
        {
            saveHistoryIndex--;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        ImGui::SameLine();
        if (ImGui::Button("Forward##save") && saveHistoryIndex + 1 < (int)saveDirHistory.size())
        {
            saveHistoryIndex++;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        ImGui::SameLine();
        if (ImGui::Button("Root##save"))
            pushSaveDir(!driveRoots.empty() ? driveRoots.front() : "/");
        ImGui::SameLine();
        if (ImGui::Button("Home##save"))
            pushSaveDir(detectHomePath());

        ImGui::Separator();

        if (ImGui::BeginChild("##savefilebrowser", ImVec2(500, 200), true))
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

        ImGui::InputText("Filename##save", saveFilename, sizeof(saveFilename));

        // Format selector
        if (ImGui::Combo("Format", &selectedSaveFormat,
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
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", saveStatusMsg);

        ImGui::Separator();

        if (ImGui::Button("Save"))
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
    // ---- end Save As dialog --------------------------------------------

    // ---- Export Image dialog -------------------------------------------
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

    ImGui::SetNextWindowSize(ImVec2(560, 500), ImGuiCond_FirstUseEver);
    bool exportImageOpen = true;
    if (ImGui::BeginPopupModal("Export Image", &exportImageOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        auto pushExportDir = [&](const std::string& dir) {
            exportDir = dir;
            pushDirectoryHistory(exportDirHistory, exportHistoryIndex, exportDir);
        };

        ImGui::Text("Current folder: %s", exportDir.c_str());
        ImGui::SameLine();
        if (ImGui::Button("..##export"))
            pushExportDir(parentPath(exportDir));

        ImGui::SameLine();
        if (ImGui::Button("Back##export") && exportHistoryIndex > 0)
        {
            exportHistoryIndex--;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        ImGui::SameLine();
        if (ImGui::Button("Forward##export") && exportHistoryIndex + 1 < (int)exportDirHistory.size())
        {
            exportHistoryIndex++;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        ImGui::SameLine();
        if (ImGui::Button("Root##export"))
            pushExportDir(!driveRoots.empty() ? driveRoots.front() : "/");
        ImGui::SameLine();
        if (ImGui::Button("Home##export"))
            pushExportDir(detectHomePath());

        ImGui::Separator();

        if (ImGui::BeginChild("##exportfilebrowser", ImVec2(500, 200), true))
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

        ImGui::InputText("Filename##export", exportFilename, sizeof(exportFilename));

        if (ImGui::Combo("Format##export", &selectedExportFormat,
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

        if (exportStatusMsg[0] != '\0')
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", exportStatusMsg);

        ImGui::Separator();

        if (ImGui::Button("Export"))
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
                requestImageExport = true;

                std::cout << "[Operation] Image export requested: " << pendingImageExport.outputPath
                          << " (format=" << selectedExt << ", background="
                          << (exportIncludeBackground ? "on" : "off") << ")" << std::endl;

                showNotification(std::string("Image export started: ") + finalName, false);
                exportStatusMsg[0] = '\0';
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
    if (!exportImageOpen)
        ImGui::CloseCurrentPopup();
    // ---- end Export Image dialog ---------------------------------------

    if (showLatticePlanesDialog)
    {
        ImGui::OpenPopup("Lattice Planes");
        showLatticePlanesDialog = false;
    }

    ImGui::SetNextWindowSize(ImVec2(780.0f, 520.0f), ImGuiCond_FirstUseEver);
    bool latticePlanesOpen = true;
    if (ImGui::BeginPopupModal("Lattice Planes", &latticePlanesOpen, ImGuiWindowFlags_NoResize))
    {
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("Current structure has no unit cell. Lattice planes are unavailable.");
            if (ImGui::Button("Close", ImVec2(120, 0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
        else
        {
            ImGui::Checkbox("Show lattice planes", &showLatticePlanes);
            ImGui::TextDisabled("Live edit enabled: HKL and offset updates are shown immediately.");
            ImGui::Separator();

            ImGui::Text("New plane:");
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("H##plane", &latticePlaneInputH, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("K##plane", &latticePlaneInputK, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("L##plane", &latticePlaneInputL, 0.2f, -50, 50);

            ImGui::SetNextItemWidth(210.0f);
            ImGui::DragFloat("Offset n##plane", &latticePlaneInputOffset, 0.01f, -50.0f, 50.0f, "%.3f");
            ImGui::SetNextItemWidth(210.0f);
            ImGui::SliderFloat("Opacity##plane", &latticePlaneInputOpacity, 0.0f, 1.0f, "%.2f");
            ImGui::SameLine();
            ImGui::ColorEdit3("Color##plane", latticePlaneInputColor, ImGuiColorEditFlags_NoInputs);

            const bool invalidMiller =
                (latticePlaneInputH == 0 && latticePlaneInputK == 0 && latticePlaneInputL == 0);
            if (invalidMiller)
                ImGui::TextDisabled("(h, k, l) cannot all be zero.");

            if (ImGui::Button("Add Plane", ImVec2(140.0f, 0.0f)) && !invalidMiller)
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
            else if (ImGui::BeginChild("##lattice-plane-list", ImVec2(0.0f, 220.0f), true))
            {
                int deleteIndex = -1;
                for (size_t i = 0; i < latticePlanes.size(); ++i)
                {
                    ImGui::PushID((int)i + 30000);
                    LatticePlane& plane = latticePlanes[i];

                    ImGui::Checkbox("##visible", &plane.visible);
                    ImGui::SameLine();
                    ImGui::Text("Plane %d", (int)i + 1);

                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("H##plane-h", &plane.h, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("K##plane-k", &plane.k, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("L##plane-l", &plane.l, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(120.0f);
                    ImGui::DragFloat("n##plane-offset", &plane.offset, 0.01f, -50.0f, 50.0f, "%.3f");
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(95.0f);
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
                    if (ImGui::Button("Delete"))
                        deleteIndex = (int)i;

                    ImGui::Separator();

                    ImGui::PopID();
                }

                if (deleteIndex >= 0)
                    latticePlanes.erase(latticePlanes.begin() + deleteIndex);

                ImGui::EndChild();
            }

            ImGui::Spacing();
            if (ImGui::Button("Close", ImVec2(120, 0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }
    if (!latticePlanesOpen)
        ImGui::CloseCurrentPopup();

    if (showManual)
    {
        ImGui::OpenPopup("Manual");
        showManual = false;
    }

    ImGui::SetNextWindowSize(ImVec2(800.0f, 700.0f), ImGuiCond_Appearing);
    bool manualOpen = true;
    if (ImGui::BeginPopupModal("Manual", &manualOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("AtomForge Manual");
        ImGui::TextDisabled("Complete reference for all features and operations");
        ImGui::Separator();

        if (ImGui::BeginChild("##manual-scroll", ImVec2(0.0f, 640.0f), false))
        {
            ImGui::Text("Getting Started");
            ImGui::BulletText("Open a structure from File -> Open or press Ctrl+O.");
            ImGui::BulletText("Use left-drag to rotate the scene and the scroll wheel to zoom.");
            ImGui::BulletText("Use View -> Reset Default View to restore the fitted isometric camera.");

            ImGui::Spacing();
            ImGui::Text("Selection and Editing");
            ImGui::BulletText("Left-click selects one atom.");
            ImGui::BulletText("Ctrl+left-click adds or removes atoms from the current selection.");
            ImGui::BulletText("Ctrl+A selects all atoms. Ctrl+D or Esc clears the selection.");
            ImGui::BulletText("Delete removes the selected atoms from the structure.");
            ImGui::BulletText("Right-click opens the context menu when atoms are selected.");
            ImGui::BulletText("Enable Edit -> Box Select Mode to select with a right-drag screen rectangle.");

            ImGui::Spacing();
            ImGui::Text("Keyboard Shortcuts");
            ImGui::BulletText("Ctrl+O: open structure file.");
            ImGui::BulletText("Ctrl+S: save structure as.");
            ImGui::BulletText("Ctrl+Shift+S: export the current rendered view.");
            ImGui::BulletText("Ctrl+Z: undo. Ctrl+Y or Ctrl+Shift+Z: redo.");

            ImGui::Spacing();
            ImGui::Text("File Menu");
            ImGui::BulletText("Open loads supported structure formats such as CIF, MOL, PDB, XYZ, SDF, VASP, MOL2, PWI, and GJF.");
            ImGui::BulletText("Save As exports the current structure to supported chemistry and crystal formats.");
            ImGui::BulletText("Export Image writes PNG, JPG, or SVG output with optional background.");
            ImGui::BulletText("Close unloads the current structure.");

            ImGui::Spacing();
            ImGui::Text("Edit Menu");
            ImGui::BulletText("Undo and Redo track structure edits and style changes.");
            ImGui::BulletText("Edit Structure modifies lattice vectors and the atom list directly.");
            ImGui::BulletText("Atomic Sizes adjusts per-element radii used for display and some builders.");
            ImGui::BulletText("Element Colors adjusts per-element color and shininess.");
            ImGui::BulletText("Transform Structure applies a 3x3 transformation matrix to periodic structures.");

            ImGui::Spacing();
            ImGui::Text("Build Menu");
            ImGui::BulletText("Bulk Crystal builds a full periodic unit cell from crystal system, space group, lattice parameters, and asymmetric-unit atoms.");
            ImGui::BulletText("CSL Grain Boundary builds cubic bicrystals from ideal sc, bcc, fcc, or diamond source lattices.");
            ImGui::BulletText("Nanocrystal carves a finite particle from a loaded reference structure using sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, or cuboctahedron shapes.");
            ImGui::BulletText("Polycrystal builds a Voronoi-based polycrystalline microstructure from a loaded reference single crystal with controllable box size, grain count, seed, and grain orientations.");

            ImGui::Spacing();
            ImGui::Text("Bulk Crystal Builder");
            ImGui::BulletText("Builds a full periodic unit cell from crystal system, space group, lattice parameters, and asymmetric-unit atoms.");
            ImGui::BulletText("Useful for creating reference single crystals and simulation-ready periodic cells.");

            ImGui::Spacing();
            ImGui::Text("CSL Grain Boundary Builder");
            ImGui::BulletText("Builds cubic bicrystals from ideal sc, bcc, fcc, or diamond source lattices.");
            ImGui::BulletText("Controls Sigma selection, grain-boundary plane, replication, overlap removal, and rigid translation.");

            ImGui::Spacing();
            ImGui::Text("Nanocrystal Builder");
            ImGui::BulletText("Uses the currently loaded structure as the reference source for carving.");
            ImGui::BulletText("You can drag and drop a supported structure file into the reference preview area while the dialog is open.");
            ImGui::BulletText("Preview controls: left-drag orbits the preview, scroll zooms the preview camera.");
            ImGui::BulletText("Options include auto-centering, manual center coordinates, auto-replication for periodic inputs, and rectangular output-cell padding.");

            ImGui::Spacing();
            ImGui::Text("Polycrystal Builder");
            ImGui::BulletText("Requires a loaded reference single crystal with unit-cell information.");
            ImGui::BulletText("The left panel accepts drag-and-drop and displays a 3D preview of the reference crystal.");
            ImGui::BulletText("The builder creates grains by Voronoi tessellation inside a user-defined simulation box.");
            ImGui::BulletText("Grain orientations can be all-random, all user-specified with Euler angles, or partially specified.");
            ImGui::BulletText("Generated structures include per-atom IPF-Z orientation colors for crystal-orientation visualization.");

            ImGui::Spacing();
            ImGui::Text("View Menu");
            ImGui::BulletText("Show Element toggles element labels.");
            ImGui::BulletText("Show Bonds toggles bond-cylinder rendering.");
            ImGui::BulletText("Isometric View and Orthographic View switch the camera projection mode.");
            ImGui::BulletText("Color Structure By switches between element-type colors and crystal-orientation IPF coloring.");
            ImGui::BulletText("Structure Info shows composition, lattice metrics, positions, and symmetry when available.");
            ImGui::BulletText("Measure Distance, Measure Angle, and Atom Info open the corresponding dialogs for the current selection.");

            ImGui::Spacing();
            ImGui::Text("Crystal Orientation / IPF");
            ImGui::BulletText("Crystal Orientation coloring uses cubic IPF-Z colors and displays an IPF triangle legend in the main view.");
            ImGui::BulletText("When available, IPF data is restored from a saved sidecar file named basename.atomforge-ipf.");
            ImGui::BulletText("If no saved IPF metadata is present, AtomForge can fall back to geometry-based orientation reconstruction.");

            ImGui::Spacing();
            ImGui::Text("Analysis Menu");
            ImGui::BulletText("Common Neighbour Analysis reports pair signatures and per-atom structural environments.");
            ImGui::BulletText("Radial Distribution Function plots RDF with configurable species filters, normalization, radius range, bin count, and smoothing.");

            ImGui::Spacing();
            ImGui::Text("Context Menu");
            ImGui::BulletText("Substitute Atom replaces the selected atoms with a chosen element.");
            ImGui::BulletText("Insert Atom at Midpoint places a new atom at the centroid of the current selection.");
            ImGui::BulletText("Measure Distance, Measure Angle, Atom Info, Delete, and Deselect are available when selection rules are satisfied.");

            ImGui::Spacing();
            ImGui::Text("Display and Measurement Features");
            ImGui::BulletText("Element labels can be shown for periodic-image atoms as well.");
            ImGui::BulletText("Bonds are inferred from covalent radii and rendered with split element colors.");
            ImGui::BulletText("Selected atoms are highlighted and helper overlays are drawn for distance and angle tools.");
            ImGui::BulletText("Periodic boundary visualization includes duplicated boundary atoms and transformed supercell views when applicable.");

            ImGui::EndChild();
        }

        if (ImGui::Button("Close", ImVec2(120, 0)))
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
    if (!manualOpen)
        ImGui::CloseCurrentPopup();

    if (showAbout)
    {
        ImGui::OpenPopup("About");
        showAbout = false;
    }

    ImGui::SetNextWindowSize(ImVec2(760.0f, 520.0f), ImGuiCond_Appearing);
    bool aboutOpen = true;
    if (ImGui::BeginPopupModal("About", &aboutOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("AtomForge");
        ImGui::TextDisabled("Molecular structure viewer and editor with periodic-cell tools");
        ImGui::Separator();

        if (ImGui::BeginChild("##about-scroll", ImVec2(0.0f, 440.0f), false))
        {
            ImGui::Text("Creator");
            ImGui::BulletText("Albert Linda");

            ImGui::Spacing();
            ImGui::Text("Core Libraries");
            ImGui::BulletText("ImGui: immediate-mode GUI framework (by Omar Cornut)");
            ImGui::BulletText("OpenGL 3.3+: GPU rendering and visualization");
            ImGui::BulletText("GLFW 3: window and input management");
            ImGui::BulletText("GLEW: OpenGL extension loader");
            ImGui::BulletText("GLM: mathematics library for transformations");

            ImGui::Spacing();
            ImGui::Text("Chemistry Libraries");
            ImGui::BulletText("Open Babel: molecular file I/O (.cif, .mol, .pdb, .xyz, .sdf, .mol2, etc.)");
            ImGui::BulletText("spglib: crystallographic space group and symmetry operations");

            ImGui::Spacing();
            ImGui::Text("Data & References");
            ImGui::BulletText("Periodic table element data (atomic numbers, masses, default colors)");
            ImGui::BulletText("Cordero et al. (Dalton Trans. 2008): Van der Waals radii");
            ImGui::BulletText("Standard crystallographic conventions for cell vectors and positions");

            ImGui::Spacing();
            ImGui::Text("Features");
            ImGui::BulletText("Interactive 3D visualization with real-time rendering");
            ImGui::BulletText("Full undo/redo support for all editing operations");
            ImGui::BulletText("Multi-format structure loading and exporting");
            ImGui::BulletText("Crystallographic analysis: space groups, symmetry, cell metrics");
            ImGui::BulletText("Customizable atom colors, sizes, and materials");
            ImGui::BulletText("Distance and angle measurement tools");
            ImGui::BulletText("Supercell transformation with periodic boundary visualization");

            ImGui::EndChild();
        }

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

    ImGui::SetNextWindowSize(ImVec2(950, 600), ImGuiCond_FirstUseEver);
    bool editElementColorsOpen = true;
    if (ImGui::BeginPopupModal("Edit Element Colors", &editElementColorsOpen, ImGuiWindowFlags_NoResize))
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

    ImDrawList* drawList = ImGui::GetForegroundDrawList();
    ImGuiIO& io = ImGui::GetIO();
    const float margin = 16.0f;
    const float boxPadX = 10.0f;
    const float boxPadY = 8.0f;
    const float gap = 10.0f;
    const float wrapWidth = 360.0f;
    float y = 48.0f;

    for (size_t i = 0; i < notifications.size(); ++i)
    {
        const ToastNotification& notification = notifications[i];
        const ImVec2 textSize = ImGui::CalcTextSize(notification.message.c_str(), nullptr, false, wrapWidth);
        const float boxW = std::min(wrapWidth, textSize.x) + boxPadX * 2.0f;
        const float boxH = textSize.y + boxPadY * 2.0f;
        const float x = io.DisplaySize.x - margin - boxW;

        const ImU32 bg = notification.isError
            ? IM_COL32(115, 34, 38, 230)
            : IM_COL32(28, 42, 57, 230);
        const ImU32 border = notification.isError
            ? IM_COL32(230, 111, 114, 255)
            : IM_COL32(120, 200, 255, 255);

        drawList->AddRectFilled(ImVec2(x, y), ImVec2(x + boxW, y + boxH), bg, 6.0f);
        drawList->AddRect(ImVec2(x, y), ImVec2(x + boxW, y + boxH), border, 6.0f, 0, 1.5f);
        drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                  ImVec2(x + boxPadX, y + boxPadY),
                  IM_COL32(245, 245, 245, 255),
                  notification.message.c_str(), nullptr, wrapWidth);
        y += boxH + gap;
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
