#include "FileBrowser.h"
#include "ElementData.h"
#include "ImGuiSetup.h"
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
// Persistent sidebar width shared across all three file-browser dialogs.
static float s_sidebarW = 130.0f;
 
FileBrowser::FileBrowser()
        : useLightTheme(false),
            showAbout(false),
            showManual(false),
            showEditColors(false),
            showElementLabels(false),
            showBonds(false),
            showLatticePlanes(false),
            showLatticePlanesDialog(false),
            showMillerDirections(false),
            showMillerDirectionsDialog(false),
            showVoronoi(false),
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
    latticePlaneInputColor[0] = 0.95f;
    latticePlaneInputColor[1] = 0.62f;
    latticePlaneInputColor[2] = 0.20f;
    millerDirInputColor[0] = 0.25f;
    millerDirInputColor[1] = 0.75f;
    millerDirInputColor[2] = 1.00f;
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
            ImGui::Separator();
            mergeStructuresDialog.drawMenuItem(true);

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Build"))
        {
            bulkCrystalDialog.drawMenuItem(true);
            substitutionalSolidSolutionDialog.drawMenuItem(true);
            cslDialog.drawMenuItem(true);
            nanoCrystalDialog.drawMenuItem(true);
            customStructureDialog.drawMenuItem(true);
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

            ImGui::MenuItem("Show Miller Directions", nullptr, &showMillerDirections, structure.hasUnitCell);
            if (ImGui::MenuItem("Miller Directions", nullptr, false, structure.hasUnitCell))
                showMillerDirectionsDialog = true;

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
            ImGui::MenuItem("Show Voronoi Volume", nullptr, &showVoronoi, !structure.atoms.empty());
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
            ImGui::Separator();
            if (ImGui::BeginMenu("Select Theme"))
            {
                if (ImGui::MenuItem("Dark", nullptr, !useLightTheme))
                {
                    if (useLightTheme) { useLightTheme = false; applyDarkTheme(); }
                }
                if (ImGui::MenuItem("Light", nullptr, useLightTheme))
                {
                    if (!useLightTheme) { useLightTheme = true; applyLightTheme(); }
                }
                ImGui::EndMenu();
            }
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

            // Rotate Crystal section
            const bool hasAtoms = !structure.atoms.empty();
            ImGui::AlignTextToFramePadding();
            ImGui::Text("Rotate:");
            ImGui::SameLine(0.0f, 8.0f);
            ImGui::SetNextItemWidth(68.0f);
            ImGui::DragFloat("##rotangle", &rotateCrystalAngle, 1.0f, -360.0f, 360.0f, "%.1f\xc2\xb0");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Rotation angle in degrees");
            ImGui::SameLine(0.0f, 6.0f);
            if (!hasAtoms) ImGui::BeginDisabled();
            if (ImGui::Button("X##rot", ImVec2(36.0f, 0.0f)))
                requestRotateCrystalX = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("Y##rot", ImVec2(36.0f, 0.0f)))
                requestRotateCrystalY = true;
            ImGui::SameLine(0.0f, 4.0f);
            if (ImGui::Button("Z##rot", ImVec2(36.0f, 0.0f)))
                requestRotateCrystalZ = true;
            if (!hasAtoms) ImGui::EndDisabled();

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

    auto updateFromBuilder = [&](Structure& updatedStructure) {
        updateBuffers(updatedStructure);
        requestResetDefaultView = true;
    };

    bulkCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors, updateFromBuilder);
    cslDialog.drawDialog(structure, editMenuDialogs.elementColors,
                         editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                         updateFromBuilder);
    nanoCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateFromBuilder);
    customStructureDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                      updateFromBuilder);
    mergeStructuresDialog.drawDialog(structure, updateFromBuilder);
    interfaceBuilderDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                      editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                      updateFromBuilder);
    polyCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                 editMenuDialogs.elementRadii, editMenuDialogs.elementShininess,
                                 updateFromBuilder);
    substitutionalSolidSolutionDialog.drawDialog(structure, editMenuDialogs.elementColors,
                                                  editMenuDialogs.elementRadii,
                                                  editMenuDialogs.elementShininess,
                                                  updateFromBuilder);
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

    ImGui::SetNextWindowSize(ImVec2(720.0f, 510.0f), ImGuiCond_Appearing);
    bool openStructureOpen = true;
    if (ImGui::BeginPopupModal("Open Structure", &openStructureOpen, 0))
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
        if (ImGui::Button("\xe2\x86\x90##openBack", ImVec2(navBtnW, 0.0f)) && historyIndex > 0)
        {
            historyIndex--;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x92##openFwd", ImVec2(navBtnW, 0.0f)) && historyIndex + 1 < (int)dirHistory.size())
        {
            historyIndex++;
            openDir = dirHistory[historyIndex];
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x91##openUp", ImVec2(navBtnW, 0.0f)))
        {
            openDir = parentPath(openDir);
            pushHistory(openDir);
            openStatusMsg[0] = '\0';
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(0.0f, 8.0f);
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
        const float listH = 300.0f;

        if (ImGui::BeginChild("##opensidebar", ImVec2(s_sidebarW, listH), true))
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
        ImGui::SameLine(0.0f, 0.0f);
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
        ImGui::SameLine(0.0f, 0.0f);

        bool openFileDoubleClicked = false;
        if (ImGui::BeginChild("##filebrowser", ImVec2(0.0f, listH), true))
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
            static const char* kOpenExtHint = ".xyz  .cif  .pdb  .sdf  .mol  .vasp  .mol2  .pwi  .gjf";
            const float extHintW = ImGui::CalcTextSize(kOpenExtHint).x;
            const float labelW   = ImGui::CalcTextSize("Filename").x
                                  + ImGui::GetStyle().ItemInnerSpacing.x * 2.0f;
            const float spacing  = ImGui::GetStyle().ItemSpacing.x;
            const float inputW   = ImGui::GetContentRegionAvail().x - labelW - extHintW - spacing * 2.0f;
            ImGui::SetNextItemWidth(inputW > 120.0f ? inputW : 120.0f);
            if (ImGui::InputText("Filename##open", openFilename, sizeof(openFilename),
                                 ImGuiInputTextFlags_EnterReturnsTrue))
                openFileDoubleClicked = true;
            ImGui::SameLine(0.0f, spacing);
            ImGui::AlignTextToFramePadding();
            ImGui::TextDisabled("%s", kOpenExtHint);
        }
        if (openStatusMsg[0] != '\0')
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", openStatusMsg);

        bool doLoad = openFileDoubleClicked;
        if (ImGui::Button("Load", ImVec2(120.0f, 0.0f)))
            doLoad = true;
        ImGui::SameLine(0.0f, 8.0f);
        if (ImGui::Button("Cancel##openCancel", ImVec2(120.0f, 0.0f)))
            ImGui::CloseCurrentPopup();

        if (doLoad)
        {
            std::string fullPath = joinPath(openDir, openFilename);
            Structure newStructure;
            std::string loadError;
            if (loadStructureFromPath(fullPath, newStructure, loadError))
            {
                structure = std::move(newStructure);
                latticePlanes.clear();
                showLatticePlanes = false;
                millerDirections.clear();
                showMillerDirections = false;
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

    ImGui::SetNextWindowSize(ImVec2(720.0f, 480.0f), ImGuiCond_Appearing);
    bool saveAsOpen = true;
    if (ImGui::BeginPopupModal("Save As", &saveAsOpen, 0))
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
        if (ImGui::Button("\xe2\x86\x90##saveBack", ImVec2(navBtnW, 0.0f)) && saveHistoryIndex > 0)
        {
            saveHistoryIndex--;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x92##saveFwd", ImVec2(navBtnW, 0.0f)) && saveHistoryIndex + 1 < (int)saveDirHistory.size())
        {
            saveHistoryIndex++;
            saveDir = saveDirHistory[saveHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x91##saveUp", ImVec2(navBtnW, 0.0f)))
            pushSaveDir(parentPath(saveDir));
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(0.0f, 8.0f);
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

        if (ImGui::BeginChild("##savesidebar", ImVec2(s_sidebarW, listH), true))
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
        ImGui::SameLine(0.0f, 0.0f);
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
        ImGui::SameLine(0.0f, 0.0f);

        if (ImGui::BeginChild("##savefilebrowser", ImVec2(0.0f, listH), true))
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

        bool saveEnterPressed = ImGui::InputText("Filename##save", saveFilename, sizeof(saveFilename),
                                                 ImGuiInputTextFlags_EnterReturnsTrue);

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

        bool doSave = saveEnterPressed;
        if (ImGui::Button("Save", ImVec2(120.0f, 0.0f)))
            doSave = true;
        ImGui::SameLine(0.0f, 8.0f);
        if (ImGui::Button("Cancel##saveCancel", ImVec2(120.0f, 0.0f)))
            ImGui::CloseCurrentPopup();

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

    ImGui::SetNextWindowSize(ImVec2(720.0f, 500.0f), ImGuiCond_Appearing);
    bool exportImageOpen = true;
    if (ImGui::BeginPopupModal("Export Image", &exportImageOpen, 0))
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
        if (ImGui::Button("\xe2\x86\x90##exportBack", ImVec2(navBtnW, 0.0f)) && exportHistoryIndex > 0)
        {
            exportHistoryIndex--;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Back");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x92##exportFwd", ImVec2(navBtnW, 0.0f)) && exportHistoryIndex + 1 < (int)exportDirHistory.size())
        {
            exportHistoryIndex++;
            exportDir = exportDirHistory[exportHistoryIndex];
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Forward");
        ImGui::SameLine(0.0f, 2.0f);
        if (ImGui::Button("\xe2\x86\x91##exportUp", ImVec2(navBtnW, 0.0f)))
            pushExportDir(parentPath(exportDir));
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Up one level");
        ImGui::SameLine(0.0f, 8.0f);
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

        // --- Sidebar + file list ---
        const float listH = 200.0f;

        if (ImGui::BeginChild("##exportsidebar", ImVec2(s_sidebarW, listH), true))
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
        ImGui::SameLine(0.0f, 0.0f);
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
        ImGui::SameLine(0.0f, 0.0f);

        if (ImGui::BeginChild("##exportfilebrowser", ImVec2(0.0f, listH), true))
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

        ImGui::SliderInt("Resolution scale", &exportResolutionScale, 1, 8, "%dx");
        if (exportResolutionScale > 1)
            ImGui::TextDisabled("Output will be %dx the current window size.", exportResolutionScale);

        if (exportStatusMsg[0] != '\0')
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", exportStatusMsg);

        ImGui::Separator();

        if (ImGui::Button("Export", ImVec2(120.0f, 0.0f)))
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
        ImGui::SameLine(0.0f, 8.0f);
        if (ImGui::Button("Cancel##exportCancel", ImVec2(120.0f, 0.0f)))
            ImGui::CloseCurrentPopup();
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

    // ---- Miller Directions dialog ----------------------------------------
    if (showMillerDirectionsDialog)
    {
        ImGui::OpenPopup("Miller Directions");
        showMillerDirectionsDialog = false;
    }

    ImGui::SetNextWindowSize(ImVec2(680.0f, 480.0f), ImGuiCond_FirstUseEver);
    bool millerDirOpen = true;
    if (ImGui::BeginPopupModal("Miller Directions", &millerDirOpen, ImGuiWindowFlags_NoResize))
    {
        if (!structure.hasUnitCell)
        {
            ImGui::TextDisabled("Current structure has no unit cell. Miller directions are unavailable.");
            if (ImGui::Button("Close", ImVec2(120, 0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
        else
        {
            ImGui::Checkbox("Show Miller directions", &showMillerDirections);
            ImGui::TextDisabled("Arrows are drawn from the cell origin along the [uvw] direction.");
            ImGui::Separator();

            ImGui::Text("New direction:");
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("U##dir", &millerDirInputU, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("V##dir", &millerDirInputV, 0.2f, -50, 50);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(85.0f);
            ImGui::DragInt("W##dir", &millerDirInputW, 0.2f, -50, 50);

            ImGui::SetNextItemWidth(210.0f);
            ImGui::DragFloat("Length (A)##dir", &millerDirInputLength, 0.1f, 0.5f, 100.0f, "%.1f");
            ImGui::SameLine();
            ImGui::ColorEdit3("Color##dir", millerDirInputColor, ImGuiColorEditFlags_NoInputs);

            const bool invalidDir =
                (millerDirInputU == 0 && millerDirInputV == 0 && millerDirInputW == 0);
            if (invalidDir)
                ImGui::TextDisabled("(u, v, w) cannot all be zero.");

            if (ImGui::Button("Add Direction", ImVec2(140.0f, 0.0f)) && !invalidDir)
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
            else if (ImGui::BeginChild("##miller-dir-list", ImVec2(0.0f, 200.0f), true))
            {
                int deleteIndex = -1;
                for (size_t i = 0; i < millerDirections.size(); ++i)
                {
                    ImGui::PushID((int)i + 40000);
                    MillerDirection& md = millerDirections[i];

                    ImGui::Checkbox("##vis", &md.visible);
                    ImGui::SameLine();
                    ImGui::Text("Dir %d", (int)i + 1);

                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("U##du", &md.u, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("V##dv", &md.v, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(80.0f);
                    ImGui::DragInt("W##dw", &md.w, 0.2f, -50, 50);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(100.0f);
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
                    if (ImGui::Button("Delete##ddel"))
                        deleteIndex = (int)i;

                    ImGui::Separator();
                    ImGui::PopID();
                }

                if (deleteIndex >= 0)
                    millerDirections.erase(millerDirections.begin() + deleteIndex);

                ImGui::EndChild();
            }

            ImGui::Spacing();
            if (ImGui::Button("Close##dclose", ImVec2(120, 0)))
                ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }
    if (!millerDirOpen)
        ImGui::CloseCurrentPopup();
    // ---- end Miller Directions dialog ------------------------------------

    if (showManual)
    {
        ImGui::OpenPopup("Manual");
        showManual = false;
    }

    auto wrappedBullet = [](const char* text) {
        ImGui::Bullet();
        ImGui::SameLine();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextWrapped("%s", text);
        ImGui::PopTextWrapPos();
    };

    ImGui::SetNextWindowSize(ImVec2(800.0f, 700.0f), ImGuiCond_Appearing);
    bool manualOpen = true;
    if (ImGui::BeginPopupModal("Manual", &manualOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("AtomForge Manual");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Complete reference for all features and operations");
        ImGui::PopTextWrapPos();
        ImGui::Separator();

        if (ImGui::BeginChild("##manual-scroll", ImVec2(0.0f, 640.0f), false))
        {
            ImGui::Text("Getting Started");
            wrappedBullet("Open a structure from File -> Open or press Ctrl+O.");
            wrappedBullet("Use left-drag to rotate the scene and the scroll wheel to zoom.");
            wrappedBullet("Use View -> Reset Default View to restore the fitted isometric camera.");
            wrappedBullet("The interface scales automatically on HiDPI / high-resolution displays so text and controls remain readable.");

            ImGui::Spacing();
            ImGui::Text("Selection and Editing");
            wrappedBullet("Left-click selects one atom.");
            wrappedBullet("Ctrl+left-click adds or removes atoms from the current selection.");
            wrappedBullet("Ctrl+A selects all atoms. Ctrl+D or Esc clears the selection.");
            wrappedBullet("Delete removes the selected atoms from the structure.");
            wrappedBullet("Right-click opens the context menu when atoms are selected.");
            wrappedBullet("Enable Edit -> Box Select Mode to select with a right-drag screen rectangle.");

            ImGui::Spacing();
            ImGui::Text("Keyboard Shortcuts");
            wrappedBullet("Ctrl+O: open structure file.");
            wrappedBullet("Ctrl+S: save structure as.");
            wrappedBullet("Ctrl+Shift+S: export the current rendered view.");
            wrappedBullet("Ctrl+Z: undo. Ctrl+Y or Ctrl+Shift+Z: redo.");

            ImGui::Spacing();
            ImGui::Text("File Menu");
            wrappedBullet("Open loads supported structure formats such as CIF, MOL, PDB, XYZ, SDF, VASP, MOL2, PWI, and GJF.");
            wrappedBullet("Save As exports the current structure to supported chemistry and crystal formats.");
            wrappedBullet("Export Image writes PNG, JPG, or SVG output with optional background and adjustable resolution scale for high-resolution output.");
            wrappedBullet("Close unloads the current structure.");

            ImGui::Spacing();
            ImGui::Text("Edit Menu");
            wrappedBullet("Undo and Redo track structure edits and style changes.");
            wrappedBullet("Edit Structure modifies lattice vectors and the atom list directly.");
            wrappedBullet("Atomic Sizes adjusts per-element radii used for display and some builders.");
            wrappedBullet("Element Colors adjusts per-element color and shininess.");
            wrappedBullet("Transform Structure applies a 3x3 transformation matrix to periodic structures.");
            wrappedBullet("Merge Structures opens a 3D arrangement dialog to load, place, rotate, and merge multiple structures.");

            ImGui::Spacing();
            ImGui::Text("Build Menu");
            wrappedBullet("Contains builders for Bulk Crystal, CSL Grain Boundary, Nanocrystal, Polycrystal, and Custom Structure.");
            wrappedBullet("See the individual builder sections below for details on each.");

            ImGui::Spacing();
            ImGui::Text("Bulk Crystal Builder");
            wrappedBullet("Builds a full periodic unit cell from crystal system, space group, lattice parameters, and asymmetric-unit atoms.");
            wrappedBullet("Useful for creating reference single crystals and simulation-ready periodic cells.");

            ImGui::Spacing();
            ImGui::Text("CSL Grain Boundary Builder");
            wrappedBullet("Builds cubic bicrystals from ideal sc, bcc, fcc, or diamond source lattices.");
            wrappedBullet("Controls Sigma selection, grain-boundary plane, replication, overlap removal, and rigid translation.");

            ImGui::Spacing();
            ImGui::Text("Nanocrystal Builder");
            wrappedBullet("Uses the currently loaded structure as the reference source for carving.");
            wrappedBullet("You can drag and drop a supported structure file into the reference preview area while the dialog is open.");
            wrappedBullet("Preview controls: left-drag orbits the preview, scroll zooms the preview camera.");
            wrappedBullet("Options include auto-centering, manual center coordinates, auto-replication for periodic inputs, and rectangular output-cell padding.");

            ImGui::Spacing();
            ImGui::Text("Polycrystal Builder");
            wrappedBullet("Requires a loaded reference single crystal with unit-cell information.");
            wrappedBullet("The left panel accepts drag-and-drop and displays a 3D preview of the reference crystal.");
            wrappedBullet("The builder creates grains by Voronoi tessellation inside a user-defined simulation box.");
            wrappedBullet("Grain orientations can be all-random, all user-specified with Euler angles, or partially specified.");
            wrappedBullet("Generated structures include per-atom IPF-Z orientation colors for crystal-orientation visualization.");

            ImGui::Spacing();
            ImGui::Text("Custom Structure Builder");
            wrappedBullet("Accepts drag-and-drop input for both the source crystal and the target 3D model.");
            wrappedBullet("Displays live 3D previews of the reference crystal and the imported model side by side.");
            wrappedBullet("The model preview is rendered as a shaded surface instead of a wireframe for easier inspection.");
            wrappedBullet("Use it to generate finite atomistic structures constrained by imported mesh geometry.");

            ImGui::Spacing();
            ImGui::Text("Merge Structures");
            wrappedBullet("Located in Edit -> Merge Structures.");
            wrappedBullet("Provides a large 3D preview with per-structure selection and translate/rotate gizmo controls.");
            wrappedBullet("Supports drag-and-drop loading and optional bounding-box display before committing merge output.");

            ImGui::Spacing();
            ImGui::Text("View Menu");
            wrappedBullet("Show Element toggles element labels.");
            wrappedBullet("Show Bonds toggles bond-cylinder rendering.");
            wrappedBullet("Isometric View and Orthographic View switch the camera projection mode.");
            wrappedBullet("Select Theme switches between the default dark theme and the light theme.");
            wrappedBullet("Color Structure By switches between element-type colors and crystal-orientation IPF coloring.");
            wrappedBullet("Structure Info shows composition, lattice metrics, positions, and symmetry when available.");
            wrappedBullet("Measure Distance, Measure Angle, and Atom Info open the corresponding dialogs for the current selection.");

            ImGui::Spacing();
            ImGui::Text("Crystal Orientation / IPF");
            wrappedBullet("Crystal Orientation coloring uses cubic IPF-Z colors and displays an IPF triangle legend in the main view.");
            wrappedBullet("When available, IPF data is restored from a saved sidecar file named basename.atomforge-ipf.");
            wrappedBullet("If no saved IPF metadata is present, AtomForge can fall back to geometry-based orientation reconstruction.");

            ImGui::Spacing();
            ImGui::Text("Analysis Menu");
            wrappedBullet("Common Neighbour Analysis reports pair signatures and per-atom structural environments.");
            wrappedBullet("Radial Distribution Function plots RDF with configurable species filters, normalization, radius range, bin count, and smoothing.");

            ImGui::Spacing();
            ImGui::Text("Context Menu");
            wrappedBullet("Substitute Atom replaces the selected atoms with a chosen element.");
            wrappedBullet("Insert Atom at Midpoint places a new atom at the centroid of the current selection.");
            wrappedBullet("Measure Distance, Measure Angle, Atom Info, Delete, and Deselect are available when selection rules are satisfied.");

            ImGui::Spacing();
            ImGui::Text("Display and Measurement Features");
            wrappedBullet("Element labels can be shown for periodic-image atoms as well.");
            wrappedBullet("Bonds are inferred from covalent radii and rendered with split element colors.");
            wrappedBullet("Selected atoms are highlighted and helper overlays are drawn for distance and angle tools.");
            wrappedBullet("Periodic boundary visualization includes duplicated boundary atoms and transformed supercell views when applicable.");
            wrappedBullet("Overlay, gizmo, and bounding-box colors adapt automatically to the selected dark or light theme.");

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

    ImGui::SetNextWindowSize(ImVec2(860.0f, 620.0f), ImGuiCond_Appearing);
    bool aboutOpen = true;
    if (ImGui::BeginPopupModal("About", &aboutOpen, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("AtomForge");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Version " ATOMFORGE_VERSION "  –  Atomistic structure builder, viewer, and editor with periodic-cell and finite-geometry workflows");
        ImGui::PopTextWrapPos();
        ImGui::Separator();

        if (ImGui::BeginChild("##about-scroll", ImVec2(0.0f, -ImGui::GetFrameHeightWithSpacing() - 6.0f), false))
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

            ImGui::EndChild();
        }

        if (ImGui::Button("Close##about", ImVec2(120.0f, 0.0f)))
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

void FileBrowser::initSubstitutionalSolidSolutionRenderResources(Renderer& renderer)
{
    substitutionalSolidSolutionDialog.initRenderResources(renderer);
}

bool FileBrowser::isSubstitutionalSolidSolutionDialogOpen() const
{
    return substitutionalSolidSolutionDialog.isOpen();
}

void FileBrowser::feedDropToSubstitutionalSolidSolutionDialog(const std::string& path)
{
    substitutionalSolidSolutionDialog.feedDroppedFile(path);
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
