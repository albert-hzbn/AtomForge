#include "FileBrowser.h"
#include "ElementData.h"
#include "io/StructureLoader.h"
#include "graphics/StructureInstanceBuilder.h"
#include "ui/PeriodicTableDialog.h"

#include "imgui.h"
#include <GLFW/glfw3.h>
#include <openbabel3/openbabel/elements.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

namespace
{
using DirectoryEntry = std::pair<std::string, bool>;

std::string normalizePathSeparators(const std::string& path)
{
    std::string out = path;
    std::replace(out.begin(), out.end(), '\\', '/');
    return out;
}

std::string joinPath(const std::string& base, const std::string& name)
{
    if (base.empty() || base == ".")
        return name;
    if (base.back() == '/' || base.back() == '\\')
        return base + name;
    return base + "/" + name;
}

bool isDriveRootPath(const std::string& path)
{
    if (path.size() < 2)
        return false;
    if (!std::isalpha((unsigned char)path[0]) || path[1] != ':')
        return false;
    return (path.size() == 2) || (path.size() == 3 && (path[2] == '/' || path[2] == '\\'));
}

std::string parentPath(const std::string& path)
{
    if (path.empty() || path == ".")
        return ".";

    std::string out = path;
    while (out.size() > 1 && (out.back() == '/' || out.back() == '\\'))
    {
        if (isDriveRootPath(out))
            return normalizePathSeparators(out);
        out.pop_back();
    }

    if (out == "/" || out == "\\" || isDriveRootPath(out))
        return normalizePathSeparators(out);

    std::size_t pos = out.find_last_of("/\\");
    if (pos == std::string::npos)
        return ".";
    if (pos == 0)
        return out.substr(0, 1);
    if (pos == 2 && std::isalpha((unsigned char)out[0]) && out[1] == ':')
        return normalizePathSeparators(out.substr(0, 3));
    return out.substr(0, pos);
}

std::string detectHomePath()
{
#ifdef _WIN32
    if (const char* userProfile = std::getenv("USERPROFILE"))
        return normalizePathSeparators(userProfile);

    const char* homeDrive = std::getenv("HOMEDRIVE");
    const char* homePath = std::getenv("HOMEPATH");
    if (homeDrive && homePath)
        return normalizePathSeparators(std::string(homeDrive) + homePath);
#endif

    if (const char* home = std::getenv("HOME"))
        return normalizePathSeparators(home);

    return ".";
}

void appendUniquePath(std::vector<std::string>& paths, const std::string& value)
{
    if (value.empty())
        return;

    std::string normalized = normalizePathSeparators(value);
    if (isDriveRootPath(normalized) && normalized.size() == 2)
        normalized += "/";

    if (std::find(paths.begin(), paths.end(), normalized) == paths.end())
        paths.push_back(normalized);
}

#ifdef _WIN32
std::string toNativePath(const std::string& path)
{
    std::string native = path;
    std::replace(native.begin(), native.end(), '/', '\\');
    return native;
}
#endif

bool loadDirectoryEntries(const std::string& directory,
                          bool filterFiles,
                          const std::function<bool(const std::string&)>& includeFile,
                          std::vector<DirectoryEntry>& entries)
{
#ifdef _WIN32
    std::string nativeDir = toNativePath(directory.empty() ? "." : directory);
    if (!nativeDir.empty() && nativeDir.back() != '\\' && nativeDir.back() != '/')
        nativeDir.push_back('\\');
    std::string searchPattern = nativeDir + "*";

    WIN32_FIND_DATAA findData;
    HANDLE handle = FindFirstFileA(searchPattern.c_str(), &findData);
    if (handle == INVALID_HANDLE_VALUE)
        return false;

    do
    {
        std::string name(findData.cFileName);
        if (name == "." || name == "..")
            continue;

        bool isDir = (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

        if (!isDir && filterFiles && !includeFile(name))
            continue;

        entries.emplace_back(name, isDir);
    } while (FindNextFileA(handle, &findData) != 0);

    FindClose(handle);
#else
    DIR* dir = opendir(directory.c_str());
    if (!dir)
        return false;

    struct dirent* de;
    while ((de = readdir(dir)) != nullptr)
    {
        std::string name(de->d_name);
        if (name == "." || name == "..")
            continue;

        std::string fullPath = joinPath(directory, name);
        struct stat st;
        bool isDir = (stat(fullPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode));

        if (!isDir && filterFiles && !includeFile(name))
            continue;

        entries.emplace_back(name, isDir);
    }
    closedir(dir);
#endif

    std::sort(entries.begin(), entries.end(),
              [](const DirectoryEntry& a, const DirectoryEntry& b) {
        if (a.second != b.second)
            return a.second > b.second;
        return a.first < b.first;
    });

    return true;
}

void drawDirectoryEntries(const std::vector<DirectoryEntry>& entries,
                          char* selectedFilename,
                          int idBase,
                          const std::function<void(const std::string&)>& onEnterDirectory)
{
    for (size_t i = 0; i < entries.size(); ++i)
    {
        const std::string& name = entries[i].first;
        bool isDir = entries[i].second;

        ImGui::PushID((int)i + idBase);
        if (isDir)
        {
            std::string label = std::string("[DIR] ") + name + "##" + name;
            if (ImGui::Selectable(label.c_str()))
                onEnterDirectory(name);
        }
        else
        {
            std::string label = name + "##" + name;
            bool selected = (std::string(selectedFilename) == name);
            if (ImGui::Selectable(label.c_str(), selected))
                std::snprintf(selectedFilename, 1024, "%s", name.c_str());
        }
        ImGui::PopID();
    }
}
}

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

FileBrowser::FileBrowser()
        : showAbout(false),
            showManual(false),
            showEditColors(false),
            showElementLabels(false),
            showBonds(false),
            viewMode(ViewMode::Isometric),
            boxSelectMode(false),
            requestMeasureDistance(false),
            requestMeasureAngle(false),
            requestAtomInfo(false),
            requestResetDefaultView(false),
            requestStructureInfo(false),
            requestUndo(false),
            requestRedo(false),
            openStructurePopup(false),
            saveStructurePopup(false),
            loadErrorPopupRequested(false),
            openDir("."),
      historyIndex(-1),
      saveDir("."),
      saveHistoryIndex(-1),
      selectedSaveFormat(0),
      selectedAtomicNumber(1)
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
    loadErrorMsg[0] = '\0';
    saveFilename[0] = '\0';
    saveStatusMsg[0] = '\0';
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
            if (ImGui::MenuItem("Open...",  "Ctrl+O"))
                openStructurePopup = true;

            if (ImGui::MenuItem("Save As...", "Ctrl+S"))
            {
                saveStructurePopup = true;
                saveDir = openDir;
                saveDirHistory = dirHistory;
                saveHistoryIndex = historyIndex;
                saveStatusMsg[0] = '\0';
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
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("View"))
        {
            ImGui::MenuItem("Show Element", nullptr, &showElementLabels);
            ImGui::MenuItem("Show Bonds", nullptr, &showBonds);
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
            if (ImGui::MenuItem("Measure Distance (2 selected)"))
                requestMeasureDistance = true;
            if (ImGui::MenuItem("Measure Angle (3 selected)"))
                requestMeasureAngle = true;
            if (ImGui::MenuItem("Atom Info (1 selected)"))
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

    editMenuDialogs.drawPopups(structure, updateBuffers);

    transformDialog.drawDialog([&]() { updateBuffers(structure); });
    bulkCrystalDialog.drawDialog(structure, editMenuDialogs.elementColors, updateBuffers);
    cslDialog.drawDialog(structure, updateBuffers);
    cnaDialog.drawDialog(structure);
    rdfDialog.drawDialog(structure);

    if (loadErrorPopupRequested)
    {
        ImGui::OpenPopup("Load Error");
        loadErrorPopupRequested = false;
    }

    ImGui::SetNextWindowSize(ImVec2(720.0f, 0.0f), ImGuiCond_Appearing);
    bool loadErrorOpen = true;
    if (ImGui::BeginPopupModal("Load Error", &loadErrorOpen, ImGuiWindowFlags_NoResize))
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
            if (loadStructureFromFile(fullPath, newStructure, loadError))
            {
                structure = std::move(newStructure);
                applyElementColorOverrides(structure);

                updateBuffers(structure);
                requestResetDefaultView = true;
                openStatusMsg[0] = '\0';
                ImGui::CloseCurrentPopup();
            }
            else
            {
                std::snprintf(
                    openStatusMsg,
                    sizeof(openStatusMsg),
                    "%s",
                    loadError.empty() ? "Failed to load file." : loadError.c_str());
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
            if (saveHistoryIndex + 1 < (int)saveDirHistory.size())
                saveDirHistory.erase(saveDirHistory.begin() + saveHistoryIndex + 1, saveDirHistory.end());
            saveDirHistory.push_back(saveDir);
            saveHistoryIndex = (int)saveDirHistory.size() - 1;
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
                    if (saveExtFilter.empty())
                        return true;
                    if (lowerName.size() < saveExtFilter.size())
                        return false;
                    return lowerName.compare(lowerName.size() - saveExtFilter.size(), saveExtFilter.size(), saveExtFilter) == 0;
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
            // Auto-update file extension when format changes
            std::string fn(saveFilename);
            auto dot = fn.find_last_of('.');
            std::string base = (dot != std::string::npos) ? fn.substr(0, dot) : fn;
            if (base.empty()) base = "structure";
            std::snprintf(saveFilename, sizeof(saveFilename), "%s%s",
                          base.c_str(), kSaveFormats[selectedSaveFormat].ext);
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
                    // When a supercell transform is active, expand to the full
                    // supercell so the saved file contains all visible atoms.
                    Structure structureToSave = (isTransformMatrixEnabled() && structure.hasUnitCell)
                        ? buildSupercell(structure, getTransformMatrix())
                        : structure;
                    bool ok = saveStructure(structureToSave, fullPath, kSaveFormats[selectedSaveFormat].fmt);
                if (ok)
                {
                    ImGui::CloseCurrentPopup();
                }
                else
                {
                    std::snprintf(saveStatusMsg, sizeof(saveStatusMsg),
                                  "Error: failed to save (format may not support this structure).");
                }
            }
        }
        ImGui::EndPopup();
    }
    if (!saveAsOpen)
        ImGui::CloseCurrentPopup();
    // ---- end Save As dialog --------------------------------------------

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
            ImGui::Text("Mouse Controls");
            ImGui::BulletText("Rotate: left-click drag");
            ImGui::BulletText("Zoom: scroll wheel");
            ImGui::BulletText("Select atom: left-click on atom");
            ImGui::BulletText("Multi-select: Ctrl+left-click to add/remove atoms");
            ImGui::BulletText("Right-click: open context menu (when atoms are selected)");
            ImGui::BulletText("Box Select Mode: right-drag rectangle to select atoms");

            ImGui::Spacing();
            ImGui::Text("Keyboard Shortcuts");
            ImGui::BulletText("Ctrl+A: select all atoms");
            ImGui::BulletText("Ctrl+D or Esc: clear selection");
            ImGui::BulletText("Delete: remove selected atoms from structure");
            ImGui::BulletText("Ctrl+Z: undo last change");
            ImGui::BulletText("Ctrl+Y or Ctrl+Shift+Z: redo");
            ImGui::BulletText("Ctrl+O: open structure file");
            ImGui::BulletText("Ctrl+S: save structure as");

            ImGui::Spacing();
            ImGui::Text("File Menu");
            ImGui::BulletText("Open...: load from .cif, .mol, .pdb, .xyz, .sdf");
            ImGui::BulletText("Save As...: export to .xyz, .cif, .vasp, .pdb, .sdf, .mol2, .pwi, .gjf");
            ImGui::BulletText("Quit: exit application");

            ImGui::Spacing();
            ImGui::Text("Edit Menu");
            ImGui::BulletText("Undo / Redo: history navigation for structure and style edits");
            ImGui::BulletText("Box Select Mode: enable right-drag rectangular atom selection");
            ImGui::BulletText("Edit Structure...: modify lattice vectors and atom list");
            ImGui::BulletText("Atomic Sizes...: adjust per-element covalent radii");
            ImGui::BulletText("Element Colors...: adjust per-element colors and material shininess");
            ImGui::BulletText("Transform Structure...: apply a 3x3 integer transformation matrix");

            ImGui::Spacing();
            ImGui::Text("Build Menu");
            ImGui::BulletText("Bulk Crystal...: build a full unit cell from crystal system, space group, lattice parameters, and asymmetric-unit atoms");
            ImGui::BulletText("CSL Grain Boundary...: cubic GB builder with Sigma-list selection and in-plane supercell replication");

            ImGui::Spacing();
            ImGui::Text("View Menu");
            ImGui::BulletText("Show Element: toggle element labels");
            ImGui::BulletText("Show Bonds: toggle bond rendering");
            ImGui::BulletText("Isometric View: fit and reset to the default angled 3D perspective view");
            ImGui::BulletText("Orthographic View: switch to orthographic projection while keeping the default fitted viewing angle");
            ImGui::BulletText("Structure Info: composition, lattice metrics, positions, and symmetry");
            ImGui::BulletText("Measure Distance (2 selected): open distance dialog and overlay");
            ImGui::BulletText("Measure Angle (3 selected): open angle dialog and overlay");
            ImGui::BulletText("Atom Info (1 selected): metadata, coordinates, and PBC-aware coordination/bond-length stats");
            ImGui::BulletText("Reset Default View: restore fitted isometric camera view");

            ImGui::Spacing();
            ImGui::Text("Analysis Menu");
            ImGui::BulletText("Common Neighbour Analysis...: run CNA and inspect pair signatures and per-atom environments");
            ImGui::BulletText("Radial Distribution Function...: plot RDF with configurable species filters, radius range, bins, smoothing, and normalization");

            ImGui::Spacing();
            ImGui::Text("Context Menu (Selection)");
            ImGui::BulletText("Substitute Atom...: replace selected atoms with a chosen element");
            ImGui::BulletText("Insert Atom at Midpoint...: add atom at centroid of selected atoms");
            ImGui::BulletText("Measure Distance / Angle / Atom Info shortcuts for selected atoms");
            ImGui::BulletText("Delete / Deselect selected atoms");

            ImGui::Spacing();
            ImGui::Text("Display Features");
            ImGui::BulletText("Element labels with periodic-image notation");
            ImGui::BulletText("Split-color bonds based on bonded element colors");
            ImGui::BulletText("Yellow highlighting for selected atoms");
            ImGui::BulletText("Distance/angle helper overlays in the 3D scene");
            ImGui::BulletText("Periodic boundary visualization and supercell transforms");

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

void FileBrowser::showLoadError(const std::string& message)
{
    std::snprintf(
        loadErrorMsg,
        sizeof(loadErrorMsg),
        "%s",
        message.empty() ? "Failed to load file." : message.c_str());
    loadErrorPopupRequested = true;
}

void FileBrowser::pushHistory(const std::string& dir)
{
    if (historyIndex + 1 < (int)dirHistory.size())
        dirHistory.erase(dirHistory.begin() + historyIndex + 1, dirHistory.end());

    dirHistory.push_back(dir);
    historyIndex = (int)dirHistory.size() - 1;
}

std::string FileBrowser::toLower(const std::string& s)
{
    std::string out = s;
    for (auto& c : out)
        c = (char)std::tolower((unsigned char)c);
    return out;
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
