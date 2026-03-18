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
#include <dirent.h>
#include <sys/stat.h>

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
            openDir("."),
      historyIndex(-1),
      saveDir("."),
      saveHistoryIndex(-1),
      selectedSaveFormat(0),
      selectedAtomicNumber(1)
{
    allowedExtensions = {".cif", ".mol", ".pdb", ".xyz", ".sdf"};

    driveRoots.push_back("/");
    if (const char* home = std::getenv("HOME"))
        driveRoots.push_back(home);
    driveRoots.push_back("/mnt");
    driveRoots.push_back("/media");

    openFilename[0] = '\0';
    saveFilename[0] = '\0';
    saveStatusMsg[0] = '\0';
}

void FileBrowser::initFromPath(const std::string& initialPath)
{
    openDir = ".";
    auto pos = initialPath.find_last_of("/\\");
    if (pos != std::string::npos)
        openDir = initialPath.substr(0, pos);

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

        if (ImGui::BeginMenu("View"))
        {
            ImGui::MenuItem("Show Element", nullptr, &showElementLabels);
            ImGui::MenuItem("Show Bonds", nullptr, &showBonds);
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

    if (openStructurePopup)
    {
        ImGui::OpenPopup("Open Structure");
        openStructurePopup = false;
    }

    bool openStructureOpen = true;
    if (ImGui::BeginPopupModal("Open Structure", &openStructureOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Current folder: %s", openDir.c_str());
        ImGui::SameLine();
        if (ImGui::Button(".."))
        {
            auto pos = openDir.find_last_of("/\\");
            if (pos != std::string::npos)
                openDir = openDir.substr(0, pos);
            else
                openDir = ".";
            pushHistory(openDir);
        }

        ImGui::SameLine();
        if (ImGui::Button("Back") && historyIndex > 0)
        {
            historyIndex--;
            openDir = dirHistory[historyIndex];
        }
        ImGui::SameLine();
        if (ImGui::Button("Forward") && historyIndex + 1 < (int)dirHistory.size())
        {
            historyIndex++;
            openDir = dirHistory[historyIndex];
        }

        ImGui::SameLine();
        if (ImGui::Button("Root"))
        {
            openDir = "/";
            pushHistory(openDir);
        }
        ImGui::SameLine();
        if (ImGui::Button("Home") && !driveRoots.empty())
        {
            openDir = driveRoots[1];
            pushHistory(openDir);
        }

        ImGui::Separator();

        if (ImGui::BeginChild("##filebrowser", ImVec2(500, 300), true))
        {
            DIR* dir = opendir(openDir.c_str());
            if (dir)
            {
                struct dirent* de;
                std::vector<std::pair<std::string, bool>> entries;

                while ((de = readdir(dir)) != nullptr)
                {
                    std::string name(de->d_name);
                    if (name == "." || name == "..")
                        continue;

                    std::string fullPath = openDir + "/" + name;
                    struct stat st;
                    bool isDir = (stat(fullPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode));

                    if (!isDir && !isAllowedFile(name))
                        continue;

                    entries.emplace_back(name, isDir);
                }
                closedir(dir);

                // directories first
                std::sort(entries.begin(), entries.end(),
                          [](const std::pair<std::string, bool>& a,
                             const std::pair<std::string, bool>& b) {
                    if (a.second != b.second)
                        return a.second > b.second;
                    return a.first < b.first;
                });

                for (size_t i = 0; i < entries.size(); ++i)
                {
                    const std::string& name = entries[i].first;
                    bool isDir = entries[i].second;

                    ImGui::PushID((int)i);
                    if (isDir)
                    {
                        std::string label = std::string("[DIR] ") + name + "##" + name;
                        if (ImGui::Selectable(label.c_str()))
                        {
                            if (openDir == ".")
                                openDir = name;
                            else
                                openDir = openDir + "/" + name;
                            pushHistory(openDir);
                        }
                    }
                    else
                    {
                        std::string label = name + "##" + name;
                        bool selected = (std::string(openFilename) == name);
                        if (ImGui::Selectable(label.c_str(), selected))
                        {
                            std::snprintf(openFilename, sizeof(openFilename), "%s", name.c_str());
                        }
                    }
                    ImGui::PopID();
                }
            }
            else
            {
                ImGui::TextDisabled("Unable to open folder");
            }

            ImGui::EndChild();
        }

        ImGui::InputText("Filename", openFilename, sizeof(openFilename));

        if (ImGui::Button("Load"))
        {
            std::string fullPath = openDir + "/" + openFilename;
            Structure newStructure = loadStructure(fullPath);
            if (!newStructure.atoms.empty())
            {
                structure = std::move(newStructure);

                // Apply any user-specified element color overrides
                for (auto& atom : structure.atoms)
                {
                    int atomicNumber = atom.atomicNumber;
                    auto it = elementColorOverrides.find(atomicNumber);
                    if (it != elementColorOverrides.end())
                    {
                        atom.r = it->second[0];
                        atom.g = it->second[1];
                        atom.b = it->second[2];
                    }
                }

                updateBuffers(structure);
                requestResetDefaultView = true;
            }
            ImGui::CloseCurrentPopup();
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
            auto pos = saveDir.find_last_of("/\\");
            pushSaveDir(pos != std::string::npos ? saveDir.substr(0, pos) : ".");
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
            pushSaveDir("/");
        ImGui::SameLine();
        if (ImGui::Button("Home##save") && !driveRoots.empty())
            pushSaveDir(driveRoots[1]);

        ImGui::Separator();

        if (ImGui::BeginChild("##savefilebrowser", ImVec2(500, 200), true))
        {
            DIR* dir = opendir(saveDir.c_str());
            if (dir)
            {
                struct dirent* de;
                std::vector<std::pair<std::string, bool>> entries;
                while ((de = readdir(dir)) != nullptr)
                {
                    std::string name(de->d_name);
                    if (name == "." || name == "..")
                        continue;
                    std::string fullPath = saveDir + "/" + name;
                    struct stat st;
                    bool isDir = (stat(fullPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
                    entries.emplace_back(name, isDir);
                }
                closedir(dir);

                std::sort(entries.begin(), entries.end(),
                          [](const std::pair<std::string, bool>& a,
                             const std::pair<std::string, bool>& b) {
                    if (a.second != b.second) return a.second > b.second;
                    return a.first < b.first;
                });

                for (size_t i = 0; i < entries.size(); ++i)
                {
                    const std::string& name = entries[i].first;
                    bool isDir = entries[i].second;
                    ImGui::PushID((int)(i + 10000));
                    if (isDir)
                    {
                        std::string lbl = std::string("[DIR] ") + name + "##" + name;
                        if (ImGui::Selectable(lbl.c_str()))
                        {
                            pushSaveDir(saveDir == "." ? name : saveDir + "/" + name);
                        }
                    }
                    else
                    {
                        std::string lbl = name + "##" + name;
                        bool sel = (std::string(saveFilename) == name);
                        if (ImGui::Selectable(lbl.c_str(), sel))
                            std::snprintf(saveFilename, sizeof(saveFilename), "%s", name.c_str());
                    }
                    ImGui::PopID();
                }
            }
            else
            {
                ImGui::TextDisabled("Unable to open folder");
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
                std::string fullPath = saveDir + "/" + saveFilename;
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
        ImGui::Text("Atoms Editor Manual");
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
            ImGui::Text("View Menu");
            ImGui::BulletText("Show Element: toggle element labels");
            ImGui::BulletText("Show Bonds: toggle bond rendering");
            ImGui::BulletText("Structure Info: composition, lattice metrics, positions, and symmetry");
            ImGui::BulletText("Measure Distance (2 selected): open distance dialog and overlay");
            ImGui::BulletText("Measure Angle (3 selected): open angle dialog and overlay");
            ImGui::BulletText("Atom Info (1 selected): show atom metadata and coordinates");
            ImGui::BulletText("Reset Default View: restore fitted isometric camera view");

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
        ImGui::Text("Atoms Editor");
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
