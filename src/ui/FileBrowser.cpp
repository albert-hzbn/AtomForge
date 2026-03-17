#include <array>

// Static array of element symbols for atomic numbers 1-118
static const std::array<const char*, 119> elementSymbols = {
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
#include "FileBrowser.h"
#include "io/StructureLoader.h"
#include "graphics/StructureInstanceBuilder.h"

#include "ui/PeriodicTableDialog.h"

#include "imgui.h"
#include <GLFW/glfw3.h>
#include <openbabel3/openbabel/elements.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>

FileBrowser::FileBrowser()
    : showAbout(false),
      showEditColors(false),
    showElementLabels(false),
    requestMeasureDistance(false),      requestMeasureAngle(false),      openStructurePopup(false),
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
                       const std::function<void(Structure&)>& updateBuffers)
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
            editMenuDialogs.drawMenuItems();
            ImGui::Separator();
            transformDialog.drawMenuItem(structure.hasUnitCell);

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("View"))
        {
            ImGui::MenuItem("Show Element", nullptr, &showElementLabels);
            if (ImGui::MenuItem("Measure Distance (2 selected)"))
                requestMeasureDistance = true;
            if (ImGui::MenuItem("Measure Angle (3 selected)"))
                requestMeasureAngle = true;
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("About"))
        {
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

    if (ImGui::BeginPopupModal("Open Structure", NULL, ImGuiWindowFlags_AlwaysAutoResize))
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
            }
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel"))
        {
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    // ---- Save As dialog ------------------------------------------------
    if (saveStructurePopup)
    {
        ImGui::OpenPopup("Save As");
        saveStructurePopup = false;
    }

    // Format table (label shown in combo, file extension, OpenBabel id)
    static const struct { const char* label; const char* ext; const char* fmt; } kFmts[] = {
        { "XYZ (.xyz)",                ".xyz",  "xyz"   },
        { "CIF (.cif)",                ".cif",  "cif"   },
        { "VASP POSCAR (.vasp)",       ".vasp", "vasp"  },
        { "PDB (.pdb)",                ".pdb",  "pdb"   },
        { "SDF (.sdf)",                ".sdf",  "sdf"   },
        { "Mol2 (.mol2)",              ".mol2", "mol2"  },
        { "Quantum ESPRESSO (.pwi)",   ".pwi",  "pwscf" },
        { "Gaussian Input (.gjf)",     ".gjf",  "gjf"   },
    };
    static const int kNumFmts = (int)(sizeof(kFmts) / sizeof(kFmts[0]));

    ImGui::SetNextWindowSize(ImVec2(560, 480), ImGuiCond_FirstUseEver);
    if (ImGui::BeginPopupModal("Save As", NULL, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Current folder: %s", saveDir.c_str());
        ImGui::SameLine();
        if (ImGui::Button("..##save"))
        {
            auto pos = saveDir.find_last_of("/\\");
            if (pos != std::string::npos)
                saveDir = saveDir.substr(0, pos);
            else
                saveDir = ".";
            if (saveHistoryIndex + 1 < (int)saveDirHistory.size())
                saveDirHistory.erase(saveDirHistory.begin() + saveHistoryIndex + 1, saveDirHistory.end());
            saveDirHistory.push_back(saveDir);
            saveHistoryIndex = (int)saveDirHistory.size() - 1;
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
        {
            saveDir = "/";
            if (saveHistoryIndex + 1 < (int)saveDirHistory.size())
                saveDirHistory.erase(saveDirHistory.begin() + saveHistoryIndex + 1, saveDirHistory.end());
            saveDirHistory.push_back(saveDir);
            saveHistoryIndex = (int)saveDirHistory.size() - 1;
        }
        ImGui::SameLine();
        if (ImGui::Button("Home##save") && !driveRoots.empty())
        {
            saveDir = driveRoots[1];
            if (saveHistoryIndex + 1 < (int)saveDirHistory.size())
                saveDirHistory.erase(saveDirHistory.begin() + saveHistoryIndex + 1, saveDirHistory.end());
            saveDirHistory.push_back(saveDir);
            saveHistoryIndex = (int)saveDirHistory.size() - 1;
        }

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
                            saveDir = (saveDir == "." ? name : saveDir + "/" + name);
                            if (saveHistoryIndex + 1 < (int)saveDirHistory.size())
                                saveDirHistory.erase(saveDirHistory.begin() + saveHistoryIndex + 1, saveDirHistory.end());
                            saveDirHistory.push_back(saveDir);
                            saveHistoryIndex = (int)saveDirHistory.size() - 1;
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
        static const char* kFmtLabels[] = {
            "XYZ (.xyz)",
            "CIF (.cif)",
            "VASP POSCAR (.vasp)",
            "PDB (.pdb)",
            "SDF (.sdf)",
            "Mol2 (.mol2)",
            "Quantum ESPRESSO (.pwi)",
            "Gaussian Input (.gjf)",
        };
        if (ImGui::Combo("Format", &selectedSaveFormat, kFmtLabels, kNumFmts))
        {
            // Auto-update file extension when format changes
            std::string fn(saveFilename);
            auto dot = fn.find_last_of('.');
            std::string base = (dot != std::string::npos) ? fn.substr(0, dot) : fn;
            if (base.empty()) base = "structure";
            std::string newName = base + kFmts[selectedSaveFormat].ext;
            std::snprintf(saveFilename, sizeof(saveFilename), "%s", newName.c_str());
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
                    bool ok = saveStructure(structureToSave, fullPath, kFmts[selectedSaveFormat].fmt);
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
        ImGui::SameLine();
        if (ImGui::Button("Cancel##save"))
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
    // ---- end Save As dialog --------------------------------------------

    if (showAbout)
    {
        ImGui::OpenPopup("About");
        showAbout = false;
    }

    if (ImGui::BeginPopupModal("About", NULL, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Atoms Editor");
        ImGui::Text("Open-source molecular structure viewer and editor.");
        ImGui::Separator();

        ImGui::Text("Camera");
        ImGui::BulletText("Left drag: rotate");
        ImGui::BulletText("Right drag: pan");
        ImGui::BulletText("Scroll: zoom");

        ImGui::Spacing();
        ImGui::Text("Selection");
        ImGui::BulletText("Left click: select atom");
        ImGui::BulletText("Ctrl + click: add/remove from selection");
        ImGui::BulletText("Ctrl+A: select all");
        ImGui::BulletText("Ctrl+D / Escape: deselect all");
        ImGui::BulletText("Delete: remove selected atoms");

        ImGui::Spacing();
        ImGui::Text("Right-click context menu");
        ImGui::BulletText("Substitute Atom: replace selected atoms with a new element");
        ImGui::BulletText("Insert Atom at Midpoint: place a new atom at the centroid");
        ImGui::BulletText("  (requires >= 2 atoms selected)");
        ImGui::BulletText("Delete / Deselect");

        ImGui::Spacing();
        ImGui::Text("Edit menu");
        ImGui::BulletText("Atomic Sizes: adjust per-element covalent radii");
        ImGui::BulletText("  (defaults: Cordero et al., Dalton Trans. 2008)");
        ImGui::BulletText("Element Colors: override CPK colours per element");
        ImGui::BulletText("Transform Atoms: apply a 3x3 matrix to all atom positions");

        ImGui::Spacing();
        ImGui::Text("File → Open  (Ctrl+O)");
        ImGui::BulletText("Supported: .cif  .mol  .pdb  .xyz  .sdf");

        ImGui::Spacing();
        ImGui::Text("File → Save As  (Ctrl+S)");
        ImGui::BulletText("XYZ (.xyz)");
        ImGui::BulletText("CIF (.cif)");
        ImGui::BulletText("VASP POSCAR (.vasp)");
        ImGui::BulletText("PDB (.pdb)");
        ImGui::BulletText("SDF (.sdf)");
        ImGui::BulletText("Mol2 (.mol2)");
        ImGui::BulletText("Quantum ESPRESSO (.pwi)");
        ImGui::BulletText("Gaussian Input (.gjf)");

        ImGui::Separator();
        if (ImGui::Button("OK"))
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }

    if (showEditColors)
    {
        ImGui::OpenPopup("Edit Element Colors");
        showEditColors = false;

        if (selectedAtomicNumber < 1 || selectedAtomicNumber > 118)
            selectedAtomicNumber = 1;
    }

    ImGui::SetNextWindowSize(ImVec2(950, 600), ImGuiCond_FirstUseEver);
    if (ImGui::BeginPopupModal("Edit Element Colors", NULL, ImGuiWindowFlags_NoResize))
    {
        ImGui::Text("Select an element to edit its color.");
        ImGui::Separator();

        drawPeriodicTableInlineSelector(selectedAtomicNumber);

        if (selectedAtomicNumber >= 1 && selectedAtomicNumber <= 118)
        {
            const char* selectedElementSymbol = (selectedAtomicNumber >= 1 && selectedAtomicNumber <= 118) ? elementSymbols[selectedAtomicNumber] : "?";

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

        ImGui::Separator();
        if (ImGui::Button("Close"))
            ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
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
