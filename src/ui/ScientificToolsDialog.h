#pragma once

#include "model/Structure.h"
#include "science/AtomProperties.h"
#include "science/Json.h"
#include "science/ResultPlots.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include <array>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

// Scientific analyses and simulations computed natively by AtomForge (src/science).
class ScientificToolsDialog
{
public:
    void drawMenuItems(const char* category);
    // Opens the dialog on a catalog tool id (for example "msd"); false if unknown.
    bool open(const std::string& toolId);
    // Trajectory currently loaded in Trajectory playback (empty when none).
    void setTrajectorySource(const std::string& path) { m_trajectory = path; }
    // Cell parameter filled together with an active-structure position input.
    static const char* companionCell(const std::string& parameter);
    // Whether a tool parameter can take the loaded trajectory directly.
    static bool acceptsTrajectory(const std::string& tool, const std::string& parameter);
    // Each family of tools has its own window: per-atom analyses beside the view,
    // plot workspaces, simulation runs and input-file generators.
    enum class Layout { Atoms, Trajectory, Properties, Plot, Simulation, Generator };
    static Layout layoutFor(const std::string& tool);
    // ImGui window id suffix ("###...") of a layout.
    static const char* windowId(Layout layout);
    // colourAtoms receives a per-atom property to show in the viewport.
    using ColourAtoms = std::function<void(const std::string&, const std::vector<double>&)>;
    void draw(const Structure& structure, const std::function<void(Structure&)>& loadResult,
              const ColourAtoms& colourAtoms = {});
    // Selected tool, inputs and last result as JSON text, for project files.
    std::string snapshot() const;
    // Restores a snapshot(); empty text keeps the current state.
    void restore(const std::string& state);
private:
    struct Field {
        std::array<char, 8192> value{};
        bool enabled = true;
        std::string path;
        bool useFile = false;
        bool selectColumn = false;
        int column = 0;
        std::array<char, 256> field{};
        // Interatomic potential choice for simulation tools.
        int potential = 0;
        double epsilon = 0.0104;
        double sigma = 3.40;
        double cutoff = 8.5;
        std::string potentialFile;
        int eamFormat = 0;
    };
    struct Result { std::filesystem::path output; std::string report; std::filesystem::path structures;
                    std::vector<atomforge::science::PlotSpec> plots;
                    std::vector<atomforge::science::AtomProperty> properties;
                    std::vector<std::pair<std::string, std::string>> files;
                    std::vector<std::pair<std::string, std::string>> summary;
                    std::vector<std::pair<std::string, std::vector<std::vector<double>>>> matrices; };
    static std::vector<std::pair<std::string, std::string>> summaryOf(const atomforge::science::Json& result);
    static std::vector<std::pair<std::string, std::vector<std::vector<double>>>> matricesOf(const atomforge::science::Json& result);
    void drawFieldGrid(const Structure& structure, int group, int skip, int columns);
    void drawTrajectoryLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult);
    void drawPropertiesLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult);
    bool isSource(std::size_t index) const;
    int activeInput() const;
    bool hasFields(int group, int skip) const;
    void useActiveStructure(std::size_t index, const Structure& structure);
    void drawField(std::size_t index, const Structure& structure);
    void drawFields(const Structure& structure, int group, int skip);
    void startCalculation(const Structure& structure);
    void drawStatus();
    bool drawRunButton(const char* label, const Structure& structure, float width);
    void drawSummary(const char* id, int pairsPerRow = 1);
    static std::string prettyLabel(const std::string& key);
    bool suppliedByActive(std::size_t index) const;
    void drawPropertyControls(const Structure& structure, const ColourAtoms& colourAtoms);
    void drawPlots(float height);
    void drawSaveButtons(const std::function<void(Structure&)>& loadResult);
    void drawDetails();
    void drawAtomLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult, const ColourAtoms& colourAtoms);
    void drawPlotLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult);
    void drawSimulationLayout(const Structure& structure, const std::function<void(Structure&)>& loadResult);
    void drawGeneratorLayout(const Structure& structure);
    void selectTool(int index);
    bool m_open = false;
    int m_tool = -1;
    int m_pickerTarget = -1;
    int m_property = 0;
    int m_plot = 0;
    int m_file = 0;
    bool m_useActive = true;
    bool m_autoColour = true;
    std::string m_trajectory;
    std::vector<Field> m_fields;
    PathPicker m_picker;
    atomforge::BackgroundTask<Result> m_task;
    std::string m_error;
    Result m_result;
};
