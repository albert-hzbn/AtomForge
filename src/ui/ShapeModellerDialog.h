#pragma once
// Build > Shape Modeller: crystals carved into shapes modelled as in CAD
// (algorithms/ShapeModel). The shape is a script of primitives, booleans,
// transforms, patterns and regions, inserted from menus or typed, with live
// checking, section views through the shape and an atom-count estimate. The
// result is the build-shape pipeline step.

#include "algorithms/ShapeModel.h"
#include "model/Structure.h"
#include "pipeline/Pipeline.h"
#include "ui/PathPicker.h"
#include "ui/StepEditing.h"

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

class ShapeModellerDialog
{
public:
    void drawMenuItem(bool enabled);
    void open() { m_edit.finish(); m_open = true; }
    void close() { m_open = false; m_edit.finish(); }
    bool isOpen() const { return m_open; }

    struct Callbacks
    {
        std::function<void(Structure&)> build;                    // a carved structure (opened in a new tab)
        std::function<void(Structure&)> update;                   // the structure edited in place (undoable)
        std::function<void(const std::vector<int>&)> select;      // atoms to select in the view
        std::function<void(const std::vector<atomforge::pipeline::Modifier>&)> addToPipeline;
    };
    void draw(Structure& structure, const Callbacks& callbacks);

    // Opens the dialog on a build-shape step; "Update Step" writes the settings back.
    bool editStep(StepEdit edit);
    atomforge::pipeline::Json stepParameters() const;
    // The script being edited (tests and the tour).
    void setScript(const std::string& script);

private:
    void reparse(const Structure& structure);
    void drawSection(const Structure& structure, float size);
    void insert(const std::string& snippet);

    bool m_open = false;
    StepEdit m_edit;
    std::array<char, 32768> m_script{};
    atomforge::pipeline::Json m_parameters = atomforge::pipeline::Json::object();
    // The parsed script, for the status line and the section view.
    std::string m_parsedScript;
    const void* m_parsedFor = nullptr;   // the structure whose cell the facets used
    std::array<double, 9> m_parsedCell{};
    std::optional<atomforge::shape::Model> m_model;
    std::string m_parseError;
    int m_errorLine = 0;
    double m_volume = 0;
    // Section view: plane (0 xy, 1 xz, 2 yz), position as a fraction of the extent.
    int m_plane = 1;
    float m_sectionAt = 0.5f;
    std::vector<unsigned int> m_raster;
    int m_rasterSize = 0;
    std::string m_rasterKey;
    std::string m_message;
    bool m_error = false;
    char m_meshPath[1024]{};
    PathPicker m_meshPicker;
};
