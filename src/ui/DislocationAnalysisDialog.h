#pragma once

#include "algorithms/DislocationFit.h"
#include "algorithms/DragPrep.h"
#include "algorithms/NyeTensor.h"
#include "algorithms/PatternMatch.h"
#include "algorithms/VitekMap.h"
#include "model/Structure.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"

#include <string>

// Post-processing analyses that compare a deformed structure against its
// undeformed reference: Nye tensor, Vitek differential-displacement map,
// crystallographic pattern matching, drag/migration-barrier preparation,
// and dislocation position/Burgers-vector fitting. Each is a native
// reimplementation of the corresponding BABEL program (see the header
// comment in each src/algorithms/*.h for the literature followed). Like
// WannierAnalysisDialog/LobsterAnalysisDialog this is a self-contained
// modal analysis dialog: load reference/deformed structures, compute,
// inspect -- no workspace persistence, and results are exported to a file
// rather than written back into the live scene (except drag preparation's
// interpolated structure, which is saved directly since it IS a structure).
struct DislocationAnalysisDialog
{
    void drawMenuItem();
    void drawDialog();

private:
    void loadReference(const std::string& path);
    void loadDeformed(const std::string& path);

    bool m_openRequested = false;
    bool m_open = false;

    PathPicker m_referencePicker;
    std::string m_referencePath;
    std::string m_referenceError;
    Structure m_reference;
    bool m_referenceLoaded = false;

    PathPicker m_deformedPicker;
    std::string m_deformedPath;
    std::string m_deformedError;
    Structure m_deformed;
    bool m_deformedLoaded = false;

    // Nye tensor
    float m_nyeCutoff = 3.0f;
    bool m_nyeNoPbc = false;
    atomforge::BackgroundTask<atomforge::NyeTensorResult> m_nyeTask;
    atomforge::NyeTensorResult m_nyeResult;
    PathPicker m_nyeExportPicker;

    // Vitek differential-displacement map
    float m_vitekLine[3] = {0.0f, 0.0f, 1.0f};
    float m_vitekBurgers = 0.0f;
    float m_vitekCutoff = 0.0f; // <=0: auto-detect
    bool m_vitekNoPbc = false;
    atomforge::BackgroundTask<atomforge::VitekResult> m_vitekTask;
    atomforge::VitekResult m_vitekResult;
    PathPicker m_vitekExportPicker;

    // Pattern matching
    float m_patternCutoff = 3.2f;
    bool m_patternNoPbc = false;
    atomforge::BackgroundTask<atomforge::CrystalPattern> m_patternBuildTask;
    atomforge::CrystalPattern m_pattern;
    float m_patternAngleThreshold = 10.0f;
    atomforge::BackgroundTask<atomforge::PatternDetectResult> m_patternDetectTask;
    atomforge::PatternDetectResult m_patternResult;
    PathPicker m_patternExportPicker;

    // Drag/migration-barrier preparation
    float m_dragZeta = 0.5f;
    bool m_dragClipDisplacement = false;
    bool m_dragNoPbc = false;
    atomforge::BackgroundTask<atomforge::DragPrepResult> m_dragTask;
    atomforge::DragPrepResult m_dragResult;
    PathPicker m_dragSavePicker;

    // Dislocation position/Burgers-vector fit
    float m_fitLine[3] = {0.0f, 0.0f, 1.0f};
    float m_fitCutoff = 3.0f;
    float m_fitArea = 1.0f;
    bool m_fitNoPbc = false;
    atomforge::BackgroundTask<atomforge::DislocationFitResult> m_fitTask;
    atomforge::DislocationFitResult m_fitResult;
};
