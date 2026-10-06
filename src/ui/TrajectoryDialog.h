#pragma once
#include "io/Trajectory.h"
#include "science/AtomProperties.h"
#include "science/GifWriter.h"
#include "science/ScienceData.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include "util/TrajectoryPlayback.h"
#include <functional>
#include <memory>

struct TrajectoryDialog
{
    void drawMenuItem();
    // showAtomProperty switches the view to atom-property colouring (per-frame analysis).
    void draw(Structure& structure, const std::function<void(Structure&)>& update,
              const std::function<void()>& showAtomProperty = {});
    // Also true while exporting a GIF, so frame changes stay out of undo history.
    bool isPlaying() const { return playback.playing || gifFrame>=0; }
    // File of the currently loaded trajectory (empty when none).
    const std::string& loadedPath() const { return loaded; }
    // GIF export: draw() shows one frame, the main loop renders it and hands the
    // pixels back here (RGBA, rows top to bottom).
    bool gifCapturePending() const { return gifFrame>=0 && gifShown; }
    void addGifFrame(const std::vector<unsigned char>& rgba,int width,int height);
    // Per-atom analyses that can colour every displayed frame (catalog tool ids).
    static const std::vector<std::pair<std::string, std::string>>& frameAnalyses();
    // Runs a per-frame analysis on one structure (its first per-atom property).
    static atomforge::science::AtomProperty analyseFrame(const Structure& frame, const std::string& tool, double cutoff);
private:
    void finishGif(const std::string& message);
    std::size_t frameCount() const { return stream ? stream->size() : 0; }
    // Loads frame `index` into `structure` (with the per-frame analysis, when enabled).
    void show(int index, Structure& structure, const std::function<void(Structure&)>& update, const std::function<void()>& showAtomProperty);
    bool open=false;
    TrajectoryPlayback playback;
    float fps=12;
    std::shared_ptr<const atomforge::science::TrajectoryStream> stream;
    std::size_t shownAtoms=0;
    PathPicker picker;
    PathPicker gifPicker;
    atomforge::BackgroundTask<std::shared_ptr<const atomforge::science::TrajectoryStream>> task;
    std::string pending;
    std::string loaded;
    std::string error;
    std::unique_ptr<atomforge::science::GifWriter> gif;
    std::string gifPath;
    std::string gifMessage;
    int gifFrame=-1;
    int gifStride=1;
    int gifRestoreFrame=0;
    bool gifShown=false;
    bool gifRestorePending=false;
    bool analyse=false;
    int analysis=0;
    float analysisCutoff=3.0f;
    std::string analysisMessage;
};
