#pragma once
#include "io/Trajectory.h"
#include "science/GifWriter.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include "util/TrajectoryPlayback.h"
#include <functional>
#include <memory>

struct TrajectoryDialog
{
    void drawMenuItem();
    void draw(Structure& structure,const std::function<void(Structure&)>& update);
    // Also true while exporting a GIF, so frame changes stay out of undo history.
    bool isPlaying() const { return playback.playing || gifFrame>=0; }
    // File of the currently loaded trajectory (empty when none).
    const std::string& loadedPath() const { return loaded; }
    // GIF export: draw() shows one frame, the main loop renders it and hands the
    // pixels back here (RGBA, rows top to bottom).
    bool gifCapturePending() const { return gifFrame>=0 && gifShown; }
    void addGifFrame(const std::vector<unsigned char>& rgba,int width,int height);
private:
    void finishGif(const std::string& message);
    bool open=false;
    TrajectoryPlayback playback;
    float fps=12;
    std::vector<Structure> frames;
    PathPicker picker;
    PathPicker gifPicker;
    atomforge::BackgroundTask<std::vector<Structure>> task;
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
};
