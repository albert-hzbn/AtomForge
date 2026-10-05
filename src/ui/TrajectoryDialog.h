#pragma once
#include "io/Trajectory.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include "util/TrajectoryPlayback.h"
#include <functional>

struct TrajectoryDialog
{
    void drawMenuItem();
    void draw(Structure& structure,const std::function<void(Structure&)>& update);
    bool isPlaying() const { return playback.playing; }
    // File of the currently loaded trajectory (empty when none).
    const std::string& loadedPath() const { return loaded; }
private:
    bool open=false;
    TrajectoryPlayback playback;
    float fps=12;
    std::vector<Structure> frames;
    PathPicker picker;
    atomforge::BackgroundTask<std::vector<Structure>> task;
    std::string pending;
    std::string loaded;
    std::string error;
};
