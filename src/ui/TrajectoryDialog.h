#pragma once
#include "io/Trajectory.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"
#include <functional>

struct TrajectoryDialog
{
    void drawMenuItem();
    void draw(Structure& structure,const std::function<void(Structure&)>& update);
    bool isPlaying() const { return playing; }
private:
    bool open=false, playing=false;
    int frame=0;
    float fps=12;
    double previous=0;
    std::vector<Structure> frames;
    PathPicker picker;
    atomforge::BackgroundTask<std::vector<Structure>> task;
    std::string error;
};
