#pragma once
// Trajectory lines (as OVITO's Generate Trajectory Lines): the continuous,
// unwrapped paths of atoms through a trajectory, drawn over the view up to the
// frame shown.

#include <glm/glm.hpp>

#include <vector>

struct ImDrawList;

struct TrajectoryLinesOverlay
{
    std::vector<std::vector<glm::vec3>> paths;  // one per traced atom, one point per sampled frame
    std::vector<int> frames;                    // trajectory frame of each point
    int currentFrame = 0;                       // draw points up to this frame
    bool visible = false;
    bool wholePath = false;                     // the whole path instead of up to the current frame
    float width = 1.5f;
};

TrajectoryLinesOverlay& trajectoryLinesOverlay();

void drawTrajectoryLinesOverlay(ImDrawList* drawList, const glm::mat4& projection, const glm::mat4& view,
                                int framebufferWidth, int framebufferHeight);
