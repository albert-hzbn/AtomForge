#pragma once

// Pure playback state machine for trajectory scrubbing, kept free of ImGui/
// GLFW so it is independently unit-testable (see tests/core_regressions.cpp).
// Callers translate this frame's UI events into the three booleans below and
// apply frames[frame] to the active structure whenever tick() returns true.
struct TrajectoryPlayback
{
    bool playing = false;
    int frame = 0;
    double previous = 0.0;

    // togglePlay: the Play/Pause button was just pressed this frame.
    // jumpToFirst: the "First frame" button was just pressed this frame.
    // sliderChanged: the frame slider was just dragged this frame.
    // frameCount: total frames available (0 if none loaded).
    // fps: playback rate; ignored unless playing.
    // now: a monotonic clock reading in seconds (e.g. ImGui::GetTime()).
    // Returns true exactly when `frame` now holds a value the caller should
    // apply -- a jump, a slider drag, a fresh (un)pause, or an autoplay step.
    bool tick(bool togglePlay, bool jumpToFirst, bool sliderChanged, int frameCount, float fps, double now)
    {
        bool changed = sliderChanged;
        if (togglePlay) { playing = !playing; previous = now; changed = true; }
        if (jumpToFirst) { frame = 0; playing = false; changed = true; }
        if (playing && frameCount > 0 && now - previous >= 1.0 / fps)
        {
            frame = (frame + 1) % frameCount;
            previous = now;
            changed = true;
        }
        return changed;
    }
};
