#pragma once

#include <algorithm>
#include <atomic>
#include <stdexcept>

namespace atomforge
{
struct TaskControl
{
    std::atomic<bool> cancelled{false};
    // Negative means no numerical progress estimate is available for this stage.
    std::atomic<float> progress{-1.0f};
};

inline thread_local TaskControl* activeTaskControl = nullptr;

class TaskControlScope
{
public:
    explicit TaskControlScope(TaskControl* control) : previous(activeTaskControl) { activeTaskControl=control; }
    ~TaskControlScope() { activeTaskControl=previous; }
private:
    TaskControl* previous;
};

inline void taskCheckpoint()
{
    if (activeTaskControl && activeTaskControl->cancelled.load(std::memory_order_relaxed))
        throw std::runtime_error("Calculation cancelled");
}

inline void taskProgress(double fraction)
{
    taskCheckpoint();
    if (activeTaskControl)
        activeTaskControl->progress.store(static_cast<float>(std::clamp(fraction,0.0,1.0)),std::memory_order_relaxed);
}
}
