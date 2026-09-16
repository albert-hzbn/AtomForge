#pragma once
#include "util/TaskControl.h"

#include <chrono>
#include <exception>
#include <future>
#include <optional>
#include <memory>
#include <string>
#include <utility>

namespace atomforge
{
// Single-owner task state: start/poll/access belong to the UI thread. Workers
// receive value captures and return a result; they never mutate the owner.
// Destruction waits for an outstanding worker through std::future's RAII contract.
template<class Result>
class BackgroundTask
{
public:
    BackgroundTask() = default;
    BackgroundTask(BackgroundTask&&) = default;
    BackgroundTask& operator=(BackgroundTask&&) = default;
    ~BackgroundTask() { cancel(); }
    [[nodiscard]] bool running() const noexcept { return m_future.valid(); }
    [[nodiscard]] const std::optional<Result>& result() const noexcept { return m_result; }
    [[nodiscard]] const std::string& error() const noexcept { return m_error; }
    void cancel() noexcept { if (m_control) m_control->cancelled=true; }
    [[nodiscard]] float progress() const noexcept { return m_control ? m_control->progress.load() : -1.0f; }

    template<class Work>
    bool start(Work&& work)
    {
        if (running()) return false;
        clearResult();
        try
        {
            m_control = std::make_shared<TaskControl>();
            m_future = std::async(std::launch::async, [control=m_control, work=std::forward<Work>(work)]() mutable {
                TaskControlScope scope(control.get());
                taskCheckpoint();
                auto result=work();
                taskCheckpoint();
                taskProgress(1.0);
                return result;
            });
            return true;
        }
        catch (const std::exception& error)
        {
            m_error = error.what();
            return false;
        }
    }

    // Returns true exactly once when a result or worker exception is consumed.
    bool poll()
    {
        if (!running() || m_future.wait_for(std::chrono::seconds(0)) != std::future_status::ready)
            return false;
        try { m_result = m_future.get(); }
        catch (const std::exception& error) { m_error = error.what(); }
        catch (...) { m_error = "Analysis failed with an unknown exception."; }
        return true;
    }

    void clearResult()
    {
        m_result.reset();
        m_error.clear();
    }

private:
    std::optional<Result> m_result;
    std::string m_error;
    std::shared_ptr<TaskControl> m_control;
    std::future<Result> m_future;
};
} // namespace atomforge
