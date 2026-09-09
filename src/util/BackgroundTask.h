#pragma once

#include <chrono>
#include <exception>
#include <future>
#include <optional>
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
    [[nodiscard]] bool running() const noexcept { return m_future.valid(); }
    [[nodiscard]] const std::optional<Result>& result() const noexcept { return m_result; }
    [[nodiscard]] const std::string& error() const noexcept { return m_error; }

    template<class Work>
    bool start(Work&& work)
    {
        if (running()) return false;
        clearResult();
        try
        {
            m_future = std::async(std::launch::async, std::forward<Work>(work));
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
    std::future<Result> m_future;
};
} // namespace atomforge
