#pragma once

// Independent per-item work spread over threads. Results are identical to a
// serial loop when each item writes only its own output.

#include "util/TaskControl.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace atomforge::science
{
// Calls body(begin, end) on chunks of [0, n). The calling thread works too and
// alone reports progress and checks cancellation (the task control is
// per-thread), stopping the others when the task is cancelled or a chunk
// throws; the first exception is rethrown after all threads have stopped.
template <class Body>
void parallelChunks(std::size_t n, Body body, std::size_t chunk = 64)
{
    const std::size_t hardware = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t threads = std::min<std::size_t>(std::min<std::size_t>(hardware, 32), (n + chunk - 1) / std::max<std::size_t>(chunk, 1));
    if (threads <= 1) {
        for (std::size_t begin = 0; begin < n; begin += chunk) {
            taskProgress(static_cast<double>(begin) / static_cast<double>(n));
            body(begin, std::min(n, begin + chunk));
        }
        return;
    }
    std::atomic<std::size_t> next{0}, done{0};
    std::atomic<bool> stop{false};
    std::exception_ptr failure;
    std::mutex failureLock;
    const auto run = [&](bool reporting) {
        while (!stop.load(std::memory_order_relaxed)) {
            const std::size_t begin = next.fetch_add(chunk);
            if (begin >= n) break;
            try {
                if (reporting) taskProgress(static_cast<double>(done.load()) / static_cast<double>(n));
                const std::size_t end = std::min(n, begin + chunk);
                body(begin, end);
                done.fetch_add(end - begin);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failureLock);
                if (!failure) failure = std::current_exception();
                stop = true;
            }
        }
    };
    // Checked here as well: the other threads may take every chunk first.
    taskCheckpoint();
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < threads; ++t) pool.emplace_back(run, false);
    run(true);
    for (auto& thread : pool) thread.join();
    if (failure) std::rethrow_exception(failure);
    taskCheckpoint();
}
}  // namespace atomforge::science
