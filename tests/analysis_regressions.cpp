#include "algorithms/CommonNeighbourAnalysis.h"
#include "algorithms/RadialDistributionAnalysis.h"
#include "util/BackgroundTask.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

namespace
{
void check(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(1); }
}

template<class T>
void awaitResult(atomforge::BackgroundTask<T>& task)
{
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (!task.poll())
    {
        check(std::chrono::steady_clock::now() < deadline, "Worker did not complete");
        std::this_thread::yield();
    }
}
} // namespace

int main()
{
    using namespace atomforge::analysis;
    Structure structure;
    check(!computeCna(structure).valid, "Empty CNA must fail gracefully");
    check(!computeRdf(structure).valid, "Empty RDF must fail gracefully");
    AtomSite iron;
    iron.symbol = "Fe";
    iron.atomicNumber = 26;
    structure.atoms.push_back(iron);
    iron.x = 1.0;
    structure.atoms.push_back(iron);

    const auto cna = computeCna(structure);
    check(cna.valid && cna.atomCount == 2 && cna.pairCount == 1,
          "CNA must find the single Fe-Fe bond");
    check(cna.atomRows.size() == 2 && cna.atomRows[0].coordination == 1
          && cna.atomRows[1].coordination == 1, "CNA coordination changed");

    RdfParams options;
    options.rMax = 2.0f;
    options.binCount = 8;
    options.normalize = false;
    options.enableDistortionAnalysis = false;
    const auto rdf = computeRdf(structure, options);
    check(rdf.valid && rdf.bins.size() == 8 && rdf.bins[4].rawCount == 2.0f,
          "RDF must count both directed Fe-Fe pairs at r=1");
    check(rdf.bins.back().cumulative == 1.0f, "RDF cumulative coordination changed");
    options.refAtomicNumberFilter = 28;
    check(!computeRdf(structure, options).valid, "Missing species must fail gracefully");

    structure.hasUnitCell = true;
    structure.cellVectors = {{{4, 0, 0}, {0, 4, 0}, {0, 0, 4}}};
    structure.atoms[0].x = 0.1;
    structure.atoms[1].x = 3.9;
    options.refAtomicNumberFilter = 0;
    options.rMax = 1.0f;
    const auto periodicRdf = computeRdf(structure, options);
    check(periodicRdf.pbcUsed && periodicRdf.bins[1].rawCount == 2.0f,
          "Periodic RDF must find the boundary-spanning pair");
    check(computeCna(structure).pairCount == 1, "Periodic CNA pair changed");
    CnaParams invalidCna;
    invalidCna.cutoffScale = std::numeric_limits<float>::quiet_NaN();
    check(!computeCna(structure, invalidCna).valid, "CNA must reject a non-finite cutoff");
    options.rMax = std::numeric_limits<float>::quiet_NaN();
    check(!computeRdf(structure, options).valid, "RDF must reject a non-finite radius");

    atomforge::BackgroundTask<int> task;
    std::promise<void> release;
    auto gate = release.get_future();
    int value = 42;
    check(task.start([snapshot = value, gate = std::move(gate)]() mutable {
        gate.wait();
        return snapshot;
    }), "Task failed to start");
    value = 99;
    check(task.running() && !task.poll(), "Polling an unfinished task must not block");
    check(!task.start([] { return -1; }), "Running task must reject replacement");
    release.set_value();
    awaitResult(task);
    check(!task.running() && task.result() == 42, "Worker must use its captured snapshot");
    check(!task.poll(), "Task result must be consumed exactly once");

    check(task.start([]() -> int { throw std::runtime_error("expected failure"); }),
          "Completed task must be reusable");
    awaitResult(task);
    check(!task.result() && task.error() == "expected failure",
          "Worker exception must reach its owner without terminating");
    check(task.start([] { return 7; }), "Failed task must be reusable");
    awaitResult(task);
    check(task.result() == 7 && task.error().empty(), "Restart must clear previous errors");

    std::atomic<bool> finished{false};
    {
        atomforge::BackgroundTask<int> owned;
        owned.start([&finished] { finished = true; return 1; });
    }
    check(finished, "Task destruction must wait for its worker");
    std::cout << "Analysis and background-task regressions passed\n";
}
