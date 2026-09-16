#include "algorithms/CommonNeighbourAnalysis.h"
#include "algorithms/AngularDistributionAnalysis.h"
#include "algorithms/RadialDistributionAnalysis.h"
#include "util/BackgroundTask.h"

#include <atomic>
#include <algorithm>
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

    // Ideal periodic crystals exercise the full signature distribution, not just pairs.
    for (bool bcc : {false, true}) {
        Structure crystal;
        crystal.hasUnitCell = true;
        crystal.cellVectors = {{{4, 0, 0}, {0, 4, 0}, {0, 0, 4}}};
        const std::vector<std::array<double, 3>> basis = bcc
            ? std::vector<std::array<double, 3>>{{0,0,0}, {.5,.5,.5}}
            : std::vector<std::array<double, 3>>{{0,0,0}, {0,.5,.5}, {.5,0,.5}, {.5,.5,0}};
        for (int x=0; x<4; ++x) for (int y=0; y<4; ++y) for (int z=0; z<4; ++z)
            for (const auto& point : basis) {
                AtomSite atom = iron;
                atom.x=x+point[0]; atom.y=y+point[1]; atom.z=z+point[2];
                crystal.atoms.push_back(atom);
            }
        CnaParams params;
        params.cutoffRadius = bcc ? 1.1f : .85f;
        const auto result = computeCna(crystal, params);
        check(result.valid, "Ideal crystal CNA must succeed");
        for (const auto& row : result.atomRows)
            check(row.environment == (bcc ? "BCC-like" : "FCC-like"),
                  "CNA must identify every atom in ideal periodic FCC/BCC");
    }

    {
        Structure hcp;
        const double root3=std::sqrt(3.0),height=std::sqrt(8.0/3);
        hcp.hasUnitCell=true;
        hcp.cellVectors={{{4,0,0},{-2,2*root3,0},{0,0,4*height}}};
        for(int x=0;x<4;++x) for(int y=0;y<4;++y) for(int z=0;z<4;++z)
            for(int basis=0;basis<2;++basis) {
                AtomSite atom=iron;
                atom.x=x-.5*y+.5*basis; atom.y=root3*y/2+root3*basis/6;
                atom.z=height*(z+.5*basis); hcp.atoms.push_back(atom);
            }
        CnaParams params; params.cutoffRadius=1.1f;
        const auto result=computeCna(hcp,params);
        check(result.valid,"HCP CNA must succeed");
        for(const auto& row:result.atomRows) check(row.environment=="HCP-like","Ideal HCP classification failed");
        Structure ico; AtomSite origin=iron; origin.x=origin.y=origin.z=0; ico.atoms.push_back(origin);
        const double phi=(1+std::sqrt(5.0))/2, scale=std::sqrt(1+phi*phi);
        for(int axis=0;axis<3;++axis) for(int a:{-1,1}) for(int b:{-1,1}) {
            const double v[]={0,a/scale,b*phi/scale}; AtomSite atom=origin;
            atom.x=v[axis]; atom.y=v[(axis+1)%3]; atom.z=v[(axis+2)%3]; ico.atoms.push_back(atom);
        }
        check(computeCna(ico,params).atomRows.front().environment=="ICO-like","Icosahedral centre classification failed");
    }

    Structure triplet;
    AtomSite centre = iron; centre.symbol="O"; centre.x=centre.y=centre.z=0;
    AtomSite first = centre; first.symbol="H"; first.x=1;
    AtomSite second = centre; second.symbol="C"; second.y=1;
    triplet.atoms={centre,first,second};
    AdfParams angles;
    angles.usePbc=false; angles.rCutoff=1.1f; angles.binCount=180;
    angles.centreMode=AdfCentreMode::ByPair; angles.centreSymbol="O";
    angles.neighSymbol1="H"; angles.neighSymbol2="C";
    angles.smoothPasses=0;
    const auto forward = computeADF(triplet,angles);
    std::swap(triplet.atoms[1],triplet.atoms[2]);
    const auto reverse = computeADF(triplet,angles);
    check(forward.valid && reverse.valid && forward.nTriplets==1 && reverse.nTriplets==1,
          "ADF element-pair selection must not depend on atom order");
    check(forward.bins[90].count==1 && reverse.bins[90].count==1,
          "Orthogonal bonds must produce a 90 degree triplet");
    angles.smoothPasses=4;
    const auto smoothed=computeADF(triplet,angles);
    check(smoothed.bins[90].count==1, "ADF must retain raw counts when smoothing");
    for (int bins : {0, -1, 100001}) {
        angles.binCount=bins;
        check(!computeADF(triplet,angles).valid, "ADF must reject invalid bin counts");
    }
    angles.binCount=180; angles.rCutoff=std::numeric_limits<float>::quiet_NaN();
    check(!computeADF(triplet,angles).valid, "ADF must reject a NaN cutoff");

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

    std::promise<void> cancelGate;
    auto cancelReady=cancelGate.get_future();
    check(task.start([gate=std::move(cancelReady)]() mutable {
        gate.wait(); atomforge::taskCheckpoint(); return 99;
    }), "Cancellable task must start");
    task.cancel(); cancelGate.set_value(); awaitResult(task);
    check(!task.result() && task.error()=="Calculation cancelled", "Cancelled work must not publish a result");
    check(task.start([] { atomforge::taskProgress(.5); return 3; }), "Cancelled task must be reusable");
    awaitResult(task);
    check(task.result()==3 && task.progress()==1, "Completed task must report full progress");

    std::atomic<bool> finished{false};
    std::promise<void> entered;
    auto started=entered.get_future();
    {
        atomforge::BackgroundTask<int> owned;
        check(owned.start([&finished,&entered] {
            entered.set_value();
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            finished=true; return 1;
        }),"Owned task must start");
        check(started.wait_for(std::chrono::seconds(5))==std::future_status::ready,"Owned worker did not start");
    }
    check(finished, "Task destruction must wait for its worker");
    std::cout << "Analysis and background-task regressions passed\n";
}
