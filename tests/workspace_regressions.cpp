#include "model/Workspace.h"
#include "io/Trajectory.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

void check(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }

int main()
{
    const auto root=std::filesystem::temp_directory_path()/
        ("atomforge-project-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    try {
        const auto xyz=root/"trajectory.xyz";
        { std::ofstream file(xyz); file << "2\nLattice=\"3 0 0 0 3 0 0 0 3\" Properties=species:S:1:pos:R:3\nCu 0 0 0\nO 1 1 1\n1\nnext\nH .5 .5 .5\n"; }
        const auto frames=loadXyzTrajectory(xyz.string());
        check(frames.size()==2 && frames[0].atoms[0].atomicNumber==29,"Trajectory elements or frame count changed");
        check(frames[0].hasUnitCell && frames[0].cellVectors[2][2]==3,"Trajectory cell lost");
        { std::ofstream file(xyz,std::ios::app); file << "2\ntruncated\nH 0 0 0\n"; }
        bool rejected=false;
        try { loadXyzTrajectory(xyz.string()); } catch (const std::exception&) { rejected=true; }
        check(rejected,"Trajectory must reject partial final frames");
        atomforge::Workspace source;
        source.title="Project with spaces, quotes \" and newline\n";
        AtomSite atom; atom.symbol="Cu"; atom.atomicNumber=29; atom.x=.123456789012345;
        source.structure.atoms={atom}; source.structure.hasUnitCell=true;
        source.structure.cellVectors={{{2,0,0},{.2,2,0},{0,0,3}}};
        source.structure.grainColors={{{.1f,.2f,.3f}}}; source.structure.grainRegionIds={7};
        atomforge::electronic::Grid grid;
        grid.values={1,2,3,4,5,6,7,8}; grid.name="charge"; grid.unit="e/A^3";
        source.volume.fields={grid}; source.reference=source.volume;
        source.surface.vertices={{0,0,0},{1,0,0},{0,1,0}}; source.surface.colors={1,2,3};
        source.table={1,2,3,4}; source.columns=2; source.heading="x,y";
        source.settings={{"opacity",.35}}; source.sliceSettings={{"normal_x",1}};
        source.camera={1,2,3,4,5,6,7}; source.history={"load", "subtract"};
        const auto path=(root/"saved.afproject").string();
        atomforge::saveWorkspace({source,source},path);
        const auto restored=atomforge::loadWorkspace(path);
        check(restored.size()==2 && restored[0].title==source.title,"Project tabs or strings changed");
        check(restored[0].structure.atoms[0].x==atom.x,"Coordinates lost precision");
        check(restored[0].structure.grainRegionIds==source.structure.grainRegionIds,"Grain metadata lost");
        check(restored[0].volume.fields[0].values==grid.values && restored[0].reference.fields[0].values==grid.values,"Density fields changed");
        check(restored[0].table==source.table && restored[0].camera==source.camera && restored[0].settings==source.settings,"Workspace settings changed");
        check(restored[0].history==source.history && restored[0].surface.colors==source.surface.colors,"History or surface changed");
        source.volume.fields[0].values[0]=std::numeric_limits<double>::quiet_NaN();
        bool failed=false;
        try { atomforge::saveWorkspace({source},path); } catch (const std::exception&) { failed=true; }
        check(failed,"Project must reject non-finite values");
        check(atomforge::loadWorkspace(path).size()==2,"Failed save replaced the existing project");
        const auto truncated=root/"truncated.afproject";
        std::filesystem::copy_file(path,truncated);
        std::filesystem::resize_file(truncated,std::filesystem::file_size(truncated)-5);
        failed=false;
        try { atomforge::loadWorkspace(truncated.string()); } catch (const std::exception&) { failed=true; }
        check(failed,"Truncated project must be rejected");
        std::ofstream corrupt(root/"corrupt.afproject",std::ios::binary);
        for (int i=0;i<8;++i) corrupt.put(static_cast<char>(255));
        corrupt.close(); failed=false;
        try { atomforge::loadWorkspace((root/"corrupt.afproject").string()); } catch (const std::exception&) { failed=true; }
        check(failed,"Unbounded string allocation must be rejected");
        std::filesystem::remove_all(root);
        std::cout << "Workspace round-trip and failure recovery passed\n";
    } catch (...) { std::filesystem::remove_all(root); throw; }
}
