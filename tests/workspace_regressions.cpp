#include "model/Workspace.h"
#include "io/Trajectory.h"
#include "util/TrajectoryPlayback.h"
#include <chrono>
#include <cmath>
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

        // LAMMPS custom dump: orthogonal box, two frames, explicit element column.
        const auto lammpsOrtho=root/"ortho.lammpstrj";
        { std::ofstream file(lammpsOrtho);
          file << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n"
                  "0.0 4.0\n0.0 5.0\n0.0 6.0\nITEM: ATOMS id type element x y z\n1 1 Cu 0.0 0.0 0.0\n2 2 O 1.0 1.0 1.0\n"
                  "ITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n"
                  "0.0 4.0\n0.0 5.0\n0.0 6.0\nITEM: ATOMS id type element x y z\n1 1 Cu 0.1 0.0 0.0\n2 2 O 1.1 1.0 1.0\n"; }
        const auto lammpsFrames=loadLammpsTrajectory(lammpsOrtho.string());
        check(lammpsFrames.size()==2,"LAMMPS trajectory frame count changed");
        check(lammpsFrames[0].atoms[0].atomicNumber==29 && lammpsFrames[0].atoms[1].atomicNumber==8,
              "LAMMPS 'element' column must map to atomic numbers, not 'type'");
        check(lammpsFrames[0].hasUnitCell && lammpsFrames[0].cellVectors[0][0]==4 &&
              lammpsFrames[0].cellVectors[1][1]==5 && lammpsFrames[0].cellVectors[2][2]==6,
              "LAMMPS orthogonal box bounds lost");
        check(lammpsFrames[1].atoms[0].x==0.1,"LAMMPS second frame not read");
        check(loadTrajectory(lammpsOrtho.string()).size()==2,"loadTrajectory must auto-detect a LAMMPS dump");
        const auto xyzFresh=root/"fresh.xyz";
        { std::ofstream file(xyzFresh); file << "1\nfixture\nCu 0 0 0\n"; }
        check(loadTrajectory(xyzFresh.string()).size()==1,"loadTrajectory must still auto-detect XYZ");

        // Triclinic box bounds (LAMMPS "xy xz yz" convention): bounding box
        // (0,4)/(0,5)/(0,6) with tilts xy=1, xz=0.5, yz=0.2 gives true cell
        // edges a=(2.5,0,0), b=(1,4.8,0), c=(0.5,0.2,6) -- verified by hand
        // against the LAMMPS manual's triclinic-box conversion formula.
        const auto lammpsTri=root/"tri.lammpstrj";
        { std::ofstream file(lammpsTri);
          file << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS xy xz yz pp pp pp\n"
                  "0.0 4.0 1.0\n0.0 5.0 0.5\n0.0 6.0 0.2\nITEM: ATOMS id type element x y z\n1 1 Cu 0.0 0.0 0.0\n"; }
        const auto triFrames=loadLammpsTrajectory(lammpsTri.string());
        check(std::abs(triFrames[0].cellVectors[0][0]-2.5)<1e-9 && triFrames[0].cellVectors[1][0]==1 &&
              std::abs(triFrames[0].cellVectors[1][1]-4.8)<1e-9 && triFrames[0].cellVectors[2][0]==.5 &&
              triFrames[0].cellVectors[2][1]==.2 && triFrames[0].cellVectors[2][2]==6,
              "LAMMPS triclinic box-bounds conversion is wrong");

        // A LAMMPS dump without an explicit "element" column must be
        // rejected outright: the numeric "type" column is an arbitrary
        // per-simulation label, never reliably an atomic number.
        const auto lammpsNoElement=root/"noelement.lammpstrj";
        { std::ofstream file(lammpsNoElement);
          file << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
                  "0.0 4.0\n0.0 5.0\n0.0 6.0\nITEM: ATOMS id type x y z\n1 1 0.0 0.0 0.0\n"; }
        bool rejectedNoElement=false;
        try { loadLammpsTrajectory(lammpsNoElement.string()); } catch (const std::exception&) { rejectedNoElement=true; }
        check(rejectedNoElement,"LAMMPS dump without an 'element' column must be rejected, not guessed");

        // TrajectoryPlayback::tick(): pure playback state machine.
        {
            TrajectoryPlayback playback;
            check(!playback.tick(false,false,false,0,10,0.0),"No input and no frames must not report a change");
            const bool startedPlaying=playback.tick(true,false,false,3,10,0.0);
            check(startedPlaying && playback.playing && playback.frame==0,
                  "Pressing Play must report a change immediately (this was the original bug: changed was computed as !playing)");
            check(!playback.tick(false,false,false,3,10,0.05),"No autoplay step before 1/fps has elapsed");
            check(playback.tick(false,false,false,3,10,0.11) && playback.frame==1,"Autoplay must advance one frame after 1/fps");
            check(playback.tick(false,false,false,3,10,0.22) && playback.frame==2,"Autoplay must keep advancing");
            check(playback.tick(false,false,false,3,10,0.33) && playback.frame==0,"Autoplay must wrap around at the frame count");
            const bool paused=playback.tick(true,false,false,3,10,0.40);
            check(paused && !playback.playing,"Pressing Pause must stop playback and report a change");
            check(playback.tick(false,true,false,3,10,0.50) && playback.frame==0 && !playback.playing,
                  "Jump to first frame must reset frame/playing and report a change");
            playback.frame=2;
            check(playback.tick(false,false,true,3,10,0.60),"A slider drag alone must report a change");
        }
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
