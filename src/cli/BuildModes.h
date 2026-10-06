#pragma once
// Builder modes of `AtomForge --build <mode>`. Each mode has a run function and
// its own help page; CLIMode.cpp registers them in one table.
//
// Adding a mode: implement printHelpX/runX in the matching family file,
// declare them here and add a row to kBuildModes in CLIMode.cpp.

namespace cli
{
// BuildCrystalCLI.cpp
void printHelpBulk();
int runBulk(int argc, char* argv[]);
void printHelpCustom();
int runCustom(int argc, char* argv[]);
void printHelpSSS();
int runSSS(int argc, char* argv[]);
void printHelpPrimitive();
int runPrimitive(int argc, char* argv[]);
void printHelpSurface();
int runSurface(int argc, char* argv[]);
void printHelpSQS();
int runSQS(int argc, char* argv[]);
void printHelpVacancy();
int runVacancy(int argc, char* argv[]);
void printHelpStrain();
int runStrain(int argc, char* argv[]);
// BuildDefectCLI.cpp
void printHelpDislocation();
int runDislocation(int argc, char* argv[]);
void printHelpGB();
int runGB(int argc, char* argv[]);
void printHelpInterface();
int runInterface(int argc, char* argv[]);
void printHelpStackingFault();
int runStackingFault(int argc, char* argv[]);
// BuildNanoCLI.cpp
void printHelpPoly();
int runPoly(int argc, char* argv[]);
void printHelpNano();
int runNano(int argc, char* argv[]);
void printHelpAmorphous();
int runAmorphous(int argc, char* argv[]);
void printHelpNanowire();
int runNanowire(int argc, char* argv[]);
void printHelpCoreShell();
int runCoreShell(int argc, char* argv[]);
}
