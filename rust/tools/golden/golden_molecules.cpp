#include "pipeline/MoleculeLibrary.h"
#include <cstdio>
int main(){ for (const auto& n : atomforge::pipeline::moleculeNames()) for (const auto& a : atomforge::pipeline::molecule(n)) std::printf("%s %s %a %a %a\n", n.c_str(), a.symbol.c_str(), a.position[0], a.position[1], a.position[2]); }
