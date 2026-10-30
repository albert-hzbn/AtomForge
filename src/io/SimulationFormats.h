#pragma once
// Input structures of simulation codes, written and read natively:
// Quantum ESPRESSO pw.x input (pwi), CASTEP .cell, DFTB+ .gen, XCrySDen .xsf
// and a CP2K &SUBSYS section (cp2k; written only).

#include "model/Structure.h"

#include <ostream>
#include <string>

namespace atomforge::io
{
// Format names: "pwi" (also "pwscf"), "cell", "gen", "xsf", "cp2k".
bool isSimulationFormat(const std::string& format);
bool canReadSimulationFormat(const std::string& format);

// Writes the structure; false with a reason (e.g. no cell where one is needed).
bool writeSimulationFormat(const Structure& structure, const std::string& format, std::ostream& out, std::string& error);
// Reads a structure (symbols, positions, cell, periodicity, fixed coordinates
// where the format has them); false with a reason.
bool readSimulationFormat(const std::string& text, const std::string& format, Structure& structure, std::string& error);
}
