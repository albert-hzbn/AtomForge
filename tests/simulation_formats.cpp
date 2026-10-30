// Simulation-code structure formats: exact round trips (positions, cell,
// periodicity, fixed coordinates) and hand-written inputs in other units.
#include "io/SimulationFormats.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
using atomforge::io::readSimulationFormat;
using atomforge::io::writeSimulationFormat;

void check(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

AtomSite site(const std::string& symbol, int z, double x, double y, double zc, std::uint8_t fixed = 0)
{
    AtomSite a;
    a.symbol = symbol;
    a.atomicNumber = z;
    a.x = x; a.y = y; a.z = zc;
    a.fixed = fixed;
    return a;
}

Structure sample()
{
    Structure s;
    s.hasUnitCell = true;
    s.cellVectors = {{{4.1, 0.0, 0.0}, {1.3, 3.9, 0.0}, {0.7, -0.4, 6.2}}};
    s.atoms = {site("Fe", 26, 0.1, 0.2, 0.3, 0b101), site("O", 8, 1.234567891, 2.5, 3.75), site("Fe", 26, 3.0, 1.0, 5.0, 0b010),
               site("Mo", 42, -0.5, 0.25, 1.125)};
    return s;
}

double gap(const Structure& a, const Structure& b)
{
    double worst = 0;
    for (std::size_t i = 0; i < a.atoms.size(); ++i)
        worst = std::max({worst, std::abs(a.atoms[i].x - b.atoms[i].x), std::abs(a.atoms[i].y - b.atoms[i].y), std::abs(a.atoms[i].z - b.atoms[i].z)});
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::abs(a.cellVectors[r][c] - b.cellVectors[r][c]));
    return worst;
}

Structure roundTrip(const Structure& s, const std::string& format)
{
    std::ostringstream out;
    std::string error;
    check(writeSimulationFormat(s, format, out, error), format + " write: " + error);
    Structure back;
    check(readSimulationFormat(out.str(), format, back, error), format + " read: " + error + "\n" + out.str());
    return back;
}
}

int main()
{
    try {
        const Structure s = sample();
        for (const std::string format : {"pwi", "cell", "gen", "xsf"}) {
            const Structure back = roundTrip(s, format);
            check(back.atoms.size() == s.atoms.size(), format + ": atom count");
            check(back.hasUnitCell, format + ": cell kept");
            check(gap(s, back) < 1e-9, format + ": positions and cell within 1e-9 A");
            for (std::size_t i = 0; i < s.atoms.size(); ++i) {
                check(back.atoms[i].symbol == s.atoms[i].symbol && back.atoms[i].atomicNumber == s.atoms[i].atomicNumber, format + ": species");
                if (format == "pwi" || format == "cell") check(back.atoms[i].fixed == s.atoms[i].fixed, format + ": fixed coordinates");
            }
        }
        // Magnetic moments: CASTEP SPIN= round trip; Quantum ESPRESSO species per moment.
        {
            Structure magnetic = sample();
            magnetic.atoms[0].magmom = 2.5;
            magnetic.atoms[2].magmom = -2.5;
            const Structure back = roundTrip(magnetic, "cell");
            check(back.atoms[0].magmom == 2.5 && back.atoms[2].magmom == -2.5 && back.atoms[1].magmom == 0, "cell SPIN= round trip");
            Structure spaced;
            std::string error;
            check(readSimulationFormat("%BLOCK LATTICE_CART\n3 0 0\n0 3 0\n0 0 3\n%ENDBLOCK LATTICE_CART\n%BLOCK POSITIONS_ABS\nNi 0 0 0 SPIN = -1.5\n%ENDBLOCK POSITIONS_ABS\n",
                                       "cell", spaced, error) && spaced.atoms[0].magmom == -1.5, "cell SPIN = x");
            std::ostringstream out;
            check(writeSimulationFormat(magnetic, "pwi", out, error), "pwi magnetic write");
            const std::string text = out.str();
            check(text.find("nspin = 2") != std::string::npos && text.find("ntyp = 4") != std::string::npos, "pwi species per moment: " + text);
            check(text.find("starting_magnetization(1) = 0.5") != std::string::npos && text.find("starting_magnetization(3) = -0.5") != std::string::npos,
                  "pwi starting magnetization");
            check(text.find("Fe1 ") != std::string::npos && text.find("Fe2 ") != std::string::npos, "pwi numbered species labels");
            Structure parsed;
            check(readSimulationFormat(text, "pwi", parsed, error) && parsed.atoms[2].symbol == "Fe", "pwi labels read as elements: " + error);
        }
        // Periodicity in XSF: slabs and polymers.
        Structure slab = sample();
        slab.pbc = {true, true, false};
        check(roundTrip(slab, "xsf").pbc == slab.pbc, "xsf slab periodicity");
        Structure wire = sample();
        wire.pbc = {true, false, false};
        check(roundTrip(wire, "xsf").pbc == wire.pbc, "xsf polymer periodicity");
        // Molecules: gen (C) and xsf (ATOMS) without a cell; pwi, cell, cp2k refuse.
        Structure molecule;
        molecule.atoms = {site("O", 8, 0, 0, 0.1197), site("H", 1, 0, 0.7572, -0.4786), site("H", 1, 0, -0.7572, -0.4786)};
        for (const std::string format : {"gen", "xsf"}) {
            const Structure back = roundTrip(molecule, format);
            check(!back.hasUnitCell && back.atoms.size() == 3 && std::abs(back.atoms[1].y - 0.7572) < 1e-9, format + ": molecule");
        }
        for (const std::string format : {"pwi", "cell", "cp2k"}) {
            std::ostringstream out;
            std::string error;
            check(!writeSimulationFormat(molecule, format, out, error) && error.find("unit cell") != std::string::npos, format + ": needs a cell");
        }
        // CP2K: cell, periodicity and coordinates.
        {
            std::ostringstream out;
            std::string error;
            check(writeSimulationFormat(slab, "cp2k", out, error), "cp2k write");
            const std::string text = out.str();
            check(text.find("PERIODIC XY\n") != std::string::npos && text.find("&KIND Mo") != std::string::npos &&
                  text.find("  &COORD") != std::string::npos, "cp2k sections");
        }
        // pw.x input: Bohr cell, crystal coordinates, celldm alat, Fortran exponents.
        {
            const std::string text =
                "&SYSTEM\n ibrav = 0, celldm(1) = 6.8306, nat = 2, ntyp = 1\n/\n"
                "CELL_PARAMETERS alat\n 0.5 0.5 0.0\n 0.0 0.5 0.5\n 0.5 0.0 0.5\n"
                "ATOMIC_POSITIONS {crystal}\n Cu 0.0 0.0 0.0\n Cu 0.25d0 0.25d0 0.25d0 0 1 1\n";
            Structure cu;
            std::string error;
            check(readSimulationFormat(text, "pwi", cu, error), "pwi alat read: " + error);
            const double alat = 6.8306 * 0.529177210903;
            check(std::abs(cu.cellVectors[0][0] - 0.5 * alat) < 1e-9, "pwi alat cell");
            check(std::abs(cu.atoms[1].x - 0.25 * (0.5 + 0.0 + 0.5) * alat) < 1e-9, "pwi crystal coordinates");
            check(cu.atoms[1].fixed == 1, "pwi if_pos 0 fixes x");
            Structure other;
            check(!readSimulationFormat("&SYSTEM\n ibrav = 2\n/\n", "pwi", other, error) && error.find("ibrav") != std::string::npos,
                  "pwi refuses other ibrav");
        }
        // CASTEP: lattice lengths and angles, fractional positions, Bohr units.
        {
            const std::string text =
                "%block lattice_abc\nang\n3 4 5\n90 90 120\n%endblock lattice_abc\n"
                "%BLOCK POSITIONS_FRAC\nSi 0 0 0\nSi:tag 0.5 0.5 0.5\n%ENDBLOCK POSITIONS_FRAC\n";
            Structure si;
            std::string error;
            check(readSimulationFormat(text, "cell", si, error), "cell abc read: " + error);
            check(std::abs(si.cellVectors[1][0] + 2.0) < 1e-9 && std::abs(si.cellVectors[1][1] - 4 * std::sqrt(3.0) / 2) < 1e-9, "cell gamma 120");
            check(std::abs(si.atoms[1].z - 2.5) < 1e-9 && si.atoms[1].symbol == "Si", "cell fractional and labels");
            const std::string bohr = "%BLOCK LATTICE_CART\nbohr\n10 0 0\n0 10 0\n0 0 10\n%ENDBLOCK LATTICE_CART\n"
                                     "%BLOCK POSITIONS_ABS\nbohr\nC 1 0 0\n%ENDBLOCK POSITIONS_ABS\n";
            check(readSimulationFormat(bohr, "cell", si, error), "cell bohr read: " + error);
            check(std::abs(si.cellVectors[0][0] - 5.29177210903) < 1e-9 && std::abs(si.atoms[0].x - 0.529177210903) < 1e-12, "cell bohr units");
        }
        // DFTB+: fractional (F) geometries and the origin.
        {
            const std::string text = "2 F\n Ga As\n 1 1 0.0 0.0 0.0\n 2 2 0.25 0.25 0.25\n 0.0 0.0 0.0\n 0.0 2.8 2.8\n 2.8 0.0 2.8\n 2.8 2.8 0.0\n";
            Structure gaas;
            std::string error;
            check(readSimulationFormat(text, "gen", gaas, error), "gen F read: " + error);
            check(std::abs(gaas.atoms[1].x - 1.4) < 1e-12 && gaas.atoms[1].symbol == "As", "gen fractional");
            check(!readSimulationFormat("1 C\n Cu\n 1 3 0 0 0\n", "gen", gaas, error), "gen species index checked");
        }
        // XSF with atomic numbers or symbols; unreadable files report.
        {
            Structure x;
            std::string error;
            check(readSimulationFormat("ATOMS\n 6 0 0 0\n O 0 0 1.128\n", "xsf", x, error), "xsf atoms read: " + error);
            check(x.atoms.size() == 2 && x.atoms[0].symbol == "C" && x.atoms[1].symbol == "O", "xsf atomic numbers and symbols");
            check(!readSimulationFormat("CRYSTAL\nPRIMVEC\n1 0 0\n", "xsf", x, error) && !error.empty(), "truncated xsf reports");
        }
    } catch (const std::exception& e) {
        std::cerr << "simulation_formats failed: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "simulation_formats: all checks passed" << std::endl;
    return 0;
}
