#include "electronic/Lobster.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
void close(double actual, double expected, double tolerance = 1e-9)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error("Expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}
template<class Work> void rejects(Work work)
{
    try { work(); } catch (const std::exception&) { return; }
    throw std::runtime_error("Invalid input was accepted");
}
std::filesystem::path writeFile(const std::string& name, const std::string& content)
{
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream(path) << content;
    return path;
}
}

int main()
{
    try
    {
        // Non-spin-polarized COHPCAR: 2 bonds (+ average), 3 energy points.
        // Column layout (verified against pymatgen.io.lobster.outputs.Cohpcar):
        // energy, avgCOHP, avgICOHP, bond1COHP, bond1ICOHP, bond2COHP, bond2ICOHP.
        {
            const std::string text =
                "created by AtomForge test fixture\n"
                "3 1 -10.0 10.0 -2.5\n"
                "No.1  COHP  ICOHP\n"
                "No.1:Fe1->Fe9(2.45)\n"
                "No.2:Fe1->O12(1.9)\n"
                "-1.0 0.1 0.01 0.2 0.02 0.3 0.03\n"
                "0.0 0.4 0.04 0.5 0.05 0.6 0.06\n"
                "1.0 0.7 0.07 0.8 0.08 0.9 0.09\n";
            const auto data = readCohpcar(writeFile("atomforge_cohpcar_fixture", text).string());
            if (data.spinPolarized) throw std::runtime_error("Expected non-spin-polarized COHPCAR");
            close(data.fermiEnergy, -2.5);
            if (data.bonds.size() != 2) throw std::runtime_error("Expected 2 bonds");
            if (data.bonds[0].atom1 != 0 || data.bonds[0].atom2 != 8) throw std::runtime_error("Bad Fe1->Fe9 site indices");
            close(data.bonds[0].length, 2.45);
            if (data.bonds[1].atom1 != 0 || data.bonds[1].atom2 != 11) throw std::runtime_error("Bad Fe1->O12 site indices");
            close(data.energies[0], -1.0); close(data.energies[2], 1.0);
            close(data.averageCohp[0][1], 0.4); close(data.averageIcohp[0][1], 0.04);
            close(data.cohp[0][0][2], 0.8); close(data.icohp[0][0][2], 0.08);
            close(data.cohp[0][0][0], 0.2); close(data.icohp[0][0][0], 0.02);
            close(data.cohp[0][1][0], 0.3); close(data.icohp[0][1][0], 0.03);
            rejects([&] { readCohpcar((std::filesystem::temp_directory_path() / "missing_cohpcar").string()); });
        }
        // Spin-polarized COHPCAR: 1 bond (+ average), 2 energy points; the
        // up-block (2*(numBonds+1) columns) is followed by the down-block.
        {
            const std::string text =
                "created by AtomForge test fixture\n"
                "2 2 -10.0 10.0 0.0\n"
                "No.1  COHP  ICOHP\n"
                "No.1:H1->H2(1.0)\n"
                "-1.0 1.0 0.1 1.5 0.15 2.0 0.2 2.5 0.25\n"
                "1.0 3.0 0.3 3.5 0.35 4.0 0.4 4.5 0.45\n";
            const auto data = readCohpcar(writeFile("atomforge_cohpcar_spin_fixture", text).string());
            if (!data.spinPolarized) throw std::runtime_error("Expected spin-polarized COHPCAR");
            close(data.averageCohp[0][0], 1.0); close(data.averageCohp[1][0], 2.0);
            close(data.cohp[0][0][1], 3.5); close(data.cohp[1][0][1], 4.5);
        }
        // Real LOBSTER output commonly interleaves each bond's total entry
        // with orbital-resolved sub-entries for the same bond (e.g.
        // "No.1:Fe1->Fe2(2.0)" immediately followed by
        // "No.1:Fe1[3d]->Fe2[3d](2.0)"). Those sub-entries must be skipped,
        // not merely ignored in the header list: this checks that the data
        // *columns* belonging to a skipped orbital slot are also correctly
        // skipped, not misread as the next bond's data.
        {
            const std::string text =
                "created by AtomForge test fixture\n5 1 -10.0 10.0 0.0\nNo.1  COHP  ICOHP\n"
                "No.1:Fe1->Fe2(2.0)\nNo.1:Fe1[3d]->Fe2[3d](2.0)\nNo.2:Fe1->Fe3(3.0)\nNo.2:Fe1[3d]->Fe3[3d](3.0)\n"
                "-1.0 0.1 0.01 0.2 0.02 0.99 0.99 0.3 0.03 0.88 0.88\n";
            const auto data = readCohpcar(writeFile("atomforge_cohpcar_orbital_fixture", text).string());
            if (data.bonds.size() != 2) throw std::runtime_error("Orbital sub-entries must be skipped, not counted as bonds");
            if (data.bonds[0].atom1 != 0 || data.bonds[0].atom2 != 1) throw std::runtime_error("Bad bond0 site indices");
            if (data.bonds[1].atom1 != 0 || data.bonds[1].atom2 != 2) throw std::runtime_error("Bad bond1 site indices");
            close(data.cohp[0][0][0], 0.2); close(data.icohp[0][0][0], 0.02);
            close(data.cohp[0][1][0], 0.3); close(data.icohp[0][1][0], 0.03);  // must skip the orbital columns, not read them
        }
        // A COHPCAR with only orbital-resolved headers (no bond totals at
        // all) has nothing usable and must be rejected outright.
        {
            const std::string text =
                "created by AtomForge test fixture\n2 1 -10.0 10.0 0.0\nNo.1  COHP  ICOHP\n"
                "No.1:Fe1[3d]->Fe2[3d](2.0)\n-1.0 0.1 0.01 0.2 0.02\n1.0 0.3 0.03 0.4 0.04\n";
            rejects([&] { readCohpcar(writeFile("atomforge_cohpcar_allorbital_fixture", text).string()); });
        }

        // ICOHPLIST.lobster, LOBSTER >=3.1.1 format (8 columns, with a translation vector).
        {
            const std::string text =
                "COHP#  atom1  atom2  distance  translation  ICOHP(eV)\n"
                "1 Fe1 Fe9 2.45000 0 0 0 -1.23456\n"
                "2 Fe1 O12 1.90000 0 0 0 -0.54321\n"
                "\n";
            const auto list = readIcohplist(writeFile("atomforge_icohplist_311_fixture", text).string());
            if (list.spinPolarized) throw std::runtime_error("Expected non-spin-polarized ICOHPLIST");
            if (list.entries.size() != 2) throw std::runtime_error("Expected 2 ICOHPLIST entries");
            if (list.entries[0].atom1 != 0 || list.entries[0].atom2 != 8) throw std::runtime_error("Bad site indices");
            close(list.entries[0].length, 2.45); close(list.entries[0].icohp[0], -1.23456);
            close(list.entries[1].length, 1.90); close(list.entries[1].icohp[0], -0.54321);
        }
        // ICOHPLIST.lobster, older <=2.2.1 format (6 columns: label atom1 atom2 distance icohp numBonds).
        {
            const std::string text = "COHP#  atom1  atom2  distance  ICOHP(eV)  numBonds\n1 H1 H2 1.00000 -0.75000 1\n\n";
            const auto list = readIcohplist(writeFile("atomforge_icohplist_221_fixture", text).string());
            close(list.entries[0].icohp[0], -0.75); if (list.entries[0].numBonds != 1) throw std::runtime_error("Bad numBonds");
        }
        // Spin-polarized ICOHPLIST: n_bonds up-entries, one repeated header, n_bonds down-entries.
        {
            const std::string text =
                "COHP#  atom1  atom2  distance  translation  ICOHP(eV)\n"
                "1 H1 H2 1.00000 0 0 0 -0.40000\n"
                "COHP#  atom1  atom2  distance  translation  ICOHP(eV)\n"
                "1 H1 H2 1.00000 0 0 0 -0.35000\n"
                "\n";
            const auto list = readIcohplist(writeFile("atomforge_icohplist_spin_fixture", text).string());
            if (!list.spinPolarized) throw std::runtime_error("Expected spin-polarized ICOHPLIST");
            if (list.entries.size() != 1) throw std::runtime_error("Expected 1 bond (up+down folded into one entry)");
            close(list.entries[0].icohp[0], -0.40); close(list.entries[0].icohp[1], -0.35);
        }
        // Real LOBSTER ICOHPLIST output interleaves each bond's total with
        // orbital-resolved sub-rows (second column containing '_', e.g.
        // "Fe1_3d"); only the totals must be kept.
        {
            const std::string text =
                "COHP#  atom1  atom2  distance  translation  ICOHP(eV)\n"
                "1 Fe1 Fe9 2.45000 0 0 0 -1.23456\n"
                "1 Fe1_3d Fe9_3d 2.45000 0 0 0 -0.99999\n"
                "2 Fe1 O12 1.90000 0 0 0 -0.54321\n"
                "2 Fe1_3d O12_2p 1.90000 0 0 0 -0.11111\n"
                "\n";
            const auto list = readIcohplist(writeFile("atomforge_icohplist_orbital_fixture", text).string());
            if (list.entries.size() != 2) throw std::runtime_error("Orbital-resolved rows must be skipped, not counted as bonds");
            close(list.entries[0].icohp[0], -1.23456); close(list.entries[1].icohp[0], -0.54321);
        }
        rejects([&] { readIcohplist((std::filesystem::temp_directory_path() / "missing_icohplist").string()); });

        std::cout << "LOBSTER analytic regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
