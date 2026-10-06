#include "science/LammpsExport.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
// Box parameters of a cell in the LAMMPS frame (a along x, b in the xy plane).
Mat3 lammpsBox(const Mat3& c)
{
    const double lx = norm(c[0]);
    const Vec3 ax = scale(c[0], 1 / lx);
    const double xy = dot(c[1], ax);
    const double ly = norm(cross(ax, c[1]));
    const double xz = dot(c[2], ax);
    const double yz = (dot(c[1], c[2]) - xy * xz) / ly;
    const double lz = std::sqrt(std::max(0.0, dot(c[2], c[2]) - xz * xz - yz * yz));
    return {{{lx, 0, 0}, {xy, ly, 0}, {xz, yz, lz}}};
}

// The same lattice with tilts reduced by lattice translations (original frame).
Mat3 reduceTilts(Mat3 cell)
{
    for (int pass = 0; pass < 16; ++pass) {
        bool changed = false;
        Mat3 b = lammpsBox(cell);
        if (std::abs(b[2][1]) > 0.5 * b[1][1] + 1e-9) { cell[2] = sub(cell[2], scale(cell[1], std::round(b[2][1] / b[1][1]))); changed = true; b = lammpsBox(cell); }
        if (std::abs(b[2][0]) > 0.5 * b[0][0] + 1e-9) { cell[2] = sub(cell[2], scale(cell[0], std::round(b[2][0] / b[0][0]))); changed = true; b = lammpsBox(cell); }
        if (std::abs(b[1][0]) > 0.5 * b[0][0] + 1e-9) { cell[1] = sub(cell[1], scale(cell[0], std::round(b[1][0] / b[0][0]))); changed = true; }
        if (!changed) break;
    }
    return cell;
}
}

Mat3 lammpsCell(const Mat3& input)
{
    if (determinant(input) <= 0) throw std::runtime_error("LAMMPS needs a right-handed cell; reorder or negate a cell vector");
    return lammpsBox(reduceTilts(input));
}

namespace
{
std::string fileName(const std::string& path) { return std::filesystem::u8path(path).filename().u8string(); }
}

ToolOutput lammpsExport(const Parameters& p)
{
    const StructureInput& input = p.structure("structure");
    const Structure& structure = input.structure;
    if (structure.atoms.empty()) throw std::runtime_error("The structure has no atoms");
    const bool periodic = structure.hasUnitCell && input.pbc[0] && input.pbc[1] && input.pbc[2];
    std::vector<std::string> species;
    for (const auto& atom : structure.atoms)
        if (std::find(species.begin(), species.end(), atom.symbol) == species.end()) species.push_back(atom.symbol);
    // Box and coordinates in the LAMMPS frame.
    Mat3 box{};
    std::vector<Vec3> positions;
    Vec3 low{0, 0, 0};
    Mat3 reduced{};  // lattice basis of the LAMMPS box, in the original frame
    if (periodic) {
        Mat3 cell{};
        for (int r = 0; r < 3; ++r) cell[r] = structure.cellVectors[r];
        box = lammpsCell(cell);
        reduced = reduceTilts(cell);
        for (const auto& atom : structure.atoms) {
            Vec3 f = fractional({atom.x, atom.y, atom.z}, reduced);
            for (double& v : f) v -= std::floor(v);
            positions.push_back(rowTimes(f, box));
        }
    } else {
        Vec3 high = {structure.atoms[0].x, structure.atoms[0].y, structure.atoms[0].z};
        low = high;
        for (const auto& atom : structure.atoms) {
            const Vec3 x{atom.x, atom.y, atom.z};
            for (int k = 0; k < 3; ++k) { low[k] = std::min(low[k], x[k]); high[k] = std::max(high[k], x[k]); }
        }
        const double pad = 10.0;
        for (int k = 0; k < 3; ++k) { low[k] -= pad; high[k] += pad; }
        box = {{{high[0] - low[0], 0, 0}, {0, high[1] - low[1], 0}, {0, 0, high[2] - low[2]}}};
        for (const auto& atom : structure.atoms) positions.push_back(sub({atom.x, atom.y, atom.z}, low));
        low = {0, 0, 0};
    }
    // Potential.
    std::string elementList;
    for (const auto& s : species) elementList += " " + s;
    std::ostringstream pair;
    std::vector<std::string> notes;
    std::map<std::string, double> charges;   // atom_style charge when non-empty
    std::map<std::string, std::string> extraFiles;
    const Json potential = p.has("calculator") ? p.json("calculator") : Json::parse(R"({"potential": "EMT"})");
    std::string name = potential.contains("potential") ? potential.at("potential").string() : "EMT";
    std::string key;
    for (char c : name) if (std::isalnum(static_cast<unsigned char>(c))) key += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (key == "eam" || key == "eamalloy" || key == "eamfs" || key == "finnissinclair") {
        const std::string file = fileName(potential.at("file").string());
        std::string format = potential.contains("format") ? potential.at("format").string() : "auto";
        if (format == "auto") {
            std::string lower;
            for (char c : file) lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            format = lower.find(".fs") != std::string::npos ? "fs" : (lower.find("alloy") != std::string::npos || lower.find("setfl") != std::string::npos ? "setfl" : "funcfl");
        }
        if (format == "funcfl") pair << "pair_style eam\npair_coeff 1 1 " << file << "\n";
        else pair << "pair_style eam/" << (format == "fs" ? "fs" : "alloy") << "\npair_coeff * * " << file << elementList << "\n";
        notes.push_back("Copy " + file + " beside in.lammps.");
    } else if (key == "lennardjones" || key == "lj") {
        const double sigma = potential.contains("sigma") ? potential.at("sigma").number() : 1.0;
        const double epsilon = potential.contains("epsilon") ? potential.at("epsilon").number() : 1.0;
        const double cutoff = potential.contains("cutoff") ? potential.at("cutoff").number() : 3 * sigma;
        pair << "pair_style lj/cut " << cutoff << "\npair_coeff * * " << epsilon << ' ' << sigma << "\npair_modify shift yes\n";
    } else if (key == "tersoff" || key == "stillingerweber" || key == "sw") {
        const bool tersoff = key == "tersoff";
        std::string file = potential.contains("file") && !potential.at("file").string().empty() ? fileName(potential.at("file").string()) : "";
        if (file.empty()) {
            // The built-in Si parameters, in LAMMPS's own file format.
            file = tersoff ? "Si.tersoff" : "Si.sw";
            extraFiles[file] = tersoff
                ? "# Si, Tersoff PRB 38 9902 (1988)\nSi Si Si 3.0 1.0 0.0 1.0039e5 16.217 -0.59825 0.78734 1.1e-6 1.7322 471.18 2.85 0.15 2.4799 1830.8\n"
                : "# Si, Stillinger and Weber PRB 31 5262 (1985)\nSi Si Si 2.1683 2.0951 1.80 21.0 1.20 -0.333333333333 7.049556277 0.6022245584 4.0 0.0 0.0\n";
        } else notes.push_back("Copy " + file + " beside in.lammps.");
        pair << "pair_style " << (tersoff ? "tersoff" : "sw") << "\npair_coeff * * " << file << elementList << "\n";
    } else if (key == "buckingham" || key == "buck") {
        const double cutoff = potential.contains("cutoff") ? potential.at("cutoff").number() : 10.0;
        if (const Json* table = potential.find("charges"); table && table->isObject())
            for (const auto& [element, q] : table->members()) charges[element] = q.number();
        pair << (charges.empty() ? "pair_style buck " : "pair_style buck/coul/long ") << cutoff << "\n";
        // Unlisted pairs interact only through their charges.
        for (std::size_t a = 0; a < species.size(); ++a)
            for (std::size_t b = a; b < species.size(); ++b) {
                double A = 0, rho = 1, C = 0;
                if (const Json* list = potential.find("pairs"); list && list->isArray())
                    for (const auto& entry : list->items()) {
                        const auto& e = entry.at("elements").items();
                        if ((e[0].string() == species[a] && e[1].string() == species[b]) || (e[0].string() == species[b] && e[1].string() == species[a])) {
                            A = entry.at("A").number();
                            rho = entry.at("rho").number();
                            C = entry.contains("C") ? entry.at("C").number() : 0.0;
                        }
                    }
                pair << "pair_coeff " << a + 1 << ' ' << b + 1 << ' ' << A << ' ' << rho << ' ' << C << "  # " << species[a] << '-' << species[b] << "\n";
            }
        if (!charges.empty()) {
            if (!periodic) throw std::runtime_error("LAMMPS export of Coulomb interactions needs a periodic structure");
            pair << "kspace_style ewald " << (potential.contains("ewald_accuracy") ? potential.at("ewald_accuracy").number() : 1e-6) << "\n";
        }
    } else {
        pair << "# EMT has no LAMMPS pair style: choose a potential for" << elementList << ", e.g.\n# pair_style eam/alloy\n# pair_coeff * * FILE.eam.alloy" << elementList << "\n";
        notes.push_back("EMT is not available in LAMMPS; the pair style is left as a commented placeholder.");
    }
    std::ostringstream data;
    data << std::setprecision(12);
    data << "LAMMPS data file written by AtomForge\n\n" << structure.atoms.size() << " atoms\n" << species.size() << " atom types\n\n";
    data << "0.0 " << box[0][0] << " xlo xhi\n0.0 " << box[1][1] << " ylo yhi\n0.0 " << box[2][2] << " zlo zhi\n";
    if (periodic && (std::abs(box[1][0]) > 1e-12 || std::abs(box[2][0]) > 1e-12 || std::abs(box[2][1]) > 1e-12))
        data << box[1][0] << ' ' << box[2][0] << ' ' << box[2][1] << " xy xz yz\n";
    data << "\nMasses\n\n";
    for (std::size_t t = 0; t < species.size(); ++t) data << t + 1 << ' ' << atomicMass(species[t]) << "  # " << species[t] << "\n";
    data << (charges.empty() ? "\nAtoms  # atomic\n\n" : "\nAtoms  # charge\n\n");
    for (std::size_t i = 0; i < structure.atoms.size(); ++i) {
        const std::size_t type = static_cast<std::size_t>(std::find(species.begin(), species.end(), structure.atoms[i].symbol) - species.begin()) + 1;
        data << i + 1 << ' ' << type << ' ';
        if (!charges.empty()) {
            const auto charge = charges.find(structure.atoms[i].symbol);
            if (charge == charges.end()) throw std::runtime_error("No charge given for " + structure.atoms[i].symbol);
            data << charge->second << ' ';
        }
        data << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << "\n";
    }

    const std::string task = p.has("task") ? p.json("task").string() : "minimize";
    const double temperature = p.number("temperature_K", 300.0);
    const double pressureBar = p.number("pressure_GPa", 0.0) * 10000.0;
    const double timestepPs = positive(p.number("timestep_fs", 1.0), "timestep_fs") / 1000.0;
    const long long steps = integer(p.number("steps", 10000), "steps");
    const long long interval = integer(p.number("sample_interval", 100), "sample_interval");
    const double thermostat = positive(p.number("thermostat_fs", 100.0), "thermostat_fs") / 1000.0;
    const double barostat = positive(p.number("barostat_fs", 1000.0), "barostat_fs") / 1000.0;
    const long long seed = integer(p.number("seed", 12345), "seed");
    std::ostringstream in;
    in << std::setprecision(10);
    in << "# LAMMPS input written by AtomForge (" << task << ")\nunits metal\natom_style " << (charges.empty() ? "atomic" : "charge") << "\nboundary " << (periodic ? "p p p" : "f f f") << "\n";
    if (task == "neb") in << "atom_modify map array\n";
    in << "read_data data.lammps\n\n" << pair.str() << "\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nthermo " << interval << "\n";
    in << "thermo_style custom step temp pe etotal press vol\n";
    // Frames with element names are read back by AtomForge's trajectory playback.
    const std::string dump = "dump traj all custom " + std::to_string(interval) + " dump.lammpstrj id element x y z\ndump_modify traj element" + elementList + " sort id\n";
    Json files = Json::object();
    if (task == "minimize") {
        in << dump << "min_style cg\n";
        if (p.boolean("relax_cell", false)) {
            if (!periodic) throw std::runtime_error("Cell relaxation needs a periodic structure");
            in << "fix relax all box/relax tri " << pressureBar << " vmax 0.001\n";
        }
        in << "minimize 0.0 " << p.number("fmax", 0.01) << ' ' << steps << ' ' << 10 * steps << "\nwrite_data relaxed.data\n";
    } else if (task == "nvt" || task == "npt") {
        in << "velocity all create " << temperature << ' ' << seed << " dist gaussian\ntimestep " << timestepPs << "\n" << dump;
        if (task == "nvt") in << "fix md all nvt temp " << temperature << ' ' << temperature << ' ' << thermostat << "\n";
        else {
            if (!periodic) throw std::runtime_error("NPT needs a periodic structure");
            in << "fix md all npt temp " << temperature << ' ' << temperature << ' ' << thermostat << " iso " << pressureBar << ' ' << pressureBar << ' ' << barostat << "\n";
        }
        in << "run " << steps << "\n";
    } else if (task == "neb") {
        if (!p.has("final")) throw std::runtime_error("NEB export needs the final structure");
        const Structure& final = p.structure("final").structure;
        if (final.atoms.size() != structure.atoms.size()) throw std::runtime_error("NEB endpoints need the same atoms");
        std::ostringstream neb;
        neb << std::setprecision(12) << final.atoms.size() << "\n";
        // Final coordinates in the same (LAMMPS) frame: apply the initial structure's mapping.
        for (std::size_t i = 0; i < final.atoms.size(); ++i) {
            const auto& a = final.atoms[i];
            const auto& s0 = structure.atoms[i];
            // Displacements rotate with the cell into the LAMMPS frame.
            const Vec3 shift = sub({a.x, a.y, a.z}, {s0.x, s0.y, s0.z});
            const Vec3 x = periodic ? add(positions[i], rowTimes(rowTimes(shift, inverse(reduced)), box)) : add(positions[i], shift);
            neb << i + 1 << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << "\n";
        }
        files["final.neb"] = neb.str();
        const long long images = integer(p.number("images", 7), "images", 3);
        in << "variable u uloop " << images << "\nfix neb all neb " << p.number("spring_eV_per_A2", 0.1) << "\nmin_style quickmin\ntimestep 0.01\n";
        in << "neb 0.0 " << p.number("fmax", 0.03) << ' ' << steps << ' ' << steps << ' ' << interval << " final final.neb\n";
        notes.push_back("Run NEB with one partition per image, e.g. mpirun -np " + std::to_string(images) + " lmp -partition " + std::to_string(images) + "x1 -in in.lammps");
    } else throw std::runtime_error("task must be minimize, nvt, npt or neb");
    for (const auto& [file, text] : extraFiles) files[file] = text;
    files["data.lammps"] = data.str();
    files["in.lammps"] = in.str();
    Json result = Json::object();
    result["task"] = task;
    result["atom_types"] = [&] { Json list = Json::array(); for (const auto& s : species) list.push(s); return list; }();
    result["box"] = toJson(box);
    result["triclinic"] = periodic && (std::abs(box[1][0]) > 1e-12 || std::abs(box[2][0]) > 1e-12 || std::abs(box[2][1]) > 1e-12);
    result["notes"] = [&] { Json list = Json::array(); for (const auto& n : notes) list.push(n); return list; }();
    result["files"] = files;
    ToolOutput output;
    output.result = result;
    return output;
}
}
