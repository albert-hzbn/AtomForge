#include "science/VaspElectronic.h"
#include "science/Analysis.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
const double kPi = std::acos(-1.0);

std::vector<std::string> tokens(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> result;
    for (std::string token; stream >> token;) result.push_back(token);
    return result;
}

double number(const std::string& token, const char* what)
{
    try { return std::stod(token); }
    catch (const std::exception&) { throw std::runtime_error(std::string("Invalid number in ") + what + ": " + token); }
}

bool isBlank(const std::string& line) { return line.find_first_not_of(" \t\r") == std::string::npos; }
}

Eigenvalues readEigenval(const std::filesystem::path& path)
{
    std::istringstream input(readTextFile(path));
    std::string line;
    std::getline(input, line);
    const auto first = tokens(line);
    if (first.size() < 4) throw std::runtime_error("Invalid EIGENVAL header");
    const int spins = static_cast<int>(number(first[3], "EIGENVAL"));
    if (spins != 1 && spins != 2) throw std::runtime_error("EIGENVAL ISPIN must be 1 or 2");
    for (int k = 0; k < 4; ++k) std::getline(input, line);
    std::getline(input, line);
    const auto sizes = tokens(line);
    if (sizes.size() < 3) throw std::runtime_error("Invalid EIGENVAL dimensions");
    const int kpoints = static_cast<int>(number(sizes[1], "EIGENVAL")), bands = static_cast<int>(number(sizes[2], "EIGENVAL"));
    Eigenvalues result;
    result.energies.assign(static_cast<std::size_t>(spins), std::vector<std::vector<double>>(static_cast<std::size_t>(kpoints), std::vector<double>(static_cast<std::size_t>(bands))));
    for (int k = 0; k < kpoints; ++k) {
        do { if (!std::getline(input, line)) throw std::runtime_error("Truncated EIGENVAL"); } while (isBlank(line));
        const auto kline = tokens(line);
        if (kline.size() < 4) throw std::runtime_error("Invalid EIGENVAL k point");
        result.kpoints.push_back({number(kline[0], "EIGENVAL"), number(kline[1], "EIGENVAL"), number(kline[2], "EIGENVAL")});
        result.weights.push_back(number(kline[3], "EIGENVAL"));
        for (int b = 0; b < bands; ++b) {
            if (!std::getline(input, line)) throw std::runtime_error("Truncated EIGENVAL");
            const auto row = tokens(line);
            if (row.size() < static_cast<std::size_t>(1 + spins)) throw std::runtime_error("Invalid EIGENVAL band row");
            for (int s = 0; s < spins; ++s) result.energies[static_cast<std::size_t>(s)][static_cast<std::size_t>(k)][static_cast<std::size_t>(b)] = number(row[static_cast<std::size_t>(1 + s)], "EIGENVAL");
        }
    }
    return result;
}

DensityOfStates readDoscar(const std::filesystem::path& path)
{
    std::istringstream input(readTextFile(path));
    std::string line;
    std::getline(input, line);
    const auto first = tokens(line);
    if (first.size() < 3) throw std::runtime_error("Invalid DOSCAR header");
    const int atoms = static_cast<int>(number(first[0], "DOSCAR"));
    const bool partial = number(first[2], "DOSCAR") != 0;
    for (int k = 0; k < 4; ++k) std::getline(input, line);
    std::getline(input, line);
    const auto header = tokens(line);
    if (header.size() < 4) throw std::runtime_error("Invalid DOSCAR energy header");
    const int points = static_cast<int>(number(header[2], "DOSCAR"));
    DensityOfStates dos;
    dos.fermi = number(header[3], "DOSCAR");
    int spins = 0;
    for (int i = 0; i < points; ++i) {
        if (!std::getline(input, line)) throw std::runtime_error("Truncated DOSCAR");
        const auto row = tokens(line);
        if (i == 0) {
            spins = row.size() >= 5 ? 2 : 1;
            dos.total.assign(static_cast<std::size_t>(spins), {});
        }
        dos.energy.push_back(number(row[0], "DOSCAR"));
        for (int s = 0; s < spins; ++s) dos.total[static_cast<std::size_t>(s)].push_back(number(row[static_cast<std::size_t>(1 + s)], "DOSCAR"));
    }
    if (partial)
        for (int a = 0; a < atoms; ++a) {
            if (!std::getline(input, line) || isBlank(line)) break;  // header of the site block
            std::vector<std::vector<double>> site(static_cast<std::size_t>(spins), std::vector<double>(static_cast<std::size_t>(points), 0.0));
            for (int i = 0; i < points; ++i) {
                if (!std::getline(input, line)) throw std::runtime_error("Truncated DOSCAR site block");
                const auto row = tokens(line);
                // Spin-polarized orbital columns alternate up/down.
                for (std::size_t c = 1; c < row.size(); ++c)
                    site[static_cast<std::size_t>(spins == 2 ? (c - 1) % 2 : 0)][static_cast<std::size_t>(i)] += number(row[c], "DOSCAR");
            }
            dos.projected.push_back(site);
        }
    return dos;
}

std::vector<std::string> readLineModeLabels(const std::filesystem::path& path, int& pointsPerSegment)
{
    std::istringstream input(readTextFile(path));
    std::string line;
    std::getline(input, line);
    std::getline(input, line);
    const auto count = tokens(line);
    if (count.empty()) throw std::runtime_error("Invalid KPOINTS file");
    pointsPerSegment = static_cast<int>(number(count[0], "KPOINTS"));
    std::getline(input, line);
    if (line.find_first_not_of(" \t") == std::string::npos || std::tolower(static_cast<unsigned char>(line[line.find_first_not_of(" \t")])) != 'l')
        throw std::runtime_error("KPOINTS is not in line mode");
    std::getline(input, line);  // Reciprocal / Cartesian
    std::vector<std::string> labels;
    while (std::getline(input, line)) {
        if (isBlank(line)) continue;
        std::string label;
        const auto bang = line.find('!');
        if (bang != std::string::npos) {
            const auto rest = tokens(line.substr(bang + 1));
            if (!rest.empty()) label = rest[0];
            line = line.substr(0, bang);
        }
        const auto row = tokens(line);
        if (row.size() < 3) continue;
        if (label.empty() && row.size() >= 4) label = row[3];
        if (label == "\\Gamma" || label == "Gamma" || label == "G" || label == "GAMMA") label = "G";
        labels.push_back(label);
    }
    return labels;
}

Projections readProcar(const std::filesystem::path& path)
{
    std::istringstream input(readTextFile(path));
    Projections result;
    std::string line;
    int spin = -1, k = -1, band = -1;
    bool inBlock = false, blockDone = false;
    while (std::getline(input, line)) {
        const auto row = tokens(line);
        if (row.empty()) continue;
        if (line.find("# of k-points") != std::string::npos) {
            ++spin;
            result.weight.emplace_back();
            k = -1;
            continue;
        }
        if (row[0] == "k-point") {
            ++k;
            result.weight[static_cast<std::size_t>(spin)].emplace_back();
            band = -1;
            continue;
        }
        if (row[0] == "band") {
            ++band;
            result.weight[static_cast<std::size_t>(spin)][static_cast<std::size_t>(k)].emplace_back();
            inBlock = false;
            blockDone = false;
            continue;
        }
        if (row[0] == "ion") {
            // The first ion block of a band holds the projections (a later one,
            // with LORBIT = 12, holds phases).
            if (!blockDone) {
                inBlock = true;
                if (result.orbitals.empty()) result.orbitals.assign(row.begin() + 1, row.end() - 1);
            }
            continue;
        }
        if (row[0] == "tot") {
            if (inBlock) blockDone = true;
            inBlock = false;
            continue;
        }
        if (inBlock && spin >= 0 && k >= 0 && band >= 0 && std::isdigit(static_cast<unsigned char>(row[0][0]))) {
            std::vector<double> values;
            for (std::size_t c = 1; c < row.size() - 1 && values.size() < result.orbitals.size(); ++c) values.push_back(number(row[c], "PROCAR"));
            result.weight[static_cast<std::size_t>(spin)][static_cast<std::size_t>(k)][static_cast<std::size_t>(band)].push_back(values);
        }
    }
    if (result.weight.empty() || result.orbitals.empty()) throw std::runtime_error("No projections found in PROCAR");
    return result;
}

ToolOutput vaspElectronic(const Parameters& p)
{
    const Eigenvalues eigen = readEigenval(std::filesystem::u8path(p.json("eigenval_file").string()));
    const std::size_t spins = eigen.energies.size(), kpoints = eigen.kpoints.size();
    const std::size_t bands = eigen.energies[0][0].size();
    Json result = Json::object();
    std::string fermiSource = "none";
    double fermi = 0;
    DensityOfStates dos;
    const bool haveDos = p.has("doscar_file");
    if (haveDos) dos = readDoscar(std::filesystem::u8path(p.json("doscar_file").string()));
    if (p.has("fermi_eV")) { fermi = p.number("fermi_eV"); fermiSource = "input"; }
    else if (haveDos) { fermi = dos.fermi; fermiSource = "DOSCAR"; }
    // Reciprocal lattice (2 pi convention) from the structure, else fractional distances.
    Mat3 reciprocal = identity();
    std::vector<std::string> symbols;
    if (p.has("structure")) {
        const auto& structure = p.structure("structure").structure;
        if (!structure.hasUnitCell) throw std::runtime_error("The structure needs a cell for k-space distances");
        Mat3 cell{};
        for (int r = 0; r < 3; ++r) cell[r] = structure.cellVectors[r];
        const Mat3 inv = inverse(cell);
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) reciprocal[i][j] = 2 * kPi * inv[j][i];
        for (const auto& atom : structure.atoms) symbols.push_back(atom.symbol);
        result["distance_unit"] = "1/Angstrom (2 pi included)";
    } else result["distance_unit"] = "fractional reciprocal coordinates";
    // Path distances, with line-mode labels and path breaks.
    std::vector<std::string> pointLabels(kpoints);
    std::vector<bool> breakBefore(kpoints, false);
    if (p.has("kpoints_file")) {
        int perSegment = 0;
        const auto labels = readLineModeLabels(std::filesystem::u8path(p.json("kpoints_file").string()), perSegment);
        const std::size_t segments = labels.size() / 2;
        if (perSegment > 1 && segments * static_cast<std::size_t>(perSegment) == kpoints)
            for (std::size_t s = 0; s < segments; ++s) {
                const std::size_t start = s * static_cast<std::size_t>(perSegment), end = start + static_cast<std::size_t>(perSegment) - 1;
                pointLabels[start] = labels[2 * s];
                pointLabels[end] = labels[2 * s + 1];
                if (s > 0 && labels[2 * s] != labels[2 * s - 1]) breakBefore[start] = true;
            }
        else result["kpoints_note"] = "KPOINTS segments do not match the EIGENVAL k points; labels omitted";
    }
    std::vector<double> distance(kpoints, 0.0);
    for (std::size_t k = 1; k < kpoints; ++k)
        distance[k] = distance[k - 1] + (breakBefore[k] ? 0.0 : norm(rowTimes(sub(eigen.kpoints[k], eigen.kpoints[k - 1]), reciprocal)));
    Json markers = Json::array();
    for (std::size_t k = 0; k < kpoints; ++k) {
        if (pointLabels[k].empty()) continue;
        std::string label = pointLabels[k];
        if (k + 1 < kpoints && breakBefore[k + 1] && !pointLabels[k + 1].empty()) label += "|" + pointLabels[k + 1];
        if (breakBefore[k]) continue;  // merged into the previous point's label
        // A continuous junction repeats its point at the end and start of two segments.
        if (markers.size() > 0) {
            const Json& last = markers.items().back();
            if (std::abs(last.at("distance").number() - distance[k]) < 1e-12 && last.at("label").string() == label) continue;
        }
        Json marker = Json::object();
        marker["distance"] = distance[k];
        marker["label"] = label;
        markers.push(marker);
    }
    Json energies = Json::array();
    for (std::size_t s = 0; s < spins; ++s) {
        Json spin = Json::array();
        for (std::size_t k = 0; k < kpoints; ++k) {
            std::vector<double> row;
            for (double e : eigen.energies[s][k]) row.push_back(e - fermi);
            spin.push(toJson(row));
        }
        energies.push(spin);
    }
    result["fermi_eV"] = fermi;
    result["fermi_source"] = fermiSource;
    result["spins"] = spins;
    result["kpoints"] = kpoints;
    result["bands"] = bands;
    result["distance"] = toJson(distance);
    result["labels"] = markers;
    result["energies_minus_fermi_eV"] = energies;
    if (fermiSource != "none") {
        try {
            Json request = Json::object();
            request["energies_eV"] = energies;
            request["fermi_eV"] = 0.0;
            result["gap"] = bandGap(Parameters("band-gap", request, "."));
        } catch (const std::exception& error) {
            result["gap"] = Json();
            result["gap_note"] = error.what();
        }
    }
    if (haveDos) {
        Json dosJson = Json::object();
        std::vector<double> shifted;
        for (double e : dos.energy) shifted.push_back(e - fermi);
        dosJson["energy_minus_fermi_eV"] = toJson(shifted);
        Json total = Json::array();
        for (const auto& spin : dos.total) total.push(toJson(spin));
        dosJson["total"] = total;
        if (!dos.projected.empty() && dos.projected.size() == symbols.size()) {
            std::map<std::string, std::vector<std::vector<double>>> byElement;
            for (std::size_t a = 0; a < symbols.size(); ++a) {
                auto& sum = byElement[symbols[a]];
                if (sum.empty()) sum = dos.projected[a];
                else for (std::size_t s = 0; s < sum.size(); ++s) for (std::size_t i = 0; i < sum[s].size(); ++i) sum[s][i] += dos.projected[a][s][i];
            }
            Json projected = Json::object();
            for (const auto& [element, spinsData] : byElement) {
                Json data = Json::array();
                for (const auto& spin : spinsData) data.push(toJson(spin));
                projected[element] = data;
            }
            dosJson["projected_by_element"] = projected;
        }
        result["dos"] = dosJson;
    }
    if (p.has("procar_file")) {
        const Projections projections = readProcar(std::filesystem::u8path(p.json("procar_file").string()));
        if (projections.weight.size() != spins || projections.weight[0].size() != kpoints)
            throw std::runtime_error("PROCAR k points or spins do not match EIGENVAL");
        const std::string orbital = p.has("projection_orbital") ? p.json("projection_orbital").string() : "all";
        std::vector<std::string> elements;
        if (p.has("projection_element")) elements = {p.json("projection_element").string()};
        else { elements = symbols; std::sort(elements.begin(), elements.end()); elements.erase(std::unique(elements.begin(), elements.end()), elements.end()); }
        if (elements.empty()) throw std::runtime_error("Fat bands need the structure (for element names) or projection_element");
        auto orbitalSelected = [&](const std::string& name) {
            return orbital == "all" || name == orbital || (orbital.size() == 1 && !name.empty() && name[0] == orbital[0] && name != "tot");
        };
        Json fat = Json::object();
        for (const auto& element : elements) {
            Json bySpin = Json::array();
            for (std::size_t s = 0; s < spins; ++s) {
                Json byK = Json::array();
                for (std::size_t k = 0; k < kpoints; ++k) {
                    std::vector<double> row;
                    for (std::size_t b = 0; b < projections.weight[s][k].size(); ++b) {
                        const auto& ions = projections.weight[s][k][b];
                        double selected = 0, total = 0;
                        for (std::size_t ion = 0; ion < ions.size(); ++ion)
                            for (std::size_t o = 0; o < ions[ion].size(); ++o) {
                                total += ions[ion][o];
                                if (ion < symbols.size() ? symbols[ion] == element : false)
                                    if (orbitalSelected(projections.orbitals[o])) selected += ions[ion][o];
                            }
                        row.push_back(total > 0 ? selected / total : 0.0);
                    }
                    byK.push(toJson(row));
                }
                bySpin.push(byK);
            }
            fat[element] = bySpin;
        }
        result["fat_band_orbital"] = orbital;
        result["fat_band_weights"] = fat;
        result["procar_orbitals"] = [&] { Json names = Json::array(); for (const auto& o : projections.orbitals) names.push(o); return names; }();
    }
    ToolOutput output;
    output.result = result;
    return output;
}
}
