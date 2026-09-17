#include "electronic/Lobster.h"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
std::vector<std::string> readAllLines(const std::string& path)
{
    std::ifstream in{std::filesystem::u8path(path)};
    if (!in) throw std::invalid_argument("Cannot open LOBSTER output file");
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(std::move(line));
    }
    return lines;
}
std::vector<std::string> words(const std::string& s)
{
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string token; in >> token;) out.push_back(token);
    return out;
}
bool blank(const std::string& s) { return words(s).empty(); }

// "Fe12" -> 11 (zero-based site index); LOBSTER numbers sites to match the
// structure/POSCAR order, so this index is directly usable against it.
int siteIndexFromLabel(const std::string& label)
{
    std::size_t digitsStart = label.size();
    while (digitsStart > 0 && std::isdigit(static_cast<unsigned char>(label[digitsStart - 1]))) --digitsStart;
    if (digitsStart == label.size()) throw std::invalid_argument("LOBSTER site label has no trailing index: " + label);
    return std::stoi(label.substr(digitsStart)) - 1;
}

std::vector<double> toDoubles(const std::vector<std::string>& tokens)
{
    std::vector<double> values(tokens.size());
    for (std::size_t i = 0; i < tokens.size(); ++i) values[i] = std::stod(tokens[i]);
    return values;
}
}

CohpData readCohpcar(const std::string& path)
{
    const auto lines = readAllLines(path);
    if (lines.size() < 4) throw std::invalid_argument("Truncated COHPCAR file");
    const auto parameters = words(lines[1]);
    if (parameters.size() < 3) throw std::invalid_argument("Invalid COHPCAR parameter line");
    const int totalBonds = std::stoi(parameters[0]);  // includes the "average" entry
    const int spinFlag = std::stoi(parameters[1]);
    if (totalBonds < 1 || (spinFlag != 1 && spinFlag != 2))
        throw std::invalid_argument("Invalid COHPCAR bond count or spin flag");
    const int numBonds = totalBonds - 1;
    const int numSpins = spinFlag;

    CohpData result;
    result.spinPolarized = numSpins == 2;
    result.fermiEnergy = std::stod(parameters.back());
    // lines[2] is a column-label line, not used; bond headers follow at lines[3..3+numBonds-1].
    // Real LOBSTER output commonly interleaves each bond's total entry with
    // its orbital-resolved contributions ("No.7:Fe1[3d_xy]->Fe2[3d_xy](...)"
    // immediately after "No.7:Fe1->Fe2(...)"); only the totals are kept
    // here, identified by the absence of '[' in the header. rawIndex records
    // each kept bond's position among *all* header lines (including
    // orbital ones), since that raw position -- not the stored bond count
    // -- is what the data-column offsets below are keyed on.
    if (static_cast<int>(lines.size()) < 3 + numBonds) throw std::invalid_argument("Truncated COHPCAR bond header block");
    std::vector<int> rawIndex;
    for (int bond = 0; bond < numBonds; ++bond)
    {
        const std::string& header = lines[3 + bond];
        if (header.find('[') != std::string::npos) continue;
        const auto openParen = header.rfind('(');
        if (openParen == std::string::npos || header.back() != ')')
            throw std::invalid_argument("Invalid COHPCAR bond header: " + header);
        CohpBond info;
        info.length = std::stod(header.substr(openParen + 1, header.size() - openParen - 2));
        std::string prefix = header.substr(0, openParen);
        for (auto pos = prefix.find("->"); pos != std::string::npos; pos = prefix.find("->"))
            prefix.replace(pos, 2, ":");
        std::vector<std::string> parts;
        std::istringstream in(prefix);
        for (std::string part; std::getline(in, part, ':');) parts.push_back(part);
        if (parts.size() < 3) throw std::invalid_argument("Invalid COHPCAR bond header: " + header);
        info.atom1 = siteIndexFromLabel(parts[1]);
        info.atom2 = siteIndexFromLabel(parts[2]);
        result.bonds.push_back(info);
        rawIndex.push_back(bond);
    }
    if (result.bonds.empty()) throw std::invalid_argument("COHPCAR file contains no bond-total headers");
    const int storedBonds = static_cast<int>(result.bonds.size());

    std::vector<std::vector<double>> rows;
    for (std::size_t i = 3 + static_cast<std::size_t>(numBonds); i < lines.size(); ++i)
    {
        if (blank(lines[i])) continue;
        rows.push_back(toDoubles(words(lines[i])));
    }
    if (rows.empty()) throw std::invalid_argument("COHPCAR file contains no data rows");
    const std::size_t expectedColumns = 1 + 2 * static_cast<std::size_t>(numBonds + 1) * static_cast<std::size_t>(numSpins);
    result.energies.resize(rows.size());
    result.averageCohp.assign(numSpins, std::vector<double>(rows.size()));
    result.averageIcohp.assign(numSpins, std::vector<double>(rows.size()));
    result.cohp.assign(numSpins, std::vector<std::vector<double>>(static_cast<std::size_t>(storedBonds), std::vector<double>(rows.size())));
    result.icohp.assign(numSpins, std::vector<std::vector<double>>(static_cast<std::size_t>(storedBonds), std::vector<double>(rows.size())));
    for (std::size_t r = 0; r < rows.size(); ++r)
    {
        const auto& row = rows[r];
        if (row.size() != expectedColumns) throw std::invalid_argument("Inconsistent COHPCAR data row width");
        result.energies[r] = row[0];
        for (int s = 0; s < numSpins; ++s)
        {
            result.averageCohp[s][r] = row[1 + 2 * s * (numBonds + 1)];
            result.averageIcohp[s][r] = row[2 + 2 * s * (numBonds + 1)];
            for (int b = 0; b < storedBonds; ++b)
            {
                result.cohp[s][static_cast<std::size_t>(b)][r] = row[2 * (rawIndex[static_cast<std::size_t>(b)] + s * (numBonds + 1)) + 3];
                result.icohp[s][static_cast<std::size_t>(b)][r] = row[2 * (rawIndex[static_cast<std::size_t>(b)] + s * (numBonds + 1)) + 4];
            }
        }
    }
    return result;
}

IcohpList readIcohplist(const std::string& path)
{
    auto lines = readAllLines(path);
    if (lines.size() < 2) throw std::invalid_argument("Truncated ICOHPLIST file");
    lines.erase(lines.begin());  // drop the column-label header line
    while (!lines.empty() && blank(lines.back())) lines.pop_back();
    if (lines.empty()) throw std::invalid_argument("ICOHPLIST file contains no data");

    const bool spinPolarized = lines[lines.size() / 2].find("distance") != std::string::npos;

    // Real LOBSTER output commonly interleaves each bond's total ICOHP with
    // its orbital-resolved contributions (second column containing '_', e.g.
    // "Na1_3s"); only the totals are kept (matching pymatgen's Icohplist
    // "data_without_orbitals" filtering). The repeated spin-down header line
    // (present when spinPolarized) has no '_' in its second token either, so
    // it naturally survives this filter, keeping the up/down split intact.
    std::vector<std::string> totals;
    for (const auto& line : lines)
    {
        const auto tokens = words(line);
        if (tokens.size() > 1 && tokens[1].find('_') != std::string::npos) continue;
        totals.push_back(line);
    }
    if (totals.empty()) throw std::invalid_argument("ICOHPLIST file contains no bond-total rows");

    const auto firstTokens = words(totals[0]);
    int version;  // 8 columns: LOBSTER >=3.1.1 (with a translation vector); 6 columns: LOBSTER <=2.2.1
    if (firstTokens.size() == 8) version = 311;
    else if (firstTokens.size() == 6) version = 221;
    else throw std::invalid_argument("Unrecognized ICOHPLIST column count");

    const bool totalsSpinPolarized = totals[totals.size() / 2].find("distance") != std::string::npos;
    const std::size_t numBonds = totalsSpinPolarized ? totals.size() / 2 : totals.size();
    if (numBonds == 0) throw std::invalid_argument("ICOHPLIST file contains no data");

    IcohpList result;
    result.spinPolarized = spinPolarized;
    for (std::size_t bond = 0; bond < numBonds; ++bond)
    {
        const auto tokens = words(totals[bond]);
        if (static_cast<int>(tokens.size()) != (version == 311 ? 8 : 6))
            throw std::invalid_argument("Inconsistent ICOHPLIST row width");
        IcohpEntry entry;
        entry.atom1 = siteIndexFromLabel(tokens[1]);
        entry.atom2 = siteIndexFromLabel(tokens[2]);
        entry.length = std::stod(tokens[3]);
        if (version == 311) { entry.icohp[0] = std::stod(tokens[7]); entry.numBonds = 1; }
        else { entry.icohp[0] = std::stod(tokens[4]); entry.numBonds = std::stoi(tokens[5]); }
        if (spinPolarized)
        {
            const auto downTokens = words(totals[bond + numBonds + 1]);
            if (static_cast<int>(downTokens.size()) != (version == 311 ? 8 : 6))
                throw std::invalid_argument("Inconsistent ICOHPLIST spin-down row width");
            entry.icohp[1] = std::stod(downTokens[version == 311 ? 7 : 4]);
        }
        result.entries.push_back(entry);
    }
    return result;
}
}
