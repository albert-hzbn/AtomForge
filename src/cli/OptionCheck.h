#pragma once
// Builders read their options by name, so an unknown one (a typo such as
// --spacegrup) would be ignored and the default used. The options a builder
// knows are those its help text names; anything else is refused.

#include <algorithm>
#include <cctype>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace cliOptions
{
// The --flags a help function prints.
inline std::set<std::string> flagsInHelp(void (*help)())
{
    std::ostringstream text;
    std::streambuf* previous = std::cout.rdbuf(text.rdbuf());
    help();
    std::cout.rdbuf(previous);
    const std::string help_text = text.str();
    std::set<std::string> flags;
    for (std::size_t at = help_text.find("--"); at != std::string::npos; at = help_text.find("--", at + 2))
    {
        std::size_t end = at + 2;
        while (end < help_text.size() && (std::isalnum(static_cast<unsigned char>(help_text[end])) || help_text[end] == '-')) ++end;
        if (end > at + 2 && std::isalpha(static_cast<unsigned char>(help_text[at + 2]))) flags.insert(help_text.substr(at, end - at));
    }
    return flags;
}

// Edit distance for "did you mean": insertions, deletions, substitutions
// and swapped neighbours ("--widht") each count one.
inline std::size_t distance(const std::string& a, const std::string& b)
{
    std::vector<std::vector<std::size_t>> d(a.size() + 1, std::vector<std::size_t>(b.size() + 1));
    for (std::size_t i = 0; i <= a.size(); ++i) d[i][0] = i;
    for (std::size_t j = 0; j <= b.size(); ++j) d[0][j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i)
        for (std::size_t j = 1; j <= b.size(); ++j)
        {
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1)});
            if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1])
                d[i][j] = std::min(d[i][j], d[i - 2][j - 2] + 1);
        }
    return d[a.size()][b.size()];
}

// "" when every --option in argv is one the mode knows, else the message
// (pointing at helpCommand, by default "AtomForge --help <mode>").
inline std::string unknownOption(const std::vector<std::string>& args, void (*help)(), const std::string& mode, std::string helpCommand = {})
{
    std::set<std::string> known = flagsInHelp(help);
    for (const char* common : {"--build", "--output", "--input", "--help", "--analyze"}) known.insert(common);
    for (const std::string& arg : args)
    {
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0 || !std::isalpha(static_cast<unsigned char>(arg[2]))) continue;
        const std::string flag = arg.substr(0, arg.find('='));
        if (known.count(flag)) continue;
        std::string suggestion;
        // Suggest only close matches: a third of the name's letters may differ.
        std::size_t best = std::max<std::size_t>(1, (flag.size() - 2) / 3) + 1;
        for (const std::string& candidate : known)
            if (const std::size_t d = distance(flag, candidate); d < best) { best = d; suggestion = candidate; }
        if (helpCommand.empty()) helpCommand = "AtomForge --help " + mode;
        return "unknown option " + flag + " for " + mode + (suggestion.empty() ? "" : " (did you mean " + suggestion + "?)") +
               "; see " + helpCommand;
    }
    return {};
}
}
