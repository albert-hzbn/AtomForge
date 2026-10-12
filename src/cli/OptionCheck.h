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

// Edit distance, for "did you mean".
inline std::size_t distance(const std::string& a, const std::string& b)
{
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) row[j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i)
    {
        std::size_t diagonal = row[0];
        row[0] = i;
        for (std::size_t j = 1; j <= b.size(); ++j)
        {
            const std::size_t above = row[j];
            row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] == b[j - 1] ? 0 : 1)});
            diagonal = above;
        }
    }
    return row[b.size()];
}

// "" when every --option in argv is one the mode knows, else the message.
inline std::string unknownOption(const std::vector<std::string>& args, void (*help)(), const std::string& mode)
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
        return "unknown option " + flag + " for " + mode + (suggestion.empty() ? "" : " (did you mean " + suggestion + "?)") +
               "; see AtomForge --help " + mode;
    }
    return {};
}
}
