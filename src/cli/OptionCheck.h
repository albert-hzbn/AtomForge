#pragma once
// Builders read their options by name, so an unknown one (a typo such as
// --spacegrup) would be ignored and the default used. The options a builder
// knows are those its help text names; anything else is refused.

#include "util/Suggest.h"

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

inline std::size_t distance(const std::string& a, const std::string& b) { return suggest::distance(a, b); }

// "" when every --option in argv is one the mode knows, else the message
// (pointing at helpCommand, by default "AtomForge --help <mode>").
inline std::string unknownOption(const std::vector<std::string>& args, void (*help)(), const std::string& mode, std::string helpCommand = {})
{
    std::set<std::string> known = flagsInHelp(help);
    for (const char* common : {"--build", "--output", "--input", "--help", "--analyze"}) known.insert(common);
    for (const std::string& arg : args)
    {
        // "-q": options start with two dashes (a value such as -0.5 is a number).
        if (arg.size() >= 2 && arg[0] == '-' && arg[1] != '-' && std::isalpha(static_cast<unsigned char>(arg[1])) && arg != "-h")
        {
            if (helpCommand.empty()) helpCommand = "AtomForge --help " + mode;
            const std::string suggestion = arg.size() > 2 ? suggest::closest("-" + arg, known) : std::string();
            return "unknown option " + arg + " for " + mode + " (options start with --" +
                   (suggestion.empty() ? std::string() : "; did you mean " + suggestion + "?") + "); see " + helpCommand;
        }
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0 || !std::isalpha(static_cast<unsigned char>(arg[2]))) continue;
        const std::string flag = arg.substr(0, arg.find('='));
        if (known.count(flag)) continue;
        // Suggest only close matches: a third of the name's letters may differ.
        const std::string suggestion = suggest::closest(flag, known);
        if (helpCommand.empty()) helpCommand = "AtomForge --help " + mode;
        return "unknown option " + flag + " for " + mode + (suggestion.empty() ? "" : " (did you mean " + suggestion + "?)") +
               "; see " + helpCommand;
    }
    return {};
}
}
