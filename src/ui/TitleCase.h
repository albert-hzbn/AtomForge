#pragma once
// Title case for names shown in menus and window titles ("Powder X-ray
// diffraction" -> "Powder X-ray Diffraction"). Names in data (the science
// catalog, shared with the command line and Python) keep their own case.
#include <cctype>
#include <cstring>
#include <string>

namespace uiText
{
inline std::string titleCase(const char* text)
{
    static const char* const kSmall[] = {"a", "an", "the", "and", "or", "to", "in", "of", "as", "for",
                                         "on", "at", "by", "from", "with", "into", "per", "vs"};
    const std::string input = text ? text : "";
    std::string out;
    out.reserve(input.size());
    std::size_t start = 0;
    bool first = true;
    while (start <= input.size()) {
        const std::size_t end = std::min(input.find(' ', start), input.size());
        std::string word = input.substr(start, end - start);
        const bool last = end >= input.size();
        // Only plain lower-case words change: symbols ("g(r)"), units and
        // acronyms keep their case; the part after a hyphen is left as is.
        bool plain = !word.empty() && std::islower(static_cast<unsigned char>(word[0]));
        for (char c : word) plain = plain && (std::isalpha(static_cast<unsigned char>(c)) || c == '-');
        bool small = false;
        for (const char* s : kSmall) small = small || word == s;
        if (plain && (!small || first || last)) word[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(word[0])));
        out += word;
        if (end < input.size()) out += ' ';
        first = false;
        start = end + 1;
    }
    return out;
}
}
