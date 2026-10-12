#pragma once
// Units for display: "e/A^3" -> "e/Å³", "cutoff (A)" -> "cutoff (Å)",
// "(Angstrom)" -> "(Å)". Unit strings in data stay ASCII (they are compared
// and written to files); only what is shown changes.
#include <cctype>
#include <cstring>
#include <string>

namespace unitText
{
inline std::string display(std::string text)
{
    struct Power { const char* ascii; const char* shown; };
    static const Power kPowers[] = {{"^-3", "\xe2\x81\xbb\xc2\xb3"}, {"^-2", "\xe2\x81\xbb\xc2\xb2"}, {"^-1", "\xe2\x81\xbb\xc2\xb9"},
                                    {"^3", "\xc2\xb3"}, {"^2", "\xc2\xb2"}};
    const auto alnum = [](char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0; };
    std::string out;
    out.reserve(text.size() + 8);
    for (std::size_t i = 0; i < text.size(); ++i) {
        const bool lone = i == 0 || !alnum(text[i - 1]);
        // The word: "Angstrom", "Angstroms", "angstrom".
        for (const char* word : {"Angstroms", "angstroms", "Angstrom", "angstrom"}) {
            const std::size_t length = std::strlen(word);
            if (lone && text.compare(i, length, word) == 0 && (i + length == text.size() || !alnum(text[i + length]))) {
                out += "\xc3\x85";
                i += length;
                break;
            }
        }
        if (i >= text.size()) break;
        // An Angstrom sign: a lone "A" followed by a power.
        if (text[i] == 'A' && lone && i + 1 < text.size() && text[i + 1] == '^') {
            bool matched = false;
            for (const Power& power : kPowers) {
                const std::string ascii = power.ascii;
                if (text.compare(i + 1, ascii.size(), ascii) == 0) {
                    out += "\xc3\x85";
                    out += power.shown;
                    i += ascii.size();
                    matched = true;
                    break;
                }
            }
            if (matched) continue;
        }
        // A unit in brackets or per Angstrom: "(A)", "(A, ...", "1/A".
        if (text[i] == 'A' && i > 0 && (text[i - 1] == '(' || text[i - 1] == '/')
            && (i + 1 == text.size() || text[i + 1] == ')' || text[i + 1] == ',' || text[i + 1] == ' ')) {
            out += "\xc3\x85";
            continue;
        }
        out += text[i];
    }
    return out;
}
}
