#pragma once
// Units for display: "e/A^3" -> "e/Å³". Unit strings in data stay ASCII
// (they are compared and written to files); only what is shown changes.
#include <cctype>
#include <string>

namespace unitText
{
inline std::string display(std::string text)
{
    struct Power { const char* ascii; const char* shown; };
    static const Power kPowers[] = {{"^-3", "\xe2\x81\xbb\xc2\xb3"}, {"^-2", "\xe2\x81\xbb\xc2\xb2"}, {"^-1", "\xe2\x81\xbb\xc2\xb9"},
                                    {"^3", "\xc2\xb3"}, {"^2", "\xc2\xb2"}};
    std::string out;
    out.reserve(text.size() + 8);
    for (std::size_t i = 0; i < text.size(); ++i) {
        // An Angstrom sign: a lone "A" followed by a power.
        const bool lone = i == 0 || !std::isalnum(static_cast<unsigned char>(text[i - 1]));
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
        out += text[i];
    }
    return out;
}
}
