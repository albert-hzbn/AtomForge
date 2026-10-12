#pragma once
// "Did you mean ...?" for a misspelt name (an option, tool, step or
// parameter): the closest known name, if close enough.

#include <algorithm>
#include <string>
#include <vector>

namespace suggest
{
// Edit distance: insertions, deletions, substitutions and swapped
// neighbours ("widht") each count one.
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

// The closest candidate when at most about a third of the name's letters
// differ (ignoring a leading "--"), else "".
template <typename Names>
std::string closest(const std::string& name, const Names& candidates)
{
    const std::size_t letters = name.compare(0, 2, "--") == 0 ? name.size() - 2 : name.size();
    std::size_t best = std::max<std::size_t>(1, letters / 3) + 1;
    std::string found;
    for (const auto& candidate : candidates)
        if (const std::size_t d = distance(name, std::string(candidate)); d < best) { best = d; found = candidate; }
    return found;
}

// " (did you mean X?)" or "".
template <typename Names>
std::string didYouMean(const std::string& name, const Names& candidates)
{
    const std::string found = closest(name, candidates);
    return found.empty() ? std::string() : " (did you mean " + found + "?)";
}

// "a, b, c".
template <typename Names>
std::string list(const Names& names)
{
    std::string text;
    for (const auto& name : names) text += (text.empty() ? "" : ", ") + std::string(name);
    return text;
}
}
