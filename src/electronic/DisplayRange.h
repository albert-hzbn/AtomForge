#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace atomforge::electronic
{
struct DisplayRange
{
    double low, high, suggested;
    bool contains(double value) const { return std::isfinite(value) && value>=low && value<=high; }
};

// A display heuristic, not a physical charge-partition threshold. Use zero for
// signed fields; otherwise prefer the 90th percentile over an outlier-led maximum.
inline DisplayRange displayRange(const std::vector<double>& values)
{
    if (values.empty()) throw std::invalid_argument("Cannot estimate an empty field");
    const auto bounds=std::minmax_element(values.begin(),values.end());
    DisplayRange result{*bounds.first,*bounds.second,*bounds.first};
    if (result.low==result.high) return result;
    if (result.low<0 && result.high>0) { result.suggested=0; return result; }
    std::vector<double> sample;
    const std::size_t stride=std::max<std::size_t>(1,(values.size()+262143)/262144);
    for (std::size_t i=0;i<values.size();i+=stride) sample.push_back(values[i]);
    const auto index=(sample.size()-1)*9/10;
    std::nth_element(sample.begin(),sample.begin()+index,sample.end());
    result.suggested=sample[index];
    if (result.suggested<=result.low || result.suggested>=result.high)
        result.suggested=result.low*.95+result.high*.05;
    return result;
}
}
