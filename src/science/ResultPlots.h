#pragma once

#include "science/Json.h"

#include <string>
#include <vector>

namespace atomforge::science
{
struct PlotSeries
{
    std::string name;
    std::vector<double> x;
    std::vector<double> y;
    bool points = false;  // draw markers instead of a line
};

struct PlotMarker
{
    double x = 0;
    std::string label;
};

struct PlotSpec
{
    std::string title;
    std::string xLabel;
    std::string yLabel;
    std::vector<PlotSeries> series;
    std::vector<PlotMarker> markers;  // labelled vertical lines (e.g. k-path points)
    std::vector<double> horizontal;   // reference horizontal lines
    bool logY = false;
};

// Plots describing a tool result; empty when the result has no curve to show.
std::vector<PlotSpec> resultPlots(const std::string& tool, const Json& result);

// Fixed-width histogram of the finite values (centres and counts).
PlotSeries histogram(const std::string& name, const std::vector<double>& values, int bins);

// CSV text of a plot's series (one x/y column pair per series).
std::string plotCsv(const PlotSpec& plot);
}
