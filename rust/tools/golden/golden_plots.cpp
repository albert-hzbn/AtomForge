// Reference plots of scientific tool results (science/ResultPlots) for
// af-science::plots: the plots, CSV text and an overlay of every result in
// tests/golden/tools.json (the C++ tool outputs), plus histogram edge cases.
// Usage (from the repository root): golden_plots <tools.json> <output dir>.

#include "science/Json.h"
#include "science/ResultPlots.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

using atomforge::science::Json;
using namespace atomforge::science;

namespace
{
Json numbers(const std::vector<double>& values)
{
    // Non-finite values as strings (JSON has no NaN).
    Json out = Json::array();
    for (double v : values) {
        if (std::isfinite(v)) out.push(v);
        else out.push(std::string(std::isnan(v) ? "nan" : v > 0 ? "inf" : "-inf"));
    }
    return out;
}

Json seriesJson(const PlotSeries& s)
{
    Json o = Json::object();
    o["name"] = s.name;
    o["x"] = numbers(s.x);
    o["y"] = numbers(s.y);
    o["points"] = s.points;
    return o;
}

Json plotJson(const PlotSpec& p)
{
    Json o = Json::object();
    o["title"] = p.title;
    o["x_label"] = p.xLabel;
    o["y_label"] = p.yLabel;
    Json series = Json::array();
    for (const auto& s : p.series) series.push(seriesJson(s));
    o["series"] = series;
    Json markers = Json::array();
    for (const auto& m : p.markers) markers.push(Json::array({m.x, m.label}));
    o["markers"] = markers;
    o["horizontal"] = numbers(p.horizontal);
    o["log_y"] = p.logY;
    o["csv"] = plotCsv(p);
    return o;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_plots <tools.json> <output dir>\n"); return 1; }
    std::ifstream in(argv[1]);
    std::stringstream text;
    text << in.rdbuf();
    const Json cases = Json::parse(text.str());
    Json out = Json::array();
    std::vector<std::pair<std::string, std::vector<PlotSpec>>> earlier;
    for (const auto& item : cases.items()) {
        Json entry = Json::object();
        const std::string tool = item.at("tool").string();
        entry["tool"] = tool;
        if (!item.contains("result")) { out.push(entry); continue; }
        try {
            const auto plots = resultPlots(tool, item.at("result"));
            Json list = Json::array();
            for (const auto& p : plots) list.push(plotJson(p));
            entry["plots"] = list;
            // The first plot over the earlier runs of the same tool.
            if (!plots.empty()) {
                std::vector<std::pair<std::string, std::vector<PlotSpec>>> same;
                for (const auto& [label, runs] : earlier) if (label.rfind(tool + "#", 0) == 0) same.push_back({label, runs});
                entry["overlay"] = plotJson(overlayPlots(plots[0], "now", same));
            }
            earlier.push_back({tool + "#" + std::to_string(earlier.size()), plots});
        } catch (const std::exception& e) {
            entry["error"] = std::string(e.what());
        }
        out.push(entry);
    }
    // Histogram edge cases.
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    Json histograms = Json::array();
    for (const auto& values : std::vector<std::vector<double>>{{}, {nan, inf}, {2.5}, {-1, -1, -1}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1e-9, 3e-9, nan, 2e-9}, {-5, 0.25, 7.75, 1e6}})
        for (int bins : {1, 3, 40})
            histograms.push(seriesJson(histogram("h", values, bins)));
    Json result = Json::object();
    result["tools"] = out;
    result["histograms"] = histograms;
    std::ofstream(std::string(argv[2]) + "/plots.json") << result.dump();
    std::printf("wrote %s/plots.json (%zu results)\n", argv[2], out.size());
    return 0;
}
