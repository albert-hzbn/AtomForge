// Reference results of the C++ science tools (tools.json) for the Rust
// af-science tests: each request of tests/fixtures/tools/requests.json run
// through runTool, with its report or the error.
// Usage (from the repository root): golden_tools <requests.json> <output dir>.

#include "science/Json.h"
#include "science/ScienceData.h"
#include "science/ScienceTools.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace atomforge::science;

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: golden_tools <requests.json> <output dir>\n");
        return 1;
    }
    std::ifstream in(argv[1]);
    std::stringstream text;
    text << in.rdbuf();
    const Json cases = Json::parse(text.str());
    Json out = Json::array();
    for (const auto& item : cases.items()) {
        const std::string tool = item.at("tool").string();
        Json entry = Json::object();
        entry["tool"] = tool;
        try {
            const ToolOutput output = runTool(tool, item.at("request"), ".");
            entry["result"] = output.result;
            entry["report"] = resultReport(tool, output.result);
            Json frames = Json::array();
            for (const auto& frame : output.frames) frames.push(structureJson(frame));
            entry["frames"] = frames;
            Json times = Json::array();
            for (double t : output.times) times.push(t);
            entry["times"] = times;
        } catch (const std::exception& e) {
            entry["error"] = std::string(e.what());
        }
        out.push(entry);
    }
    std::ofstream(std::string(argv[2]) + "/tools.json") << out.dump();
    std::printf("wrote %s/tools.json\n", argv[2]);
    return 0;
}
