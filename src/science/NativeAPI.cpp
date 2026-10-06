// C interface of the native scientific tools for the Python package
// (atomforge.native, via ctypes). Requests and results are JSON text, the same
// documents used by `AtomForge --science`.
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>

#ifdef _WIN32
#define AFS_API extern "C" __declspec(dllexport)
#else
#define AFS_API extern "C" __attribute__((visibility("default")))
#endif

namespace
{
using atomforge::science::Json;

char* copyOut(const std::string& text)
{
    char* out = static_cast<char*>(std::malloc(text.size() + 1));
    if (out) std::memcpy(out, text.c_str(), text.size() + 1);
    return out;
}

Json catalogJson()
{
    Json tools = Json::array();
    for (const auto& tool : scienceToolCatalog()) {
        Json entry = Json::object();
        entry["id"] = tool.id;
        entry["title"] = tool.title;
        entry["category"] = tool.category;
        entry["help"] = tool.help;
        Json parameters = Json::array();
        for (const auto& parameter : tool.parameters) {
            Json item = Json::object();
            item["name"] = parameter.name;
            item["label"] = parameter.label;
            item["kind"] = parameter.kind;
            item["required"] = parameter.required;
            item["default"] = parameter.value[0] ? Json::parse(parameter.value) : Json();
            parameters.push(item);
        }
        entry["parameters"] = parameters;
        tools.push(entry);
    }
    return tools;
}
}

AFS_API int afs_version() { return 1; }

AFS_API void afs_free(char* text) { std::free(text); }

// JSON array of {id, title, category, help, parameters[]}.
AFS_API char* afs_catalog()
{
    try { return copyOut(catalogJson().dump()); }
    catch (...) { return nullptr; }
}

// Runs a tool. Returns {"result": ..., "frames": [structures], "velocities": ...,
// "times_fs": ...} with *status = 0, or an error message with *status = 1.
// Relative file references resolve against base_dir.
AFS_API char* afs_run(const char* tool, const char* request, const char* baseDir, int* status)
{
    *status = 1;
    try {
        const auto output = atomforge::science::runTool(tool, Json::parse(request ? request : "{}"),
                                                         std::filesystem::u8path(baseDir && *baseDir ? baseDir : "."));
        Json document = Json::object();
        document["result"] = output.result;
        Json frames = Json::array();
        for (const auto& frame : output.frames) frames.push(atomforge::science::structureJson(frame));
        document["frames"] = frames;
        Json velocities = Json::array();
        for (const auto& frame : output.velocities) {
            Json rows = Json::array();
            for (const auto& v : frame) rows.push(Json::array({v[0], v[1], v[2]}));
            velocities.push(rows);
        }
        document["velocities"] = velocities;
        Json times = Json::array();
        for (double t : output.times) times.push(t);
        document["times_fs"] = times;
        *status = 0;
        return copyOut(document.dump());
    } catch (const std::exception& error) {
        return copyOut(error.what());
    } catch (...) {
        return copyOut("Unknown native error");
    }
}
