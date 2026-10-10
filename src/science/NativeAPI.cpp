// C interface of the native scientific tools for the Python package
// (atomforge.native, via ctypes). Requests and results are JSON text, the same
// documents used by `AtomForge --science`.
#include "science/ScienceCatalog.h"
#include "science/ScienceTools.h"
#include "pipeline/Pipeline.h"

#include <cmath>

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

// JSON array of the structure-pipeline modifiers:
// {id, title, category, help, parameters: [{name, label, kind, default, options}]}.
AFS_API char* afs_pipeline_modifiers()
{
    try {
        Json types = Json::array();
        for (const auto& type : atomforge::pipeline::modifierTypes()) {
            Json entry = Json::object();
            entry["id"] = type.id;
            entry["title"] = type.title;
            entry["category"] = type.category;
            entry["help"] = type.help;
            const auto defaults = atomforge::pipeline::makeModifier(type.id).parameters;
            Json parameters = Json::array();
            for (const auto& p : type.parameters) {
                Json item = Json::object();
                item["name"] = p.name;
                item["label"] = p.label;
                item["kind"] = p.kind;
                item["default"] = defaults.at(p.name);
                item["options"] = p.options;
                parameters.push(item);
            }
            entry["parameters"] = parameters;
            types.push(entry);
        }
        return copyOut(types.dump());
    } catch (...) { return nullptr; }
}

// Runs a structure pipeline: request {"pipeline": JSON object or pipe text,
// "structure": structure JSON, "property": [per-atom values] (optional)}.
// Returns {"structure", "selected", "property", "property_name", "stages",
// "pipeline"} with *status = 0, or an error message with *status = 1.
AFS_API char* afs_pipeline_run(const char* request, int* status)
{
    *status = 1;
    try {
        const Json input = Json::parse(request ? request : "{}");
        const Json& definition = input.at("pipeline");
        const auto pipeline = definition.isString() ? atomforge::pipeline::Pipeline::parse(definition.string())
                                                    : atomforge::pipeline::Pipeline::fromJson(definition);
        Structure structure = atomforge::science::structureFromJson(input.at("structure"));
        if (const Json* property = input.find("property"); property && property->isArray() && property->size()) {
            if (property->size() != structure.atoms.size()) throw std::runtime_error("property needs one value per atom");
            for (const auto& value : property->items()) structure.atomProperty.push_back(value.isNumber() ? value.number() : std::nan(""));
            structure.atomPropertyName = input.contains("property_name") ? input.at("property_name").string() : "property";
        }
        std::vector<atomforge::pipeline::StageResult> stages;
        const auto data = pipeline.evaluate(structure, &stages);
        Json document = Json::object();
        document["structure"] = atomforge::science::structureJson(data.structure);
        Json selected = Json::array();
        for (char s : data.selected) selected.push(static_cast<bool>(s));
        document["selected"] = selected;
        Json property = Json::array();
        if (data.structure.atomProperty.size() == data.structure.atoms.size())
            for (double v : data.structure.atomProperty) property.push(std::isfinite(v) ? Json(v) : Json());
        document["property"] = property;
        document["property_name"] = data.structure.atomPropertyName;
        Json stageList = Json::array();
        for (std::size_t i = 0; i < stages.size(); ++i) {
            Json stage = Json::object();
            stage["type"] = pipeline.modifiers[i].type;
            stage["atoms"] = stages[i].atoms;
            stage["selected"] = stages[i].selected;
            stage["milliseconds"] = stages[i].milliseconds;
            stage["skipped"] = stages[i].skipped;
            stage["error"] = stages[i].error;
            Json notes = Json::array();
            for (const auto& note : stages[i].notes) notes.push(note);
            stage["notes"] = notes;
            stageList.push(stage);
        }
        document["stages"] = stageList;
        document["pipeline"] = pipeline.toJson();
        document["text"] = pipeline.toText();
        *status = 0;
        return copyOut(document.dump());
    } catch (const std::exception& error) {
        return copyOut(error.what());
    }
}

// Parses pipe text into pipeline JSON (validating every step), or *status = 1 with the message.
AFS_API char* afs_pipeline_parse(const char* text, int* status)
{
    *status = 1;
    try {
        const auto pipeline = atomforge::pipeline::Pipeline::parse(text ? text : "");
        *status = 0;
        return copyOut(pipeline.toJson().dump());
    } catch (const std::exception& error) {
        return copyOut(error.what());
    }
}

// Pipe text for pipeline JSON (validating it), or *status = 1 with the message.
AFS_API char* afs_pipeline_text(const char* json, int* status)
{
    *status = 1;
    try {
        const auto pipeline = atomforge::pipeline::Pipeline::fromJson(Json::parse(json ? json : "{}"));
        *status = 0;
        return copyOut(pipeline.toText());
    } catch (const std::exception& error) {
        return copyOut(error.what());
    }
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
