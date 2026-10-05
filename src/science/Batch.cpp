#include "science/Batch.h"
#include "science/ScienceTools.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
bool matches(const std::string& pattern, const std::string& text, std::size_t p = 0, std::size_t t = 0)
{
    while (p < pattern.size()) {
        if (pattern[p] == '*') {
            for (std::size_t skip = t; skip <= text.size(); ++skip)
                if (matches(pattern, text, p + 1, skip)) return true;
            return false;
        }
        if (t >= text.size() || (pattern[p] != '?' && pattern[p] != text[t])) return false;
        ++p; ++t;
    }
    return t == text.size();
}

std::vector<std::filesystem::path> expand(const std::string& pattern, const std::filesystem::path& base)
{
    const auto full = (base / std::filesystem::u8path(pattern)).lexically_normal();
    const std::string name = full.filename().u8string();
    if (name.find_first_of("*?") == std::string::npos) return {full};
    std::vector<std::filesystem::path> result;
    const auto folder = full.parent_path();
    if (std::filesystem::is_directory(folder))
        for (const auto& entry : std::filesystem::directory_iterator(folder))
            if (entry.is_regular_file() && matches(name, entry.path().filename().u8string())) result.push_back(entry.path());
    std::sort(result.begin(), result.end());
    if (result.empty()) throw std::runtime_error("No files match " + full.u8string());
    return result;
}

std::string cell(const Json& value)
{
    if (value.isNull()) return "";
    if (value.isString()) {
        std::string text = value.string();
        if (text.find_first_of(",\"\n") != std::string::npos) {
            std::string quoted = "\"";
            for (char c : text) { if (c == '"') quoted += '"'; quoted += c; }
            return quoted + "\"";
        }
        return text;
    }
    if (value.isNumber()) {
        char buffer[40];
        std::snprintf(buffer, sizeof(buffer), "%.12g", value.number());
        return buffer;
    }
    if (value.isBool()) return value.boolean() ? "true" : "false";
    return cell(Json(value.dump()));
}
}

Json valueAt(const Json& value, const std::string& path)
{
    const Json* current = &value;
    std::size_t start = 0;
    while (start <= path.size()) {
        const std::size_t end = path.find('.', start);
        const std::string key = path.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (current->isArray()) {
            std::size_t index = 0;
            try { index = static_cast<std::size_t>(std::stoul(key)); } catch (const std::exception&) { return Json(); }
            if (index >= current->size()) return Json();
            current = &current->items()[index];
        } else if (current->isObject()) {
            current = current->find(key);
            if (!current) return Json();
        } else return Json();
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return *current;
}

Json runBatch(const Json& batch, const std::filesystem::path& base, const StructureReader& reader)
{
    if (!batch.isObject() || !batch.contains("tool")) throw std::runtime_error("A batch needs a \"tool\"");
    const std::string tool = batch.at("tool").string();
    const Json baseParameters = batch.contains("base") ? batch.at("base") : Json::object();
    if (!baseParameters.isObject()) throw std::runtime_error("\"base\" must be an object");
    const bool continueOnError = !batch.contains("continue_on_error") || batch.at("continue_on_error").boolean();
    // Axes of the run grid: each is (parameter, list of values).
    std::vector<std::pair<std::string, std::vector<Json>>> axes;
    if (const Json* sweep = batch.find("sweep")) {
        for (const auto& [name, values] : sweep->members()) {
            if (!values.isArray() || values.size() == 0) throw std::runtime_error("sweep." + name + " must be a nonempty list");
            axes.push_back({name, values.items()});
        }
    }
    if (const Json* files = batch.find("files")) {
        for (const auto& [name, pattern] : files->members()) {
            std::vector<Json> references;
            for (const auto& path : expand(pattern.string(), base)) {
                // *_file parameters take a path; data and structure inputs take a reference.
                if (name.size() > 5 && name.compare(name.size() - 5, 5, "_file") == 0) references.push_back(path.u8string());
                else { Json reference = Json::object(); reference["file"] = path.u8string(); references.push_back(reference); }
            }
            axes.push_back({name, references});
        }
    }
    std::vector<std::string> collect;
    if (const Json* fields = batch.find("collect"))
        for (const auto& field : fields->items()) collect.push_back(field.string());
    std::size_t total = 1;
    for (const auto& axis : axes) total *= axis.second.size();
    if (total > 100000) throw std::runtime_error("The batch has more than 100000 runs");
    Json runs = Json::array(), rows = Json::array(), columns = Json::array();
    for (const auto& axis : axes) columns.push(axis.first);
    for (const auto& field : collect) columns.push(field);
    columns.push("error");
    std::size_t failures = 0;
    std::vector<std::size_t> index(axes.size(), 0);
    for (std::size_t run = 0; run < total; ++run) {
        taskProgress(static_cast<double>(run) / static_cast<double>(total));
        Json parameters = baseParameters;
        Json row = Json::array();
        for (std::size_t a = 0; a < axes.size(); ++a) {
            const Json& value = axes[a].second[index[a]];
            parameters[axes[a].first] = value;
            row.push(value.isObject() && value.contains("file") ? Json(std::filesystem::u8path(value.at("file").string()).filename().u8string()) : value);
        }
        Json record = Json::object();
        record["parameters"] = parameters;
        try {
            const ToolOutput output = runTool(tool, parameters, base, reader);
            record["result"] = output.result;
            for (const auto& field : collect) row.push(valueAt(output.result, field));
            row.push(Json());
        } catch (const std::exception& error) {
            if (std::string(error.what()) == "Calculation cancelled") throw;
            if (!continueOnError) throw std::runtime_error("Run " + std::to_string(run + 1) + ": " + error.what());
            ++failures;
            record["error"] = error.what();
            for (std::size_t k = 0; k < collect.size(); ++k) row.push(Json());
            row.push(error.what());
        }
        runs.push(record);
        rows.push(row);
        for (std::size_t a = axes.size(); a-- > 0;) {
            if (++index[a] < axes[a].second.size()) break;
            index[a] = 0;
        }
    }
    Json table = Json::object();
    table["columns"] = columns;
    table["rows"] = rows;
    Json result = Json::object();
    result["tool"] = tool;
    result["run_count"] = total;
    result["failures"] = failures;
    result["table"] = table;
    result["runs"] = runs;
    return result;
}

std::string batchCsv(const Json& batchResult)
{
    const Json& table = batchResult.at("table");
    std::string text;
    const auto& columns = table.at("columns").items();
    for (std::size_t c = 0; c < columns.size(); ++c) text += (c ? "," : "") + cell(columns[c]);
    text += "\n";
    for (const auto& row : table.at("rows").items()) {
        for (std::size_t c = 0; c < row.size(); ++c) text += (c ? "," : "") + cell(row.items()[c]);
        text += "\n";
    }
    return text;
}
}
