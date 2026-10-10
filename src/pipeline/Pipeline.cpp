#include "pipeline/Pipeline.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace atomforge::pipeline
{
std::size_t PipelineData::selectedCount() const
{
    return static_cast<std::size_t>(std::count(selected.begin(), selected.end(), 1));
}

const ModifierType* findModifierType(const std::string& id)
{
    for (const auto& type : modifierTypes())
        if (id == type.id) return &type;
    return nullptr;
}

namespace
{
const ModifierType& requireType(const std::string& id)
{
    if (const ModifierType* type = findModifierType(id)) return *type;
    throw std::runtime_error("Unknown modifier '" + id + "'");
}

const ModifierParameter* findParameter(const ModifierType& type, const std::string& name)
{
    for (const auto& p : type.parameters)
        if (name == p.name) return &p;
    return nullptr;
}

std::vector<std::string> words(const std::string& text)
{
    std::istringstream in(text);
    std::vector<std::string> out;
    for (std::string w; in >> w;) out.push_back(w);
    return out;
}

// Parses one value in the text syntax for a parameter kind.
Json valueFor(const ModifierParameter& p, const std::string& text)
{
    const std::string kind = p.kind;
    const auto number = [&](const std::string& token) {
        try {
            std::size_t used = 0;
            const double value = std::stod(token, &used);
            if (used != token.size()) throw std::invalid_argument(token);
            return value;
        } catch (const std::exception&) {
            throw std::runtime_error(std::string("Parameter ") + p.name + ": '" + token + "' is not a number");
        }
    };
    if (kind == "int" || kind == "float") {
        const double value = number(text);
        if (kind == "int" && value != std::floor(value)) throw std::runtime_error(std::string("Parameter ") + p.name + " must be an integer");
        return Json(value);
    }
    if (kind == "bool") {
        if (text == "true" || text == "1" || text == "yes" || text == "on") return Json(true);
        if (text == "false" || text == "0" || text == "no" || text == "off") return Json(false);
        throw std::runtime_error(std::string("Parameter ") + p.name + " must be true or false");
    }
    if (kind == "vector" || kind == "ivec3") {
        std::string spaced = text;
        std::replace(spaced.begin(), spaced.end(), ',', ' ');
        const auto parts = words(spaced);
        if (parts.size() != 3) throw std::runtime_error(std::string("Parameter ") + p.name + " needs three numbers");
        Json list = Json::array();
        for (const auto& part : parts) {
            const double value = number(part);
            if (kind == "ivec3" && value != std::floor(value)) throw std::runtime_error(std::string("Parameter ") + p.name + " needs integers");
            list.push(value);
        }
        return list;
    }
    if (kind == "choice") {
        const std::string options = std::string("|") + p.options + "|";
        if (options.find("|" + text + "|") == std::string::npos)
            throw std::runtime_error(std::string("Parameter ") + p.name + " must be one of " + p.options);
    }
    return Json(text);
}

std::string textFor(const ModifierParameter& p, const Json& value)
{
    const std::string kind = p.kind;
    std::ostringstream out;
    out.precision(12);
    if (kind == "bool") return value.boolean() ? "true" : "false";
    if (kind == "int" || kind == "float") { out << value.number(); return out.str(); }
    if (kind == "vector" || kind == "ivec3") {
        const auto& items = value.items();
        out << items[0].number() << ',' << items[1].number() << ',' << items[2].number();
        return out.str();
    }
    return value.string();
}

// Splits on '|' outside quotes; "||" is the logical operator, not a pipe.
std::vector<std::string> steps(const std::string& text)
{
    std::vector<std::string> out(1);
    char quote = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (quote) { if (c == quote) quote = 0; out.back() += c; continue; }
        if (c == '"' || c == '\'') { quote = c; out.back() += c; continue; }
        if (c == '|') {
            if (i + 1 < text.size() && text[i + 1] == '|') { out.back() += "||"; ++i; continue; }
            out.emplace_back();
            continue;
        }
        out.back() += c;
    }
    if (quote) throw std::runtime_error("Unterminated quote in the pipeline text");
    return out;
}

// Tokens of one step: whitespace-separated, quotes group (and are removed).
std::vector<std::string> tokens(const std::string& step)
{
    std::vector<std::string> out;
    std::string current;
    bool inToken = false;
    char quote = 0;
    for (char c : step) {
        if (quote) { if (c == quote) quote = 0; else current += c; continue; }
        if (c == '"' || c == '\'') { quote = c; inToken = true; continue; }
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (inToken) { out.push_back(current); current.clear(); inToken = false; }
            continue;
        }
        current += c;
        inToken = true;
    }
    if (inToken) out.push_back(current);
    return out;
}
}

Modifier makeModifier(const std::string& id)
{
    const ModifierType& type = requireType(id);
    Modifier modifier;
    modifier.type = id;
    for (const auto& p : type.parameters) modifier.parameters[p.name] = valueFor(p, p.value);
    return modifier;
}

Json parameter(const Json& parameters, const ModifierType& type, const std::string& name)
{
    if (const Json* value = parameters.find(name)) return *value;
    const ModifierParameter* p = findParameter(type, name);
    if (!p) throw std::runtime_error(std::string(type.id) + " has no parameter " + name);
    return valueFor(*p, p->value);
}

PipelineData Pipeline::evaluate(const Structure& input, std::vector<StageResult>* stages) const
{
    PipelineCache cache;
    PipelineData out = evaluate(input, cache);
    if (stages) *stages = cache.stages;
    return out;
}

PipelineData Pipeline::evaluate(const Structure& input, PipelineCache& cache) const
{
    cache.outputs.resize(modifiers.size());
    cache.stages.resize(modifiers.size());
    cache.valid = std::min(cache.valid, modifiers.size());
    PipelineData data;
    if (cache.valid > 0) data = cache.outputs[cache.valid - 1];
    else {
        data.structure = input;
        data.clearSelection();
    }
    bool failed = false;
    for (std::size_t i = 0; i < cache.valid; ++i) failed = failed || !cache.stages[i].error.empty();
    for (std::size_t i = cache.valid; i < modifiers.size(); ++i) {
        taskCheckpoint();
        StageResult stage;
        const Modifier& modifier = modifiers[i];
        if (failed || !modifier.enabled) {
            stage.skipped = true;
        } else {
            const auto start = std::chrono::steady_clock::now();
            try {
                const ModifierType& type = requireType(modifier.type);
                PipelineData next = data;
                next.notes.clear();
                type.apply(modifier.parameters, next);
                if (next.selected.size() != next.structure.atoms.size()) next.selected.resize(next.structure.atoms.size(), 0);
                data = std::move(next);
                stage.notes = data.notes;
            } catch (const std::exception& error) {
                stage.error = error.what();
                failed = true;
            }
            stage.milliseconds = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        }
        stage.atoms = data.structure.atoms.size();
        stage.selected = data.selectedCount();
        cache.outputs[i] = data;
        cache.stages[i] = stage;
    }
    cache.valid = modifiers.size();
    return modifiers.empty() ? data : cache.outputs.back();
}

Json Pipeline::toJson() const
{
    Json list = Json::array();
    for (const auto& m : modifiers) {
        Json item = Json::object();
        item["type"] = m.type;
        item["enabled"] = m.enabled;
        if (!m.label.empty()) item["label"] = m.label;
        item["parameters"] = m.parameters;
        list.push(item);
    }
    Json json = Json::object();
    json["format"] = "atomforge-pipeline";
    json["version"] = 1;
    json["modifiers"] = list;
    return json;
}

Pipeline Pipeline::fromJson(const Json& json)
{
    if (!json.isObject() || !json.contains("modifiers")) throw std::runtime_error("A pipeline needs a \"modifiers\" list");
    Pipeline pipeline;
    for (const auto& item : json.at("modifiers").items()) {
        Modifier modifier = makeModifier(item.at("type").string());
        const ModifierType& type = requireType(modifier.type);
        if (const Json* parameters = item.find("parameters"); parameters && parameters->isObject())
            for (const auto& [name, value] : parameters->members()) {
                const ModifierParameter* p = findParameter(type, name);
                if (!p) throw std::runtime_error(modifier.type + " has no parameter " + name);
                // Values are validated through the text form, which checks kinds and choices.
                modifier.parameters[name] = valueFor(*p, value.isString() ? value.string() : textFor(*p, value));
            }
        if (const Json* enabled = item.find("enabled")) modifier.enabled = enabled->boolean();
        if (const Json* label = item.find("label")) modifier.label = label->string();
        pipeline.modifiers.push_back(std::move(modifier));
    }
    return pipeline;
}

Pipeline Pipeline::parse(const std::string& text)
{
    Pipeline pipeline;
    for (const auto& step : steps(text)) {
        const auto parts = tokens(step);
        if (parts.empty()) {
            if (step.find_first_not_of(" \t\r\n") == std::string::npos && steps(text).size() == 1) break;
            throw std::runtime_error("Empty step in the pipeline text");
        }
        Modifier modifier = makeModifier(parts[0]);
        const ModifierType& type = requireType(modifier.type);
        std::size_t positional = 0;
        std::vector<std::string> pending;  // positional tokens
        for (std::size_t t = 1; t < parts.size(); ++t) {
            const std::string& token = parts[t];
            const auto equals = token.find('=');
            // key=value (but not comparisons such as x==1 or x<=1 inside expressions)
            if (equals != std::string::npos && equals > 0 && token[equals - 1] != '<' && token[equals - 1] != '>' && token[equals - 1] != '!' &&
                (equals + 1 >= token.size() || token[equals + 1] != '=') && findParameter(type, token.substr(0, equals))) {
                const ModifierParameter& p = *findParameter(type, token.substr(0, equals));
                modifier.parameters[p.name] = valueFor(p, token.substr(equals + 1));
            } else pending.push_back(token);
        }
        for (std::size_t t = 0; t < pending.size();) {
            if (positional >= type.parameters.size()) throw std::runtime_error(std::string(type.id) + ": too many values (" + pending[t] + ")");
            const ModifierParameter& p = type.parameters[positional++];
            const std::string kind = p.kind;
            std::string value;
            if (kind == "expression" || kind == "elements" || kind == "options") {  // the rest of the step
                // Words that were quoted (they contain spaces) keep their quotes.
                for (; t < pending.size(); ++t) {
                    const std::string& word = pending[t];
                    // (Command-line options only: there quotes group a flag's value.)
                    const bool spaced = kind == "options" && word.find_first_of(" \t") != std::string::npos;
                    const char mark = word.find('"') != std::string::npos ? '\'' : '"';
                    value += (value.empty() ? "" : " ") + (spaced ? mark + word + mark : word);
                }
            } else if ((kind == "vector" || kind == "ivec3") && pending[t].find(',') == std::string::npos) {
                if (t + 3 > pending.size()) throw std::runtime_error(std::string(type.id) + ": " + p.name + " needs three numbers");
                value = pending[t] + " " + pending[t + 1] + " " + pending[t + 2];
                t += 3;
            } else value = pending[t++];
            modifier.parameters[p.name] = valueFor(p, value);
        }
        pipeline.modifiers.push_back(std::move(modifier));
    }
    return pipeline;
}

std::string Pipeline::toText() const
{
    std::string out;
    for (const auto& m : modifiers) {
        if (!m.enabled) continue;
        const ModifierType& type = requireType(m.type);
        if (!out.empty()) out += " | ";
        out += m.type;
        for (const auto& p : type.parameters) {
            const Json value = parameter(m.parameters, type, p.name);
            const std::string text = textFor(p, value);
            if (text == textFor(p, valueFor(p, p.value))) continue;  // defaults are implied
            const bool quote = text.find_first_of(" \t|\"'") != std::string::npos || text.empty();
            const char mark = text.find('"') != std::string::npos ? '\'' : '"';
            out += std::string(" ") + p.name + "=" + (quote ? mark + text + mark : text);
        }
    }
    return out;
}
}
