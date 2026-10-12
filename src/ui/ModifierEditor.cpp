#include "ui/UnitText.h"
#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>

using atomforge::pipeline::Json;

namespace modifierFlags
{
namespace
{
// Splits a command line at spaces, keeping "quoted words" together (unquoted).
std::vector<std::string> tokens(const std::string& text)
{
    std::vector<std::string> out;
    std::string word;
    bool quoted = false, any = false;
    for (char c : text) {
        if (c == '"') { quoted = !quoted; any = true; }
        else if (!quoted && (c == ' ' || c == '\t' || c == '\n')) {
            if (any) out.push_back(word);
            word.clear();
            any = false;
        } else { word += c; any = true; }
    }
    if (any) out.push_back(word);
    return out;
}

std::string quote(const std::string& value)
{
    return value.find_first_of(" \t") == std::string::npos ? value : "\"" + value + "\"";
}

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t");
    if (first == std::string::npos) return {};
    return text.substr(first, text.find_last_not_of(" \t") - first + 1);
}
}

Usage parse(const std::string& help)
{
    // [preamble], "TITLE (--build x)", usage lines (starting with "--" or "["),
    // description. Without a title line there is no usage.
    Usage usage;
    const std::size_t title = help.find("(--");
    if (title == std::string::npos) return usage;
    std::stringstream stream(help.substr(help.find('\n', title) == std::string::npos ? help.size() : help.find('\n', title) + 1));
    std::string line, usageText;
    bool inUsage = true;
    while (std::getline(stream, line)) {
        const std::string text = trim(line);
        if (text.empty() && usage.description.empty()) continue;
        if (inUsage && (text.rfind("--", 0) == 0 || text.rfind("[", 0) == 0)) { usageText += text + " "; continue; }
        inUsage = false;
        usage.description += (usage.description.empty() ? "" : " ") + text;
    }
    int depth = 0;
    for (std::size_t i = 0; i < usageText.size();) {
        const char c = usageText[i];
        if (c == '[') { ++depth; ++i; continue; }
        if (c == ']') { --depth; ++i; continue; }
        if (usageText.compare(i, 2, "--") != 0) { ++i; continue; }
        std::size_t j = i + 2;
        while (j < usageText.size() && (std::isalnum(static_cast<unsigned char>(usageText[j])) || usageText[j] == '-')) ++j;
        Flag flag;
        flag.name = usageText.substr(i + 2, j - i - 2);
        flag.optional = depth > 0;
        std::size_t k = j;
        while (k < usageText.size() && usageText[k] == ' ') ++k;
        if (k < usageText.size() && usageText[k] == '"') {
            const std::size_t end = usageText.find('"', k + 1);
            flag.value = usageText.substr(k + 1, end - k - 1);
            k = end + 1;
        } else if (k < usageText.size() && usageText.compare(k, 2, "--") != 0 && !std::strchr("[]|(", usageText[k])) {
            std::size_t end = k;
            while (end < usageText.size() && !std::strchr(" ]|", usageText[end])) ++end;
            flag.value = usageText.substr(k, end - k);
            k = end;
        }
        const std::size_t next = usageText.find("--", k);
        flag.repeatable = usageText.substr(k, next == std::string::npos ? std::string::npos : next - k).find("repeatable") != std::string::npos;
        if (flag.name != "input" && flag.name != "output") usage.flags.push_back(flag);
        i = k;
    }
    return usage;
}

std::vector<std::string> values(const Usage& usage, const std::string& options, std::string& other)
{
    std::vector<std::string> out(usage.flags.size());
    other.clear();
    const std::vector<std::string> words = tokens(options);
    for (std::size_t i = 0; i < words.size(); ++i) {
        std::size_t f = usage.flags.size();
        if (words[i].rfind("--", 0) == 0)
            for (std::size_t k = 0; k < usage.flags.size(); ++k)
                if (words[i].substr(2) == usage.flags[k].name) f = k;
        if (f == usage.flags.size()) { other += (other.empty() ? "" : " ") + quote(words[i]); continue; }
        if (usage.flags[f].value.empty()) { out[f] = "1"; continue; }
        const std::string value = i + 1 < words.size() ? words[++i] : std::string();
        out[f] = usage.flags[f].repeatable && !out[f].empty() ? out[f] + "; " + value : value;
    }
    return out;
}

std::string compose(const Usage& usage, const std::vector<std::string>& values, const std::string& other)
{
    std::string out;
    const auto add = [&](const std::string& word) { out += (out.empty() ? "" : " ") + word; };
    for (std::size_t k = 0; k < usage.flags.size() && k < values.size(); ++k) {
        const Flag& flag = usage.flags[k];
        if (flag.value.empty()) { if (values[k] == "1") add("--" + flag.name); continue; }
        std::vector<std::string> entries;
        if (flag.repeatable) {
            std::stringstream stream(values[k]);
            for (std::string entry; std::getline(stream, entry, ';');) entries.push_back(trim(entry));
        } else entries.push_back(trim(values[k]));
        for (const auto& entry : entries)
            if (!entry.empty()) add("--" + flag.name + " " + quote(entry));
    }
    if (!trim(other).empty()) add(trim(other));
    return out;
}
}

std::string modifierDescription(const atomforge::pipeline::ModifierType& type)
{
    for (const auto& p : type.parameters)
        if (std::string(p.kind) == "options") {
            modifierFlags::Usage usage = modifierFlags::parse(type.help);
            if (!usage.flags.empty()) return usage.description;
        }
    return type.help;
}

namespace
{
// "min-separation" -> "Min separation"; single letters stay as they are.
std::string flagLabel(const std::string& name)
{
    if (name == "symprec") return "Symmetry tolerance";
    if (name == "nmax") return "Max Miller index";
    if (name == "no-primitive") return "Keep conventional cell";
    std::string label = name;
    std::replace(label.begin(), label.end(), '-', ' ');
    if (label.size() > 1) label[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(label[0])));
    return label;
}

// One table row per flag; true when `options` changed.
bool drawOptionFlags(const modifierFlags::Usage& usage, std::string& options)
{
    std::string other;
    std::vector<std::string> values = modifierFlags::values(usage, options, other);
    bool edited = false;
    for (std::size_t k = 0; k < usage.flags.size(); ++k) {
        const modifierFlags::Flag& flag = usage.flags[k];
        ImGui::PushID(static_cast<int>(k));
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(flagLabel(flag.name).c_str());
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("--%s %s", flag.name.c_str(), flag.value.c_str());
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (flag.value.empty()) {
            bool on = values[k] == "1";
            if (ImGui::Checkbox("##flag", &on)) { values[k] = on ? "1" : ""; edited = true; }
        } else {
            char buffer[512];
            std::snprintf(buffer, sizeof(buffer), "%s", values[k].c_str());
            char* end = nullptr;
            std::strtod(flag.value.c_str(), &end);
            // Defaults read "default 4"; one-letter placeholders (P, N, V) "number".
            std::string hint = end && *end == '\0' ? "default " + flag.value : flag.value.size() == 1 ? "number" : flag.value;
            std::transform(hint.begin(), hint.end(), hint.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (flag.repeatable) hint += "; ...";
            if (!flag.optional) hint += " (required)";
            if (ImGui::InputTextWithHint("##flag", hint.c_str(), buffer, sizeof(buffer))) { values[k] = buffer; edited = true; }
        }
        ImGui::PopID();
    }
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextDisabled("Other flags");
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
    char buffer[1024];
    std::snprintf(buffer, sizeof(buffer), "%s", other.c_str());
    if (ImGui::InputTextWithHint("##other", "optional", buffer, sizeof(buffer))) { other = buffer; edited = true; }
    if (edited) options = modifierFlags::compose(usage, values, other);
    return edited;
}
}

void drawModifierParameters(const atomforge::pipeline::ModifierType& type, const Json& parameters,
                            const std::function<void(const std::string&, const Json&)>& changed)
{
    // A two-column form: labels on the left, fields filling the right.
    if (type.parameters.empty() || !ImGui::BeginTable("##parameters", 2, ImGuiTableFlags_SizingStretchProp)) return;
    ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.42f);
    ImGui::TableSetupColumn("field", ImGuiTableColumnFlags_WidthStretch, 0.58f);
    for (const auto& p : type.parameters) {
        ImGui::PushID(p.name);
        const Json value = atomforge::pipeline::parameter(parameters, type, p.name);
        const std::string kind = p.kind;
        if (kind == "options") {
            const modifierFlags::Usage usage = modifierFlags::parse(type.help);
            if (!usage.flags.empty()) {
                std::string options = value.string();
                if (drawOptionFlags(usage, options)) changed(p.name, Json(options));
                ImGui::PopID();
                continue;
            }
        }
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::AlignTextToFramePadding();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextUnformatted(unitText::display(p.label).c_str());
        ImGui::PopTextWrapPos();
        ImGui::TableSetColumnIndex(1);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (kind == "bool") {
            bool v = value.boolean();
            if (ImGui::Checkbox("##value", &v)) changed(p.name, Json(v));
        } else if (kind == "int") {
            int v = static_cast<int>(value.number());
            if (ImGui::InputInt("##value", &v)) changed(p.name, Json(v));
        } else if (kind == "float") {
            double v = value.number();
            if (ImGui::InputDouble("##value", &v, 0, 0, "%.6g")) changed(p.name, Json(v));
        } else if (kind == "vector" || kind == "ivec3") {
            double v[3] = {value.items()[0].number(), value.items()[1].number(), value.items()[2].number()};
            bool edited = false;
            if (kind == "ivec3") {
                int n[3] = {static_cast<int>(v[0]), static_cast<int>(v[1]), static_cast<int>(v[2])};
                edited = ImGui::InputInt3("##value", n);
                for (int k = 0; k < 3; ++k) v[k] = n[k];
            } else edited = ImGui::InputScalarN("##value", ImGuiDataType_Double, v, 3, nullptr, nullptr, "%.6g");
            if (edited) changed(p.name, Json::array({v[0], v[1], v[2]}));
        } else if (kind == "choice") {
            std::vector<std::string> options;
            std::stringstream stream(p.options);
            for (std::string option; std::getline(stream, option, '|');) options.push_back(option);
            const std::string current = value.string();
            if (ImGui::BeginCombo("##value", current.c_str())) {
                for (const auto& option : options)
                    if (ImGui::Selectable(option.c_str(), option == current)) changed(p.name, Json(option));
                ImGui::EndCombo();
            }
        } else {
            // string, element(s), expression or command-line options
            char buffer[2048];
            std::snprintf(buffer, sizeof(buffer), "%s", value.string().c_str());
            const char* hint = kind == "expression" ? "e.g. fz > 0.5 && element == O"
                             : kind == "options" ? "command-line flags, e.g. --count 4"
                             : kind == "elements" ? "e.g. Cu Ni" : "";
            if (ImGui::InputTextWithHint("##value", hint, buffer, sizeof(buffer))) changed(p.name, Json(std::string(buffer)));
        }
        ImGui::PopID();
    }
    ImGui::EndTable();
}
