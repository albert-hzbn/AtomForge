#pragma once
// Command-line option text of Build steps ("--count 4 --frac \"Cu=0.7,Zn=0.3\""):
// splitting it into tokens, reading flags from it and writing it back. The
// Build and Edit dialogs use these to edit a pipeline step's options.
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace atomforge::pipeline::options
{
// Shell-style tokens (quotes group words and are removed).
inline std::vector<std::string> tokens(const std::string& text)
{
    std::vector<std::string> out;
    std::string current;
    bool inToken = false;
    char quote = 0;
    for (char c : text) {
        if (quote) { if (c == quote) quote = 0; else current += c; continue; }
        if (c == '"' || c == '\'') { quote = c; inToken = true; continue; }
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (inToken) { out.push_back(current); current.clear(); inToken = false; }
            continue;
        }
        current += c;
        inToken = true;
    }
    if (quote) throw std::runtime_error("Unterminated quote in the options");
    if (inToken) out.push_back(current);
    return out;
}

// Flags of an option text. A flag's value is the token after it unless that
// token is another flag (a negative number is a value, not a flag).
class Reader
{
public:
    explicit Reader(const std::string& text) : m_tokens(tokens(text)) {}

    bool has(const std::string& flag) const
    {
        for (const auto& t : m_tokens) if (t == flag) return true;
        return false;
    }
    // Every value given to `flag`, in order (flags such as --atom repeat).
    std::vector<std::string> values(const std::string& flag) const
    {
        std::vector<std::string> out;
        for (std::size_t i = 0; i + 1 < m_tokens.size(); ++i)
            if (m_tokens[i] == flag && !isFlag(m_tokens[i + 1])) out.push_back(m_tokens[i + 1]);
        return out;
    }
    std::string text(const std::string& flag, const std::string& fallback = "") const
    {
        const auto all = values(flag);
        return all.empty() ? fallback : all.back();
    }
    double number(const std::string& flag, double fallback) const
    {
        const std::string value = text(flag);
        if (value.empty()) return fallback;
        try { return std::stod(value); } catch (const std::exception&) { return fallback; }
    }
    int integer(const std::string& flag, int fallback) const
    {
        const std::string value = text(flag);
        if (value.empty()) return fallback;
        try { return std::stoi(value); } catch (const std::exception&) { return fallback; }
    }
    // Numbers in a value such as "0 0 1" or "1,1,0".
    std::vector<double> numbers(const std::string& flag) const
    {
        std::string value = text(flag);
        for (char& c : value) if (c == ',') c = ' ';
        std::istringstream in(value);
        std::vector<double> out;
        for (double x; in >> x;) out.push_back(x);
        return out;
    }

    static bool isFlag(const std::string& token)
    {
        return token.size() > 1 && token[0] == '-' && !std::isdigit(static_cast<unsigned char>(token[1])) && token[1] != '.';
    }

private:
    std::vector<std::string> m_tokens;
};

// Builds an option text; values with spaces or special characters are quoted.
class Writer
{
public:
    Writer& flag(const std::string& name) { append(name); return *this; }
    Writer& add(const std::string& name, const std::string& value) { append(name); append(quoted(value)); return *this; }
    Writer& add(const std::string& name, const char* value) { return add(name, std::string(value)); }
    Writer& add(const std::string& name, double value) { return add(name, number(value)); }
    Writer& add(const std::string& name, float value) { return add(name, number(value)); }
    Writer& add(const std::string& name, int value) { return add(name, std::to_string(value)); }
    // Several numbers as one value ("0 0 1").
    Writer& add(const std::string& name, const std::vector<double>& values)
    {
        std::string text;
        for (double v : values) text += (text.empty() ? "" : " ") + number(v);
        return add(name, text);
    }
    const std::string& str() const { return m_text; }

    static std::string number(double value)
    {
        std::ostringstream out;
        out.precision(10);
        out << value;
        return out.str();
    }
    static std::string quoted(const std::string& value)
    {
        if (!value.empty() && value.find_first_of(" \t\"'|") == std::string::npos) return value;
        const char quote = value.find('"') == std::string::npos ? '"' : '\'';
        return quote + value + quote;
    }

private:
    void append(const std::string& token) { m_text += (m_text.empty() ? "" : " ") + token; }
    std::string m_text;
};
}
