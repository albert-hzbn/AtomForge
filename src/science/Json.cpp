#include "science/Json.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace atomforge::science
{
namespace
{
class Parser
{
public:
    explicit Parser(const std::string& text) : m_text(text) {}

    Json document()
    {
        Json value = parseValue(0);
        skipSpace();
        if (m_position != m_text.size()) fail("Unexpected trailing characters");
        return value;
    }

private:
    [[noreturn]] void fail(const std::string& message) const
    {
        throw std::runtime_error("Invalid JSON at character " + std::to_string(m_position) + ": " + message);
    }

    void skipSpace()
    {
        while (m_position < m_text.size() &&
               (m_text[m_position] == ' ' || m_text[m_position] == '\t' ||
                m_text[m_position] == '\n' || m_text[m_position] == '\r'))
            ++m_position;
    }

    bool consume(const char* literal)
    {
        std::size_t length = 0;
        while (literal[length]) ++length;
        if (m_text.compare(m_position, length, literal) != 0) return false;
        m_position += length;
        return true;
    }

    Json parseValue(int depth)
    {
        if (depth > 256) fail("Nesting is too deep");
        skipSpace();
        if (m_position >= m_text.size()) fail("Unexpected end of input");
        const char c = m_text[m_position];
        if (c == '{') return parseObject(depth);
        if (c == '[') return parseArray(depth);
        if (c == '"') return Json(parseString());
        if (consume("true")) return Json(true);
        if (consume("false")) return Json(false);
        if (consume("null")) return Json();
        if (consume("NaN")) return Json(std::nan(""));
        if (consume("Infinity")) return Json(HUGE_VAL);
        if (consume("-Infinity")) return Json(-HUGE_VAL);
        return parseNumber();
    }

    Json parseNumber()
    {
        const char* begin = m_text.c_str() + m_position;
        char* end = nullptr;
        const double value = std::strtod(begin, &end);
        if (end == begin) fail("Expected a value");
        m_position += static_cast<std::size_t>(end - begin);
        return Json(value);
    }

    static void appendUtf8(std::string& out, unsigned code)
    {
        if (code < 0x80) out += static_cast<char>(code);
        else if (code < 0x800) {
            out += static_cast<char>(0xC0 | (code >> 6));
            out += static_cast<char>(0x80 | (code & 0x3F));
        } else if (code < 0x10000) {
            out += static_cast<char>(0xE0 | (code >> 12));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (code & 0x3F));
        } else {
            out += static_cast<char>(0xF0 | (code >> 18));
            out += static_cast<char>(0x80 | ((code >> 12) & 0x3F));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (code & 0x3F));
        }
    }

    unsigned parseHex()
    {
        if (m_position + 4 > m_text.size()) fail("Truncated unicode escape");
        unsigned code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = m_text[m_position++];
            code <<= 4;
            if (c >= '0' && c <= '9') code |= static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f') code |= static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') code |= static_cast<unsigned>(c - 'A' + 10);
            else fail("Invalid unicode escape");
        }
        return code;
    }

    std::string parseString()
    {
        ++m_position;
        std::string out;
        while (true) {
            if (m_position >= m_text.size()) fail("Unterminated string");
            const char c = m_text[m_position++];
            if (c == '"') return out;
            if (c != '\\') { out += c; continue; }
            if (m_position >= m_text.size()) fail("Unterminated escape");
            const char escape = m_text[m_position++];
            switch (escape) {
            case '"': out += '"'; break;
            case '\\': out += '\\'; break;
            case '/': out += '/'; break;
            case 'b': out += '\b'; break;
            case 'f': out += '\f'; break;
            case 'n': out += '\n'; break;
            case 'r': out += '\r'; break;
            case 't': out += '\t'; break;
            case 'u': {
                unsigned code = parseHex();
                if (code >= 0xD800 && code < 0xDC00 && m_text.compare(m_position, 2, "\\u") == 0) {
                    m_position += 2;
                    const unsigned low = parseHex();
                    code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                }
                appendUtf8(out, code);
                break;
            }
            default: fail("Invalid escape");
            }
        }
    }

    Json parseArray(int depth)
    {
        ++m_position;
        Json result = Json::array();
        skipSpace();
        if (m_position < m_text.size() && m_text[m_position] == ']') { ++m_position; return result; }
        while (true) {
            result.push(parseValue(depth + 1));
            skipSpace();
            if (m_position >= m_text.size()) fail("Unterminated array");
            const char c = m_text[m_position++];
            if (c == ']') return result;
            if (c != ',') fail("Expected ',' or ']'");
        }
    }

    Json parseObject(int depth)
    {
        ++m_position;
        Json result = Json::object();
        skipSpace();
        if (m_position < m_text.size() && m_text[m_position] == '}') { ++m_position; return result; }
        while (true) {
            skipSpace();
            if (m_position >= m_text.size() || m_text[m_position] != '"') fail("Expected a member name");
            const std::string key = parseString();
            skipSpace();
            if (m_position >= m_text.size() || m_text[m_position] != ':') fail("Expected ':'");
            ++m_position;
            result[key] = parseValue(depth + 1);
            skipSpace();
            if (m_position >= m_text.size()) fail("Unterminated object");
            const char c = m_text[m_position++];
            if (c == '}') return result;
            if (c != ',') fail("Expected ',' or '}'");
        }
    }

    const std::string& m_text;
    std::size_t m_position = 0;
};

void writeString(std::string& out, const std::string& value)
{
    out += '"';
    for (unsigned char c : value) {
        switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
            if (c < 0x20) {
                char buffer[8];
                std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
                out += buffer;
            } else out += static_cast<char>(c);
        }
    }
    out += '"';
}

void newline(std::string& out, int indent, int depth)
{
    if (indent < 0) return;
    out += '\n';
    out.append(static_cast<std::size_t>(indent * depth), ' ');
}
}

Json Json::array(std::vector<Json> items)
{
    Json value;
    value.m_type = Type::Array;
    value.m_items = std::move(items);
    return value;
}

Json Json::object()
{
    Json value;
    value.m_type = Type::Object;
    return value;
}

Json Json::parse(const std::string& text)
{
    return Parser(text).document();
}

bool Json::boolean() const
{
    if (m_type != Type::Bool) throw std::runtime_error("Expected a boolean, not a string or number");
    return m_bool;
}

double Json::number() const
{
    if (m_type != Type::Number) throw std::runtime_error("Expected a number");
    return m_number;
}

const std::string& Json::string() const
{
    if (m_type != Type::String) throw std::runtime_error("Expected a string");
    return m_string;
}

const std::vector<Json>& Json::items() const
{
    if (m_type != Type::Array) throw std::runtime_error("Expected an array");
    return m_items;
}

std::vector<Json>& Json::items()
{
    if (m_type != Type::Array) throw std::runtime_error("Expected an array");
    return m_items;
}

const std::vector<std::pair<std::string, Json>>& Json::members() const
{
    if (m_type != Type::Object) throw std::runtime_error("Expected an object");
    return m_members;
}

const Json* Json::find(const std::string& key) const
{
    if (m_type != Type::Object) return nullptr;
    for (const auto& member : m_members)
        if (member.first == key) return &member.second;
    return nullptr;
}

const Json& Json::at(const std::string& key) const
{
    if (const Json* value = find(key)) return *value;
    throw std::runtime_error("Missing JSON field: " + key);
}

Json& Json::operator[](const std::string& key)
{
    if (m_type == Type::Null) m_type = Type::Object;
    if (m_type != Type::Object) throw std::runtime_error("Expected an object");
    for (auto& member : m_members)
        if (member.first == key) return member.second;
    m_members.emplace_back(key, Json());
    return m_members.back().second;
}

void Json::push(Json value)
{
    if (m_type == Type::Null) m_type = Type::Array;
    if (m_type != Type::Array) throw std::runtime_error("Expected an array");
    m_items.push_back(std::move(value));
}

std::size_t Json::size() const
{
    if (m_type == Type::Array) return m_items.size();
    if (m_type == Type::Object) return m_members.size();
    return 0;
}

std::string Json::dump(int indent) const
{
    std::string out;
    write(out, indent, 0);
    return out;
}

void Json::write(std::string& out, int indent, int depth) const
{
    switch (m_type) {
    case Type::Null: out += "null"; break;
    case Type::Bool: out += m_bool ? "true" : "false"; break;
    case Type::Number: {
        if (!std::isfinite(m_number)) { out += "null"; break; }
        char buffer[32];
        if (m_number == std::floor(m_number) && std::abs(m_number) < 1e15)
            std::snprintf(buffer, sizeof(buffer), "%.0f", m_number);
        else
            std::snprintf(buffer, sizeof(buffer), "%.17g", m_number);
        out += buffer;
        break;
    }
    case Type::String: writeString(out, m_string); break;
    case Type::Array: {
        out += '[';
        // Keep numeric rows on one line so large arrays stay readable.
        bool scalar = true;
        for (const auto& item : m_items) scalar = scalar && !item.isArray() && !item.isObject();
        for (std::size_t i = 0; i < m_items.size(); ++i) {
            if (i) out += scalar && indent >= 0 ? ", " : ",";
            if (!scalar) newline(out, indent, depth + 1);
            m_items[i].write(out, indent, depth + 1);
        }
        if (!scalar && !m_items.empty()) newline(out, indent, depth);
        out += ']';
        break;
    }
    case Type::Object: {
        out += '{';
        for (std::size_t i = 0; i < m_members.size(); ++i) {
            if (i) out += ',';
            newline(out, indent, depth + 1);
            writeString(out, m_members[i].first);
            out += indent >= 0 ? ": " : ":";
            m_members[i].second.write(out, indent, depth + 1);
        }
        if (!m_members.empty()) newline(out, indent, depth);
        out += '}';
        break;
    }
    }
}
}
