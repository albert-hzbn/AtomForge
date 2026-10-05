#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace atomforge::science
{
// Small insertion-ordered JSON value for scientific requests and results.
// Non-finite numbers are written as null, matching strict JSON output.
class Json
{
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Json() = default;
    Json(std::nullptr_t) {}
    Json(bool value) : m_type(Type::Bool), m_bool(value) {}
    Json(int value) : m_type(Type::Number), m_number(value) {}
    Json(long long value) : m_type(Type::Number), m_number(static_cast<double>(value)) {}
    Json(std::size_t value) : m_type(Type::Number), m_number(static_cast<double>(value)) {}
    Json(double value) : m_type(Type::Number), m_number(value) {}
    Json(const char* value) : m_type(Type::String), m_string(value) {}
    Json(std::string value) : m_type(Type::String), m_string(std::move(value)) {}

    static Json array(std::vector<Json> items = {});
    static Json object();
    static Json parse(const std::string& text);

    Type type() const { return m_type; }
    bool isNull() const { return m_type == Type::Null; }
    bool isBool() const { return m_type == Type::Bool; }
    bool isNumber() const { return m_type == Type::Number; }
    bool isString() const { return m_type == Type::String; }
    bool isArray() const { return m_type == Type::Array; }
    bool isObject() const { return m_type == Type::Object; }

    bool boolean() const;
    double number() const;
    const std::string& string() const;
    const std::vector<Json>& items() const;
    std::vector<Json>& items();
    const std::vector<std::pair<std::string, Json>>& members() const;

    const Json* find(const std::string& key) const;
    bool contains(const std::string& key) const { return find(key) != nullptr; }
    const Json& at(const std::string& key) const;
    // Inserts a null member when the key is absent.
    Json& operator[](const std::string& key);
    void push(Json value);
    std::size_t size() const;

    std::string dump(int indent = -1) const;

private:
    void write(std::string& out, int indent, int depth) const;

    Type m_type = Type::Null;
    bool m_bool = false;
    double m_number = 0.0;
    std::string m_string;
    std::vector<Json> m_items;
    std::vector<std::pair<std::string, Json>> m_members;
};
}
