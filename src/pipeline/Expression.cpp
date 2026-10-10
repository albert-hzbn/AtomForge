#include "pipeline/Expression.h"

#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace atomforge::pipeline
{
// A value is a number or, for element comparisons, a string.
struct Value
{
    double number = 0;
    std::string text;
    bool isText = false;
    Value() = default;
    Value(double value) : number(value) {}
    static Value of(std::string name) { Value v; v.text = std::move(name); v.isText = true; return v; }
};

struct Expression::Node
{
    std::function<Value(const AtomVariables&)> evaluate;
};

namespace
{
using NodePtr = std::unique_ptr<Expression::Node>;

NodePtr node(std::function<Value(const AtomVariables&)> f)
{
    auto n = std::make_unique<Expression::Node>();
    n->evaluate = std::move(f);
    return n;
}

double numberOf(const Value& v, const char* context)
{
    if (v.isText) throw std::runtime_error(std::string("'") + v.text + "' is not a number (" + context + ")");
    return v.number;
}

class Parser
{
public:
    explicit Parser(const std::string& text) : m_text(text) {}

    NodePtr parse()
    {
        NodePtr result = parseOr();
        skip();
        if (m_pos != m_text.size()) fail("unexpected '" + m_text.substr(m_pos, 12) + "'");
        return result;
    }

private:
    void skip() { while (m_pos < m_text.size() && std::isspace(static_cast<unsigned char>(m_text[m_pos]))) ++m_pos; }

    bool accept(const std::string& token)
    {
        skip();
        if (m_text.compare(m_pos, token.size(), token) != 0) return false;
        // Words must not run into identifier characters (e.g. "and" vs "andx").
        if (std::isalpha(static_cast<unsigned char>(token[0])) && m_pos + token.size() < m_text.size() &&
            (std::isalnum(static_cast<unsigned char>(m_text[m_pos + token.size()])) || m_text[m_pos + token.size()] == '_'))
            return false;
        m_pos += token.size();
        return true;
    }

    [[noreturn]] void fail(const std::string& message) const
    {
        throw std::runtime_error("Expression error at position " + std::to_string(m_pos + 1) + ": " + message);
    }

    NodePtr parseOr()
    {
        NodePtr left = parseAnd();
        while (accept("||") || accept("or")) {
            std::shared_ptr<Expression::Node> a(left.release()), b(parseAnd().release());
            left = node([a, b](const AtomVariables& v) { return Value{numberOf(a->evaluate(v), "or") != 0 || numberOf(b->evaluate(v), "or") != 0 ? 1.0 : 0.0}; });
        }
        return left;
    }

    NodePtr parseAnd()
    {
        NodePtr left = parseNot();
        while (accept("&&") || accept("and")) {
            std::shared_ptr<Expression::Node> a(left.release()), b(parseNot().release());
            left = node([a, b](const AtomVariables& v) { return Value{numberOf(a->evaluate(v), "and") != 0 && numberOf(b->evaluate(v), "and") != 0 ? 1.0 : 0.0}; });
        }
        return left;
    }

    NodePtr parseNot()
    {
        skip();
        if ((m_pos < m_text.size() && m_text[m_pos] == '!' && (m_pos + 1 >= m_text.size() || m_text[m_pos + 1] != '=')) || accept("not")) {
            if (m_text[m_pos] == '!') ++m_pos;
            std::shared_ptr<Expression::Node> a(parseNot().release());
            return node([a](const AtomVariables& v) { return Value{numberOf(a->evaluate(v), "not") == 0 ? 1.0 : 0.0}; });
        }
        return parseComparison();
    }

    NodePtr parseComparison()
    {
        NodePtr left = parseSum();
        for (const char* op : {"<=", ">=", "==", "!=", "<", ">"}) {
            if (!accept(op)) continue;
            std::shared_ptr<Expression::Node> a(left.release()), b(parseSum().release());
            const std::string o = op;
            return node([a, b, o](const AtomVariables& v) {
                const Value x = a->evaluate(v), y = b->evaluate(v);
                if (x.isText || y.isText) {
                    if (o != "==" && o != "!=") throw std::runtime_error("Element names can only be compared with == or !=");
                    const bool equal = x.isText == y.isText && x.text == y.text;
                    return Value{(o == "==") == equal ? 1.0 : 0.0};
                }
                bool r = false;
                if (o == "<") r = x.number < y.number;
                else if (o == "<=") r = x.number <= y.number;
                else if (o == ">") r = x.number > y.number;
                else if (o == ">=") r = x.number >= y.number;
                else if (o == "==") r = x.number == y.number;
                else r = x.number != y.number;
                return Value{r ? 1.0 : 0.0};
            });
        }
        return left;
    }

    NodePtr parseSum()
    {
        NodePtr left = parseProduct();
        while (true) {
            skip();
            if (m_pos >= m_text.size() || (m_text[m_pos] != '+' && m_text[m_pos] != '-')) return left;
            const char op = m_text[m_pos++];
            std::shared_ptr<Expression::Node> a(left.release()), b(parseProduct().release());
            left = node([a, b, op](const AtomVariables& v) {
                const double x = numberOf(a->evaluate(v), "+/-"), y = numberOf(b->evaluate(v), "+/-");
                return Value{op == '+' ? x + y : x - y};
            });
        }
    }

    NodePtr parseProduct()
    {
        NodePtr left = parseUnary();
        while (true) {
            skip();
            if (m_pos >= m_text.size() || (m_text[m_pos] != '*' && m_text[m_pos] != '/')) return left;
            const char op = m_text[m_pos++];
            std::shared_ptr<Expression::Node> a(left.release()), b(parseUnary().release());
            left = node([a, b, op](const AtomVariables& v) {
                const double x = numberOf(a->evaluate(v), "* /"), y = numberOf(b->evaluate(v), "* /");
                return Value{op == '*' ? x * y : x / y};
            });
        }
    }

    NodePtr parseUnary()
    {
        skip();
        if (m_pos < m_text.size() && (m_text[m_pos] == '-' || m_text[m_pos] == '+')) {
            const char op = m_text[m_pos++];
            std::shared_ptr<Expression::Node> a(parseUnary().release());
            return node([a, op](const AtomVariables& v) { const double x = numberOf(a->evaluate(v), "sign"); return Value{op == '-' ? -x : x}; });
        }
        return parsePower();
    }

    NodePtr parsePower()
    {
        NodePtr base = parsePrimary();
        skip();
        if (m_pos < m_text.size() && m_text[m_pos] == '^') {
            ++m_pos;
            std::shared_ptr<Expression::Node> a(base.release()), b(parseUnary().release());
            return node([a, b](const AtomVariables& v) { return Value{std::pow(numberOf(a->evaluate(v), "^"), numberOf(b->evaluate(v), "^"))}; });
        }
        return base;
    }

    NodePtr parsePrimary()
    {
        skip();
        if (m_pos >= m_text.size()) fail("expression ends early");
        const char c = m_text[m_pos];
        if (c == '(') {
            ++m_pos;
            NodePtr inner = parseOr();
            if (!accept(")")) fail("missing ')'");
            return inner;
        }
        if (c == '"' || c == '\'') {
            const std::size_t end = m_text.find(c, m_pos + 1);
            if (end == std::string::npos) fail("unterminated string");
            const Value text = Value::of(m_text.substr(m_pos + 1, end - m_pos - 1));
            m_pos = end + 1;
            return node([text](const AtomVariables&) { return text; });
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            std::size_t used = 0;
            double value = 0;
            try { value = std::stod(m_text.substr(m_pos), &used); } catch (const std::exception&) { fail("invalid number"); }
            m_pos += used;
            return node([value](const AtomVariables&) { return Value{value}; });
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            std::size_t end = m_pos;
            while (end < m_text.size() && (std::isalnum(static_cast<unsigned char>(m_text[end])) || m_text[end] == '_')) ++end;
            const std::string name = m_text.substr(m_pos, end - m_pos);
            m_pos = end;
            skip();
            if (m_pos < m_text.size() && m_text[m_pos] == '(') return parseFunction(name);
            return variable(name);
        }
        fail(std::string("unexpected '") + c + "'");
    }

    NodePtr parseFunction(const std::string& name)
    {
        ++m_pos;  // '('
        std::vector<std::shared_ptr<Expression::Node>> arguments;
        if (!accept(")")) {
            do { arguments.emplace_back(parseOr().release()); } while (accept(","));
            if (!accept(")")) fail("missing ')' after the arguments of " + name);
        }
        static const std::vector<std::pair<std::string, double (*)(double)>> unary = {
            {"abs", [](double x) { return std::abs(x); }}, {"sqrt", [](double x) { return std::sqrt(x); }},
            {"exp", [](double x) { return std::exp(x); }}, {"log", [](double x) { return std::log(x); }},
            {"sin", [](double x) { return std::sin(x); }}, {"cos", [](double x) { return std::cos(x); }},
            {"tan", [](double x) { return std::tan(x); }}, {"floor", [](double x) { return std::floor(x); }},
            {"ceil", [](double x) { return std::ceil(x); }}};
        for (const auto& [function, f] : unary)
            if (name == function) {
                if (arguments.size() != 1) fail(name + " takes one argument");
                auto a = arguments[0];
                return node([a, f](const AtomVariables& v) { return Value{f(numberOf(a->evaluate(v), "function"))}; });
            }
        if (name == "min" || name == "max") {
            if (arguments.size() != 2) fail(name + " takes two arguments");
            auto a = arguments[0], b = arguments[1];
            const bool isMin = name == "min";
            return node([a, b, isMin](const AtomVariables& v) {
                const double x = numberOf(a->evaluate(v), "min/max"), y = numberOf(b->evaluate(v), "min/max");
                return Value{isMin ? std::min(x, y) : std::max(x, y)};
            });
        }
        fail("unknown function " + name);
    }

    NodePtr variable(const std::string& name)
    {
        using V = AtomVariables;
        static const std::vector<std::pair<std::string, double V::*>> numbers = {
            {"x", &V::x}, {"y", &V::y}, {"z", &V::z}, {"fx", &V::fx}, {"fy", &V::fy}, {"fz", &V::fz},
            {"index", &V::index}, {"count", &V::count}, {"property", &V::property}, {"selected", &V::selected},
            {"Z", &V::atomicNumber}, {"atomic_number", &V::atomicNumber}};
        for (const auto& [n, member] : numbers)
            if (name == n) return node([member](const AtomVariables& v) { return Value{v.*member}; });
        if (name == "element" || name == "type") return node([](const AtomVariables& v) { return Value::of(v.element); });
        if (name == "pi") return node([](const AtomVariables&) { return Value{std::acos(-1.0)}; });
        // Anything else is an element name, as in element == Cu.
        const Value text = Value::of(name);
        return node([text](const AtomVariables&) { return text; });
    }

    const std::string& m_text;
    std::size_t m_pos = 0;
};
}

Expression::Expression(const std::string& text)
{
    if (text.find_first_not_of(" \t\r\n") == std::string::npos) throw std::runtime_error("The selection expression is empty");
    m_root = Parser(text).parse();
}

Expression::~Expression() = default;
Expression::Expression(Expression&&) noexcept = default;
Expression& Expression::operator=(Expression&&) noexcept = default;

double Expression::evaluate(const AtomVariables& variables) const
{
    const Value value = m_root->evaluate(variables);
    if (value.isText) throw std::runtime_error("The expression gives the text '" + value.text + "', not a condition");
    return value.number;
}
}
