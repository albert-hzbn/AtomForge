#include "algorithms/ShapeModel.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <set>

namespace atomforge::shape
{
using glm::dmat3;
using glm::dvec2;
using glm::dvec3;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kDegree = kPi / 180.0;
constexpr double kHuge = 1e30;

// ---------------------------------------------------------------- values
struct Value
{
    enum class Kind { Number, Vector, Text, Boolean } kind = Kind::Number;
    double number = 0;
    std::vector<Value> items;
    std::string text;

    static Value of(double x) { Value v; v.number = x; return v; }
    static Value list(std::vector<Value> items)
    {
        Value v;
        v.kind = Kind::Vector;
        v.items = std::move(items);
        return v;
    }
    static Value string(std::string s)
    {
        Value v;
        v.kind = Kind::Text;
        v.text = std::move(s);
        return v;
    }
    static Value boolean(bool b)
    {
        Value v;
        v.kind = Kind::Boolean;
        v.number = b ? 1 : 0;
        return v;
    }
};

std::string describe(const Value& v)
{
    switch (v.kind) {
    case Value::Kind::Number: return "a number";
    case Value::Kind::Vector: return "a vector";
    case Value::Kind::Text: return "text";
    case Value::Kind::Boolean: return "true/false";
    }
    return "a value";
}

// ---------------------------------------------------------------- lexer
struct Token
{
    enum class Kind { Identifier, Number, Text, Symbol, End } kind = Kind::End;
    std::string text;
    double number = 0;
    int line = 1;
};

std::vector<Token> tokenize(const std::string& source)
{
    std::vector<Token> out;
    int line = 1;
    std::size_t i = 0;
    while (i < source.size()) {
        const char c = source[i];
        if (c == '\n') { ++line; ++i; continue; }
        if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
        if (c == '#' || (c == '/' && i + 1 < source.size() && source[i + 1] == '/')) {
            while (i < source.size() && source[i] != '\n') ++i;
            continue;
        }
        if (c == '/' && i + 1 < source.size() && source[i + 1] == '*') {
            i += 2;
            while (i + 1 < source.size() && !(source[i] == '*' && source[i + 1] == '/')) {
                if (source[i] == '\n') ++line;
                ++i;
            }
            i += 2;
            continue;
        }
        Token t;
        t.line = line;
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            std::size_t j = i;
            while (j < source.size() && (std::isalnum(static_cast<unsigned char>(source[j])) || source[j] == '_')) ++j;
            t.kind = Token::Kind::Identifier;
            t.text = source.substr(i, j - i);
            i = j;
        } else if (std::isdigit(static_cast<unsigned char>(c)) || (c == '.' && i + 1 < source.size() && std::isdigit(static_cast<unsigned char>(source[i + 1])))) {
            char* end = nullptr;
            t.kind = Token::Kind::Number;
            t.number = std::strtod(source.c_str() + i, &end);
            t.text = source.substr(i, static_cast<std::size_t>(end - (source.c_str() + i)));
            i = static_cast<std::size_t>(end - source.c_str());
        } else if (c == '"' || c == '\'') {
            std::size_t j = i + 1;
            while (j < source.size() && source[j] != c && source[j] != '\n') ++j;
            if (j >= source.size() || source[j] != c) throw ShapeError(line, "unterminated text");
            t.kind = Token::Kind::Text;
            t.text = source.substr(i + 1, j - i - 1);
            i = j + 1;
        } else if (std::string("(){}[],;=+-*/^").find(c) != std::string::npos) {
            t.kind = Token::Kind::Symbol;
            t.text = std::string(1, c);
            ++i;
        } else
            throw ShapeError(line, std::string("unexpected character '") + c + "'");
        out.push_back(t);
    }
    Token end;
    end.line = line;
    out.push_back(end);
    return out;
}
}  // namespace

// ---------------------------------------------------------------- shape tree
struct Node
{
    enum class Kind {
        Sphere, Box, Cylinder, Cone, Torus, Capsule, Ellipsoid, Polygon, HalfSpace, Slab, Facets,
        Union, Intersection, Difference, SmoothUnion, SmoothIntersection, SmoothDifference,
        Transform, Shell, Offset, LinearArray, PolarArray, GridArray, Region
    };
    Kind kind = Kind::Union;
    int line = 0;
    dvec3 size{0.0};                 // box size, ellipsoid radii
    double r = 0, r2 = 0, h = 0, k = 0;
    std::vector<dvec2> polygon;      // extruded polygon (prisms too)
    std::vector<dvec3> normals;      // half-spaces, slabs, facets (unit)
    std::vector<double> offsets;     // their distances from the origin
    dmat3 forward{1.0}, inverse{1.0};
    dvec3 translation{0.0};
    double distanceScale = 1;        // transform: lower bound of the stretch
    dvec3 step{0.0};
    glm::ivec3 counts{1, 1, 1};
    dvec3 axis{0, 0, 1};
    double angle = 360;
    std::string element;
    std::vector<Node> children;
};

namespace
{
double maxComponent(const dvec3& v) { return std::max(v.x, std::max(v.y, v.z)); }

// The right-handed turn by `angle` (radians) about a unit axis (glm: m[column][row]).
dmat3 rotation(const dvec3& axis, double angle)
{
    const double c = std::cos(angle), s = std::sin(angle);
    dmat3 r(c);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r[j][i] += (1 - c) * axis[i] * axis[j];
    r[1][0] -= s * axis.z; r[0][1] += s * axis.z;
    r[2][0] += s * axis.y; r[0][2] -= s * axis.y;
    r[2][1] -= s * axis.x; r[1][2] += s * axis.x;
    return r;
}

double polygonDistance(const std::vector<dvec2>& v, const dvec2& p)
{
    double d = glm::dot(p - v[0], p - v[0]);
    double s = 1;
    for (std::size_t i = 0, j = v.size() - 1; i < v.size(); j = i, ++i) {
        const dvec2 e = v[j] - v[i], w = p - v[i];
        const dvec2 b = w - e * std::clamp(glm::dot(w, e) / glm::dot(e, e), 0.0, 1.0);
        d = std::min(d, glm::dot(b, b));
        const bool c1 = p.y >= v[i].y, c2 = p.y < v[j].y, c3 = e.x * w.y > e.y * w.x;
        if ((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s = -s;
    }
    return s * std::sqrt(d);
}

double smoothMin(double a, double b, double k)
{
    if (k <= 0) return std::min(a, b);
    const double h = std::clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return b + (a - b) * h - k * h * (1 - h);
}

double unionOf(const std::vector<Node>& children, const dvec3& p);

// Copies of an array: the offsets (grid arrays: only those near p).
template <typename F>
void forEachCopy(const Node& n, const dvec3& p, F&& visit)
{
    if (n.kind == Node::Kind::LinearArray) {
        for (int i = 0; i < n.counts.x; ++i) visit(p - n.step * double(i));
    } else if (n.kind == Node::Kind::PolarArray) {
        const int copies = n.counts.x;
        const bool full = std::abs(std::abs(n.angle) - 360) < 1e-9;
        for (int i = 0; i < copies; ++i) {
            const double turn = copies <= 1 ? 0 : n.angle * kDegree * i / (full ? copies : copies - 1);
            // Rodrigues: rotate p by -turn about the axis.
            const dvec3 a = n.axis;
            const double c = std::cos(-turn), s = std::sin(-turn);
            visit(p * c + glm::cross(a, p) * s + a * glm::dot(a, p) * (1 - c));
        }
    } else {
        // Grid array: the copies whose cells are next to p (the child should fit in one cell).
        int low[3], high[3];
        for (int a = 0; a < 3; ++a) {
            if (n.counts[a] <= 1 || std::abs(n.step[a]) < 1e-12) { low[a] = 0; high[a] = n.counts[a] - 1; continue; }
            const int nearest = static_cast<int>(std::lround(p[a] / n.step[a]));
            low[a] = std::clamp(nearest - 1, 0, n.counts[a] - 1);
            high[a] = std::clamp(nearest + 1, 0, n.counts[a] - 1);
        }
        for (int i = low[0]; i <= high[0]; ++i)
            for (int j = low[1]; j <= high[1]; ++j)
                for (int k = low[2]; k <= high[2]; ++k) visit(p - dvec3(i * n.step.x, j * n.step.y, k * n.step.z));
    }
}

double distanceOf(const Node& n, const dvec3& p)
{
    switch (n.kind) {
    case Node::Kind::Sphere: return glm::length(p) - n.r;
    case Node::Kind::Box: {
        const dvec3 q = glm::abs(p) - n.size * 0.5;
        return glm::length(glm::max(q, dvec3(0.0))) + std::min(maxComponent(q), 0.0);
    }
    case Node::Kind::Cylinder: {
        const dvec2 d(glm::length(dvec2(p.x, p.y)) - n.r, std::abs(p.z) - n.h / 2);
        return std::min(std::max(d.x, d.y), 0.0) + glm::length(glm::max(d, dvec2(0.0)));
    }
    case Node::Kind::Cone: {
        // Capped cone: radius r at z = -h/2, r2 at z = +h/2.
        const double hh = n.h / 2, r1 = n.r, r2 = n.r2;
        const dvec2 q(glm::length(dvec2(p.x, p.y)), p.z);
        const dvec2 k1(r2, hh), k2(r2 - r1, 2 * hh);
        const dvec2 ca(q.x - std::min(q.x, q.y < 0 ? r1 : r2), std::abs(q.y) - hh);
        const dvec2 cb = q - k1 + k2 * std::clamp(glm::dot(k1 - q, k2) / glm::dot(k2, k2), 0.0, 1.0);
        const double s = (cb.x < 0 && ca.y < 0) ? -1.0 : 1.0;
        return s * std::sqrt(std::min(glm::dot(ca, ca), glm::dot(cb, cb)));
    }
    case Node::Kind::Torus: return glm::length(dvec2(glm::length(dvec2(p.x, p.y)) - n.r, p.z)) - n.r2;
    case Node::Kind::Capsule: {
        dvec3 q = p;
        q.z -= std::clamp(q.z, -n.h / 2, n.h / 2);
        return glm::length(q) - n.r;
    }
    case Node::Kind::Ellipsoid: {
        const double k0 = glm::length(p / n.size);
        if (k0 < 1e-12) return -std::min(n.size.x, std::min(n.size.y, n.size.z));
        const double k1 = glm::length(p / (n.size * n.size));
        return k0 * (k0 - 1) / k1;
    }
    case Node::Kind::Polygon: {
        const dvec2 w(polygonDistance(n.polygon, dvec2(p.x, p.y)), std::abs(p.z) - n.h / 2);
        return std::min(std::max(w.x, w.y), 0.0) + glm::length(glm::max(w, dvec2(0.0)));
    }
    case Node::Kind::HalfSpace: return glm::dot(p, n.normals[0]) - n.offsets[0];
    case Node::Kind::Slab: return std::abs(glm::dot(p, n.normals[0]) - n.offsets[0]) - n.h / 2;
    case Node::Kind::Facets: {
        double d = -kHuge;
        for (std::size_t i = 0; i < n.normals.size(); ++i) d = std::max(d, glm::dot(p, n.normals[i]) - n.offsets[i]);
        return d;
    }
    case Node::Kind::Union: return unionOf(n.children, p);
    case Node::Kind::Region: return kHuge;  // regions label atoms of the shape; they add no volume
    case Node::Kind::Intersection: {
        double d = -kHuge;
        for (const Node& c : n.children) d = std::max(d, distanceOf(c, p));
        return n.children.empty() ? kHuge : d;
    }
    case Node::Kind::Difference: {
        if (n.children.empty()) return kHuge;
        double d = distanceOf(n.children[0], p);
        for (std::size_t i = 1; i < n.children.size(); ++i) d = std::max(d, -distanceOf(n.children[i], p));
        return d;
    }
    case Node::Kind::SmoothUnion: {
        double d = kHuge;
        for (const Node& c : n.children) d = d >= kHuge ? distanceOf(c, p) : smoothMin(d, distanceOf(c, p), n.k);
        return d;
    }
    case Node::Kind::SmoothIntersection: {
        double d = kHuge;
        for (const Node& c : n.children) d = d >= kHuge ? -distanceOf(c, p) : smoothMin(d, -distanceOf(c, p), n.k);
        return n.children.empty() ? kHuge : -d;
    }
    case Node::Kind::SmoothDifference: {
        if (n.children.empty()) return kHuge;
        double d = -distanceOf(n.children[0], p);
        for (std::size_t i = 1; i < n.children.size(); ++i) d = smoothMin(d, distanceOf(n.children[i], p), n.k);
        return -d;
    }
    case Node::Kind::Transform: return unionOf(n.children, n.inverse * (p - n.translation)) * n.distanceScale;
    case Node::Kind::Shell: {
        const double d = unionOf(n.children, p);
        return std::max(d, -(d + n.h));
    }
    case Node::Kind::Offset: return unionOf(n.children, p) - n.r;
    case Node::Kind::LinearArray:
    case Node::Kind::PolarArray:
    case Node::Kind::GridArray: {
        double d = kHuge;
        forEachCopy(n, p, [&](const dvec3& q) { d = std::min(d, unionOf(n.children, q)); });
        return d;
    }
    }
    return kHuge;
}

double unionOf(const std::vector<Node>& children, const dvec3& p)
{
    double d = kHuge;
    for (const Node& c : children) d = std::min(d, distanceOf(c, p));
    return d;
}

void materialOf(const Node& n, const dvec3& p, const std::string*& out)
{
    switch (n.kind) {
    case Node::Kind::Region:
        if (unionOf(n.children, p) <= 0) out = &n.element;
        for (const Node& c : n.children) materialOf(c, p, out);
        return;
    case Node::Kind::Transform:
        for (const Node& c : n.children) materialOf(c, n.inverse * (p - n.translation), out);
        return;
    case Node::Kind::LinearArray:
    case Node::Kind::PolarArray:
    case Node::Kind::GridArray:
        forEachCopy(n, p, [&](const dvec3& q) { for (const Node& c : n.children) materialOf(c, q, out); });
        return;
    default:
        for (const Node& c : n.children) materialOf(c, p, out);
    }
}

// ---------------------------------------------------------------- bounds
Bounds merged(Bounds a, const Bounds& b)
{
    if (a.empty) return b;
    if (b.empty) return a;
    if (!a.finite || !b.finite) return Bounds{};
    a.low = glm::min(a.low, b.low);
    a.high = glm::max(a.high, b.high);
    return a;
}

Bounds box(const dvec3& half)
{
    Bounds b;
    b.low = -half;
    b.high = half;
    b.finite = true;
    return b;
}

Bounds mapped(const Bounds& b, const std::function<dvec3(const dvec3&)>& f)
{
    if (!b.finite || b.empty) return b;
    Bounds out;
    out.finite = true;
    out.low = dvec3(kHuge);
    out.high = dvec3(-kHuge);
    for (int c = 0; c < 8; ++c) {
        const dvec3 corner((c & 1) ? b.high.x : b.low.x, (c & 2) ? b.high.y : b.low.y, (c & 4) ? b.high.z : b.low.z);
        const dvec3 q = f(corner);
        out.low = glm::min(out.low, q);
        out.high = glm::max(out.high, q);
    }
    return out;
}

Bounds boundsOfUnion(const std::vector<Node>& children);

Bounds facetBounds(const Node& n)
{
    // Bounded only if the normals surround the origin: no direction u with n_i . u <= 0 for all i.
    const std::size_t m = n.normals.size();
    for (int t = 0; t < 2000; ++t) {
        // Fibonacci sphere directions.
        const double z = 1 - 2 * (t + 0.5) / 2000, rr = std::sqrt(std::max(0.0, 1 - z * z)), phi = t * 2.399963229728653;
        const dvec3 u(rr * std::cos(phi), rr * std::sin(phi), z);
        bool open = true;
        for (std::size_t i = 0; i < m && open; ++i) open = glm::dot(n.normals[i], u) <= 1e-9;
        if (open) return Bounds{};
    }
    Bounds b;
    b.low = dvec3(kHuge);
    b.high = dvec3(-kHuge);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = i + 1; j < m; ++j)
            for (std::size_t k = j + 1; k < m; ++k) {
                const dmat3 a = glm::transpose(dmat3(n.normals[i], n.normals[j], n.normals[k]));
                if (std::abs(glm::determinant(a)) < 1e-10) continue;
                const dvec3 x = glm::inverse(a) * dvec3(n.offsets[i], n.offsets[j], n.offsets[k]);
                bool inside = true;
                for (std::size_t q = 0; q < m && inside; ++q) inside = glm::dot(n.normals[q], x) <= n.offsets[q] + 1e-7;
                if (!inside) continue;
                b.low = glm::min(b.low, x);
                b.high = glm::max(b.high, x);
                b.finite = true;
            }
    return b;
}

Bounds boundsOf(const Node& n)
{
    switch (n.kind) {
    case Node::Kind::Sphere: return box(dvec3(n.r));
    case Node::Kind::Box: return box(n.size * 0.5);
    case Node::Kind::Cylinder: return box(dvec3(n.r, n.r, n.h / 2));
    case Node::Kind::Cone: return box(dvec3(std::max(n.r, n.r2), std::max(n.r, n.r2), n.h / 2));
    case Node::Kind::Torus: return box(dvec3(n.r + n.r2, n.r + n.r2, n.r2));
    case Node::Kind::Capsule: return box(dvec3(n.r, n.r, n.h / 2 + n.r));
    case Node::Kind::Ellipsoid: return box(n.size);
    case Node::Kind::Polygon: {
        Bounds b = box(dvec3(0.0, 0.0, n.h / 2));
        b.low.x = b.low.y = kHuge;
        b.high.x = b.high.y = -kHuge;
        for (const dvec2& v : n.polygon) {
            b.low.x = std::min(b.low.x, v.x); b.low.y = std::min(b.low.y, v.y);
            b.high.x = std::max(b.high.x, v.x); b.high.y = std::max(b.high.y, v.y);
        }
        return b;
    }
    case Node::Kind::HalfSpace:
    case Node::Kind::Slab: return Bounds{};
    case Node::Kind::Facets: return facetBounds(n);
    case Node::Kind::Union: return boundsOfUnion(n.children);
    case Node::Kind::Region: {
        Bounds none;
        none.finite = none.empty = true;
        return none;
    }
    case Node::Kind::Intersection:
    case Node::Kind::SmoothIntersection: {
        Bounds b;
        for (const Node& c : n.children) {
            const Bounds cb = boundsOf(c);
            if (!cb.finite || cb.empty) continue;
            if (!b.finite) b = cb;
            else {
                b.low = glm::max(b.low, cb.low);
                b.high = glm::min(b.high, cb.high);
            }
        }
        if (b.finite) b.high = glm::max(b.high, b.low);
        return b;
    }
    case Node::Kind::Difference:
    case Node::Kind::SmoothDifference: return n.children.empty() ? Bounds{} : boundsOf(n.children[0]);
    case Node::Kind::SmoothUnion: {
        Bounds b = boundsOfUnion(n.children);
        if (b.finite) { b.low -= dvec3(n.k); b.high += dvec3(n.k); }
        return b;
    }
    case Node::Kind::Transform:
        return mapped(boundsOfUnion(n.children), [&](const dvec3& q) { return n.forward * q + n.translation; });
    case Node::Kind::Shell: return boundsOfUnion(n.children);
    case Node::Kind::Offset: {
        Bounds b = boundsOfUnion(n.children);
        if (b.finite && n.r > 0) { b.low -= dvec3(n.r); b.high += dvec3(n.r); }
        return b;
    }
    case Node::Kind::LinearArray: {
        const Bounds b = boundsOfUnion(n.children);
        return merged(b, mapped(b, [&](const dvec3& q) { return q + n.step * double(n.counts.x - 1); }));
    }
    case Node::Kind::GridArray: {
        const Bounds b = boundsOfUnion(n.children);
        const dvec3 last(n.step.x * (n.counts.x - 1), n.step.y * (n.counts.y - 1), n.step.z * (n.counts.z - 1));
        return merged(b, mapped(b, [&](const dvec3& q) { return q + last; }));
    }
    case Node::Kind::PolarArray: {
        const Bounds child = boundsOfUnion(n.children);
        if (!child.finite) return child;
        Bounds b = child;
        const bool full = std::abs(std::abs(n.angle) - 360) < 1e-9;
        for (int i = 1; i < n.counts.x; ++i) {
            const double turn = n.angle * kDegree * i / (full ? n.counts.x : n.counts.x - 1);
            const double c = std::cos(turn), s = std::sin(turn);
            const dvec3 a = n.axis;
            b = merged(b, mapped(child, [&](const dvec3& p) { return p * c + glm::cross(a, p) * s + a * glm::dot(a, p) * (1 - c); }));
        }
        return b;
    }
    }
    return Bounds{};
}

Bounds boundsOfUnion(const std::vector<Node>& children)
{
    Bounds b;
    b.finite = b.empty = true;
    for (const Node& c : children) b = merged(b, boundsOf(c));
    if (children.empty()) b = Bounds{};
    return b;
}

// Whether a node adds volume (anything but regions).
bool hasVolume(const Node& n)
{
    if (n.kind == Node::Kind::Region) return false;
    if (n.kind == Node::Kind::Union || n.kind == Node::Kind::Transform || n.kind == Node::Kind::LinearArray ||
        n.kind == Node::Kind::PolarArray || n.kind == Node::Kind::GridArray)
        return std::any_of(n.children.begin(), n.children.end(), hasVolume);
    return true;
}

void collectRegions(const Node& n, std::vector<std::string>& out)
{
    if (n.kind == Node::Kind::Region && std::find(out.begin(), out.end(), n.element) == out.end()) out.push_back(n.element);
    for (const Node& c : n.children) collectRegions(c, out);
}

// ---------------------------------------------------------------- parser
class Parser
{
public:
    Parser(const std::string& source, const dmat3* cell) : m_tokens(tokenize(source)), m_cell(cell)
    {
        m_variables["pi"] = Value::of(kPi);
        m_variables["PI"] = Value::of(kPi);
        m_variables["true"] = Value::boolean(true);
        m_variables["false"] = Value::boolean(false);
    }

    Node program(int& primitives)
    {
        Node root;
        root.kind = Node::Kind::Union;
        while (peek().kind != Token::Kind::End) statement(root.children);
        primitives = m_primitives;
        return root;
    }

private:
    std::vector<Token> m_tokens;
    std::size_t m_at = 0;
    const dmat3* m_cell;
    std::map<std::string, Value> m_variables;
    int m_primitives = 0;

    const Token& peek(std::size_t ahead = 0) const { return m_tokens[std::min(m_at + ahead, m_tokens.size() - 1)]; }
    Token next() { return m_tokens[std::min(m_at++, m_tokens.size() - 1)]; }
    bool symbol(const char* s) const { return peek().kind == Token::Kind::Symbol && peek().text == s; }
    void expect(const char* s)
    {
        if (!symbol(s)) throw ShapeError(peek().line, std::string("expected '") + s + "'" + found());
        ++m_at;
    }
    std::string found() const
    {
        const Token& t = peek();
        if (t.kind == Token::Kind::End) return " at the end of the script";
        return " before '" + t.text + "'";
    }

    // statement: name = expr ;  |  command ( args ) ( { statements } | statement | ; )
    void statement(std::vector<Node>& out)
    {
        if (symbol(";")) { ++m_at; return; }
        const Token name = next();
        if (name.kind != Token::Kind::Identifier) throw ShapeError(name.line, "expected a command or a variable, found '" + name.text + "'");
        if (symbol("=")) {
            ++m_at;
            m_variables[name.text] = expression();
            if (symbol(";")) ++m_at;
            return;
        }
        expect("(");
        std::vector<std::pair<std::string, Value>> arguments;
        while (!symbol(")")) {
            if (peek().kind == Token::Kind::End) throw ShapeError(name.line, "'" + name.text + "(' is not closed with ')'");
            std::string key;
            if (peek().kind == Token::Kind::Identifier && peek(1).kind == Token::Kind::Symbol && peek(1).text == "=") {
                key = next().text;
                ++m_at;
            }
            arguments.push_back({key, expression()});
            if (peek().kind == Token::Kind::End) throw ShapeError(name.line, "'" + name.text + "(' is not closed with ')'");
            if (!symbol(")") && !symbol(",")) throw ShapeError(peek().line, "expected ',' or ')' in " + name.text + "()" + found());
            if (symbol(",")) ++m_at;
        }
        expect(")");
        Node node = make(name, arguments);
        if (isOperation(node.kind)) {
            if (symbol("{")) {
                ++m_at;
                while (!symbol("}")) {
                    if (peek().kind == Token::Kind::End) throw ShapeError(name.line, "'" + name.text + "' block is not closed with '}'");
                    statement(node.children);
                }
                ++m_at;
            } else if (!symbol(";") && peek().kind != Token::Kind::End)
                statement(node.children);  // one child, OpenSCAD-style: translate([..]) sphere(1);
            else if (symbol(";"))
                ++m_at;
            if (node.children.empty()) throw ShapeError(name.line, "'" + name.text + "' needs shapes to act on, in { } or after it");
            finish(node, name);
        } else if (symbol(";"))
            ++m_at;
        out.push_back(std::move(node));
    }

    static bool isOperation(Node::Kind k) { return k >= Node::Kind::Union; }

    // ---------------------------------------------------------------- expressions
    Value expression()
    {
        Value v = term();
        while (symbol("+") || symbol("-")) {
            const Token op = next();
            v = arithmetic(v, term(), op.text[0], op.line);
        }
        return v;
    }
    Value term()
    {
        Value v = unary();
        while (symbol("*") || symbol("/")) {
            const Token op = next();
            v = arithmetic(v, unary(), op.text[0], op.line);
        }
        return v;
    }
    Value unary()
    {
        if (symbol("-")) {
            const int line = next().line;
            return arithmetic(Value::of(-1), unary(), '*', line);
        }
        if (symbol("+")) { ++m_at; return unary(); }
        Value base = primary();
        if (symbol("^")) {
            const int line = next().line;
            const Value exponent = unary();
            return Value::of(std::pow(number(base, line), number(exponent, line)));
        }
        return base;
    }
    Value primary()
    {
        const Token t = next();
        if (t.kind == Token::Kind::Number) return Value::of(t.number);
        if (t.kind == Token::Kind::Text) return Value::string(t.text);
        if (t.kind == Token::Kind::Symbol && t.text == "(") {
            Value v = expression();
            expect(")");
            return v;
        }
        if (t.kind == Token::Kind::Symbol && t.text == "[") {
            std::vector<Value> items;
            while (!symbol("]")) {
                items.push_back(expression());
                if (!symbol("]")) expect(",");
            }
            ++m_at;
            return Value::list(std::move(items));
        }
        if (t.kind == Token::Kind::Identifier) {
            if (symbol("(")) {
                ++m_at;
                std::vector<double> args;
                while (!symbol(")")) {
                    args.push_back(number(expression(), t.line));
                    if (!symbol(")")) expect(",");
                }
                ++m_at;
                return Value::of(function(t, args));
            }
            const auto it = m_variables.find(t.text);
            if (it == m_variables.end()) throw ShapeError(t.line, "unknown variable '" + t.text + "' (set it first: " + t.text + " = 10;)");
            return it->second;
        }
        throw ShapeError(t.line, "expected a value, found '" + t.text + "'");
    }
    static double function(const Token& t, const std::vector<double>& a)
    {
        const auto need = [&](std::size_t n) {
            if (a.size() != n) throw ShapeError(t.line, t.text + "() takes " + std::to_string(n) + " argument" + (n == 1 ? "" : "s"));
        };
        const std::string& f = t.text;
        if (f == "sqrt") { need(1); return std::sqrt(a[0]); }
        if (f == "sin") { need(1); return std::sin(a[0] * kDegree); }
        if (f == "cos") { need(1); return std::cos(a[0] * kDegree); }
        if (f == "tan") { need(1); return std::tan(a[0] * kDegree); }
        if (f == "asin") { need(1); return std::asin(a[0]) / kDegree; }
        if (f == "acos") { need(1); return std::acos(a[0]) / kDegree; }
        if (f == "atan") { need(1); return std::atan(a[0]) / kDegree; }
        if (f == "atan2") { need(2); return std::atan2(a[0], a[1]) / kDegree; }
        if (f == "abs") { need(1); return std::abs(a[0]); }
        if (f == "floor") { need(1); return std::floor(a[0]); }
        if (f == "ceil") { need(1); return std::ceil(a[0]); }
        if (f == "round") { need(1); return std::round(a[0]); }
        if (f == "exp") { need(1); return std::exp(a[0]); }
        if (f == "ln" || f == "log") { need(1); return std::log(a[0]); }
        if (f == "pow") { need(2); return std::pow(a[0], a[1]); }
        if (f == "min" || f == "max") {
            if (a.empty()) throw ShapeError(t.line, f + "() needs values");
            return f == "min" ? *std::min_element(a.begin(), a.end()) : *std::max_element(a.begin(), a.end());
        }
        throw ShapeError(t.line, "unknown function '" + f + "'");
    }
    static double number(const Value& v, int line)
    {
        if (v.kind != Value::Kind::Number && v.kind != Value::Kind::Boolean) throw ShapeError(line, "expected a number, found " + describe(v));
        return v.number;
    }
    static Value arithmetic(const Value& a, const Value& b, char op, int line)
    {
        const bool av = a.kind == Value::Kind::Vector, bv = b.kind == Value::Kind::Vector;
        if (!av && !bv) {
            const double x = number(a, line), y = number(b, line);
            switch (op) {
            case '+': return Value::of(x + y);
            case '-': return Value::of(x - y);
            case '*': return Value::of(x * y);
            default:
                if (y == 0) throw ShapeError(line, "division by zero");
                return Value::of(x / y);
            }
        }
        if (av && bv) {
            if (op != '+' && op != '-') throw ShapeError(line, "vectors can only be added or subtracted");
            if (a.items.size() != b.items.size()) throw ShapeError(line, "vectors of different lengths");
            std::vector<Value> items;
            for (std::size_t i = 0; i < a.items.size(); ++i) items.push_back(arithmetic(a.items[i], b.items[i], op, line));
            return Value::list(std::move(items));
        }
        const Value& vec = av ? a : b;
        const Value& scalar = av ? b : a;
        if (op == '+' || op == '-' || (op == '/' && !av)) throw ShapeError(line, "a vector and a number can only be multiplied or divided");
        std::vector<Value> items;
        for (const Value& item : vec.items) items.push_back(av ? arithmetic(item, scalar, op, line) : arithmetic(scalar, item, op, line));
        return Value::list(std::move(items));
    }

    // ---------------------------------------------------------------- arguments
    struct Arguments
    {
        const Token& command;
        std::vector<std::pair<std::string, Value>> values;
        std::vector<std::string> positional;  // names of positional arguments, in order
        std::set<std::string> used;

        const Value* find(const std::string& name)
        {
            for (std::size_t i = 0, p = 0; i < values.size(); ++i) {
                const std::string& key = values[i].first.empty() ? (p < positional.size() ? positional[p] : std::string()) : values[i].first;
                if (values[i].first.empty()) ++p;
                if (key == name) { used.insert(name); return &values[i].second; }
            }
            return nullptr;
        }
        bool has(const std::string& name) { return find(name) != nullptr; }
        double number(const std::string& name, double fallback)
        {
            const Value* v = find(name);
            if (!v) return fallback;
            if (v->kind != Value::Kind::Number && v->kind != Value::Kind::Boolean)
                throw ShapeError(command.line, command.text + "(): " + name + " must be a number, not " + describe(*v));
            return v->number;
        }
        double positiveNumber(const std::string& name, double fallback)
        {
            const double x = number(name, fallback);
            if (!(x > 0)) throw ShapeError(command.line, command.text + "(): " + name + " must be positive");
            return x;
        }
        dvec3 vector(const std::string& name, dvec3 fallback, bool scalarFills = false)
        {
            const Value* v = find(name);
            if (!v) return fallback;
            if (v->kind == Value::Kind::Number && scalarFills) return dvec3(v->number);
            if (v->kind != Value::Kind::Vector || v->items.size() < 2 || v->items.size() > 3)
                throw ShapeError(command.line, command.text + "(): " + name + " must be a vector [x, y, z]");
            dvec3 out(0.0);
            for (std::size_t i = 0; i < v->items.size(); ++i) out[static_cast<int>(i)] = Parser::number(v->items[i], command.line);
            return out;
        }
        std::string text(const std::string& name, const std::string& fallback)
        {
            const Value* v = find(name);
            if (!v) return fallback;
            if (v->kind != Value::Kind::Text) throw ShapeError(command.line, command.text + "(): " + name + " must be text in quotes");
            return v->text;
        }
        void finish()
        {
            for (std::size_t i = 0, p = 0; i < values.size(); ++i) {
                if (values[i].first.empty()) {
                    if (p >= positional.size())
                        throw ShapeError(command.line, command.text + "() takes at most " + std::to_string(positional.size()) + " unnamed argument(s)");
                    ++p;
                    continue;
                }
                if (!used.count(values[i].first)) throw ShapeError(command.line, command.text + "() has no argument '" + values[i].first + "'");
            }
        }
    };

    // The normal of {hkl} (reciprocal lattice) or a Cartesian direction.
    dvec3 planeNormal(const dvec3& v, bool miller) const
    {
        dvec3 n = v;
        // The columns of the cell are a, b, c; the rows of its inverse are the reciprocal vectors,
        // so the normal h a* + k b* + l c* is inverse^T (h, k, l).
        if (miller && m_cell) n = glm::transpose(glm::inverse(*m_cell)) * v;
        return n;
    }

    Node make(const Token& command, const std::vector<std::pair<std::string, Value>>& values)
    {
        Node n;
        n.line = command.line;
        Arguments a{command, values, {}, {}};
        const std::string& c = command.text;
        const auto radius = [&](const char* name, const char* diameter, double fallback) {
            if (a.has(diameter)) return a.positiveNumber(diameter, 2 * fallback) / 2;
            return a.positiveNumber(name, fallback);
        };
        if (c == "sphere") {
            a.positional = {"r"};
            n.kind = Node::Kind::Sphere;
            n.r = radius("r", "d", 10);
        } else if (c == "box" || c == "cube") {
            a.positional = {"size", "center"};
            n.kind = Node::Kind::Box;
            n.size = a.vector("size", dvec3(10.0), true);
            if (!(n.size.x > 0 && n.size.y > 0 && n.size.z > 0)) throw ShapeError(command.line, c + "(): size must be positive");
            if (!a.number("center", 1)) {
                const dvec3 shift = n.size * 0.5;
                return corner(std::move(n), shift, a);
            }
        } else if (c == "cylinder") {
            a.positional = {"h", "r", "r2", "center"};
            n.h = a.positiveNumber("h", 10);
            if (a.has("r1") || a.has("r2") || a.has("d1") || a.has("d2")) {
                n.kind = Node::Kind::Cone;
                n.r = a.has("d1") ? a.number("d1", 0) / 2 : a.number("r1", a.number("r", 5));
                n.r2 = a.has("d2") ? a.number("d2", 0) / 2 : a.number("r2", a.number("r", 5));
                if (n.r < 0 || n.r2 < 0 || (n.r == 0 && n.r2 == 0)) throw ShapeError(command.line, "cylinder(): radii must not be negative or both zero");
            } else {
                n.kind = Node::Kind::Cylinder;
                n.r = radius("r", "d", 5);
            }
            if (!a.number("center", 1)) {
                const dvec3 shift(0, 0, n.h / 2);
                return corner(std::move(n), shift, a);
            }
        } else if (c == "cone") {
            a.positional = {"r1", "r2", "h"};
            n.kind = Node::Kind::Cone;
            n.r = a.number("r1", 5);
            n.r2 = a.number("r2", 0);
            n.h = a.positiveNumber("h", 10);
            if (n.r < 0 || n.r2 < 0 || (n.r == 0 && n.r2 == 0)) throw ShapeError(command.line, "cone(): radii must not be negative or both zero");
        } else if (c == "torus") {
            a.positional = {"R", "r"};
            n.kind = Node::Kind::Torus;
            n.r = a.positiveNumber("R", 10);
            n.r2 = a.positiveNumber("r", 3);
        } else if (c == "capsule") {
            a.positional = {"r", "h"};
            n.kind = Node::Kind::Capsule;
            n.r = a.positiveNumber("r", 5);
            n.h = a.number("h", 10);
            if (n.h < 0) throw ShapeError(command.line, "capsule(): h must not be negative");
        } else if (c == "ellipsoid") {
            a.positional = {"r"};
            n.kind = Node::Kind::Ellipsoid;
            n.size = a.vector("r", dvec3(10, 6, 4), true);
            if (!(n.size.x > 0 && n.size.y > 0 && n.size.z > 0)) throw ShapeError(command.line, "ellipsoid(): radii must be positive");
        } else if (c == "prism") {
            // A regular polygon of n sides (circumradius r, or apothem a) extruded along z.
            a.positional = {"n", "r", "h"};
            n.kind = Node::Kind::Polygon;
            const int sides = static_cast<int>(a.number("n", 6));
            if (sides < 3 || sides > 1000) throw ShapeError(command.line, "prism(): n must be 3 or more sides");
            double r = a.has("apothem") ? a.positiveNumber("apothem", 1) / std::cos(kPi / sides) : a.positiveNumber("r", 10);
            n.h = a.positiveNumber("h", 20);
            const double start = a.number("rotation", 0) * kDegree;
            for (int i = 0; i < sides; ++i) n.polygon.push_back(r * dvec2(std::cos(start + 2 * kPi * i / sides), std::sin(start + 2 * kPi * i / sides)));
        } else if (c == "polygon" || c == "extrude") {
            a.positional = {"points", "h"};
            n.kind = Node::Kind::Polygon;
            const Value* points = a.find("points");
            if (!points || points->kind != Value::Kind::Vector || points->items.size() < 3)
                throw ShapeError(command.line, c + "(): points must list at least three [x, y] corners");
            for (const Value& p : points->items) {
                if (p.kind != Value::Kind::Vector || p.items.size() != 2) throw ShapeError(command.line, c + "(): each point is [x, y]");
                n.polygon.push_back(dvec2(number(p.items[0], command.line), number(p.items[1], command.line)));
            }
            n.h = a.positiveNumber("h", 10);
        } else if (c == "halfspace" || c == "plane") {
            // Inside: n . p <= d (below the plane, the normal pointing out).
            a.positional = {"normal", "d"};
            n.kind = Node::Kind::HalfSpace;
            const bool miller = a.has("hkl");
            const dvec3 direction = planeNormal(miller ? a.vector("hkl", dvec3(0, 0, 1)) : a.vector("normal", dvec3(0, 0, 1)), miller);
            if (glm::length(direction) < 1e-12) throw ShapeError(command.line, c + "(): the normal must not be zero");
            n.normals = {glm::normalize(direction)};
            n.offsets = {a.number("d", 0)};
        } else if (c == "slab") {
            a.positional = {"normal", "thickness", "d"};
            n.kind = Node::Kind::Slab;
            const bool miller = a.has("hkl");
            const dvec3 direction = planeNormal(miller ? a.vector("hkl", dvec3(0, 0, 1)) : a.vector("normal", dvec3(0, 0, 1)), miller);
            if (glm::length(direction) < 1e-12) throw ShapeError(command.line, "slab(): the normal must not be zero");
            n.normals = {glm::normalize(direction)};
            n.h = a.positiveNumber("thickness", 10);
            n.offsets = {a.number("d", 0)};
        } else if (c == "facets" || c == "polyhedron") {
            // Crystallographic facets: planes {hkl} at distances d (Angstrom) from the centre; with
            // family, every permutation and sign of (h, k, l) (cubic families).
            a.positional = {"hkl", "d"};
            n.kind = Node::Kind::Facets;
            const Value* hkl = a.find("hkl");
            if (!hkl || hkl->kind != Value::Kind::Vector || hkl->items.empty()) throw ShapeError(command.line, c + "(): hkl must be [h, k, l] or a list of them");
            std::vector<dvec3> planes;
            if (hkl->items[0].kind == Value::Kind::Vector)
                for (const Value& p : hkl->items) {
                    if (p.kind != Value::Kind::Vector || p.items.size() != 3) throw ShapeError(command.line, c + "(): each plane is [h, k, l]");
                    planes.push_back(dvec3(number(p.items[0], n.line), number(p.items[1], n.line), number(p.items[2], n.line)));
                }
            else {
                if (hkl->items.size() != 3) throw ShapeError(command.line, c + "(): hkl is [h, k, l]");
                planes.push_back(dvec3(number(hkl->items[0], n.line), number(hkl->items[1], n.line), number(hkl->items[2], n.line)));
            }
            std::vector<double> distances;
            const Value* d = a.find("d");
            if (!d) throw ShapeError(command.line, c + "(): give the facet distances d (Angstrom)");
            if (d->kind == Value::Kind::Vector)
                for (const Value& x : d->items) distances.push_back(number(x, n.line));
            else
                distances.assign(planes.size(), number(*d, n.line));
            if (distances.size() != planes.size()) throw ShapeError(command.line, c + "(): one distance per plane");
            const bool family = a.number("family", 1) != 0;
            for (std::size_t i = 0; i < planes.size(); ++i) {
                if (!(distances[i] > 0)) throw ShapeError(command.line, c + "(): distances must be positive");
                std::vector<dvec3> members{planes[i]};
                if (family) {
                    members.clear();
                    int order[3] = {0, 1, 2};
                    do {
                        for (int signs = 0; signs < 8; ++signs) {
                            dvec3 m(planes[i][order[0]], planes[i][order[1]], planes[i][order[2]]);
                            for (int s = 0; s < 3; ++s)
                                if (signs & (1 << s)) m[s] = -m[s];
                            bool seen = false;
                            for (const dvec3& q : members) seen = seen || glm::length(q - m) < 1e-9;
                            if (!seen) members.push_back(m);
                        }
                    } while (std::next_permutation(order, order + 3));
                }
                for (const dvec3& m : members) {
                    const dvec3 direction = planeNormal(m, true);
                    if (glm::length(direction) < 1e-12) throw ShapeError(command.line, c + "(): (0 0 0) is not a plane");
                    n.normals.push_back(glm::normalize(direction));
                    n.offsets.push_back(distances[i]);
                }
            }
        } else if (c == "octahedron") {
            a.positional = {"d"};
            n.kind = Node::Kind::Facets;
            const double d = a.positiveNumber("d", 10);
            for (int s = 0; s < 8; ++s) {
                n.normals.push_back(glm::normalize(dvec3(s & 1 ? -1 : 1, s & 2 ? -1 : 1, s & 4 ? -1 : 1)));
                n.offsets.push_back(d);
            }
        } else if (c == "union" || c == "group") {
            n.kind = Node::Kind::Union;
        } else if (c == "intersection" || c == "intersect") {
            n.kind = Node::Kind::Intersection;
        } else if (c == "difference" || c == "subtract") {
            n.kind = Node::Kind::Difference;
        } else if (c == "smooth_union" || c == "blend") {
            a.positional = {"k"};
            n.kind = Node::Kind::SmoothUnion;
            n.k = a.positiveNumber("k", 2);
        } else if (c == "smooth_intersection") {
            a.positional = {"k"};
            n.kind = Node::Kind::SmoothIntersection;
            n.k = a.positiveNumber("k", 2);
        } else if (c == "smooth_difference") {
            a.positional = {"k"};
            n.kind = Node::Kind::SmoothDifference;
            n.k = a.positiveNumber("k", 2);
        } else if (c == "translate" || c == "move") {
            a.positional = {"v"};
            n.kind = Node::Kind::Transform;
            n.translation = a.vector("v", dvec3(0.0));
        } else if (c == "rotate") {
            // rotate([ax, ay, az]): about x, then y, then z (degrees, as OpenSCAD); rotate(a, v): about an axis.
            a.positional = {"a", "v"};
            n.kind = Node::Kind::Transform;
            const Value* angle = a.find("a");
            if (!angle) throw ShapeError(command.line, "rotate(): give the angles [ax, ay, az] or an angle a and an axis v");
            if (angle->kind == Value::Kind::Vector) {
                const dvec3 e = a.vector("a", dvec3(0.0)) * kDegree;
                n.forward = rotation(dvec3(0, 0, 1), e.z) * rotation(dvec3(0, 1, 0), e.y) * rotation(dvec3(1, 0, 0), e.x);
            } else {
                const double t = number(*angle, command.line) * kDegree;
                dvec3 axis = a.vector("v", dvec3(0, 0, 1));
                if (glm::length(axis) < 1e-12) throw ShapeError(command.line, "rotate(): the axis must not be zero");
                n.forward = rotation(glm::normalize(axis), t);
            }
            n.inverse = glm::transpose(n.forward);
        } else if (c == "scale") {
            a.positional = {"v"};
            n.kind = Node::Kind::Transform;
            const dvec3 s = a.vector("v", dvec3(1.0), true);
            if (s.x == 0 || s.y == 0 || s.z == 0) throw ShapeError(command.line, "scale(): factors must not be zero");
            n.forward = dmat3(s.x, 0, 0, 0, s.y, 0, 0, 0, s.z);
            n.inverse = glm::inverse(n.forward);
            n.distanceScale = std::min(std::abs(s.x), std::min(std::abs(s.y), std::abs(s.z)));
        } else if (c == "mirror") {
            a.positional = {"v"};
            n.kind = Node::Kind::Transform;
            dvec3 normal = a.vector("v", dvec3(1, 0, 0));
            if (glm::length(normal) < 1e-12) throw ShapeError(command.line, "mirror(): the normal must not be zero");
            normal = glm::normalize(normal);
            n.forward = dmat3(1.0) - 2.0 * glm::outerProduct(normal, normal);
            n.inverse = n.forward;
        } else if (c == "shell" || c == "hollow") {
            a.positional = {"t"};
            n.kind = Node::Kind::Shell;
            n.h = a.positiveNumber("t", 2);
        } else if (c == "offset" || c == "grow") {
            a.positional = {"r"};
            n.kind = Node::Kind::Offset;
            n.r = a.number("r", 1);
        } else if (c == "linear_array") {
            a.positional = {"n", "step"};
            n.kind = Node::Kind::LinearArray;
            n.counts.x = count(a, "n", 2);
            n.step = a.vector("step", dvec3(10, 0, 0));
        } else if (c == "polar_array") {
            a.positional = {"n", "angle", "axis"};
            n.kind = Node::Kind::PolarArray;
            n.counts.x = count(a, "n", 6);
            n.angle = a.number("angle", 360);
            n.axis = a.vector("axis", dvec3(0, 0, 1));
            if (glm::length(n.axis) < 1e-12) throw ShapeError(command.line, "polar_array(): the axis must not be zero");
            n.axis = glm::normalize(n.axis);
        } else if (c == "grid_array") {
            a.positional = {"n", "step"};
            n.kind = Node::Kind::GridArray;
            const dvec3 counts = a.vector("n", dvec3(2, 2, 1));
            for (int i = 0; i < 3; ++i) {
                if (counts[i] < 1 || counts[i] != std::floor(counts[i]) || counts[i] > 10000)
                    throw ShapeError(command.line, "grid_array(): n must be whole numbers from 1");
                n.counts[i] = static_cast<int>(counts[i]);
            }
            n.step = a.vector("step", dvec3(10, 10, 10));
        } else if (c == "region" || c == "material") {
            a.positional = {"element"};
            n.kind = Node::Kind::Region;
            n.element = a.text("element", "");
            if (n.element.empty()) throw ShapeError(command.line, c + "(): give the element, e.g. region(\"Au\") { ... }");
        } else
            throw ShapeError(command.line, "unknown command '" + c + "'");
        a.finish();
        if (!isOperation(n.kind)) ++m_primitives;
        return n;
    }
    // center=false: the primitive's corner (or base) at the origin.
    Node corner(Node n, const dvec3& shift, Arguments& a)
    {
        a.finish();
        ++m_primitives;
        Node moved;
        moved.kind = Node::Kind::Transform;
        moved.line = n.line;
        moved.translation = shift;
        moved.children.push_back(std::move(n));
        return moved;
    }
    static int count(Arguments& a, const char* name, double fallback)
    {
        const double v = a.number(name, fallback);
        if (v < 1 || v != std::floor(v) || v > 100000) throw ShapeError(a.command.line, a.command.text + "(): " + name + " must be a whole number from 1");
        return static_cast<int>(v);
    }
    void finish(Node& n, const Token& command)
    {
        if (n.kind == Node::Kind::Transform && n.distanceScale <= 0) throw ShapeError(command.line, "invalid transform");
    }
};
}  // namespace

Model Model::parse(const std::string& script, const glm::dmat3* cell)
{
    Parser parser(script, cell);
    Model model;
    int primitives = 0;
    Node root = parser.program(primitives);
    if (root.children.empty()) throw ShapeError(1, "the script has no shape (try: sphere(r = 15);)");
    if (!hasVolume(root)) throw ShapeError(1, "the script has only regions: give the shape itself too (regions only assign elements inside it)");
    model.m_root = std::make_shared<Node>(std::move(root));
    model.m_primitives = primitives;
    return model;
}

double Model::distance(const glm::dvec3& p) const { return m_root ? distanceOf(*m_root, p) : kHuge; }

const std::string* Model::material(const glm::dvec3& p) const
{
    const std::string* out = nullptr;
    if (m_root) materialOf(*m_root, p, out);
    return out;
}

Bounds Model::bounds() const
{
    Bounds b = m_root ? boundsOf(*m_root) : Bounds{};
    if (b.empty) b = Bounds{};
    return b;
}

std::vector<std::string> Model::regionElements() const
{
    std::vector<std::string> out;
    if (m_root) collectRegions(*m_root, out);
    return out;
}

double Model::volume(double spacing) const
{
    const Bounds b = bounds();
    if (!b.finite || !(spacing > 0)) return 0;
    long long inside = 0;
    for (double x = b.low.x + spacing / 2; x < b.high.x; x += spacing)
        for (double y = b.low.y + spacing / 2; y < b.high.y; y += spacing)
            for (double z = b.low.z + spacing / 2; z < b.high.z; z += spacing)
                if (contains(dvec3(x, y, z))) ++inside;
    return static_cast<double>(inside) * spacing * spacing * spacing;
}

electronic::Mesh Model::mesh(double spacing) const
{
    const Bounds b = bounds();
    if (!b.finite) throw std::runtime_error("The shape is unbounded: intersect it with a box to mesh it");
    if (!(spacing > 0)) throw std::runtime_error("The mesh spacing must be positive");
    electronic::Grid grid;
    const dvec3 low = b.low - dvec3(2 * spacing), extent = (b.high - b.low) + dvec3(4 * spacing);
    for (int a = 0; a < 3; ++a) grid.shape[static_cast<std::size_t>(a)] = std::max(2, static_cast<int>(std::ceil(extent[a] / spacing)) + 1);
    const double samples = double(grid.shape[0]) * grid.shape[1] * grid.shape[2];
    if (samples > 6e7) throw std::runtime_error("The mesh grid is too fine for this shape; use a larger spacing");
    grid.cell = dmat3(extent.x, 0, 0, 0, extent.y, 0, 0, 0, extent.z);
    grid.origin = low;
    grid.periodic = false;
    grid.values.resize(static_cast<std::size_t>(samples));
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x) grid.values[grid.index(x, y, z)] = distance(grid.position(x, y, z));
    return electronic::isosurface(grid, 0.0);
}

void saveObj(const electronic::Mesh& mesh, const std::string& path)
{
    std::ofstream out{std::filesystem::u8path(path)};
    if (!out) throw std::runtime_error("Cannot write " + path);
    out << std::setprecision(10) << "# AtomForge shape: " << mesh.vertices.size() / 3 << " triangles (Angstrom)\n";
    // Shared vertices: triangle corners at the same place become one vertex.
    std::map<std::array<long long, 3>, std::size_t> index;
    std::vector<std::size_t> faces;
    std::vector<dvec3> vertices;
    for (const dvec3& v : mesh.vertices) {
        const std::array<long long, 3> key{std::llround(v.x * 1e6), std::llround(v.y * 1e6), std::llround(v.z * 1e6)};
        auto it = index.find(key);
        if (it == index.end()) {
            it = index.emplace(key, vertices.size()).first;
            vertices.push_back(v);
        }
        faces.push_back(it->second);
    }
    for (const dvec3& v : vertices) out << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';
    for (std::size_t i = 0; i + 2 < faces.size(); i += 3) out << "f " << faces[i] + 1 << ' ' << faces[i + 1] + 1 << ' ' << faces[i + 2] + 1 << '\n';
    if (!out) throw std::runtime_error("Failed while writing " + path);
}

const std::vector<std::pair<std::string, std::string>>& templates()
{
    static const std::vector<std::pair<std::string, std::string>> list = {
        {"Sphere", "sphere(r = 15);\n"},
        {"Core-shell particle", "// A particle of the crystal with a core of another element.\nR = 20;\nsphere(r = R);\nregion(\"Au\") sphere(r = 0.6 * R);\n"},
        {"Hollow sphere (nanoshell)", "shell(t = 5) sphere(r = 20);\n"},
        {"Hollow cylinder (tube)", "difference() {\n    cylinder(r = 15, h = 60);\n    cylinder(r = 10, h = 62);\n}\n"},
        {"Truncated octahedron {111} + {100}", "// Distances of the facets from the centre (Angstrom); family: all equivalent planes (cubic).\nfacets(hkl = [[1, 1, 1], [1, 0, 0]], d = [18, 20], family = true);\n"},
        {"Faceted rod {110} x {001}", "intersection() {\n    facets(hkl = [[1, 1, 0], [1, -1, 0]], d = [15, 15], family = false) ;\n    facets(hkl = [[-1, -1, 0], [-1, 1, 0]], d = [15, 15], family = false);\n    slab(normal = [0, 0, 1], thickness = 60);\n}\n"},
        {"Hexagonal nanowire segment", "prism(n = 6, r = 15, h = 60);\n"},
        {"Cube with a spherical void", "difference() {\n    box(size = [40, 40, 40]);\n    sphere(r = 10);\n}\n"},
        {"Nanoporous block (pore lattice)", "difference() {\n    box(size = [60, 60, 30]);\n    translate([-22.5, -22.5, 0])\n        grid_array(n = [4, 4, 1], step = [15, 15, 0]) cylinder(r = 4, h = 32);\n}\n"},
        {"Dumbbell (sintered particles)", "smooth_union(k = 4) {\n    translate([-12, 0, 0]) sphere(r = 10);\n    translate([12, 0, 0]) sphere(r = 10);\n}\n"},
        {"Nanoring (torus)", "torus(R = 25, r = 6);\n"},
        {"Janus particle", "R = 20;\nsphere(r = R);\n// The upper half (z >= 0) becomes silver.\nregion(\"Ag\") halfspace(normal = [0, 0, -1], d = 0);\n"},
        {"Probe tip (cone on a shank)", "union() {\n    translate([0, 0, 15]) cone(r1 = 15, r2 = 2, h = 30);\n    translate([0, 0, -10]) cylinder(r = 15, h = 20);\n}\n"},
        {"Matrix with an ellipsoidal precipitate", "box(size = [60, 60, 60]);\nregion(\"Ni\") rotate([0, 0, 30]) ellipsoid(r = [15, 8, 8]);\n"},
        {"Notched beam (fracture specimen)", "difference() {\n    box(size = [80, 20, 30]);\n    // A V-notch from the top surface, along y.\n    translate([0, 0, 15]) rotate([90, 0, 0]) prism(n = 4, r = 8, h = 22);\n}\n"},
        {"Gear-like disc (polar array)", "difference() {\n    cylinder(r = 25, h = 10);\n    polar_array(n = 8) translate([25, 0, 0]) cylinder(r = 5, h = 12);\n}\n"},
        {"Embedded nanoparticle array", "box(size = [80, 80, 30]);\nregion(\"Au\") translate([-25, -25, 0]) grid_array(n = [3, 3, 1], step = [25, 25, 0]) sphere(r = 6);\n"},
    };
    return list;
}

const std::vector<Command>& commands()
{
    static const std::vector<Command> list = {
        {"Primitives", "sphere", "sphere(r = 10);", "Sphere of radius r (or diameter d) about the origin"},
        {"Primitives", "box", "box(size = [20, 20, 20]);", "Box of the given edge lengths, centred (center = false: corner at the origin)"},
        {"Primitives", "cylinder", "cylinder(r = 8, h = 30);", "Cylinder along z, centred; r1 and r2 make a frustum"},
        {"Primitives", "cone", "cone(r1 = 10, r2 = 0, h = 20);", "Cone or frustum along z: radius r1 at the bottom, r2 at the top"},
        {"Primitives", "torus", "torus(R = 20, r = 5);", "Torus about z: ring radius R, tube radius r"},
        {"Primitives", "capsule", "capsule(r = 6, h = 20);", "Cylinder of length h with hemispherical ends"},
        {"Primitives", "ellipsoid", "ellipsoid(r = [15, 10, 6]);", "Ellipsoid with semi-axes along x, y, z"},
        {"Primitives", "prism", "prism(n = 6, r = 12, h = 30);", "Regular n-sided prism along z (circumradius r, or apothem)"},
        {"Primitives", "polygon", "polygon(points = [[0, 0], [20, 0], [10, 15]], h = 10);", "Any polygon in the xy plane extruded along z"},
        {"Primitives", "facets", "facets(hkl = [[1, 1, 1], [1, 0, 0]], d = [15, 17]);", "Crystal facets {hkl} at distances d (Angstrom); family: all equivalent planes"},
        {"Primitives", "octahedron", "octahedron(d = 12);", "Regular octahedron, faces at distance d from the centre"},
        {"Primitives", "halfspace", "halfspace(normal = [0, 0, 1], d = 0);", "Everything below a plane (n . p <= d); also hkl = [h, k, l]"},
        {"Primitives", "slab", "slab(normal = [0, 0, 1], thickness = 10);", "A slab of a thickness about a plane (unbounded: intersect it)"},
        {"Booleans", "union", "union() {\n    \n}", "All of the shapes inside"},
        {"Booleans", "intersection", "intersection() {\n    \n}", "Where all of the shapes overlap"},
        {"Booleans", "difference", "difference() {\n    \n}", "The first shape minus the others"},
        {"Booleans", "smooth_union", "smooth_union(k = 3) {\n    \n}", "Union blended over about k Angstrom (necks, fillets)"},
        {"Booleans", "smooth_difference", "smooth_difference(k = 3) {\n    \n}", "Difference with rounded edges"},
        {"Transforms", "translate", "translate([0, 0, 10])", "Move the shapes that follow"},
        {"Transforms", "rotate", "rotate([0, 0, 45])", "Turn about x, then y, then z (degrees); rotate(a = 30, v = [1, 1, 1]) about an axis"},
        {"Transforms", "scale", "scale([1, 1, 2])", "Stretch along x, y, z"},
        {"Transforms", "mirror", "mirror([1, 0, 0])", "Reflect through the plane with this normal"},
        {"Modifiers", "shell", "shell(t = 3)", "Hollow: keep a skin of thickness t inside the surface"},
        {"Modifiers", "offset", "offset(r = 2)", "Grow (r > 0) or shrink (r < 0) the shape"},
        {"Patterns", "linear_array", "linear_array(n = 3, step = [15, 0, 0])", "n copies, each moved by step"},
        {"Patterns", "polar_array", "polar_array(n = 6, angle = 360, axis = [0, 0, 1])", "n copies turned about an axis through the origin"},
        {"Patterns", "grid_array", "grid_array(n = [3, 3, 1], step = [15, 15, 0])", "A 3D grid of copies"},
        {"Regions", "region", "region(\"Au\")", "Atoms of the shape inside get this element (\"none\": removed); adds no volume; later regions win"},
    };
    return list;
}
}
