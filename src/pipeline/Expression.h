#pragma once

#include <memory>
#include <string>

namespace atomforge::pipeline
{
// Values seen by an expression for one atom.
struct AtomVariables
{
    double x = 0, y = 0, z = 0;     // Cartesian position (Angstrom)
    double fx = 0, fy = 0, fz = 0;  // fractional position (0 without a cell)
    double index = 0;               // zero-based atom index
    double count = 0;               // number of atoms
    double property = 0;            // per-atom property (NaN when absent)
    double selected = 0;            // 1 when selected
    double atomicNumber = 0;
    std::string element;
};

// A boolean/arithmetic expression over atom variables, compiled once:
//   x > 10 && element == "O",  fz < 0.5 or Z == 29,  abs(property) >= 0.2
// Operators: + - * / ^ (power), comparisons < <= > >= == !=, logic && || !
// (also and, or, not), parentheses, and the functions abs sqrt exp log sin cos
// tan floor ceil min max. Element symbols may be written with or without quotes.
// A nonzero result selects the atom.
class Expression
{
public:
    explicit Expression(const std::string& text);
    ~Expression();
    Expression(Expression&&) noexcept;
    Expression& operator=(Expression&&) noexcept;
    double evaluate(const AtomVariables& variables) const;
    struct Node;

private:
    std::unique_ptr<Node> m_root;
};
}
