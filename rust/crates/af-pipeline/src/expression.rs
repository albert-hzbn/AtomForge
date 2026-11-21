//! Boolean/arithmetic expressions over atom variables (the C++
//! pipeline/Expression), compiled once:
//!   x > 10 && element == "O",  fz < 0.5 or Z == 29,  abs(property) >= 0.2
//! Operators: + - * / % ^, comparisons < <= > >= == !=, logic && || ! (also
//! and, or, not), parentheses and the functions abs sqrt exp log sin cos tan
//! floor ceil min max. Element symbols may be written with or without quotes.

use af_core::elements;
use af_core::format::stod;
use af_core::suggest;

/// Values seen by an expression for one atom.
#[derive(Clone, Debug, Default)]
pub struct AtomVariables {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
    pub index: f64,
    pub magmom: f64,
    pub charge: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub count: f64,
    pub property: f64,
    pub frozen: f64,
    pub selected: f64,
    pub atomic_number: f64,
    pub fixed: f64,
    pub element: String,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Variable {
    X,
    Y,
    Z,
    Fx,
    Fy,
    Fz,
    Index,
    Count,
    Property,
    Frozen,
    Selected,
    AtomicNumber,
    Fixed,
    Magmom,
    Charge,
    Vx,
    Vy,
    Vz,
}

const VARIABLES: [(&str, Variable); 19] = [
    ("x", Variable::X),
    ("y", Variable::Y),
    ("z", Variable::Z),
    ("fx", Variable::Fx),
    ("fy", Variable::Fy),
    ("fz", Variable::Fz),
    ("index", Variable::Index),
    ("count", Variable::Count),
    ("property", Variable::Property),
    ("frozen", Variable::Frozen),
    ("selected", Variable::Selected),
    ("Z", Variable::AtomicNumber),
    ("atomic_number", Variable::AtomicNumber),
    ("fixed", Variable::Fixed),
    ("magmom", Variable::Magmom),
    ("charge", Variable::Charge),
    ("vx", Variable::Vx),
    ("vy", Variable::Vy),
    ("vz", Variable::Vz),
];

impl Variable {
    fn of(self, v: &AtomVariables) -> f64 {
        match self {
            Variable::X => v.x,
            Variable::Y => v.y,
            Variable::Z => v.z,
            Variable::Fx => v.fx,
            Variable::Fy => v.fy,
            Variable::Fz => v.fz,
            Variable::Index => v.index,
            Variable::Count => v.count,
            Variable::Property => v.property,
            Variable::Frozen => v.frozen,
            Variable::Selected => v.selected,
            Variable::AtomicNumber => v.atomic_number,
            Variable::Fixed => v.fixed,
            Variable::Magmom => v.magmom,
            Variable::Charge => v.charge,
            Variable::Vx => v.vx,
            Variable::Vy => v.vy,
            Variable::Vz => v.vz,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Unary {
    Abs,
    Sqrt,
    Exp,
    Log,
    Sin,
    Cos,
    Tan,
    Floor,
    Ceil,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Compare {
    Le,
    Ge,
    Eq,
    Ne,
    Lt,
    Gt,
}

#[derive(Clone, Debug)]
enum Node {
    Number(f64),
    Text(String),
    Variable(Variable),
    Element,
    Or(Box<Node>, Box<Node>),
    And(Box<Node>, Box<Node>),
    Not(Box<Node>),
    Compare(Compare, Box<Node>, Box<Node>),
    Arithmetic(u8, Box<Node>, Box<Node>),
    Negate(u8, Box<Node>),
    Power(Box<Node>, Box<Node>),
    Function(Unary, Box<Node>),
    MinMax(bool, Box<Node>, Box<Node>),
}

/// A number or, for element comparisons, a string.
enum Value {
    Number(f64),
    Text(String),
}

fn number_of(v: Value, context: &str) -> Result<f64, String> {
    match v {
        Value::Number(n) => Ok(n),
        Value::Text(t) => Err(format!("'{t}' is not a number ({context})")),
    }
}

impl Node {
    fn evaluate(&self, v: &AtomVariables) -> Result<Value, String> {
        let truth = |b: bool| Value::Number(if b { 1.0 } else { 0.0 });
        Ok(match self {
            Node::Number(n) => Value::Number(*n),
            Node::Text(t) => Value::Text(t.clone()),
            Node::Variable(var) => Value::Number(var.of(v)),
            Node::Element => Value::Text(v.element.clone()),
            // Short-circuit, as the C++ || and &&.
            Node::Or(a, b) => truth(number_of(a.evaluate(v)?, "or")? != 0.0 || number_of(b.evaluate(v)?, "or")? != 0.0),
            Node::And(a, b) => truth(number_of(a.evaluate(v)?, "and")? != 0.0 && number_of(b.evaluate(v)?, "and")? != 0.0),
            Node::Not(a) => truth(number_of(a.evaluate(v)?, "not")? == 0.0),
            Node::Compare(op, a, b) => {
                let (x, y) = (a.evaluate(v)?, b.evaluate(v)?);
                match (&x, &y) {
                    (Value::Number(x), Value::Number(y)) => truth(match op {
                        Compare::Lt => x < y,
                        Compare::Le => x <= y,
                        Compare::Gt => x > y,
                        Compare::Ge => x >= y,
                        Compare::Eq => x == y,
                        Compare::Ne => x != y,
                    }),
                    _ => {
                        if *op != Compare::Eq && *op != Compare::Ne {
                            return Err("Element names can only be compared with == or !=".into());
                        }
                        let equal = match (&x, &y) {
                            (Value::Text(a), Value::Text(b)) => a == b,
                            _ => false,
                        };
                        truth((*op == Compare::Eq) == equal)
                    }
                }
            }
            Node::Arithmetic(op, a, b) => {
                let context = if *op == b'+' || *op == b'-' { "+/-" } else { "* / %" };
                let x = number_of(a.evaluate(v)?, context)?;
                let y = number_of(b.evaluate(v)?, context)?;
                Value::Number(match op {
                    b'+' => x + y,
                    b'-' => x - y,
                    b'*' => x * y,
                    b'/' => x / y,
                    // The floored remainder (index % 2 is 0 or 1, also for negative values).
                    _ => x - y * (x / y).floor(),
                })
            }
            Node::Negate(op, a) => {
                let x = number_of(a.evaluate(v)?, "sign")?;
                Value::Number(if *op == b'-' { -x } else { x })
            }
            Node::Power(a, b) => {
                let x = number_of(a.evaluate(v)?, "^")?;
                Value::Number(x.powf(number_of(b.evaluate(v)?, "^")?))
            }
            Node::Function(f, a) => {
                let x = number_of(a.evaluate(v)?, "function")?;
                Value::Number(match f {
                    Unary::Abs => x.abs(),
                    Unary::Sqrt => x.sqrt(),
                    Unary::Exp => x.exp(),
                    Unary::Log => x.ln(),
                    Unary::Sin => x.sin(),
                    Unary::Cos => x.cos(),
                    Unary::Tan => x.tan(),
                    Unary::Floor => x.floor(),
                    Unary::Ceil => x.ceil(),
                })
            }
            Node::MinMax(is_min, a, b) => {
                let x = number_of(a.evaluate(v)?, "min/max")?;
                let y = number_of(b.evaluate(v)?, "min/max")?;
                // std::min(x, y) = y < x ? y : x; std::max(x, y) = x < y ? y : x.
                Value::Number(if *is_min {
                    if y < x {
                        y
                    } else {
                        x
                    }
                } else if x < y {
                    y
                } else {
                    x
                })
            }
        })
    }
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

struct Parser<'a> {
    text: &'a [u8],
    pos: usize,
}

type Parse = Result<Box<Node>, String>;

impl Parser<'_> {
    fn skip(&mut self) {
        while self.pos < self.text.len() && is_space(self.text[self.pos]) {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.pos).copied()
    }

    fn accept(&mut self, token: &str) -> bool {
        self.skip();
        let t = token.as_bytes();
        if !self.text[self.pos..].starts_with(t) {
            return false;
        }
        // Words must not run into identifier characters (e.g. "and" vs "andx").
        if t[0].is_ascii_alphabetic() {
            if let Some(&next) = self.text.get(self.pos + t.len()) {
                if next.is_ascii_alphanumeric() || next == b'_' {
                    return false;
                }
            }
        }
        self.pos += t.len();
        true
    }

    fn fail<T>(&self, message: &str) -> Result<T, String> {
        Err(format!("Expression error at position {}: {message}", self.pos + 1))
    }

    fn parse(&mut self) -> Parse {
        let result = self.parse_or()?;
        self.skip();
        if self.pos != self.text.len() {
            let end = (self.pos + 12).min(self.text.len());
            let rest = String::from_utf8_lossy(&self.text[self.pos..end]).into_owned();
            return self.fail(&format!("unexpected '{rest}'"));
        }
        Ok(result)
    }

    fn parse_or(&mut self) -> Parse {
        let mut left = self.parse_and()?;
        while self.accept("||") || self.accept("or") {
            left = Box::new(Node::Or(left, self.parse_and()?));
        }
        Ok(left)
    }

    fn parse_and(&mut self) -> Parse {
        let mut left = self.parse_not()?;
        while self.accept("&&") || self.accept("and") {
            left = Box::new(Node::And(left, self.parse_not()?));
        }
        Ok(left)
    }

    fn parse_not(&mut self) -> Parse {
        self.skip();
        let bang = self.peek() == Some(b'!') && self.text.get(self.pos + 1) != Some(&b'=');
        if bang || self.accept("not") {
            if self.peek() == Some(b'!') {
                self.pos += 1;
            }
            return Ok(Box::new(Node::Not(self.parse_not()?)));
        }
        self.parse_comparison()
    }

    fn parse_comparison(&mut self) -> Parse {
        let left = self.parse_sum()?;
        for (token, op) in [("<=", Compare::Le), (">=", Compare::Ge), ("==", Compare::Eq), ("!=", Compare::Ne), ("<", Compare::Lt), (">", Compare::Gt)] {
            if self.accept(token) {
                return Ok(Box::new(Node::Compare(op, left, self.parse_sum()?)));
            }
        }
        Ok(left)
    }

    fn parse_sum(&mut self) -> Parse {
        let mut left = self.parse_product()?;
        loop {
            self.skip();
            match self.peek() {
                Some(op @ (b'+' | b'-')) => {
                    self.pos += 1;
                    left = Box::new(Node::Arithmetic(op, left, self.parse_product()?));
                }
                _ => return Ok(left),
            }
        }
    }

    fn parse_product(&mut self) -> Parse {
        let mut left = self.parse_unary()?;
        loop {
            self.skip();
            match self.peek() {
                Some(op @ (b'*' | b'/' | b'%')) => {
                    self.pos += 1;
                    left = Box::new(Node::Arithmetic(op, left, self.parse_unary()?));
                }
                _ => return Ok(left),
            }
        }
    }

    fn parse_unary(&mut self) -> Parse {
        self.skip();
        if let Some(op @ (b'-' | b'+')) = self.peek() {
            self.pos += 1;
            return Ok(Box::new(Node::Negate(op, self.parse_unary()?)));
        }
        self.parse_power()
    }

    fn parse_power(&mut self) -> Parse {
        let base = self.parse_primary()?;
        self.skip();
        if self.peek() == Some(b'^') {
            self.pos += 1;
            return Ok(Box::new(Node::Power(base, self.parse_unary()?)));
        }
        Ok(base)
    }

    fn parse_primary(&mut self) -> Parse {
        self.skip();
        let Some(c) = self.peek() else { return self.fail("expression ends early") };
        if c == b'(' {
            self.pos += 1;
            let inner = self.parse_or()?;
            if !self.accept(")") {
                return self.fail("missing ')'");
            }
            return Ok(inner);
        }
        if c == b'"' || c == b'\'' {
            let Some(end) = self.text[self.pos + 1..].iter().position(|&b| b == c).map(|e| e + self.pos + 1) else {
                return self.fail("unterminated string");
            };
            let text = String::from_utf8_lossy(&self.text[self.pos + 1..end]).into_owned();
            self.pos = end + 1;
            return Ok(Box::new(Node::Text(text)));
        }
        if c.is_ascii_digit() || c == b'.' {
            let rest = String::from_utf8_lossy(&self.text[self.pos..]).into_owned();
            return match stod(&rest) {
                Ok((value, used)) => {
                    self.pos += used;
                    Ok(Box::new(Node::Number(value)))
                }
                Err(_) => self.fail("invalid number"),
            };
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let mut end = self.pos;
            while end < self.text.len() && (self.text[end].is_ascii_alphanumeric() || self.text[end] == b'_') {
                end += 1;
            }
            let name = String::from_utf8_lossy(&self.text[self.pos..end]).into_owned();
            self.pos = end;
            self.skip();
            if self.peek() == Some(b'(') {
                return self.parse_function(&name);
            }
            return self.variable(&name);
        }
        let shown = String::from_utf8_lossy(&[c]).into_owned();
        self.fail(&format!("unexpected '{shown}'"))
    }

    fn parse_function(&mut self, name: &str) -> Parse {
        self.pos += 1; // '('
        let mut arguments = Vec::new();
        if !self.accept(")") {
            loop {
                arguments.push(self.parse_or()?);
                if !self.accept(",") {
                    break;
                }
            }
            if !self.accept(")") {
                return self.fail(&format!("missing ')' after the arguments of {name}"));
            }
        }
        let unary = [
            ("abs", Unary::Abs),
            ("sqrt", Unary::Sqrt),
            ("exp", Unary::Exp),
            ("log", Unary::Log),
            ("sin", Unary::Sin),
            ("cos", Unary::Cos),
            ("tan", Unary::Tan),
            ("floor", Unary::Floor),
            ("ceil", Unary::Ceil),
        ];
        if let Some(&(_, f)) = unary.iter().find(|(n, _)| *n == name) {
            if arguments.len() != 1 {
                return self.fail(&format!("{name} takes one argument"));
            }
            return Ok(Box::new(Node::Function(f, arguments.pop().unwrap())));
        }
        if name == "min" || name == "max" {
            if arguments.len() != 2 {
                return self.fail(&format!("{name} takes two arguments"));
            }
            let b = arguments.pop().unwrap();
            let a = arguments.pop().unwrap();
            return Ok(Box::new(Node::MinMax(name == "min", a, b)));
        }
        self.fail(&format!("unknown function {name}"))
    }

    fn variable(&mut self, name: &str) -> Parse {
        if let Some(&(_, var)) = VARIABLES.iter().find(|(n, _)| *n == name) {
            return Ok(Box::new(Node::Variable(var)));
        }
        if name == "element" || name == "type" {
            return Ok(Box::new(Node::Element));
        }
        if name == "pi" {
            return Ok(Box::new(Node::Number((-1.0f64).acos())));
        }
        // An element symbol stands for itself, in any case (cu and CU mean Cu).
        let lower = name.to_ascii_lowercase();
        for z in 1..=118 {
            if elements::symbol(z).to_ascii_lowercase() == lower {
                return Ok(Box::new(Node::Text(elements::symbol(z).to_string())));
            }
        }
        let mut names = vec!["element", "Z"];
        names.extend(VARIABLES.iter().map(|(n, _)| *n));
        self.fail(&format!(
            "unknown name '{name}'{}; the variables are x y z fx fy fz index count Z element property frozen selected fixed, and other text goes in quotes",
            suggest::did_you_mean(name, names)
        ))
    }
}

/// A compiled expression.
#[derive(Clone, Debug)]
pub struct Expression {
    root: Box<Node>,
}

impl Expression {
    pub fn new(text: &str) -> Result<Self, String> {
        if text.chars().all(|c| matches!(c, ' ' | '\t' | '\r' | '\n')) {
            return Err("The selection expression is empty".into());
        }
        let root = Parser { text: text.as_bytes(), pos: 0 }.parse()?;
        Ok(Self { root })
    }

    pub fn evaluate(&self, variables: &AtomVariables) -> Result<f64, String> {
        match self.root.evaluate(variables)? {
            Value::Number(n) => Ok(n),
            Value::Text(t) => Err(format!("The expression gives the text '{t}', not a condition")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval(text: &str, v: &AtomVariables) -> Result<f64, String> {
        Expression::new(text)?.evaluate(v)
    }

    #[test]
    fn grammar() {
        let v = AtomVariables { x: 2.0, z: -1.0, index: 5.0, element: "Cu".into(), ..Default::default() };
        assert_eq!(eval("x > 1 && element == Cu", &v), Ok(1.0));
        assert_eq!(eval("element == 'O' or not x < 1", &v), Ok(1.0));
        assert_eq!(eval("-index % 2", &v), Ok(1.0));
        assert_eq!(eval("2 ^ 3 ^ 2", &v), Ok(512.0));
        assert_eq!(eval("-2 ^ 2", &v), Ok(-4.0));
        assert_eq!(eval("min(x, z) + max(1, 3)", &v), Ok(2.0));
        assert_eq!(eval("0x10 + .5", &v), Ok(16.5));
        assert!(eval("elemnt == Cu", &v).unwrap_err().contains("did you mean element?"));
        assert!(eval("element < Cu", &v).unwrap_err().contains("only be compared"));
        assert!(eval("x >", &v).unwrap_err().contains("position 4: expression ends early"));
        assert!(eval("Cu", &v).unwrap_err().contains("gives the text 'Cu'"));
        assert_eq!(eval("1 || 'a'", &v), Ok(1.0));
        assert_eq!(Expression::new("  ").unwrap_err(), "The selection expression is empty");
    }
}
