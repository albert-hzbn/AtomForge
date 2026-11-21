//! Command-line option text of Build steps (`--count 4 --frac "Cu=0.7,Zn=0.3"`):
//! splitting it into tokens, reading flags from it and writing it back
//! (the C++ pipeline/Options.h). The Build and Edit dialogs use these to edit
//! a pipeline step's options.

/// Shell-style tokens (quotes group words and are removed).
pub fn tokens(text: &str) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut current = String::new();
    let mut in_token = false;
    let mut quote: Option<char> = None;
    for c in text.chars() {
        if let Some(q) = quote {
            if c == q {
                quote = None;
            } else {
                current.push(c);
            }
            continue;
        }
        if c == '"' || c == '\'' {
            quote = Some(c);
            in_token = true;
            continue;
        }
        if c.is_ascii_whitespace() || c == '\u{0b}' {
            if in_token {
                out.push(std::mem::take(&mut current));
                in_token = false;
            }
            continue;
        }
        current.push(c);
        in_token = true;
    }
    if quote.is_some() {
        return Err("Unterminated quote in the options".into());
    }
    if in_token {
        out.push(current);
    }
    Ok(out)
}

/// The longest leading number `std::stod` accepts (after leading spaces).
fn stod(text: &str) -> Option<f64> {
    let t = text.trim_start();
    let lower = t.to_ascii_lowercase();
    for special in ["+infinity", "-infinity", "infinity", "+inf", "-inf", "inf", "+nan", "-nan", "nan"] {
        if lower.starts_with(special) {
            let v = if special.contains("nan") { f64::NAN } else { f64::INFINITY };
            return Some(if special.starts_with('-') { -v } else { v });
        }
    }
    let b = t.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    let mut digits = i - digits_start;
    if i < b.len() && b[i] == b'.' {
        i += 1;
        let f = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        digits += i - f;
    }
    if digits == 0 {
        return None;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let e = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > e {
            i = j;
        }
    }
    t[..i].parse::<f64>().ok()
}

/// The leading integer `std::stoi` accepts; None when absent or out of range.
fn stoi(text: &str) -> Option<i32> {
    let t = text.trim_start();
    let b = t.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let s = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == s {
        return None;
    }
    t[..i].parse::<i32>().ok()
}

/// Flags of an option text. A flag's value is the token after it unless that
/// token is another flag (a negative number is a value, not a flag).
#[derive(Clone, Debug, Default)]
pub struct Reader {
    tokens: Vec<String>,
}

impl Reader {
    pub fn new(text: &str) -> Result<Self, String> {
        Ok(Self { tokens: tokens(text)? })
    }

    pub fn has(&self, flag: &str) -> bool {
        self.tokens.iter().any(|t| t == flag)
    }

    /// Every value given to `flag`, in order (flags such as --atom repeat).
    pub fn values(&self, flag: &str) -> Vec<String> {
        self.tokens.windows(2).filter(|w| w[0] == flag && !Self::is_flag(&w[1])).map(|w| w[1].clone()).collect()
    }

    pub fn text(&self, flag: &str, fallback: &str) -> String {
        self.values(flag).pop().unwrap_or_else(|| fallback.to_string())
    }

    pub fn number(&self, flag: &str, fallback: f64) -> f64 {
        let value = self.text(flag, "");
        if value.is_empty() {
            return fallback;
        }
        stod(&value).unwrap_or(fallback)
    }

    pub fn integer(&self, flag: &str, fallback: i32) -> i32 {
        let value = self.text(flag, "");
        if value.is_empty() {
            return fallback;
        }
        stoi(&value).unwrap_or(fallback)
    }

    /// Numbers in a value such as "0 0 1" or "1,1,0" (stops at the first non-number).
    pub fn numbers(&self, flag: &str) -> Vec<f64> {
        let value = self.text(flag, "").replace(',', " ");
        let mut out = Vec::new();
        for word in value.split_ascii_whitespace() {
            match word.parse::<f64>() {
                Ok(x) => out.push(x),
                Err(_) => {
                    // istream >> double reads a numeric prefix, then fails.
                    if let Some(x) = stod(word) {
                        out.push(x);
                    }
                    break;
                }
            }
        }
        out
    }

    pub fn is_flag(token: &str) -> bool {
        let b = token.as_bytes();
        b.len() > 1 && b[0] == b'-' && !b[1].is_ascii_digit() && b[1] != b'.'
    }
}

/// Builds an option text; values with spaces or special characters are quoted.
#[derive(Clone, Debug, Default)]
pub struct Writer {
    text: String,
}

impl Writer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn flag(&mut self, name: &str) -> &mut Self {
        self.append(name);
        self
    }

    pub fn add(&mut self, name: &str, value: &str) -> &mut Self {
        self.append(name);
        let q = Self::quoted(value);
        self.append(&q);
        self
    }

    pub fn add_number(&mut self, name: &str, value: f64) -> &mut Self {
        self.add(name, &Self::number(value))
    }

    pub fn add_integer(&mut self, name: &str, value: i64) -> &mut Self {
        self.add(name, &value.to_string())
    }

    /// Several numbers as one value ("0 0 1").
    pub fn add_numbers(&mut self, name: &str, values: &[f64]) -> &mut Self {
        let text = values.iter().map(|&v| Self::number(v)).collect::<Vec<_>>().join(" ");
        self.add(name, &text)
    }

    pub fn str(&self) -> &str {
        &self.text
    }

    /// `ostream << value` with precision 10 (%.10g).
    pub fn number(value: f64) -> String {
        af_core::format::format_g(value, 10)
    }

    pub fn quoted(value: &str) -> String {
        if !value.is_empty() && !value.contains([' ', '\t', '"', '\'', '|']) {
            return value.to_string();
        }
        let quote = if value.contains('"') { '\'' } else { '"' };
        format!("{quote}{value}{quote}")
    }

    fn append(&mut self, token: &str) {
        if !self.text.is_empty() {
            self.text.push(' ');
        }
        self.text.push_str(token);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_and_quotes() {
        assert_eq!(tokens(r#"--count 4 --frac "Cu=0.7, Zn=0.3" --x 'a"b'"#).unwrap(), vec!["--count", "4", "--frac", "Cu=0.7, Zn=0.3", "--x", "a\"b"]);
        assert_eq!(tokens(r#"a""b"#).unwrap(), vec!["ab"]);
        assert_eq!(tokens(r#""""#).unwrap(), vec![""]);
        assert!(tokens(r#"--a "open"#).is_err());
    }

    #[test]
    fn reader() {
        let r = Reader::new("--atom \"Na 0 0 0\" --atom \"Cl .5 .5 .5\" --shift -1.5 --n 12x --hkl 1,1,0 --flag --last").unwrap();
        assert_eq!(r.values("--atom"), vec!["Na 0 0 0", "Cl .5 .5 .5"]);
        assert_eq!(r.number("--shift", 0.0), -1.5);
        assert_eq!(r.integer("--n", 3), 12);
        assert_eq!(r.numbers("--hkl"), vec![1.0, 1.0, 0.0]);
        assert!(r.has("--flag"));
        assert_eq!(r.text("--flag", "none"), "none");
        assert_eq!(r.number("--missing", 2.5), 2.5);
        assert!(Reader::is_flag("--x") && !Reader::is_flag("-1") && !Reader::is_flag("-.5") && !Reader::is_flag("-"));
    }

    #[test]
    fn writer() {
        let mut w = Writer::new();
        w.add_number("--a", 3.61).add_integer("--n", 4).add_numbers("--hkl", &[1.0, -1.0, 0.0]).add("--frac", "Cu=0.7").add("--name", "it's").add("--e", "").flag("--conventional");
        assert_eq!(w.str(), "--a 3.61 --n 4 --hkl \"1 -1 0\" --frac Cu=0.7 --name \"it's\" --e \"\" --conventional");
        assert_eq!(Writer::number(1.0 / 3.0), "0.3333333333");
        assert_eq!(Writer::number(1e-12), "1e-12");
        assert_eq!(Writer::quoted("say \"hi\""), "'say \"hi\"'");
    }
}
