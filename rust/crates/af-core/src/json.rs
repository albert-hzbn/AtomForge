//! The insertion-ordered JSON value of the C++ science/Json: the same parser
//! (NaN, Infinity, depth limit, error messages) and the same output (whole
//! numbers below 1e15 without a fraction, others with 17 digits, non-finite
//! as null, numeric arrays on one line), so pipelines, projects and the
//! native API exchange identical text.

use crate::format::{format_g, parse_prefix_f64};

#[derive(Clone, Debug, Default, PartialEq)]
pub enum Json {
    #[default]
    Null,
    Bool(bool),
    Number(f64),
    String(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl From<bool> for Json {
    fn from(v: bool) -> Self {
        Json::Bool(v)
    }
}
impl From<f64> for Json {
    fn from(v: f64) -> Self {
        Json::Number(v)
    }
}
impl From<f32> for Json {
    fn from(v: f32) -> Self {
        Json::Number(v as f64)
    }
}
impl From<i32> for Json {
    fn from(v: i32) -> Self {
        Json::Number(v as f64)
    }
}
impl From<i64> for Json {
    fn from(v: i64) -> Self {
        Json::Number(v as f64)
    }
}
impl From<usize> for Json {
    fn from(v: usize) -> Self {
        Json::Number(v as f64)
    }
}
impl From<&str> for Json {
    fn from(v: &str) -> Self {
        Json::String(v.to_string())
    }
}
impl From<String> for Json {
    fn from(v: String) -> Self {
        Json::String(v)
    }
}
impl<T: Into<Json>> From<Vec<T>> for Json {
    fn from(v: Vec<T>) -> Self {
        Json::Array(v.into_iter().map(Into::into).collect())
    }
}

pub type JsonResult<T> = Result<T, String>;

impl Json {
    pub fn object() -> Json {
        Json::Object(Vec::new())
    }
    pub fn array(items: Vec<Json>) -> Json {
        Json::Array(items)
    }

    pub fn is_null(&self) -> bool {
        matches!(self, Json::Null)
    }
    pub fn is_bool(&self) -> bool {
        matches!(self, Json::Bool(_))
    }
    pub fn is_number(&self) -> bool {
        matches!(self, Json::Number(_))
    }
    pub fn is_string(&self) -> bool {
        matches!(self, Json::String(_))
    }
    pub fn is_array(&self) -> bool {
        matches!(self, Json::Array(_))
    }
    pub fn is_object(&self) -> bool {
        matches!(self, Json::Object(_))
    }

    pub fn boolean(&self) -> JsonResult<bool> {
        match self {
            Json::Bool(b) => Ok(*b),
            _ => Err("Expected a boolean, not a string or number".into()),
        }
    }
    pub fn number(&self) -> JsonResult<f64> {
        match self {
            Json::Number(n) => Ok(*n),
            _ => Err("Expected a number".into()),
        }
    }
    pub fn string(&self) -> JsonResult<&str> {
        match self {
            Json::String(s) => Ok(s),
            _ => Err("Expected a string".into()),
        }
    }
    pub fn items(&self) -> JsonResult<&Vec<Json>> {
        match self {
            Json::Array(v) => Ok(v),
            _ => Err("Expected an array".into()),
        }
    }
    pub fn items_mut(&mut self) -> JsonResult<&mut Vec<Json>> {
        match self {
            Json::Array(v) => Ok(v),
            _ => Err("Expected an array".into()),
        }
    }
    pub fn members(&self) -> JsonResult<&Vec<(String, Json)>> {
        match self {
            Json::Object(m) => Ok(m),
            _ => Err("Expected an object".into()),
        }
    }

    pub fn find(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(m) => m.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }
    pub fn contains(&self, key: &str) -> bool {
        self.find(key).is_some()
    }
    pub fn at(&self, key: &str) -> JsonResult<&Json> {
        self.find(key).ok_or_else(|| format!("Missing JSON field: {key}"))
    }

    /// The member `key`, inserted as null when absent (a null value becomes an object).
    pub fn entry(&mut self, key: &str) -> &mut Json {
        if self.is_null() {
            *self = Json::object();
        }
        let Json::Object(members) = self else { panic!("Expected an object") };
        if let Some(i) = members.iter().position(|(k, _)| k == key) {
            return &mut members[i].1;
        }
        members.push((key.to_string(), Json::Null));
        &mut members.last_mut().unwrap().1
    }

    /// Sets a member (insertion order kept for existing keys).
    pub fn set(&mut self, key: &str, value: impl Into<Json>) {
        *self.entry(key) = value.into();
    }

    pub fn push(&mut self, value: impl Into<Json>) {
        if self.is_null() {
            *self = Json::Array(Vec::new());
        }
        match self {
            Json::Array(v) => v.push(value.into()),
            _ => panic!("Expected an array"),
        }
    }

    pub fn size(&self) -> usize {
        match self {
            Json::Array(v) => v.len(),
            Json::Object(m) => m.len(),
            _ => 0,
        }
    }

    pub fn parse(text: &str) -> JsonResult<Json> {
        let mut p = Parser { text: text.as_bytes(), source: text, position: 0 };
        let value = p.value(0)?;
        p.skip_space();
        if p.position != text.len() {
            return Err(p.fail("Unexpected trailing characters"));
        }
        Ok(value)
    }

    /// Text; `indent` < 0 for one line.
    pub fn dump(&self, indent: i32) -> String {
        let mut out = String::new();
        self.write(&mut out, indent, 0);
        out
    }

    fn write(&self, out: &mut String, indent: i32, depth: i32) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Number(n) => {
                if !n.is_finite() {
                    out.push_str("null");
                } else if *n == n.floor() && n.abs() < 1e15 {
                    out.push_str(&format!("{n:.0}"));
                } else {
                    out.push_str(&format_g(*n, 17));
                }
            }
            Json::String(s) => write_string(out, s),
            Json::Array(items) => {
                out.push('[');
                let scalar = items.iter().all(|i| !i.is_array() && !i.is_object());
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        out.push_str(if scalar && indent >= 0 { ", " } else { "," });
                    }
                    if !scalar {
                        newline(out, indent, depth + 1);
                    }
                    item.write(out, indent, depth + 1);
                }
                if !scalar && !items.is_empty() {
                    newline(out, indent, depth);
                }
                out.push(']');
            }
            Json::Object(members) => {
                out.push('{');
                for (i, (key, value)) in members.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    newline(out, indent, depth + 1);
                    write_string(out, key);
                    out.push_str(if indent >= 0 { ": " } else { ":" });
                    value.write(out, indent, depth + 1);
                }
                if !members.is_empty() {
                    newline(out, indent, depth);
                }
                out.push('}');
            }
        }
    }
}

impl std::ops::Index<&str> for Json {
    type Output = Json;
    fn index(&self, key: &str) -> &Json {
        static NULL: Json = Json::Null;
        self.find(key).unwrap_or(&NULL)
    }
}

fn write_string(out: &mut String, value: &str) {
    out.push('"');
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

fn newline(out: &mut String, indent: i32, depth: i32) {
    if indent < 0 {
        return;
    }
    out.push('\n');
    out.push_str(&" ".repeat((indent * depth) as usize));
}

struct Parser<'a> {
    text: &'a [u8],
    source: &'a str,
    position: usize,
}

impl Parser<'_> {
    fn fail(&self, message: &str) -> String {
        format!("Invalid JSON at character {}: {message}", self.position)
    }

    fn skip_space(&mut self) {
        while self.position < self.text.len() && matches!(self.text[self.position], b' ' | b'\t' | b'\n' | b'\r') {
            self.position += 1;
        }
    }

    fn consume(&mut self, literal: &str) -> bool {
        if self.source[self.position..].starts_with(literal) {
            self.position += literal.len();
            true
        } else {
            false
        }
    }

    fn value(&mut self, depth: i32) -> JsonResult<Json> {
        if depth > 256 {
            return Err(self.fail("Nesting is too deep"));
        }
        self.skip_space();
        if self.position >= self.text.len() {
            return Err(self.fail("Unexpected end of input"));
        }
        match self.text[self.position] {
            b'{' => return self.object(depth),
            b'[' => return self.array(depth),
            b'"' => return Ok(Json::String(self.string()?)),
            _ => {}
        }
        if self.consume("true") {
            return Ok(Json::Bool(true));
        }
        if self.consume("false") {
            return Ok(Json::Bool(false));
        }
        if self.consume("null") {
            return Ok(Json::Null);
        }
        if self.consume("NaN") {
            return Ok(Json::Number(f64::NAN));
        }
        if self.consume("Infinity") {
            return Ok(Json::Number(f64::INFINITY));
        }
        if self.consume("-Infinity") {
            return Ok(Json::Number(f64::NEG_INFINITY));
        }
        // strtod, which also skips leading spaces and reads "inf"/"nan".
        match parse_prefix_f64(&self.source[self.position..]) {
            Some((v, used)) => {
                self.position += used;
                Ok(Json::Number(v))
            }
            None => Err(self.fail("Expected a value")),
        }
    }

    fn hex(&mut self) -> JsonResult<u32> {
        if self.position + 4 > self.text.len() {
            return Err(self.fail("Truncated unicode escape"));
        }
        let mut code = 0u32;
        for _ in 0..4 {
            let c = self.text[self.position];
            self.position += 1;
            code <<= 4;
            code |= match c {
                b'0'..=b'9' => (c - b'0') as u32,
                b'a'..=b'f' => (c - b'a' + 10) as u32,
                b'A'..=b'F' => (c - b'A' + 10) as u32,
                _ => return Err(self.fail("Invalid unicode escape")),
            };
        }
        Ok(code)
    }

    fn string(&mut self) -> JsonResult<String> {
        self.position += 1;
        let mut out: Vec<u8> = Vec::new();
        loop {
            if self.position >= self.text.len() {
                return Err(self.fail("Unterminated string"));
            }
            let c = self.text[self.position];
            self.position += 1;
            if c == b'"' {
                return Ok(String::from_utf8_lossy(&out).into_owned());
            }
            if c != b'\\' {
                out.push(c);
                continue;
            }
            if self.position >= self.text.len() {
                return Err(self.fail("Unterminated escape"));
            }
            let escape = self.text[self.position];
            self.position += 1;
            match escape {
                b'"' => out.push(b'"'),
                b'\\' => out.push(b'\\'),
                b'/' => out.push(b'/'),
                b'b' => out.push(8),
                b'f' => out.push(12),
                b'n' => out.push(b'\n'),
                b'r' => out.push(b'\r'),
                b't' => out.push(b'\t'),
                b'u' => {
                    let mut code = self.hex()?;
                    if (0xD800..0xDC00).contains(&code) && self.source[self.position..].starts_with("\\u") {
                        self.position += 2;
                        let low = self.hex()?;
                        code = 0x10000 + ((code - 0xD800) << 10) + low.wrapping_sub(0xDC00);
                    }
                    let mut buffer = [0u8; 4];
                    match char::from_u32(code) {
                        Some(ch) => out.extend_from_slice(ch.encode_utf8(&mut buffer).as_bytes()),
                        None => out.extend_from_slice("\u{FFFD}".as_bytes()),
                    }
                }
                _ => return Err(self.fail("Invalid escape")),
            }
        }
    }

    fn array(&mut self, depth: i32) -> JsonResult<Json> {
        self.position += 1;
        let mut items = Vec::new();
        self.skip_space();
        if self.position < self.text.len() && self.text[self.position] == b']' {
            self.position += 1;
            return Ok(Json::Array(items));
        }
        loop {
            items.push(self.value(depth + 1)?);
            self.skip_space();
            if self.position >= self.text.len() {
                return Err(self.fail("Unterminated array"));
            }
            let c = self.text[self.position];
            self.position += 1;
            if c == b']' {
                return Ok(Json::Array(items));
            }
            if c != b',' {
                return Err(self.fail("Expected ',' or ']'"));
            }
        }
    }

    fn object(&mut self, depth: i32) -> JsonResult<Json> {
        self.position += 1;
        let mut result = Json::object();
        self.skip_space();
        if self.position < self.text.len() && self.text[self.position] == b'}' {
            self.position += 1;
            return Ok(result);
        }
        loop {
            self.skip_space();
            if self.position >= self.text.len() || self.text[self.position] != b'"' {
                return Err(self.fail("Expected a member name"));
            }
            let key = self.string()?;
            self.skip_space();
            if self.position >= self.text.len() || self.text[self.position] != b':' {
                return Err(self.fail("Expected ':'"));
            }
            self.position += 1;
            *result.entry(&key) = self.value(depth + 1)?;
            self.skip_space();
            if self.position >= self.text.len() {
                return Err(self.fail("Unterminated object"));
            }
            let c = self.text[self.position];
            self.position += 1;
            if c == b'}' {
                return Ok(result);
            }
            if c != b',' {
                return Err(self.fail("Expected ',' or '}'"));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_and_format() {
        let text = r#"{"a": [1, 2.5, -0, 1e300, NaN], "b": {"c": "x\nyé"}, "a": true}"#;
        let v = Json::parse(text).unwrap();
        // A repeated key replaces the value in place.
        assert_eq!(v.dump(-1), r#"{"a":true,"b":{"c":"x\nyé"}}"#);
        let n = Json::parse("[1, 2.5, -0, 1e300, 0.1, 1e15, 123456789012345]").unwrap();
        assert_eq!(n.dump(-1), "[1,2.5,-0,1.0000000000000001e+300,0.10000000000000001,1000000000000000,123456789012345]");
        assert_eq!(n.dump(2), "[1, 2.5, -0, 1.0000000000000001e+300, 0.10000000000000001, 1000000000000000, 123456789012345]");
        let nested = Json::parse(r#"{"x":[[1,2],{"y":null}]}"#).unwrap();
        assert_eq!(nested.dump(2), "{\n  \"x\": [\n    [1, 2],\n    {\n      \"y\": null\n    }\n  ]\n}");
        assert!(Json::parse("[1,]").unwrap_err().contains("Expected a value"));
    }
}
