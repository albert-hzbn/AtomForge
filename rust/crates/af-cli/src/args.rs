//! Command-line helpers shared by the modes (the C++ cli/CliArgs and
//! cli/OptionCheck), C++ stream-style parsing and the text-mode output of the
//! C++ streams (CRLF line ends on Windows).

use std::collections::BTreeSet;
use std::io::Write;

use af_core::format::{format_g, parse_prefix_f32, parse_prefix_f64, stod};
use af_core::suggest;

/// Text as a C++ text-mode stream writes it: "\n" becomes "\r\n" on Windows.
pub fn text_mode(text: &str) -> String {
    if cfg!(windows) {
        text.replace('\n', "\r\n")
    } else {
        text.to_string()
    }
}

thread_local! {
    /// (stdout, stderr) text while a pipeline step runs a mode.
    static CAPTURE: std::cell::RefCell<Option<(String, String)>> = const { std::cell::RefCell::new(None) };
}

/// Runs `work` with std::cout and std::cerr redirected to strings, as the
/// C++ rdbuf swap; returns the result and the captured (stdout, stderr).
pub fn captured<T>(work: impl FnOnce() -> T) -> (T, String, String) {
    CAPTURE.with(|c| *c.borrow_mut() = Some((String::new(), String::new())));
    let value = work();
    let (o, e) = CAPTURE.with(|c| c.borrow_mut().take()).unwrap_or_default();
    (value, o, e)
}

/// std::cout << text.
pub fn out(text: &str) {
    if CAPTURE.with(|c| c.borrow_mut().as_mut().map(|(o, _)| o.push_str(text)).is_some()) {
        return;
    }
    let mut stdout = std::io::stdout().lock();
    let _ = stdout.write_all(text_mode(text).as_bytes());
    let _ = stdout.flush();
}

/// std::cerr << text.
pub fn err(text: &str) {
    if CAPTURE.with(|c| c.borrow_mut().as_mut().map(|(_, e)| e.push_str(text)).is_some()) {
        return;
    }
    let _ = std::io::stderr().lock().write_all(text_mode(text).as_bytes());
}

/// A double as `stream << value` prints it (precision 6, %g).
pub fn g(value: f64) -> String {
    format_g(value, 6)
}

/// A double with precision 17 (%.17g), as `setprecision(17)` streams.
pub fn g17(value: f64) -> String {
    format_g(value, 17)
}

/// Writes a file in text mode (CRLF on Windows), as a C++ std::ofstream.
pub fn write_text_file(path: &str, text: &str) -> std::io::Result<()> {
    std::fs::write(path, text_mode(text))
}

/// argv without the program name.
pub struct Args(pub Vec<String>);

impl Args {
    /// Value following `flag` (not the last argument), or None.
    pub fn find(&self, flag: &str) -> Option<&str> {
        let n = self.0.len();
        (0..n.saturating_sub(1)).find(|&i| self.0[i] == flag).map(|i| self.0[i + 1].as_str())
    }

    pub fn has(&self, flag: &str) -> bool {
        self.0.iter().any(|a| a == flag)
    }

    /// Values after every occurrence of a repeatable `flag`.
    pub fn find_all(&self, flag: &str) -> Vec<String> {
        let n = self.0.len();
        (0..n.saturating_sub(1)).filter(|&i| self.0[i] == flag).map(|i| self.0[i + 1].clone()).collect()
    }

    /// A finite number; `default` when the flag is absent.
    pub fn double(&self, flag: &str, default: f64) -> Result<f64, String> {
        let Some(v) = self.find(flag) else { return Ok(default) };
        match stod(v) {
            Ok((value, used)) if used == v.len() && value.is_finite() => Ok(value),
            _ => Err(format!("invalid numeric value '{v}' for {flag}")),
        }
    }

    /// An int; `default` when the flag is absent.
    pub fn int(&self, flag: &str, default: i32) -> Result<i32, String> {
        let Some(v) = self.find(flag) else { return Ok(default) };
        match stoi(v) {
            Some((value, used)) if used == v.len() => Ok(value),
            _ => Err(format!("invalid integer value '{v}' for {flag}")),
        }
    }
}

/// std::stoi: leading spaces, sign, digits; the value and the bytes used.
pub fn stoi(text: &str) -> Option<(i32, usize)> {
    let b = text.as_bytes();
    let mut i = 0;
    while i < b.len() && is_space(b[i]) {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits {
        return None;
    }
    let value: i64 = text[start..i].parse().ok()?;
    i32::try_from(value).ok().map(|v| (v, i))
}

/// %.12g (setprecision(12)).
pub fn format_g12(value: f64) -> String {
    format_g(value, 12)
}

/// atomicNumberFromSymbol: any letter case; -1 when empty or unknown.
pub fn symbol_z(symbol: &str) -> i32 {
    if symbol.is_empty() {
        return -1;
    }
    (1..=118).find(|&z| af_core::elements::symbol(z).eq_ignore_ascii_case(symbol)).unwrap_or(-1)
}

/// A structure file, as loadStructureFromFile.
pub fn load(path: &str) -> Result<af_core::structure::Structure, String> {
    af_io::loader::load_structure_from_file(path)
}

/// Saves with the format of the file name; false on failure.
pub fn save(structure: &af_core::structure::Structure, path: &str, format: &str) -> bool {
    af_io::loader::save_structure(structure, path, format).is_ok()
}

pub fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
}

/// An std::istringstream over a string: `>>` words and numbers, with the
/// stream's fail state.
pub struct IStream<'a> {
    text: &'a str,
    at: usize,
    pub failed: bool,
}

impl<'a> IStream<'a> {
    pub fn new(text: &'a str) -> Self {
        Self { text, at: 0, failed: false }
    }

    fn skip_space(&mut self) {
        let b = self.text.as_bytes();
        while self.at < b.len() && is_space(b[self.at]) {
            self.at += 1;
        }
    }

    /// >> std::string.
    pub fn word(&mut self) -> Option<String> {
        if self.failed {
            return None;
        }
        self.skip_space();
        let b = self.text.as_bytes();
        let start = self.at;
        while self.at < b.len() && !is_space(b[self.at]) {
            self.at += 1;
        }
        if self.at == start {
            self.failed = true;
            return None;
        }
        Some(self.text[start..self.at].to_string())
    }

    /// >> int.
    pub fn int(&mut self) -> Option<i32> {
        if self.failed {
            return None;
        }
        self.skip_space();
        match stoi(&self.text[self.at..]) {
            Some((v, used)) => {
                self.at += used;
                Some(v)
            }
            None => {
                self.failed = true;
                None
            }
        }
    }

    /// >> double.
    pub fn double(&mut self) -> Option<f64> {
        if self.failed {
            return None;
        }
        self.skip_space();
        match parse_prefix_f64(&self.text[self.at..]) {
            Some((v, used)) if used > 0 => {
                self.at += used;
                Some(v)
            }
            _ => {
                self.failed = true;
                None
            }
        }
    }

    /// >> float.
    pub fn float(&mut self) -> Option<f32> {
        if self.failed {
            return None;
        }
        self.skip_space();
        match parse_prefix_f32(&self.text[self.at..]) {
            Some((v, used)) if used > 0 => {
                self.at += used;
                Some(v)
            }
            _ => {
                self.failed = true;
                None
            }
        }
    }

    /// >> char (one non-space character).
    pub fn char(&mut self) -> Option<char> {
        if self.failed {
            return None;
        }
        self.skip_space();
        match self.text[self.at..].chars().next() {
            Some(c) => {
                self.at += c.len_utf8();
                Some(c)
            }
            None => {
                self.failed = true;
                None
            }
        }
    }
}

/// "x y z" integers (parseIvec3).
pub fn parse_ivec3(text: &str) -> Option<[i32; 3]> {
    let mut s = IStream::new(text);
    Some([s.int()?, s.int()?, s.int()?])
}

/// "x y z" floats (parseVec3).
pub fn parse_vec3(text: &str) -> Option<[f32; 3]> {
    let mut s = IStream::new(text);
    Some([s.float()?, s.float()?, s.float()?])
}

/// Output format from a file name's extension (cif when there is none).
pub fn detect_format(path: &str) -> String {
    let name = match path.rfind(['/', '\\']) {
        Some(i) => &path[i + 1..],
        None => path,
    };
    let ext = name.rfind('.').map(|d| name[d + 1..].to_ascii_lowercase()).unwrap_or_default();
    let upper = name.to_ascii_uppercase();
    let vasp_name = upper.contains("POSCAR") || upper.contains("CONTCAR");
    let other_format = !ext.is_empty()
        && ext != "vasp"
        && ext != "poscar"
        && ext != "contcar"
        && ext.bytes().any(|c| !c.is_ascii_digit())
        && (af_io::loader::is_supported_structure_file(&format!("file.{ext}")) || matches!(ext.as_str(), "pwi" | "gjf" | "lammps" | "json"));
    if vasp_name && !other_format {
        return "vasp".into();
    }
    match ext.as_str() {
        "" => "cif".into(),
        "vasp" | "poscar" | "contcar" => "vasp".into(),
        "lmp" | "lammps" => "lammps".into(),
        _ => ext,
    }
}

/// The --flags a help text names.
pub fn flags_in_help(help: &str) -> BTreeSet<String> {
    let b = help.as_bytes();
    let mut flags = BTreeSet::new();
    let mut search = 0;
    while let Some(found) = help[search..].find("--") {
        let at = search + found;
        let mut end = at + 2;
        while end < b.len() && (b[end].is_ascii_alphanumeric() || b[end] == b'-') {
            end += 1;
        }
        if end > at + 2 && b[at + 2].is_ascii_alphabetic() {
            flags.insert(help[at..end].to_string());
        }
        search = at + 2;
    }
    flags
}

/// "" when every --option is one the mode's help names, else the message
/// (pointing at `help_command`, by default "AtomForge --help <mode>").
pub fn unknown_option(args: &[String], help: &str, mode: &str, help_command: &str) -> String {
    let mut known = flags_in_help(help);
    for common in ["--build", "--output", "--input", "--help", "--analyze"] {
        known.insert(common.to_string());
    }
    let command = if help_command.is_empty() { format!("AtomForge --help {mode}") } else { help_command.to_string() };
    for arg in args {
        let b = arg.as_bytes();
        if b.len() >= 2 && b[0] == b'-' && b[1] != b'-' && b[1].is_ascii_alphabetic() && arg != "-h" {
            let suggestion = if b.len() > 2 { suggest::closest(&format!("-{arg}"), &known) } else { String::new() };
            return format!(
                "unknown option {arg} for {mode} (options start with --{}); see {command}",
                if suggestion.is_empty() { String::new() } else { format!("; did you mean {suggestion}?") }
            );
        }
        if b.len() < 3 || !arg.starts_with("--") || !b[2].is_ascii_alphabetic() {
            continue;
        }
        let flag = arg.split('=').next().unwrap_or(arg);
        if known.contains(flag) {
            continue;
        }
        let suggestion = suggest::closest(flag, &known);
        return format!(
            "unknown option {flag} for {mode}{}; see {command}",
            if suggestion.is_empty() { String::new() } else { format!(" (did you mean {suggestion}?)") }
        );
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(detect_format("out/POSCAR"), "vasp");
        assert_eq!(detect_format("a.b/cu"), "cif");
        assert_eq!(detect_format("x.XYZ"), "xyz");
        assert_eq!(detect_format("POSCAR.xyz"), "xyz");
        assert_eq!(detect_format("POSCAR-001"), "vasp");
        let mut s = IStream::new(" Cu 0.5 1/3");
        assert_eq!(s.word().as_deref(), Some("Cu"));
        assert_eq!(s.double(), Some(0.5));
        assert_eq!(s.double(), Some(1.0));
        assert!(s.double().is_none());
        assert_eq!(parse_ivec3("1 -1 0"), Some([1, -1, 0]));
        let help = "  --a <x>\n  --spacegroup N\n";
        assert_eq!(unknown_option(&["--spacegrup".into()], help, "bulk", ""), "unknown option --spacegrup for bulk (did you mean --spacegroup?); see AtomForge --help bulk");
    }
}
