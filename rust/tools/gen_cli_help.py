"""Generates crates/af-cli/src/help.rs: the help texts of the C++ CLI
(src/cli/*.cpp), byte for byte. A help function qualifies when its body is
`std::cout << "literal" "literal" ... [<< std::endl];`. Run from the
repository root: python rust/tools/gen_cli_help.py"""
import glob
import os
import re

ESCAPES = {"n": "\n", "t": "\t", '"': '"', "\\": "\\", "'": "'", "0": "\0", "r": "\r"}

# (file stem, function) -> Rust constant.
NAMES = {
    ("CLIMode", "printHelp"): "MAIN",
    ("PipelineCLI", "printHelp"): "PIPE",
    ("RenderCLI", "printHelp"): "RENDER",
    ("CloudCLI", "printHelp"): "CLOUD",
    ("BuildModifiers", "helpSculpt"): "SCULPT",
}


def decode(literal):
    out, i = [], 0
    while i < len(literal):
        c = literal[i]
        if c == "\\":
            out.append(ESCAPES[literal[i + 1]])
            i += 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


def escape(line):
    return '"' + line.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\t", "\\t") + '"'


def rust_string(text):
    # concat!() of one literal per line (a trailing backslash continuation
    # would swallow the next line's indentation).
    lines = text.splitlines(keepends=True) or [""]
    return "concat!(\n    " + ",\n    ".join(escape(line) for line in lines) + ",\n)"


TOKEN = re.compile(r'\s*(?:"((?:[^"\\]|\\.)*)"|(<<)|(std::endl)|(;))', re.S)

helps = {}
for path in sorted(glob.glob("src/cli/*.cpp")):
    stem = os.path.splitext(os.path.basename(path))[0]
    source = open(path, encoding="utf-8").read()
    for match in re.finditer(r"void (\w+)\(\)\s*\{\s*std::cout\s*<<", source):
        name, at, parts, pure = match.group(1), match.end(), [], True
        while True:
            token = TOKEN.match(source, at)
            if not token:
                pure = False
                break
            at = token.end()
            if token.group(1) is not None:
                parts.append(decode(token.group(1)))
            elif token.group(3):
                parts.append("\n")
            elif token.group(4):
                break
        if pure and source[at:at + 40].lstrip().startswith("}"):
            constant = NAMES.get((stem, name)) or name.replace("printHelp", "").upper()
            helps[constant] = "".join(parts)
        else:
            print("skipped (not a plain literal):", path, name)

out = ["//! Help texts of the command line, generated from the C++ src/cli by",
       "//! rust/tools/gen_cli_help.py (do not edit).", ""]
for constant in sorted(helps):
    out.append("pub const %s: &str = %s;" % (constant, rust_string(helps[constant])))
    out.append("")
with open("rust/crates/af-cli/src/help.rs", "w", newline="\n", encoding="utf-8") as f:
    f.write("\n".join(out))
print(len(helps), "help texts:", ", ".join(sorted(helps)))
