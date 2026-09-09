"""Check source coverage and regenerate the manual's API/CLI reference.

Refresh only after reviewing changed features and updating their prose.
This detects documentation drift; it does not grade semantic completeness.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess


MANUAL = Path(__file__).resolve().parents[1]
ROOT = MANUAL.parents[1]
SNAPSHOT = MANUAL / "coverage.json"
REFERENCE = MANUAL / "chapters/12-api-signatures.tex"
WIDGET = re.compile(
    r'(?:ImGui::(?:Button|Checkbox|RadioButton|Combo|BeginCombo|MenuItem|'
    r'Input\w*|Drag\w*|Slider\w*|ColorEdit\w*|CollapsingHeader|BeginTabItem)|'
    r'dialogLayout::primaryButton|inputFloat|inputInt|inputVector|combo|slider|'
    r'lightRow|matRow|inputField)\s*\(\s*"((?:\\.|[^"\\])*)"'
)
HELP_NAMES = {
    "bulk": "Bulk", "gb": "GB", "poly": "Poly", "nano": "Nano",
    "amorphous": "Amorphous", "sss": "SSS", "dislocation": "Dislocation",
    "custom": "Custom",
}
PUBLIC_FILES = {
    "python/atomforge/_structure.py": "Atom and Structure",
    "python/atomforge/_io.py": "Structure file functions",
    "python/atomforge/_viewer.py": "Desktop viewer function",
    "python/atomforge/electronic/_grid.py": "Electronic data types and operations",
    "python/atomforge/electronic/_models.py": "Electronic model helpers",
}
FALLBACK = {
    "Atom.copy": "Return a separate atom with the same coordinates and colour.",
    "Structure.__init__": "Create an empty structure with atoms=[] and cell=None.",
    "Structure.copy": "Return a separate structure with copied atoms and cell vectors.",
    "Grid.__init__": "Create a scalar grid. shape contains three integers >=2; values has exactly their product in x-fastest order. cell has three vectors in Angstrom. origin is Cartesian. Periodic grids omit duplicate endpoints; finite grids include endpoints. The grid is limited to 100 million samples.",
    "Grid.shape": "Read-only (nx, ny, nz) tuple.",
    "Grid.cell": "Read-only tuple of the three span/lattice vectors in Angstrom.",
    "Grid.origin": "Read-only Cartesian grid origin in Angstrom.",
    "Grid.periodic": "Read-only Boolean periodic-grid convention.",
    "Grid.unit": "Read-only scalar unit label; do not infer it from a plotted colour.",
    "Grid.name": "Read-only field name, including imported/derived channel identification.",
    "Grid.with_sites": "Return a grid with supplied (atomic_number, x, y, z) site rows; positions are Cartesian Angstrom.",
    "Grid.__add__": "Return pointwise addition of aligned grids.",
    "Grid.__sub__": "Return pointwise subtraction of aligned grids, self minus other.",
    "Grid.__mul__": "Return pointwise grid multiplication or scalar scaling.",
    "Grid.__truediv__": "Return pointwise grid division or scalar division. Scalar zero is rejected.",
    "Grid.resample": "Interpolate this field on target grid geometry. This does not register displaced atoms or correct incompatible units.",
    "Grid.laplacian": "Return the Cartesian Laplacian grid, including cell geometry; units are field units per Angstrom squared.",
    "Grid.smooth": "Return Gaussian-smoothed grid. sigma is in Angstrom; radius bounds the kernel in grid steps.",
    "Grid.line_profile": "Return distance/value rows on the Cartesian segment start to end (Angstrom), using count samples.",
    "Grid.section_values": "Return the first z layer as rows indexed by y, with x-fastest entries. Use on a section result.",
    "Grid.integrate": "Return the whole-grid scalar integral with volume quadrature; density gives electrons.",
    "Grid.integrate_sphere": "Return the scalar integral within a sphere with Cartesian center and radius in Angstrom.",
    "Grid.patterson": "Return the periodic autocorrelation/Patterson field. See the operation reference for normalisation and interpretation.",
    "Grid.isosurface": "Return a Surface at level; optional aligned color grid supplies per-vertex scalar colour values.",
    "Grid.save": "Write this field. Explicit format is xsf, cube, or vasp. Without format, .cube/.cub and .xsf are recognised; other suffixes default to VASP.",
    "Volume.fields": "Read-only tuple of Grid channels. Select with a zero-based index after checking names and units.",
    "Volume.save": "Write the imported volume using an explicit supported format. Check the selected format's channel capabilities.",
    "Surface.vertices": "Return the native triangle vertex rows with scalar data; triangles are not welded into a shared-vertex topology.",
}


def plain(value):
    return value.translate(str.maketrans({"Å": "Angstrom", "×": "x", "—": "--", "–": "-", "≤": "<=", "≥": ">="}))


def tex(value):
    value = plain(value)
    value = re.sub(r"`+([^`]+)`+", r"\1", value)
    value = re.sub(r"\*\*([^*]+)\*\*", r"\1", value)
    value = re.sub(r"\*([A-Za-z_][^*\n]*?)\*", r"\1", value)
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
                    "#": r"\#", "_": r"\_\allowbreak{}", "/": r"/\allowbreak{}", "{": r"\{", "}": r"\}",
                    "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in value)


def source_files():
    paths = set((ROOT / "src/ui").glob("*.*"))
    paths.update((ROOT / "src/algorithms").glob("*.h"))
    paths.update((ROOT / "src/cli").glob("*.*"))
    paths.update((ROOT / "python/atomforge").rglob("*.py"))
    paths.update(ROOT / p for p in ("src/app/EditorApplication.cpp", "src/app/EditorOps.cpp", "CMakeLists.txt"))
    return sorted(p for p in paths if p.is_file() and p.suffix in {".cpp", ".h", ".py", ".txt"})


def documentation_for(path):
    name = path.name
    if "electronic" in path.parts or name.startswith("Electronic"):
        return ["05-electronic.tex", "06-charge.tex", "07-reference.tex", "08-python.tex", "04-display-controls.tex"]
    if "python" in path.parts:
        return ["08-python.tex", "12-api-signatures.tex"]
    if "cli" in path.parts:
        return ["08-python.tex", "12-api-signatures.tex"]
    if any(word in name for word in ("InterfaceBuilder", "CommonNeighbour", "AngularDistribution")):
        return ["11-coverage.tex"]
    if any(word in name for word in ("DislocationBuilder", "Interstitial", "CellSculptor", "MergeStructures", "TransformAtoms", "AtomContext", "EditMenu")):
        return ["03-editing.tex", "03-edit-controls.tex", "04-display-controls.tex"]
    if any(word in name for word in ("Bulk", "CSL", "Nano", "Custom", "Mesh", "PolyCrystal", "Amorphous", "Stacking", "Substitutional")):
        return ["02-builders.tex", "02-builder-controls.tex"]
    return ["01-start.tex", "04-analysis.tex", "04-display-controls.tex", "09-installation.tex"]


def snapshot():
    rows = []
    for path in source_files():
        content = path.read_text(encoding="utf-8-sig")
        controls = []
        for match in WIDGET.finditer(content):
            controls.append({"line": content.count("\n", 0, match.start()) + 1,
                             "declaration_label": match.group(1)})
        rows.append({"source": path.relative_to(ROOT).as_posix(),
                     "sha256": hashlib.sha256(content.encode()).hexdigest(),
                     "documentation": documentation_for(path), "literal_controls": controls})
    return {"scope": "Source drift and feature-family mapping; literal controls include internal IDs. Dynamic labels and runtime visibility require manual review.", "files": rows}


def cli_help():
    source = (ROOT / "src/cli/CLIMode.cpp").read_text(encoding="utf-8-sig")
    result = {}
    for mode, name in HELP_NAMES.items():
        match = re.search(r"static void printHelp" + name + r"\(\)\s*\{(.*?)\n\}", source, re.S)
        if not match:
            raise ValueError(f"Missing CLI help for {mode}")
        literals = re.findall(r'"(?:\\.|[^"\\])*"', match[1])
        result[mode] = "".join(ast.literal_eval(value) for value in literals).strip()
    return result


def api_reference():
    output = [r"\chapter{Complete Python signatures and native CLI options}\label{ch:api-signatures}",
              r"\section{Using this reference}",
              "The signatures below are generated from the public Python implementation. Property entries have no call parentheses. Private helpers and native-handle constructors are implementation details and are excluded. Atom's dataclass constructor is described explicitly. The workflow chapters explain physical interpretation, examples, and validation. Python indices are zero-based unless a function explicitly states otherwise.",
              r"\begin{lstlisting}[language=Python]", "Atom(symbol, x, y, z, r=1.0, g=1.0, b=1.0)", r"\end{lstlisting}",
              tex("Coordinates are Cartesian Angstrom. r/g/b here are colour components, not atomic radii. Structure.atoms is an editable list; Structure.cell is three lattice-vector rows or None. Prefer add_atom for automatic element colours. Structure translation/scaling edits in place; repeat and filter_species return new structures. Grid operations return new data.")]
    count = 1
    for relative, title in PUBLIC_FILES.items():
        tree = ast.parse((ROOT / relative).read_text(encoding="utf-8-sig"))
        output += [r"\section{" + tex(title) + "}"]
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                class_text = ast.get_docstring(node) or ""
                if node.name == "Structure":
                    class_text = "A collection of atoms with an optional periodic unit cell. atoms is a list of Atom objects with Cartesian positions in Angstrom. cell is three lattice-vector rows [[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]], in Angstrom, or None for a nonperiodic cluster."
                output += [r"\subsection{" + tex(node.name) + "}", tex(class_text)]
                entries = [(node.name + ".", child) for child in node.body if isinstance(child, ast.FunctionDef)]
            elif isinstance(node, ast.FunctionDef):
                entries = [("", node)]
            else:
                continue
            for prefix, function in entries:
                name = function.name
                key = prefix + name
                if name.startswith("_") and key not in FALLBACK:
                    continue
                if key in {"Volume.__init__", "Surface.__init__"}:
                    continue
                arguments = ast.unparse(function.args)
                arguments = re.sub(r"^self(?:, )?", "", arguments)
                is_property = any(isinstance(d, ast.Name) and d.id == "property" for d in function.decorator_list)
                signature = key if is_property else f"{key}({arguments})"
                description = ast.get_docstring(function) or FALLBACK.get(key)
                if not description:
                    raise ValueError(f"Add a description for public API {key}")
                listing = "\n".join([r"\begin{lstlisting}[language=Python]", plain(signature), r"\end{lstlisting}"])
                output += [r"\begin{samepage}", listing, tex(description), r"\end{samepage}", ""]
                count += 1
    output += [r"\section{Native builder options}",
               r"Use \code{AtomForge --help MODE} for detailed help. \code{AtomForge --build MODE --help} prints the general overview in this version. These option listings are extracted from the source help and checked against the normal executable. CLI defaults may differ from a dialog's remembered settings. The caret line continuations shown in help are for Windows cmd; use a single line or a subprocess argument list in PowerShell/Python."]
    for mode, help_text in cli_help().items():
        output += [r"\clearpage", r"\subsection{" + tex(mode) + "}",
                   "\n".join([r"\begin{lstlisting}", plain(help_text), r"\end{lstlisting}"])]
    return "\n\n".join(output) + "\n", count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--verify-executable", type=Path)
    args = parser.parse_args()
    reference, api_count = api_reference()
    current = snapshot()
    if args.refresh:
        SNAPSHOT.write_text(json.dumps(current, indent=2) + "\n", encoding="utf-8")
        REFERENCE.write_text(reference, encoding="utf-8")
    else:
        if not SNAPSHOT.exists() or json.loads(SNAPSHOT.read_text(encoding="utf-8")) != current:
            raise SystemExit("Source coverage changed. Review features, update prose, then refresh the coverage snapshot.")
        if not REFERENCE.exists() or REFERENCE.read_text(encoding="utf-8") != reference:
            raise SystemExit("Generated API/CLI appendix is stale; review then refresh.")
    for row in current["files"]:
        for name in row["documentation"]:
            if not (MANUAL / "chapters" / name).is_file():
                raise SystemExit(f"Missing documentation file: {name}")
    if args.verify_executable:
        for mode, expected in cli_help().items():
            result = subprocess.run([str(args.verify_executable.resolve()), "--help", mode],
                                    capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
            if result.stdout.strip() != expected:
                raise SystemExit(f"Executable/source help mismatch: {mode}")
    controls = sum(len(row["literal_controls"]) for row in current["files"])
    print(f"Coverage: {len(current['files'])} source files, {controls} literal control declarations, "
          f"{api_count} public API entries, {len(HELP_NAMES)} builder CLI modes.")


if __name__ == "__main__":
    main()
