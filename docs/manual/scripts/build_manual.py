"""Compile the checked-in LaTeX manual; no original CHGCAR is required."""
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    manual = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(manual / "scripts/check_coverage.py")], check=True)
    subprocess.run([sys.executable, str(manual / "scripts/render_tutorials.py")], check=True)
    build = manual / "build"
    build.mkdir(exist_ok=True)
    executable = shutil.which("pdflatex")
    if not executable:
        raise SystemExit("pdfLaTeX was not found. Install TeX Live or MiKTeX and add it to PATH.")
    command = [executable, "-interaction=nonstopmode", "-halt-on-error",
               "-file-line-error", "-output-directory=build", "manual.tex"]
    for iteration in range(1, 4):
        result = subprocess.run(command, cwd=manual, capture_output=True, text=True,
                                encoding="utf-8", errors="replace")
        (build / f"pass-{iteration}.log").write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise SystemExit(f"LaTeX pass {iteration} failed. See {build / f'pass-{iteration}.log'}")
    log = (build / "manual.log").read_text(encoding="utf-8", errors="replace")
    problems = [line for line in log.splitlines()
                if any(marker in line for marker in (
                    "Overfull", "undefined", "Missing character:",
                    "LaTeX Warning:", "pdfTeX warning", "Package enumitem Warning:"))]
    if problems:
        raise SystemExit("Resolve LaTeX diagnostics before publishing:\n" + "\n".join(problems))
    output = manual / "AtomForge-manual.pdf"
    # Publish the validated file without leaving a second manual in build/.
    (build / "manual.pdf").replace(output)
    print(f"Built {output} ({output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
