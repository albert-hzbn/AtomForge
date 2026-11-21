"""Runs the desktop's UI test scripts (rust/tests/ui/*.afscript) with
ATOMFORGE_UI_SCRIPT: each script drives the app (open, menu operations,
picking, undo, tools, saving) and checks the result; the app exits with code
3 when a check fails. Run from the repository root:

    python rust/tools/run_ui_tests.py [APP_EXE] [-k FILTER]
"""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-k")]
    filters = [a[2:] for a in sys.argv[1:] if a.startswith("-k")]
    app = Path(args[0]) if args else ROOT / "rust/target/release" / ("AtomForge.exe" if sys.platform == "win32" else "AtomForge")
    (ROOT / "rust/target/ui").mkdir(parents=True, exist_ok=True)
    failed = 0
    for script in sorted((ROOT / "rust/tests/ui").glob("*.afscript")):
        if filters and not any(f in script.stem for f in filters):
            continue
        env = dict(os.environ, ATOMFORGE_UI_SCRIPT=str(script))
        done = subprocess.run([str(app)], cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
        ok = done.returncode == 0
        failed += not ok
        print("%-4s %s%s" % ("ok" if ok else "FAIL", script.stem, "" if ok else " (exit %d)\n     %s" % (done.returncode, done.stderr.strip())))
    print("%d scripts failed" % failed)
    sys.exit(1 if failed else 0)


main()
