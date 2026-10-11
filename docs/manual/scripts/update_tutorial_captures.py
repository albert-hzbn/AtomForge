"""Applies a tutorial recapture to the manual.

    ATOMFORGE_TUTORIAL_CAPTURE=FOLDER AtomForge      # recaptures the screenshots
    python docs/manual/scripts/update_tutorial_captures.py FOLDER

Copies every captured image into docs/manual/figures and writes its crop
bounds (FOLDER/bounds-*.json, recorded by the application from the window the
image shows) into docs/manual/tutorials.json. Viewport captures (results in the
main view) are cropped to the centre of the view. PNG figures are cropped to
their bounds because chapter 5 includes them directly. Then run
render_tutorials.py --write and build_manual.py.
"""
import datetime
import json
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
CATALOGUE = ROOT / "tutorials.json"


def centre_crop(bounds):
    left, top, right, bottom = bounds
    width, height = right - left, bottom - top
    if width <= height * 1.3:
        return bounds
    keep = int(height * 1.05)
    left += (width - keep) // 2
    return [left, top, left + keep, bottom]


def set_jpeg_dpi(path, dpi=96):
    """Marks a JFIF image as `dpi` (render_tutorials.py converts crop bounds
    to points at 96 dpi; the application writes no density)."""
    data = bytearray(path.read_bytes())
    if data[2:4] != b"\xff\xe0" or data[6:11] != b"JFIF\x00":
        raise ValueError(f"Not a JFIF image: {path}")
    data[13] = 1  # dots per inch
    data[14:16] = dpi.to_bytes(2, "big")
    data[16:18] = dpi.to_bytes(2, "big")
    path.write_bytes(bytes(data))


def main(folder):
    folder = Path(folder)
    bounds = {}
    for part in sorted(folder.glob("bounds*.json")):
        bounds.update(json.loads(part.read_text(encoding="utf-8")))
    catalogue = json.loads(CATALOGUE.read_text(encoding="utf-8"))
    used = set()
    updated = 0

    def visit(node):
        nonlocal updated
        if isinstance(node, dict):
            name = node.get("file")
            if isinstance(name, str) and name in bounds:
                box = bounds[name]
                if name.startswith("tutorials/") and name.endswith("-result.jpg"):
                    box = centre_crop(box)
                if name.endswith(".png"):
                    node.pop("bounds", None)
                else:
                    node["bounds"] = box
                used.add(name)
                updated += 1
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)

    visit(catalogue)
    for name, box in bounds.items():
        source, target = folder / name, FIGURES / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".png"):
            with Image.open(source) as image:
                image.crop(tuple(box)).save(target, dpi=(96, 96))
        else:
            shutil.copyfile(source, target)
            set_jpeg_dpi(target)
    unused = sorted(set(bounds) - used)
    catalogue["capture_date"] = datetime.date.today().isoformat()
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    if head.returncode == 0 and head.stdout.strip():
        catalogue["application_source_baseline"] = head.stdout.strip()
    CATALOGUE.write_text(json.dumps(catalogue, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{len(bounds)} images copied, {updated} catalogue entries updated")
    if unused:
        print("not in the catalogue:", ", ".join(unused))


if __name__ == "__main__":
    main(sys.argv[1])
