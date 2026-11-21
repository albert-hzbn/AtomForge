"""Side-by-side images of the C++ and Rust UI tours, to compare the two
interfaces: tmp/cpp-tour (ATOMFORGE_UI_TOUR of the C++ desktop) on the left,
tmp/rust-tour (rust/tools/ui_tour.afscript) on the right, written to
tmp/tour-compare. A missing Rust shot shows as magenta.

    python rust/tools/compare_tour.py [NAME-FILTER]
"""
import os
import sys

from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
cpp, rust, out = (os.path.join(ROOT, "tmp", d) for d in ("cpp-tour", "rust-tour", "tour-compare"))
os.makedirs(out, exist_ok=True)
filters = sys.argv[1:]
for name in sorted(os.listdir(cpp)):
    if not name.endswith(".png") or (filters and not any(f in name for f in filters)):
        continue
    a = Image.open(os.path.join(cpp, name)).convert("RGB")
    path = os.path.join(rust, name)
    b = Image.open(path).convert("RGB") if os.path.exists(path) else Image.new("RGB", a.size, (255, 0, 255))
    w = 1400
    ha, hb = int(a.height * w / a.width), int(b.height * w / b.width)
    image = Image.new("RGB", (w * 2 + 10, max(ha, hb)), (255, 255, 255))
    image.paste(a.resize((w, ha)), (0, 0))
    image.paste(b.resize((w, hb)), (w + 10, 0))
    image.save(os.path.join(out, name))
print(out)
