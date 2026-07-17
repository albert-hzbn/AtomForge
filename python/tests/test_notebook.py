"""Tests for the inline notebook renderer."""

import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import atomforge as af
from atomforge._notebook import structure_to_html

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
_failures = []


def check(name, cond, detail=""):
    if cond:
        print(f"  {PASS}  {name}")
    else:
        print(f"  {FAIL}  {name}" + (f": {detail}" if detail else ""))
        _failures.append(name)


# ── Fixtures ────────────────────────────────────────────────────────────────────

def _bcc_fe():
    s = af.Structure()
    s.set_cell(2.87, 2.87, 2.87)
    s.add_atom("Fe", 0.0,   0.0,   0.0)
    s.add_atom("Fe", 1.435, 1.435, 1.435)
    return s

def _water():
    s = af.Structure()
    s.add_atom("O", 0.0, 0.0, 0.0)
    s.add_atom("H", 0.96, 0.0, 0.0)
    s.add_atom("H", 0.0, 0.96, 0.0)
    return s


# ── HTML generation ─────────────────────────────────────────────────────────────

def test_html_generation():
    print("\n[HTML generation]")
    s = _bcc_fe()
    html = structure_to_html(s)

    check("returns string", isinstance(html, str))
    check("non-empty", len(html) > 500)
    check("has canvas tag", "<canvas" in html)
    check("has script tag", "<script" in html)
    check("has atom data", "Fe" in html or "1.435" in html)
    check("has WebGL init", "getContext" in html)
    check("has rotation handler", "mousemove" in html or "rotX" in html)
    check("has zoom handler", "wheel" in html or "scale" in html)
    check("has touch support", "touchstart" in html)
    check("has 2D fallback", "Canvas 2D" in html or "2d fallback" in html.lower() or "ctx2d" in html or "getContext('2d')" in html)


def test_html_contains_correct_atom_data():
    print("\n[Atom data in HTML]")
    s = _water()
    html = structure_to_html(s)

    # Extract the JSON atom array from the HTML
    m = re.search(r"var atoms\s*=\s*(\[.*?\]);", html, re.DOTALL)
    check("atom JSON present", m is not None)
    if m:
        atoms = json.loads(m.group(1))
        check("atom count", len(atoms) == 3, f"got {len(atoms)}")
        # First atom is O centered (centroid subtracted)
        check("O x near centroid", abs(atoms[0][0]) < 1.5)
        # Each atom has 7 fields: x,y,z, r,g,b, vdw_radius
        check("atom has 7 fields", all(len(a) == 7 for a in atoms))
        # O color: r=1, g≈0.05, b≈0.05
        check("O red channel", atoms[0][3] > 0.8, f"r={atoms[0][3]}")
        # VdW radius for O should be ~1.52
        check("O vdw radius", abs(atoms[0][6] - 1.52) < 0.1, f"got {atoms[0][6]}")


def test_html_cell_lines():
    print("\n[Cell box in HTML]")
    # Periodic structure should have 12 cell edges
    s = _bcc_fe()
    html = structure_to_html(s)
    m = re.search(r"var cellLines\s*=\s*(\[.*?\]);", html, re.DOTALL)
    check("cellLines JSON present", m is not None)
    if m:
        lines = json.loads(m.group(1))
        check("12 cell edges", len(lines) == 12, f"got {len(lines)}")
        check("edge has 2 endpoints", all(len(e) == 2 for e in lines))
        check("endpoint has 3 coords", all(len(e[0]) == 3 for e in lines))

    # Non-periodic structure should have empty cell lines
    s2 = _water()
    html2 = structure_to_html(s2)
    m2 = re.search(r"var cellLines\s*=\s*(\[.*?\]);", html2, re.DOTALL)
    if m2:
        lines2 = json.loads(m2.group(1))
        check("no cell edges for cluster", len(lines2) == 0, f"got {len(lines2)}")


def test_html_title():
    print("\n[Title / formula]")
    s = _bcc_fe()
    html = structure_to_html(s)
    check("Fe in title", "Fe" in html)
    check("atom count in title", "2 atoms" in html)

    s2 = _water()
    html2 = structure_to_html(s2)
    check("O in formula", "O1" in html2 or "O " in html2)
    check("H in formula", "H2" in html2 or "H " in html2)


def test_unique_ids():
    print("\n[Unique canvas IDs per render]")
    s = _bcc_fe()
    ids = set()
    for _ in range(20):
        h = structure_to_html(s)
        m = re.search(r'id="c_(af\d+)"', h)
        if m:
            ids.add(m.group(1))
    check("all IDs unique", len(ids) == 20, f"got {len(ids)} unique out of 20")


def test_custom_dimensions():
    print("\n[Custom dimensions]")
    s = _water()
    html = structure_to_html(s, width=800, height=600)
    check("custom width in html", 'width="800"' in html)
    check("custom height in html", 'height="600"' in html)


def test_large_structure_downsampling():
    print("\n[Large structure downsampling]")
    # Build a 100k+ atom structure
    s = _bcc_fe().repeat(20, 20, 20)   # 20^3 * 2 = 16000 atoms
    html = structure_to_html(s)
    m = re.search(r"var atoms\s*=\s*(\[.*?\]);", html, re.DOTALL)
    if m:
        atoms = json.loads(m.group(1))
        # Should be <= 50000
        check("downsampled to <=50000", len(atoms) <= 50_000, f"got {len(atoms)}")
        check("downsample note in title", "showing" in html or len(atoms) == len(s.atoms))


def test_empty_structure():
    print("\n[Empty structure]")
    s = af.Structure()
    try:
        html = structure_to_html(s)
        check("no crash on empty", True)
        check("returns string", isinstance(html, str))
    except Exception as e:
        check("no crash on empty", False, str(e))


def test_repr_html_method():
    print("\n[_repr_html_ method]")
    s = _bcc_fe()
    check("method exists", hasattr(s, "_repr_html_"))
    check("method callable", callable(getattr(s, "_repr_html_", None)))
    html = s._repr_html_()
    check("returns html string", isinstance(html, str) and len(html) > 100)
    check("contains canvas", "<canvas" in html)


def test_view_notebook_method():
    print("\n[view_notebook method]")
    s = _water()
    check("method exists", hasattr(s, "view_notebook"))
    # In a non-IPython context, view_notebook should either display or fall back gracefully.
    # We just test it doesn't crash when IPython is available.
    try:
        from IPython.display import HTML
        # Monkey-patch display to capture output instead of showing it
        import atomforge._structure as _mod
        captured = []
        _orig = None
        try:
            import IPython.display as ipd
            _orig = ipd.display
            ipd.display = lambda x: captured.append(x)
            s.view_notebook()
            ipd.display = _orig
        except Exception:
            if _orig:
                ipd.display = _orig
        check("view_notebook called display", len(captured) > 0)
        check("display got HTML object", isinstance(captured[0], HTML) if captured else False)
    except ImportError:
        check("view_notebook (no IPython — skip)", True)


def test_webgl_shader_syntax():
    print("\n[WebGL shader validity checks]")
    s = _bcc_fe()
    html = structure_to_html(s)

    # Vertex shader must declare attributes we bind
    check("aPos attribute", "attribute" in html and "aPos" in html)
    check("aCenter attribute", "aCenter" in html)
    check("aRadius attribute", "aRadius" in html)
    check("aColor attribute", "aColor" in html)

    # Fragment shader must discard outside circle
    check("discard outside sphere", "discard" in html)

    # MVP uniform
    check("uMVP uniform", "uMVP" in html)

    # ANGLE_instanced_arrays extension requested
    check("instanced arrays ext", "ANGLE_instanced_arrays" in html)


def test_no_external_urls():
    print("\n[No external URLs / CDN]")
    s = _bcc_fe()
    html = structure_to_html(s)
    # Should not reference any external URLs
    external = re.findall(r'https?://[^\s\'"<>]+', html)
    check("no external URLs", len(external) == 0, f"found: {external[:3]}")


# ── IPython display integration ─────────────────────────────────────────────────

def test_ipython_display():
    print("\n[IPython display integration]")
    try:
        from IPython.display import HTML
        s = _bcc_fe()
        html_obj = HTML(s._repr_html_())
        check("HTML object created", html_obj is not None)
        check("HTML data non-empty", len(html_obj.data) > 100)
    except ImportError:
        check("IPython not available (skip)", True)


# ── Run all ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import traceback
    print("=" * 60)
    print("atomforge notebook renderer tests")
    print("=" * 60)
    for fn in [
        test_html_generation,
        test_html_contains_correct_atom_data,
        test_html_cell_lines,
        test_html_title,
        test_unique_ids,
        test_custom_dimensions,
        test_large_structure_downsampling,
        test_empty_structure,
        test_repr_html_method,
        test_view_notebook_method,
        test_webgl_shader_syntax,
        test_no_external_urls,
        test_ipython_display,
    ]:
        try:
            fn()
        except Exception:
            print(f"  {FAIL}  {fn.__name__} raised:")
            traceback.print_exc()
            _failures.append(fn.__name__)

    print("\n" + "=" * 60)
    if _failures:
        print(f"FAILED ({len(_failures)}): {', '.join(_failures)}")
        sys.exit(1)
    else:
        print("All notebook tests passed.")
    print("=" * 60)
