"""atomforge.render with the desktop's display settings (display=)."""

import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, os.path.dirname(__file__))
import atomforge as af
from atomforge.render import _display_file, _display_text, render
from _optional import requires


def copper():
    s = af.Structure()
    a = 3.615
    for f in ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)):
        s.add_atom("Cu", f[0] * a, f[1] * a, f[2] * a)
    s.cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    return s


class DisplayTextTests(unittest.TestCase):
    def test_mapping_becomes_saved_default_text(self):
        text = _display_text({"ambient": 0.5, "light_follows_camera": True, "background": (0, 0, 1),
                              "colour.Cu": [0.9, 0.5, 0.2], "radius.O": 0.8})
        self.assertEqual(text.splitlines(), ["ambient=0.5", "light_follows_camera=true", "background=0.0 0.0 1.0",
                                             "colour.Cu=0.9 0.5 0.2", "radius.O=0.8"])
        with self.assertRaises(ValueError):
            _display_text({"background": (0, 0)})
        with self.assertRaises(ValueError):
            _display_text({"ambient": "1\nbackground=0 0 0"})

    def test_files_and_saved(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertEqual(_display_file("saved", folder), "saved")
            written = Path(_display_file({"atom_scale": 0.5}, folder))
            self.assertEqual(written.read_text(encoding="utf-8"), "atom_scale=0.5\n")
            self.assertEqual(_display_file(written, folder), str(written))
            with self.assertRaises(FileNotFoundError):
                _display_file(Path(folder) / "missing.txt", folder)


@unittest.skipUnless(os.environ.get("ATOMFORGE_PATH"), "needs the AtomForge executable")
@requires("PIL")
class RenderWithDisplayTests(unittest.TestCase):
    def test_background_and_element_colour(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "look.png"
            try:
                render(copper(), out, width=80, height=80, yaw=20, pitch=15,
                       display={"background": (0, 0, 0), "colour.Cu": (0, 1, 0), "ambient": 1, "shadow_strength": 0})
            except Exception as error:  # no GPU/display on this machine
                if "OpenGL context" in str(getattr(error, "stderr", "")):
                    self.skipTest("no OpenGL context")
                raise
            image = Image.open(out).convert("RGB")
            self.assertEqual(image.getpixel((0, 0)), (0, 0, 0))
            self.assertTrue(any(g > 120 and r < 80 and b < 80 for r, g, b in zip(*[iter(image.tobytes())] * 3)),
                            "copper drawn in the display settings' green")
            # background= still applies on top of display=.
            render(copper(), out, width=40, height=40, display={"background": (0, 0, 0)}, background=(1, 1, 1))
            self.assertEqual(Image.open(out).convert("RGB").getpixel((0, 0)), (255, 255, 255))


if __name__ == "__main__":
    unittest.main()
