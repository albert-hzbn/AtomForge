"""Headless --render CLI test: python tests/render_cli_smoke.py path/to/AtomForge.

Exits 77 (ctest SKIP_RETURN_CODE) if no GL context/display is available,
mirroring electronic_viewport's GPU-less skip convention.
"""
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import zlib


def png_pixels(path):
    """(width, height, rows of (r, g, b)) of an 8-bit RGB/RGBA PNG, without dependencies."""
    data = Path(path).read_bytes()
    position, chunks, header = 8, [], None
    while position < len(data):
        length, kind = struct.unpack(">I4s", data[position:position + 8])
        body = data[position + 8:position + 8 + length]
        if kind == b"IHDR":
            header = struct.unpack(">IIBBBBB", body)
        elif kind == b"IDAT":
            chunks.append(body)
        position += 12 + length
    width, height, depth, colour, _, _, interlace = header
    if depth != 8 or colour not in (2, 6) or interlace:
        raise AssertionError(f"unexpected PNG format {header}")
    channels = 3 if colour == 2 else 4
    raw, stride, rows, previous = zlib.decompress(b"".join(chunks)), width * channels, [], bytearray(width * channels)
    for y in range(height):
        kind, line = raw[y * (stride + 1)], bytearray(raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)])
        for x in range(stride):
            left = line[x - channels] if x >= channels else 0
            up, corner = previous[x], previous[x - channels] if x >= channels else 0
            if kind == 1: line[x] = (line[x] + left) & 255
            elif kind == 2: line[x] = (line[x] + up) & 255
            elif kind == 3: line[x] = (line[x] + (left + up) // 2) & 255
            elif kind == 4:
                guess = left + up - corner
                best = min((abs(guess - left), 0, left), (abs(guess - up), 1, up), (abs(guess - corner), 2, corner))[2]
                line[x] = (line[x] + best) & 255
        rows.append([tuple(line[i:i + 3]) for i in range(0, stride, channels)])
        previous = line
    return width, height, rows

exe = str(Path(sys.argv[1]).resolve())

with tempfile.TemporaryDirectory(prefix="atomforge_render_") as folder:
    root = Path(folder)
    source = root / "cu.cif"
    subprocess.run(
        [exe, "--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", str(source)],
        check=True, capture_output=True, text=True, timeout=60)

    help_result = subprocess.run([exe, "--render", "--help"], capture_output=True, text=True, timeout=30)
    if help_result.returncode != 0 or "--render" not in help_result.stdout:
        raise AssertionError(f"--render --help failed:\n{help_result.stdout}\n{help_result.stderr}")

    output = root / "cu.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(output),
         "--width", "64", "--height", "64", "--color", "Cu 0.9 0.5 0.1"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        if "OpenGL context" in result.stderr:
            print("Skipping: no GPU/display available for --render:", result.stderr.strip())
            sys.exit(77)
        raise AssertionError(f"--render failed: exit {result.returncode}\n{result.stdout}\n{result.stderr}")

    data = output.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError("Rendered file is not a valid PNG")
    if len(data) < 100:
        raise AssertionError("Rendered PNG is suspiciously small")

    # The help's form: --color SYMBOL R G B and --background R G B as separate words.
    words = root / "words.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(words), "--width", "64", "--height", "64",
         "--color", "cu", "0.2", "0.5", "0.9", "--background", "0", "0", "0", "--radius", "Cu", "1.0"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0 or not words.exists():
        raise AssertionError(f"--render rejected separate-word values:\n{result.stdout}\n{result.stderr}")

    # Turntable mode writes numbered frames instead of a single file.
    turntable_base = root / "turn.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(turntable_base),
         "--width", "64", "--height", "64", "--frames", "2", "--yaw-step", "45"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise AssertionError(f"--render turntable failed: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
    for index in range(2):
        frame_path = root / f"turn-{index:03d}.png"
        if not frame_path.is_file() or frame_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise AssertionError(f"Missing or invalid turntable frame {frame_path}")

    dpi_output = root / "cu_dpi.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(dpi_output),
         "--width", "32", "--height", "32", "--dpi", "300"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise AssertionError(f"--render --dpi failed: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
    dpi_data = dpi_output.read_bytes()
    if b"pHYs" not in dpi_data:
        raise AssertionError("--dpi did not embed a pHYs chunk")
    if b"pHYs" in data:
        raise AssertionError("PNG written without --dpi should not contain a pHYs chunk")

    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(root / "bad_dpi.png"), "--dpi", "-5"],
        capture_output=True, text=True, timeout=30)
    if result.returncode == 0:
        raise AssertionError("--render accepted a non-positive --dpi value")

    # Invalid element symbols must be rejected, not silently ignored.
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(root / "bad.png"),
         "--color", "Xx 0.5 0.5 0.5"],
        capture_output=True, text=True, timeout=30)
    if result.returncode == 0:
        raise AssertionError("--render accepted an unknown element symbol in --color")

    # --display: a saved Display Settings file sets the background, cell colour,
    # element colours and sizes; --background still applies on top of it.
    look = root / "look.txt"
    look.write_text("background=0 0 0\ncell_colour=1 0 0\ncolour.Cu=0 1 0\nambient=1\nshadow_strength=0\n", encoding="utf-8")
    def draw(name, *extra):
        path = root / name
        run = subprocess.run([exe, "--render", "--input", str(source), "--output", str(path), "--width", "96", "--height", "96",
                              "--yaw", "20", "--pitch", "15", *extra], capture_output=True, text=True, timeout=60)
        if run.returncode != 0:
            raise AssertionError(f"--render {' '.join(extra)} failed:\n{run.stdout}\n{run.stderr}")
        return png_pixels(path)[2]
    pixels = [p for row in draw("look.png", "--display", str(look)) for p in row]
    if pixels[0] != (0, 0, 0):
        raise AssertionError(f"--display background not black: {pixels[0]}")
    if not any(r > 150 and g < 60 and b < 60 for r, g, b in pixels):
        raise AssertionError("--display cell colour (red) not drawn")
    if not any(g > 120 and r < 80 and b < 80 for r, g, b in pixels):
        raise AssertionError("--display element colour (green copper) not drawn")
    covered = lambda image: sum(1 for row in image for p in row if p != (0, 0, 0) and not (p[0] > 150 and p[1] < 60))
    small = root / "small.txt"
    small.write_text(look.read_text(encoding="utf-8") + "atom_scale=0.4\n", encoding="utf-8")
    if not covered(draw("small.png", "--display", str(small))) < 0.6 * covered(draw("look2.png", "--display", str(look))):
        raise AssertionError("--display atom_scale did not shrink the atoms")
    over = draw("over.png", "--display", str(look), "--background", "1 1 1")
    if over[0][0] != (255, 255, 255):
        raise AssertionError(f"--background did not apply on top of --display: {over[0][0]}")
    plain = draw("plain.png")
    if plain[0][0] != (255, 255, 255):
        raise AssertionError(f"default background not white: {plain[0][0]}")
    # Effects on a 5 x 5 x 5 block of copper (ball and stick is the CLI's style).
    block = root / "block.xyz"
    sites = [(2.55 * x, 2.55 * y, 2.55 * z) for x in range(5) for y in range(5) for z in range(5)]
    block.write_text(f"{len(sites)}\nblock\n" + "".join(f"Cu {x} {y} {z}\n" for x, y, z in sites), encoding="utf-8")
    def look_of(name, text):
        path = root / name
        path.write_text("background=0 0 0\nshadow_strength=0\n" + text, encoding="utf-8")
        image = root / (name + ".png")
        run = subprocess.run([exe, "--render", "--input", str(block), "--output", str(image), "--width", "128", "--height", "128",
                              "--yaw", "30", "--pitch", "25", "--display", str(path)], capture_output=True, text=True, timeout=60)
        if run.returncode != 0:
            raise AssertionError(f"--render with {text!r} failed:\n{run.stdout}\n{run.stderr}")
        return [p for row in png_pixels(image)[2] for p in row]
    lit = lambda pixels: [p for p in pixels if p != (0, 0, 0)]
    brightness = lambda pixels: sum(sum(p) for p in lit(pixels)) / max(1, len(lit(pixels)))
    plain = look_of("plain", "")
    fogged = look_of("fog", "depth_cueing=1\n")
    if not brightness(fogged) < 0.9 * brightness(plain):
        raise AssertionError(f"depth cueing did not fade far atoms: {brightness(fogged):.1f} vs {brightness(plain):.1f}")
    occluded = look_of("ao", "ambient_occlusion=1\n")
    if not brightness(occluded) < 0.9 * brightness(plain):
        raise AssertionError(f"ambient occlusion did not darken the block: {brightness(occluded):.1f} vs {brightness(plain):.1f}")
    outlined = look_of("outline", "outline_width=0.4\noutline_colour=1 0 1\n")
    magenta = lambda pixels: sum(1 for r, g, b in pixels if r > 150 and b > 150 and g < 90)
    if magenta(outlined) < 50 or magenta(plain) > 0:
        raise AssertionError(f"outlines not drawn in their colour ({magenta(outlined)} magenta pixels)")
    blue = lambda pixels: sum(1 for r, g, b in pixels if b > 120 and r < 60 and g < 60)
    no_metal = look_of("nometal", "bond_colour=0 0 1\nambient=1\n")
    metal = look_of("metal", "bond_colour=0 0 1\nambient=1\nmetal_bonds=true\n")
    if blue(no_metal) != 0 or blue(metal) < 20:
        raise AssertionError(f"metal_bonds did not add Cu-Cu bonds ({blue(no_metal)} -> {blue(metal)} bond pixels)")
    tight = look_of("tight", "bond_colour=0 0 1\nambient=1\nmetal_bonds=true\nbond_tolerance=0.8\n")
    limited = look_of("limited", "bond_colour=0 0 1\nambient=1\nmetal_bonds=true\nbond_atom_limit=100\n")
    if blue(tight) != 0 or blue(limited) != 0:
        raise AssertionError("bond_tolerance or bond_atom_limit did not remove the bonds")
    red = lambda pixels: sum(1 for r, g, b in pixels if r > 150 and g < 60 and b < 60)
    thin = look_of("thin", "cell_colour=1 0 0\ncell_line_width=1\n")
    thick = look_of("thick", "cell_colour=1 0 0\ncell_line_width=8\n")
    if red(thick) < red(thin) or red(thin) == 0:
        raise AssertionError(f"cell line width: {red(thin)} -> {red(thick)} cell pixels")
    print("wide lines:", "drawn" if red(thick) > red(thin) else "not supported by this driver (1 px)")

    result = subprocess.run([exe, "--render", "--input", str(source), "--output", str(root / "x.png"),
                             "--display", str(root / "missing.txt")], capture_output=True, text=True, timeout=30)
    if result.returncode == 0 or "not found" not in result.stderr:
        raise AssertionError("--render accepted a missing --display file")

print("Render CLI smoke tests passed (help, single PNG, turntable frames, DPI metadata, display settings, invalid input)")
