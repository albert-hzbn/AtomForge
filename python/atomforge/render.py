"""Headless PNG rendering of a Structure via the native --render CLI mode."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Mapping, Optional, Sequence, Tuple, Union

from .builders import _executable
from ._io import save

if TYPE_CHECKING:
    from ._structure import Structure


def _display_text(settings: Mapping[str, object]) -> str:
    """Display settings as the key=value text of a saved default."""
    lines = []
    for key, value in settings.items():
        if isinstance(value, bool):
            text = "true" if value else "false"
        elif isinstance(value, (tuple, list)):
            if len(value) != 3:
                raise ValueError(f"display[{key!r}] needs three colour values (r, g, b)")
            text = " ".join(repr(float(v)) for v in value)
        else:
            text = str(value)
        if "\n" in key or "=" in key or "\n" in text:
            raise ValueError(f"display[{key!r}] is not a single key=value line")
        lines.append(f"{key}={text}")
    return "\n".join(lines) + "\n"


def _display_file(display, directory) -> str:
    """The --display argument: "saved", an existing file, or a written mapping."""
    if isinstance(display, Mapping):
        path = Path(directory) / "display.txt"
        path.write_text(_display_text(display), encoding="utf-8")
        return str(path)
    display = os.fspath(display)
    if display != "saved" and not Path(display).is_file():
        raise FileNotFoundError(f"display settings file not found: {display}")
    return display


def render(
    structure: "Structure",
    output,
    *,
    width: int = 1600,
    height: int = 1200,
    yaw: float = 0.0,
    pitch: float = 0.0,
    roll: float = 0.0,
    distance: Optional[float] = None,
    orthographic: bool = False,
    background: Optional[Tuple[float, float, float]] = None,
    show_bonds: bool = True,
    show_box: bool = True,
    colors: Optional[Dict[str, Tuple[float, float, float]]] = None,
    radii: Optional[Dict[str, float]] = None,
    radius_scale: float = 1.0,
    dpi: Optional[float] = None,
    frames: Optional[int] = None,
    yaw_step: Optional[float] = None,
    display: Union[None, str, os.PathLike, Mapping[str, object]] = None,
    executable: Optional[str] = None,
    timeout: float = 120,
) -> str:
    """Render *structure* to a PNG at *output*, without launching the GUI.

    yaw/pitch/roll orbit the camera in degrees; distance defaults to an
    auto-fit view. colors/radii override one element's appearance at a time,
    e.g. ``colors={"Fe": (0.8, 0.4, 0.1)}``, ``radii={"Fe": 1.4}`` (Angstrom).
    Any color already set on the structure's atoms (via
    ``Structure.set_element_color`` or direct ``atom.r/g/b`` edits) is used
    for that element unless overridden by *colors*, since structure files
    written to disk for the native renderer don't carry per-atom color.

    dpi optionally embeds a physical resolution (dots per inch) in the saved
    PNG's metadata (a pHYs chunk); it does not change the pixel dimensions,
    only how image viewers/editors interpret the print size.

    display applies the desktop's Display Settings: a file written by Save as
    My Default, ``"saved"`` for that default, or a mapping of its keys, e.g.
    ``{"ambient": 0.5, "light_azimuth": 200, "atom_scale": 0.6,
    "bond_colour": (0.6, 0.6, 0.6), "background": (0, 0, 0),
    "colour.Cu": (0.9, 0.5, 0.2), "radius.O": 0.8}`` (colours as r, g, b in
    0-1; see the manual for every key). background, colors, radii and
    radius_scale apply on top of it. With display, element colours come from
    it rather than from the atoms, unless given in colors. background defaults
    to white (or the display settings' background).

    frames + yaw_step render a turntable sequence instead of a single image;
    output is then used as a base name, and each frame is written to
    ``<output>-000.png``, ``<output>-001.png``, etc. (returned as the first
    frame's path with the numeric suffix, or a formatted description).
    """
    exe = _executable(executable)
    output = os.fspath(output)

    resolved_colors: Dict[str, Tuple[float, float, float]] = {}
    seen_symbols = set()
    for atom in structure.atoms if display is None else ():
        if atom.symbol not in seen_symbols:
            seen_symbols.add(atom.symbol)
            resolved_colors[atom.symbol] = (atom.r, atom.g, atom.b)
    if colors:
        resolved_colors.update(colors)

    with tempfile.TemporaryDirectory(prefix="atomforge_render_") as directory:
        input_path = Path(directory) / ("structure.cif" if structure.cell else "structure.xyz")
        save(structure, str(input_path))

        command = [
            exe, "--render",
            "--input", str(input_path),
            "--output", output,
            "--width", str(width),
            "--height", str(height),
            "--yaw", str(yaw),
            "--pitch", str(pitch),
            "--roll", str(roll),
            "--radius-scale", str(radius_scale),
        ]
        if background is not None:
            command += ["--background", "{} {} {}".format(*background)]
        if display is not None:
            command += ["--display", _display_file(display, directory)]
        if dpi is not None:
            command += ["--dpi", str(dpi)]
        if distance is not None:
            command += ["--distance", str(distance)]
        if orthographic:
            command.append("--orthographic")
        if not show_bonds:
            command.append("--no-bonds")
        if not show_box:
            command.append("--no-box")
        for symbol, (r, g, b) in resolved_colors.items():
            command += ["--color", "{} {} {} {}".format(symbol, r, g, b)]
        if radii:
            for symbol, value in radii.items():
                command += ["--radius", "{} {}".format(symbol, value)]
        if frames is not None:
            command += ["--frames", str(frames)]
            if yaw_step is not None:
                command += ["--yaw-step", str(yaw_step)]

        subprocess.run(
            command, check=True, capture_output=True, text=True, timeout=timeout,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        )

    return output
