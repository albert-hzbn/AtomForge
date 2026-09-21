"""Headless PNG rendering of a Structure via the native --render CLI mode."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Sequence, Tuple

from .builders import _executable
from ._io import save

if TYPE_CHECKING:
    from ._structure import Structure


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
    background: Tuple[float, float, float] = (1.0, 1.0, 1.0),
    show_bonds: bool = True,
    show_box: bool = True,
    colors: Optional[Dict[str, Tuple[float, float, float]]] = None,
    radii: Optional[Dict[str, float]] = None,
    radius_scale: float = 1.0,
    dpi: Optional[float] = None,
    frames: Optional[int] = None,
    yaw_step: Optional[float] = None,
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

    frames + yaw_step render a turntable sequence instead of a single image;
    output is then used as a base name, and each frame is written to
    ``<output>-000.png``, ``<output>-001.png``, etc. (returned as the first
    frame's path with the numeric suffix, or a formatted description).
    """
    exe = _executable(executable)
    output = os.fspath(output)

    resolved_colors: Dict[str, Tuple[float, float, float]] = {}
    seen_symbols = set()
    for atom in structure.atoms:
        if atom.symbol not in seen_symbols:
            seen_symbols.add(atom.symbol)
            resolved_colors[atom.symbol] = (atom.r, atom.g, atom.b)
    if colors:
        resolved_colors.update(colors)

    with tempfile.TemporaryDirectory(prefix="atomforge_render_") as directory:
        input_path = Path(directory) / "structure.cif"
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
            "--background", "{} {} {}".format(*background),
            "--radius-scale", str(radius_scale),
        ]
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
