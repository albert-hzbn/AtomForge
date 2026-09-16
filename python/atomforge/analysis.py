"""Structural analyses using the same numerical engine as the desktop.

Coordinates and distances are Angstrom. Results are lists of column dictionaries;
native errors propagate as CalledProcessError with diagnostic stderr.
"""

import csv
import os
from pathlib import Path
import subprocess
import tempfile

from ._structure import Structure
from .builders import _executable


def analyze(source, method, *, executable=None, timeout=300, **parameters):
    """Run cna, rdf, adf, sro or interstitial and return named CSV rows.

    Parameter names match CLI flags with underscores replacing hyphens.
    ``no_pbc=True`` and ``raw=True`` are switches. SRO currently uses the
    native shell convention and does not accept a no_pbc override.
    """
    if method not in ("cna", "rdf", "adf", "sro", "interstitial"):
        raise ValueError("Unknown analysis method")
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    allowed = {
        "cna": {"cutoff", "scale", "no_pbc"},
        "rdf": {"cutoff", "minimum", "bins", "centre_z", "neighbor_z", "smooth", "no_pbc", "raw"},
        "adf": {"cutoff", "bins", "smooth", "centre", "neighbor_a", "neighbor_b", "no_pbc", "raw"},
        "sro": {"shells", "tolerance"},
        "interstitial": {"resolution", "limit", "clearance", "separation"},
    }
    if set(parameters) - allowed[method]:
        raise ValueError("Unsupported parameters: " + ", ".join(sorted(set(parameters) - allowed[method])))
    with tempfile.TemporaryDirectory(prefix="atomforge_analysis_") as directory:
        root = Path(directory)
        if isinstance(source, Structure):
            source_path = root / ("source.vasp" if source.cell else "source.xyz")
            source.save(source_path)
        else:
            source_path = Path(source).expanduser().resolve()
            if not source_path.is_file():
                raise FileNotFoundError(str(source_path))
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", method, "--input", str(source_path),
                   "--output", str(result)]
        for key, value in parameters.items():
            if key in ("no_pbc", "raw"):
                if not isinstance(value, bool):
                    raise TypeError(key + " must be boolean")
                if value:
                    command.append("--" + key.replace("_", "-"))
            else:
                command.extend(("--" + key.replace("_", "-"), str(value)))
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            for key, value in row.items():
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
        return rows


def cna(source, **parameters):
    """Common-neighbour indices and full-environment classification per atom."""
    return analyze(source, "cna", **parameters)


def rdf(source, **parameters):
    """Radial distribution, directed raw pair counts and coordination."""
    return analyze(source, "rdf", **parameters)


def adf(source, **parameters):
    """Angular distribution of unordered neighbour pairs, in degrees."""
    return analyze(source, "adf", **parameters)


def sro(source, **parameters):
    """Warren-Cowley and Rao-Curtin short-range order by neighbour shell."""
    return analyze(source, "sro", **parameters)


def interstitial_sites(source, **parameters):
    """Interstitial candidates, clearance, volume and coordination."""
    return analyze(source, "interstitial", **parameters)


def sculpt(source, slabs, *, repeat=(1,1,1), periodic=False, executable=None, timeout=300):
    """Apply intersecting Miller slabs with the desktop sculptor engine.

    Each slab is (h,k,l,lower,upper). Manual bounds are A relative to the
    supercell centroid; periodic bounds mean integer start plane and period
    count, using the native structural-period convention.
    """
    from ._io import load
    if len(repeat)!=3 or any(not isinstance(n,int) or n<1 for n in repeat):
        raise ValueError("Repeat must contain three positive integers")
    rows=[]
    for slab in slabs:
        if len(slab)!=5:
            raise ValueError("Each slab needs h,k,l,lower,upper")
        rows.append(" ".join(str(v) for v in (*slab,int(periodic))))
    if not rows:
        raise ValueError("Provide at least one slab")
    with tempfile.TemporaryDirectory(prefix="atomforge_sculpt_") as directory:
        root=Path(directory)
        if isinstance(source,Structure):
            path=root/"source.vasp"; source.save(path)
        else:
            path=Path(source).resolve()
        output=root/"result.vasp"
        command=[_executable(executable),"--analyze","sculpt","--input",str(path),
                 "--output",str(output),"--slabs",";".join(rows)]
        for flag,value in zip(("--nx","--ny","--nz"),repeat):
            command.extend((flag,str(value)))
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name=="nt" else 0)
        return load(output)
