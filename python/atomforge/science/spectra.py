"""Electronic bands, DOS/PDOS and orbital projections from VASP output."""

import math
from pathlib import Path


def read_bands(path):
    """Read EIGENVAL with validated dimensions; return kpoints, weights, spin energies/occupancies.

    Energies are absolute VASP eigenvalues in eV. Fermi energy is not present in
    EIGENVAL and must be supplied separately when plotting.
    """
    import numpy as np
    import gzip
    import bz2
    import lzma
    opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}.get(Path(path).suffix.lower(), open)
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    if len(lines) < 6:
        raise ValueError("Truncated EIGENVAL header")
    spin_count = int(lines[0].split()[-1])
    _, nkpoints, nbands = map(int, lines[5].split())
    if spin_count not in (1,2) or nkpoints < 1 or nbands < 1 or nkpoints*nbands > 100000000:
        raise ValueError("Invalid EIGENVAL dimensions")
    rows = iter(line.split() for line in lines[6:] if line.strip())
    labels = ("up", "down") if spin_count == 2 else ("total",)
    energies = {label: np.empty((nkpoints,nbands)) for label in labels}
    occupancies = {label: np.empty((nkpoints,nbands)) for label in labels}
    kpoints, weights = [], []
    try:
        for k in range(nkpoints):
            point = list(map(float,next(rows)))
            if len(point) != 4 or not all(math.isfinite(v) for v in point):
                raise ValueError("Invalid EIGENVAL k point")
            kpoints.append(point[:3]); weights.append(point[3])
            for band in range(nbands):
                values = list(map(float,next(rows)))
                if len(values) != 1+2*spin_count or values[0] != band+1 or not all(math.isfinite(v) for v in values):
                    raise ValueError("Invalid EIGENVAL band row")
                for spin, label in enumerate(labels):
                    energies[label][k,band] = values[1+spin]
                    occupancies[label][k,band] = values[1+spin_count+spin]
    except StopIteration as error:
        raise ValueError("Truncated EIGENVAL") from error
    if next(rows,None) is not None:
        raise ValueError("Unexpected EIGENVAL trailing data")
    return {"kpoints":kpoints, "weights":weights, "energies_eV":energies, "occupancies":occupancies}



def read_projections(path):
    """Read PROCAR, including spin and SOC projections via pymatgen.

    Returns the native Procar object: data[spin][kpoint,band,ion,orbital],
    eigenvalues, occupancies, orbital names and SOC xyz_data when present.
    """
    from pymatgen.io.vasp.outputs import Procar
    return Procar(str(path))


def read_dos(path, *, projected_labels=None):
    """Read DOSCAR total DOS and all projected columns without guessing LORBIT.

    Total DOS columns are labelled by spin; projected columns require caller
    labels or remain column_1, column_2, ... because DOSCAR does not identify
    the orbital/SOC convention. Integrated total DOS is retained separately.
    """
    import numpy as np
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < 7:
        raise ValueError("Truncated DOSCAR header")
    header = lines[5].split()
    points = int(header[2]); fermi = float(header[3])
    if points < 1 or points > 10000000 or not math.isfinite(fermi):
        raise ValueError("Invalid DOSCAR grid")
    def block(start):
        rows = [list(map(float,line.split())) for line in lines[start:start+points]]
        if len(rows) != points or not rows or len({len(row) for row in rows}) != 1:
            raise ValueError("Truncated or inconsistent DOSCAR block")
        values = np.asarray(rows)
        if values.shape[1] < 2 or not np.isfinite(values).all() or np.any(np.diff(values[:,0]) <= 0):
            raise ValueError("Invalid DOSCAR energy/data values")
        return values
    total = block(6)
    if total.shape[1] not in (3,5):
        raise ValueError("Total DOS must have 3 or 5 columns")
    polarized = total.shape[1] == 5
    result = {"energy_eV":total[:,0], "fermi_eV":fermi,
              "dos": {"up":total[:,1], "down":total[:,2]} if polarized else {"total":total[:,1]},
              "integrated":total[:,3:5] if polarized else total[:,2:3], "projected":[]}
    offset = 6 + points
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1; continue
        site_header = lines[offset].split()
        if len(site_header) < 4 or int(site_header[2]) != points:
            raise ValueError("Invalid projected DOS header")
        values = block(offset+1)
        if not np.allclose(values[:,0],total[:,0],rtol=0,atol=1e-7):
            raise ValueError("Projected DOS energy grids differ")
        result["projected"].append(values[:,1:])
        offset += points + 1
    if result["projected"]:
        columns = result["projected"][0].shape[1]
        if any(values.shape[1] != columns for values in result["projected"]):
            raise ValueError("Projected DOS column counts differ")
        labels = list(projected_labels) if projected_labels is not None else ["column_{}".format(i+1) for i in range(columns)]
        if len(labels) != columns:
            raise ValueError("Supply one label per projected DOS column")
        result["projected_labels"] = labels
    return result


def _figure():
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    figure = Figure(figsize=(7,5), layout="constrained")
    FigureCanvasAgg(figure)
    return figure, figure.add_subplot(111)


def plot_bands(data, output, *, fermi_eV=0, distances=None, projections=None):
    """Plot bands or fat bands; projections map spin to nonnegative (k,band) weights.

    Without reciprocal-space distances the horizontal axis is k-point index,
    avoiding an incorrect metric for non-orthogonal reciprocal cells.
    """
    import numpy as np
    figure, axes = _figure()
    for spin, energies in data["energies_eV"].items():
        energies = np.asarray(energies)
        x = np.arange(energies.shape[0]) if distances is None else np.asarray(distances)
        if x.shape != (energies.shape[0],):
            raise ValueError("One distance per k point is required")
        for band in range(energies.shape[1]):
            axes.plot(x, energies[:,band]-fermi_eV, linewidth=.8, label=spin if band==0 else None)
            if projections is not None:
                weights = np.asarray(projections[spin])
                if weights.shape != energies.shape or not np.isfinite(weights).all() or (weights<0).any():
                    raise ValueError("Projection weights must be finite, nonnegative and match the bands")
                axes.scatter(x,energies[:,band]-fermi_eV,s=weights[:,band]*50,alpha=.5)
    axes.axhline(0,color="gray",linewidth=.6)
    axes.set_xlabel("k-point index" if distances is None else "Reciprocal path distance (1/A)")
    axes.set_ylabel("Energy relative to reference (eV)"); axes.legend()
    figure.savefig(str(output),dpi=180)
    return figure


def plot_dos(data, output, *, sites=(), columns=None):
    """Save total and optionally site-projected DOS relative to the file's Fermi energy."""
    figure, axes = _figure()
    x = data["energy_eV"] - data["fermi_eV"]
    for name, values in data["dos"].items():
        axes.plot(x, values, label=name)
    for site in sites:
        values = data["projected"][site]
        selected = range(values.shape[1]) if columns is None else columns
        for column in selected:
            axes.plot(x,values[:,column],label="site {} {}".format(site,data["projected_labels"][column]))
    axes.set_xlabel("Energy - Fermi energy (eV)"); axes.set_ylabel("DOS (states/eV)"); axes.legend()
    figure.savefig(str(output),dpi=180)
    return figure
