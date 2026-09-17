"""Declarative, repeatable electronic processing without executing recipe code."""

import csv
import json
from pathlib import Path

from ._grid import Grid, Surface, load_volume


_METHODS = {
    "smooth", "gradient", "laplacian", "energy_density", "line_profile",
    "planar_average", "macroscopic_average", "section", "contours", "integrate",
    "integrate_sphere", "voronoi_integrate", "peaks", "structure_factors",
    "fourier_synthesis", "patterson", "ewald", "isosurface", "threshold_mask",
    "invert_mask", "split_density", "charge_summary", "cumulative_charge", "as_periodic",
    "reduced_density_gradient", "signed_density", "dori", "betti_curve",
}
_BINARY = {"add", "subtract", "multiply", "divide", "resample", "density_difference",
           "boolean", "apply_mask"}


def apply_operation(grid, operation, parameters=None, reference=None):
    """Evaluate one allowlisted operation; return a grid, fields, surface or table."""
    parameters = dict(parameters or {})
    if operation in _BINARY:
        if reference is None:
            raise ValueError(operation + " requires a reference grid")
        if operation == "density_difference":
            return grid.density_difference(reference, **parameters)
        if operation == "boolean":
            return grid.boolean(reference, **parameters)
        if parameters:
            raise ValueError("Unexpected parameters for " + operation)
        return {
            "add": lambda: grid + reference, "subtract": lambda: grid - reference,
            "multiply": lambda: grid * reference, "divide": lambda: grid / reference,
            "resample": lambda: grid.resample(reference),
            "apply_mask": lambda: grid.apply_mask(reference),
        }[operation]()
    if reference is not None:
        raise ValueError("Unexpected reference for " + operation)
    if operation == "scale":
        if set(parameters) != {"factor"}:
            raise ValueError("scale requires only factor")
        return grid * parameters["factor"]
    if operation not in _METHODS:
        raise ValueError("Unknown electronic operation: " + operation)
    return getattr(grid, operation)(**parameters)


def run_pipeline(grid, steps, *, base_directory="."):
    """Run JSON-compatible steps and retain all named intermediate results.

    Each step contains ``operation``, optional ``parameters``, ``name``,
    ``input`` (an earlier grid name), and ``reference`` (earlier name or file).
    ``field`` selects a component when an operation returns multiple fields.
    File references resolve relative to base_directory. Recipes cannot call
    arbitrary Python functions or execute shell commands.
    """
    results = {"input": grid}
    current = grid
    for index, step in enumerate(steps):
        if not isinstance(step, dict) or set(step) - {"operation", "parameters", "name", "input", "reference", "field"}:
            raise ValueError("Invalid pipeline step")
        name = step.get("name", "step_{}".format(index + 1))
        if not isinstance(name, str) or not name or name in results:
            raise ValueError("Step names must be unique nonempty strings")
        source = results[step["input"]] if "input" in step else current
        if not isinstance(source, Grid):
            raise ValueError("Step input must be a grid")
        reference = None
        if "reference" in step:
            key = step["reference"]
            reference = results.get(key)
            if reference is None:
                reference = load_volume(Path(base_directory) / key).fields[0]
            if not isinstance(reference, Grid):
                raise ValueError("Reference must be a grid")
        result = apply_operation(source, step["operation"], step.get("parameters"), reference)
        if "field" in step:
            field = step["field"]
            if not isinstance(field, int) or field < 0 or not isinstance(result, (list, tuple)):
                raise ValueError("field must select a returned grid component")
            result = result[field]
            if not isinstance(result, Grid):
                raise ValueError("Selected component is not a grid")
        results[name] = result
        if isinstance(result, Grid):
            current = result
    return results


def save_result(result, path, format=None):
    """Write a grid/mesh, numeric CSV, JSON summary, or numbered grid components."""
    path = Path(path)
    if isinstance(result, Grid):
        result.save(path, format)
    elif isinstance(result, Surface):
        result.save(path)
    elif isinstance(result, (list, tuple)) and result and isinstance(result[0], Grid):
        for index, field in enumerate(result):
            field.save(path.with_name(path.stem + "-{}".format(index) + path.suffix), format)
    elif isinstance(result, dict):
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    else:
        rows = [[result]] if isinstance(result, (int, float)) else result
        with path.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(rows)


def batch_process(paths, steps, output_directory, *, quantity="auto", channel=0,
                  format="xsf", base_directory=".", overwrite=False):
    """Apply one recipe to every input and write the last result plus provenance.

    Refuses name collisions and existing outputs by default. Failed inputs are
    reported individually; successful results remain available.
    """
    if format not in ("xsf", "cube", "vasp"):
        raise ValueError("Unknown grid format")
    paths = [Path(path) for path in paths]
    names = [path.stem for path in paths]
    if len(names) != len(set(names)):
        raise ValueError("Input names collide; process them in separate batches")
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    report = []
    for path in paths:
        folder = destination / path.stem
        try:
            folder.mkdir(exist_ok=overwrite)
            fields = load_volume(path, quantity).fields
            if not isinstance(channel, int) or not 0 <= channel < len(fields):
                raise ValueError("Channel index out of range")
            results = run_pipeline(fields[channel], steps, base_directory=base_directory)
            last = next(reversed(results.values()))
            suffix = format if isinstance(last, Grid) or (isinstance(last, (tuple, list)) and last and isinstance(last[0], Grid)) else "obj" if isinstance(last, Surface) else "json" if isinstance(last, dict) else "csv"
            save_result(last, folder / ("result." + suffix), format)
            (folder / "recipe.json").write_text(json.dumps({"input": str(path.resolve()),
                "quantity": quantity, "channel": channel, "steps": steps}, indent=2) + "\n", encoding="utf-8")
            report.append({"input": str(path), "output": str(folder), "success": True})
        except (OSError, ValueError, RuntimeError, TypeError, KeyError, IndexError) as error:
            report.append({"input": str(path), "success": False, "error": str(error)})
    return report
