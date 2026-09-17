"""Command-line volume conversion and electronic post-processing."""

import argparse
import csv
import json
from pathlib import Path
import sys

from . import load_volume
from .pipeline import _METHODS, _BINARY, apply_operation, run_pipeline, save_result, batch_process


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--quantity", choices=("auto", "raw", "density", "potential", "elf"), default="auto")
    parser.add_argument("--cube-coordinates", choices=("bohr", "angstrom"), default="bohr")
    parser.add_argument("--channel", type=int, default=0)
    parser.add_argument("--periodic", action="store_true", help="Verify and drop duplicate endpoint planes")
    parser.add_argument("--operation", choices=("info", "convert", "integrate", "laplacian", "patterson", "kinetic", "potential-energy", "total-energy", "planar", "peaks", "isosurface", "bader_partition", *sorted(_METHODS | _BINARY | {"scale"})), default="info")
    parser.add_argument("--axis", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--level", type=float, default=0.1)
    parser.add_argument("--output")
    parser.add_argument("--format", choices=("cube", "xsf", "vasp"))
    parser.add_argument("--reference", help="Reference field for binary operations")
    parser.add_argument("--parameters", default="{}", help="JSON object of operation parameters")
    parser.add_argument("--recipe", help="JSON file containing an array of processing steps")
    parser.add_argument("--batch", nargs="*", help="Additional inputs; --output is the batch directory")
    args = parser.parse_args(argv)
    try:
        parameters = json.loads(args.parameters)
        if not isinstance(parameters, dict):
            raise ValueError("Parameters must be a JSON object")
        steps = json.loads(Path(args.recipe).read_text(encoding="utf-8")) if args.recipe else None
        if steps is not None and not isinstance(steps, list):
            raise ValueError("Recipe must be an array of steps")
        if args.batch is not None:
            if steps is None or not args.output:
                raise ValueError("Batch processing requires --recipe and --output")
            report = batch_process([args.input, *args.batch], steps, args.output,
                quantity=args.quantity, channel=args.channel, format=args.format or "xsf",
                base_directory=Path(args.recipe).resolve().parent)
            print(json.dumps(report, indent=2))
            return 0 if all(item["success"] for item in report) else 1
        fields = load_volume(args.input, args.quantity, args.cube_coordinates).fields
        if not 0 <= args.channel < len(fields):
            raise ValueError("Channel index out of range")
        grid = fields[args.channel]
        if args.periodic:
            grid = grid.as_periodic()
        op = args.operation
        if steps is not None:
            if not args.output:
                raise ValueError("Recipe requires --output")
            results = run_pipeline(grid, steps, base_directory=Path(args.recipe).resolve().parent)
            save_result(next(reversed(results.values())), args.output, args.format)
            return 0
        legacy = {"info", "convert", "integrate", "laplacian", "patterson", "kinetic", "potential-energy", "total-energy", "planar", "peaks", "isosurface", "bader_partition"}
        if op not in legacy or parameters or args.reference:
            aliases = {"planar": "planar_average"}
            reference = load_volume(args.reference, args.quantity, args.cube_coordinates).fields[0] if args.reference else None
            result = apply_operation(grid, aliases.get(op, op), parameters, reference)
            if args.output:
                save_result(result, args.output, args.format)
            elif isinstance(result, (dict, int, float)):
                print(json.dumps(result, allow_nan=False))
            else:
                raise ValueError("This operation requires --output")
            return 0
        if op == "info":
            print("channels={}, shape={}, unit={}, periodic={}".format(len(fields), grid.shape, grid.unit, grid.periodic))
        elif op == "integrate":
            print("{:.17g}".format(grid.integrate()))
        else:
            if not args.output:
                raise ValueError("This operation requires --output")
            if op in ("planar", "peaks"):
                rows = grid.planar_average(args.axis) if op == "planar" else grid.peaks()
                header = ("distance_A", "value") if op == "planar" else ("x_A", "y_A", "z_A", "value")
                with open(args.output, "w", newline="", encoding="utf-8") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(header)
                    writer.writerows(rows)
            elif op == "bader_partition":
                partition = grid.bader_partition()
                with open(args.output, "w", newline="", encoding="utf-8") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(("basin", "charge_e", "volume_A3", "max_x_A", "max_y_A", "max_z_A"))
                    for basin_index in range(partition.num_basins):
                        basin = partition.basin(basin_index)
                        writer.writerow((basin_index, basin["charge"], basin["volume"], *basin["maximum"]))
            elif op == "isosurface":
                grid.isosurface(args.level).save(args.output)
            else:
                if op == "laplacian":
                    grid = grid.laplacian()
                elif op == "patterson":
                    grid = grid.patterson()
                elif op in ("kinetic", "potential-energy", "total-energy"):
                    grid = grid.energy_density()[("kinetic", "potential-energy", "total-energy").index(op)]
                grid.save(args.output, args.format)
    except (OSError, ValueError, RuntimeError, TypeError, KeyError, IndexError) as error:
        print("Electronic post-processing: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
