"""Command-line volume conversion and electronic post-processing."""

import argparse
import csv
import sys

from . import load_volume


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--quantity", choices=("auto", "raw", "density", "potential", "elf"), default="auto")
    parser.add_argument("--cube-coordinates", choices=("bohr", "angstrom"), default="bohr")
    parser.add_argument("--channel", type=int, default=0)
    parser.add_argument("--periodic", action="store_true", help="Verify and drop duplicate endpoint planes")
    parser.add_argument("--operation", choices=("info", "convert", "integrate", "laplacian", "patterson", "kinetic", "potential-energy", "total-energy", "planar", "peaks", "isosurface"), default="info")
    parser.add_argument("--axis", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--level", type=float, default=0.1)
    parser.add_argument("--output")
    parser.add_argument("--format", choices=("cube", "xsf", "vasp"))
    args = parser.parse_args(argv)
    try:
        fields = load_volume(args.input, args.quantity, args.cube_coordinates).fields
        if not 0 <= args.channel < len(fields):
            raise ValueError("Channel index out of range")
        grid = fields[args.channel]
        if args.periodic:
            grid = grid.as_periodic()
        op = args.operation
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
    except (OSError, ValueError, RuntimeError) as error:
        print("Electronic post-processing: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
