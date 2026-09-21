"""Scientific command-line workflows: python -m atomforge.science --help."""

import argparse
import json
from pathlib import Path
import os
import tempfile
import sys

from .workflows import TOOLS, json_result, run_tool, tool_catalog


def main(argv=None):
    parser = argparse.ArgumentParser(description="AtomForge condensed-matter analysis and simulation")
    parser.add_argument("tool", choices=list(TOOLS), nargs="?")
    parser.add_argument("--catalog", action="store_true", help="Print tools, assumptions and parameter defaults as JSON")
    parser.add_argument("--input", type=Path, help="JSON object containing tool parameters; data file paths are relative to it")
    parser.add_argument("--output", type=Path, help="Result JSON; existing files require --overwrite")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--report", type=Path, help="Optional human-readable result summary")
    parser.add_argument("--structures", type=Path, help="Optional extXYZ output for simulation images/frames or a primitive structure")
    args = parser.parse_args(argv)
    try:
        if args.catalog:
            print(json.dumps(tool_catalog(), indent=2, allow_nan=False))
            return 0
        if not args.tool or not args.input or not args.output:
            parser.error("tool, --input and --output are required unless --catalog is used")
        if args.output.exists() and not args.overwrite:
            raise FileExistsError("Output exists; choose another file or use --overwrite")
        destinations = [path.resolve() for path in (args.output, args.report, args.structures) if path is not None]
        if len(set(destinations)) != len(destinations) or args.input.resolve() in destinations:
            raise ValueError("Input, result, report and structure files must have distinct paths")
        if args.report and args.report.exists() and not args.overwrite:
            raise FileExistsError("Report exists; choose another file or use --overwrite")
        parameters = json.loads(args.input.read_text(encoding="utf-8"))
        result = run_tool(args.tool, parameters, base_directory=args.input.resolve().parent)
        if args.structures:
            from .workflows import export_structures
            export_structures(result, args.structures, overwrite=args.overwrite)
        payload = {"tool": args.tool, "parameters": parameters, "result": json_result(result)}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # Do not truncate an existing result if calculation or serialization fails.
        descriptor, temporary = tempfile.mkstemp(prefix=".atomforge-science-", dir=args.output.parent, suffix=".json")
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                json.dump(payload, stream, indent=2, allow_nan=False)
                stream.write("\n")
            os.replace(temporary, args.output)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        print("Saved " + str(args.output))
        if args.report:
            from .workflows import result_report
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(result_report(args.tool, result), encoding="utf-8")
        return 0
    except (ValueError, TypeError, KeyError, OSError, ImportError, RuntimeError) as error:
        print("AtomForge science: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
