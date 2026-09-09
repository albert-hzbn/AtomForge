# AtomForge manual

Read [AtomForge-manual.pdf](AtomForge-manual.pdf) for detailed desktop workflows, all 33 electronic operation choices, charge studies, Python/CLI usage, examples, and troubleshooting. Edit [manual.tex](manual.tex) and the modular files in `chapters/` to update it.

The expanded reference covers every Build and Edit tool, their individual controls and application behavior, all View/Settings families, file handling, structural analysis, and optional/source-only features. The generated appendix documents 71 public Python API entries and all eight native builder CLI modes. The coverage checklist distinguishes enabled workflows from source-only implementations.

From the repository root:

```sh
python docs/manual/scripts/build_manual.py
```

Requires pdfLaTeX with the packages listed in the manual. The checked-in figures are sufficient to compile; the original CHGCAR is not needed or included.

To regenerate builder examples and scientific plots, build AtomForge, install its local Python package plus NumPy and Matplotlib, then run:

```sh
python docs/manual/scripts/generate_examples.py --chgcar /path/to/CHGCAR
```

This overwrites generated files in `examples/` and the scientific figures in `figures/`. GUI screenshots must be recaptured from the real application. `build/` contains ignored compilation/render intermediates. The input fingerprint and checked numerical values are in `examples/chgcar-summary.json`; builder arguments and counts are in `examples/builder-results.json`.

Verify the published structure counts and executable Python recipes with:

```sh
python docs/manual/scripts/check_examples.py
```

Check source drift and the generated API/CLI reference:

```sh
python docs/manual/scripts/check_coverage.py --verify-executable build/AtomForge.exe
```

`coverage.json` maps 106 source files to manual chapters and records 580 literal UI control declarations, including internal IDs. It is a source-drift check, not an automated proof of semantic completeness. After a code change, review the affected features and update the prose before running `check_coverage.py --refresh`. Normal PDF builds run the coverage check automatically. On other platforms, pass the corresponding executable path or omit `--verify-executable` for a source-only check.

The dislocation example is a geometric starting configuration; its generation reports a lattice-family change warning. It preserves 500 host atoms and is included to demonstrate the operation, not to claim a relaxed defect structure.
