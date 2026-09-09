# AtomForge manual

Read [AtomForge-manual.pdf](AtomForge-manual.pdf) for detailed desktop workflows, all 33 electronic operation choices, charge studies, Python/CLI usage, examples, and troubleshooting. Edit [manual.tex](manual.tex) and the modular files in `chapters/` to update it.

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
