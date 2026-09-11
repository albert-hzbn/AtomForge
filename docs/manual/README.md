# AtomForge manual

Read [AtomForge-manual.pdf](AtomForge-manual.pdf) for detailed desktop workflows, all 33 electronic operation choices, charge studies, Python/CLI usage, examples, and troubleshooting. Edit [manual.tex](manual.tex) and the modular files in `chapters/` to update it.

The expanded reference covers every Build and Edit tool, their individual controls and application behavior, all View/Settings families, file handling, structural analysis, and optional/source-only features. The generated appendix documents 74 public Python API entries and all eight native builder CLI modes. The coverage checklist distinguishes enabled workflows from source-only implementations.

Chapter 9 adds 62 illustrated quick tutorials with starting files, numbered steps, result checks, and options to try. The screenshot collection includes populated reference/mesh previews and completed structures, alongside seven nanocrystal plots and actual AtomForge cover renders at default atomic radii. Input previews are distinguished from finished results. Each of the seven shape cuts has a reusable XYZ example in `examples/shapes/`. Screenshots show the states described in their captions; they are not all completed exercises. The SRO tutorial records an observed shell-distance discrepancy and explains how to detect it using the FCC reference and RDF.

Edit `tutorials.json`, then regenerate its LaTeX chapter and validate its image references and crop bounds:

```sh
python docs/manual/scripts/render_tutorials.py --write
```

Normal manual builds reject a stale generated tutorial chapter. Original GUI captures remain at full resolution; LaTeX applies only the documented crops. Extra captures preserve alternative option panels for reuse. The electronic tutorials also reuse the earlier CHGCAR screenshots.

The six cover images in `figures/cover/` are raw AtomForge viewport captures. Run `scripts/generate_cover_examples.py` after generating the ordinary examples to recreate their inputs. Bulk and Voronoi use `cover-bulk.cif`; the grain boundary uses `cover-gb.cif`; amorphous packing uses `sio2_pack.xyz`; Wulff uses `cover-wulff.xyz`. The polyhedral view uses the illustrative `cover-tio6.xyz` cluster, center atom ID 1, face opacity 0.34, and default edge opacity. All captures retain default atomic radii. Only camera zoom and representation settings vary. The cover crops the surrounding interface in LaTeX without editing screenshot pixels.

From the repository root:

```sh
python docs/manual/scripts/build_manual.py
```

Requires pdfLaTeX with the packages listed in the manual. The checked-in figures are sufficient to compile; the original CHGCAR is not needed or included.

`docs/manual/AtomForge-manual.pdf` is the single published manual in a source checkout. The build script moves its validated PDF into that location, leaving no intermediate PDF behind. Use this script instead of compiling `manual.tex` directly. Help > Manual and About > Manual resolve this same file when running a development build; no PDF copy is needed beside `build/AtomForge.exe`. Installed packages each include their own single manual at the location below.

Application installation and CPack archives include the checked-in PDF automatically. Portable archives place it at `AtomForge/AtomForge-manual.pdf`; standard installations place it at `share/doc/atomforge/AtomForge-manual.pdf`. Building or packaging the application does not require LaTeX. Rebuild and check in the PDF after editing its sources so packages contain the updated manual.

To regenerate builder examples and scientific plots, build AtomForge, install its local Python package plus NumPy and Matplotlib, then run:

```sh
python docs/manual/scripts/generate_examples.py --chgcar /path/to/CHGCAR
python docs/manual/scripts/generate_shape_gallery.py
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

`coverage.json` maps 109 source files to manual chapters and records 590 literal UI control declarations, including internal IDs. It is a source-drift check, not an automated proof of semantic completeness. After a code change, review the affected features and update the prose before running `check_coverage.py --refresh`. Normal PDF builds run the coverage check automatically. On other platforms, pass the corresponding executable path or omit `--verify-executable` for a source-only check.

The dislocation example is a geometric starting configuration; its generation reports a lattice-family change warning. It preserves 500 host atoms and is included to demonstrate the operation, not to claim a relaxed defect structure.
