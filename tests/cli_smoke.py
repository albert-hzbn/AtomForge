"""Headless integration tests: python tests/cli_smoke.py path/to/AtomForge."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())


def run(*args, success=True):
    result = subprocess.run([exe, *map(str, args)], capture_output=True, text=True, timeout=60)
    if success is not None and (result.returncode == 0) != success:
        raise AssertionError(f"{args}: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
    return result


with tempfile.TemporaryDirectory(prefix="atomforge_cli_") as folder:
    root = Path(folder)
    run("--version")
    for mode in ("bulk", "gb", "poly", "nano", "amorphous", "sss", "dislocation", "custom"):
        run("--help", mode)
    source = root / "cu.cif"
    run("--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", source)
    assert source.stat().st_size > 0
    for mode, options in (
        ("nano", ["--radius", "5"]),
        ("poly", ["--grains", "2", "--sizex", "10", "--sizey", "10", "--sizez", "10"]),
        ("sss", ["--frac", "Cu=0.5,Ni=0.5"]),
        ("gb", ["--sigma", "5", "--axis", "0 0 1"]),
    ):
        output = root / f"{mode}.cif"
        run("--build", mode, "--input", source, *options, "--output", output)
        assert output.stat().st_size > 0
    # Defect, lattice and nanostructure builders on a conventional fcc cell.
    conventional = root / "cu_conv.vasp"
    run("--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", conventional)
    for mode in ("vacancy", "strain", "primitive", "surface", "sqs", "nanowire", "core-shell"):
        assert mode.upper().replace("-", "") in run("--help", mode).stdout.upper().replace("-", "")
    for mode, options in (
        ("vacancy", ["--count", "1"]),
        ("strain", ["--exx", "0.01"]),
        ("primitive", []),
        ("surface", ["--h", "1", "--k", "1", "--l", "1", "--layers", "3"]),
        ("sqs", ["--element", "Cu 0.5", "--element", "Ni 0.5", "--steps", "200"]),
        ("nanowire", ["--radius", "6"]),
        ("core-shell", ["--core-radius", "2", "--core-element", "Au", "--shell-element", "Ag"]),
    ):
        output = root / f"{mode}.xyz"
        run("--build", mode, "--input", conventional, *options, "--output", output)
        lines = output.read_text().splitlines()
        assert int(lines[0]) > 0 and len(lines) >= int(lines[0]) + 2, mode
    run("--build", "vacancy", "--input", conventional, "--count", "99", "--output", root / "too_many.xyz", success=False)
    batch = root / "batch.json"
    batch.write_text(json.dumps({"tool": "harmonic-thermodynamics", "base": {"energies_eV": [[0.01, 0.02]]},
                                 "sweep": {"temperatures_K": [[100], [300]]}, "collect": ["heat_capacity_eV_per_K.0"]}))
    run("--science-batch", batch, "--output", root / "batch_result.json", "--csv", root / "batch.csv")
    table = (root / "batch.csv").read_text().splitlines()
    assert table[0] == "temperatures_K,heat_capacity_eV_per_K.0,error" and len(table) == 3, table
    # Structure pipeline: text steps, files, and a real stdin/stdout chain.
    piped = root / "piped.xyz"
    run("--pipe", "replicate 2 2 2 | select-expression fz > 0.5 | delete-selected", "--input", source, "--output", piped, "--quiet")
    count = int(piped.read_text().splitlines()[0])
    first = subprocess.run([exe, "--pipe", "replicate 2 2 2 | select-element Cu | invert-selection", "--input", str(source), "--output", "-", "--quiet"],
                           capture_output=True, text=True, timeout=60)
    assert first.returncode == 0 and "selection:I:1" in first.stdout.splitlines()[1], first.stderr
    second = subprocess.run([exe, "--pipe", "select-expression fz > 0.5 | delete-selected", "--input", "-", "--output", str(root / "chain.xyz")],
                            input=first.stdout, capture_output=True, text=True, timeout=60)
    assert second.returncode == 0 and "2 delete-selected" in second.stderr, second.stderr
    assert int((root / "chain.xyz").read_text().splitlines()[0]) == count, "a shell pipe gives the same result as one pipeline"
    saved = root / "steps.json"
    run("--pipe", "replicate 2 1 1 | wrap", "--input", source, "--output", root / "wrapped.xyz", "--save-pipeline", saved, "--quiet")
    assert json.loads(saved.read_text())["modifiers"][0]["type"] == "replicate"
    run("--pipeline", saved, "--input", source, "--output", root / "again.xyz", "--quiet")
    # Atom clouds for very large structures: generate, convert (dump, XYZ, other formats), describe.
    cloud = root / "cu.afcloud"
    run("--cloud", "generate", "--lattice", "fcc", "--a", "3.615", "--cells", "10", "10", "10", "--element", "Cu", "--output", cloud, "--quiet")
    info = run("--cloud", "info", cloud).stdout
    assert "atoms    4000" in info and "Cu" in info, info
    dump = root / "md.dump"
    dump.write_text("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                    "ITEM: ATOMS id type xs ys zs\n1 1 0.1 0.1 0.1\n2 2 0.5 0.5 0.5\n3 1 0.9 0.2 0.3\n")
    run("--cloud", "build", "--input", dump, "--types", "Fe,Cr", "--output", root / "md.afcloud", "--quiet")
    info = run("--cloud", "info", root / "md.afcloud").stdout
    assert "atoms    3" in info and "Fe" in info and "Cr" in info, info
    run("--cloud", "build", "--input", source, "--output", root / "cif.afcloud", "--quiet")
    assert "atoms    4" in run("--cloud", "info", root / "cif.afcloud").stdout
    run("--cloud", "info", root / "missing.afcloud", success=False)
    run("--cloud", "generate", "--lattice", "nope", "--a", "3", "--element", "Cu", "--output", root / "bad.afcloud", success=False)
    # Every pipeline step, Build and Edit operations included, is also in the menus.
    steps = json.loads(run("--pipe", "--list-json").stdout)
    missing = [step["id"] for step in steps if not step.get("menu")]
    assert not missing, f"pipeline steps without a menu entry: {missing}"
    assert {"build-bulk", "build-vacancy", "build-surface", "insert-interstitials", "select-expression"} <= {step["id"] for step in steps}
    listing = run("--pipe", "--list").stdout
    assert "select-expression" in listing and "replicate" in listing, listing[:200]
    bad = run("--pipe", "replicate 2 2", "--input", source, "--output", root / "bad.xyz", success=False)
    assert "counts needs three numbers" in bad.stderr, bad.stderr
    failing = run("--pipe", "select-property 0 1", "--input", source, "--output", root / "bad.xyz", success=False)
    assert "step 1" in failing.stderr, failing.stderr
    # Native science tools: request file in, result JSON and generated files out.
    request = root / "dft.json"
    request.write_text(json.dumps({"structure": {"file": conventional.name}, "points_per_segment": 20}))
    dft = run("--science", "dft-inputs", "--input", request, "--output", root / "dft_result.json", "--files", root / "dft_files", success=None)
    if dft.returncode == 0:
        assert (root / "dft_files" / "KPOINTS").read_text().splitlines()[2] == "Line-mode"
        assert "K_POINTS crystal_b" in (root / "dft_files" / "qe_band_cards.in").read_text()
    else:
        assert "spglib" in dft.stderr, dft.stderr  # builds without spglib cannot find symmetry paths
    # Animated GIF from a trajectory; skipped where no OpenGL context exists (headless CI).
    trajectory = root / "traj.xyz"
    trajectory.write_text("".join(f'2\nLattice="6 0 0 0 6 0 0 0 6" Properties=species:S:1:pos:R:3\nCu 1 1 1\nCu {2 + 0.3 * k} 3 3\n' for k in range(6)))
    gif = root / "traj.gif"
    animation = run("--render", "--input", trajectory, "--output", gif, "--width", 64, "--height", 48, "--every", 2, success=None)
    if animation.returncode == 0:
        data = gif.read_bytes()
        assert data[:6] == b"GIF89a" and data.endswith(b";") and b"NETSCAPE2.0" in data, data[:16]
        assert "Saved 3 frames" in animation.stdout, animation.stdout
    else:
        assert "OpenGL" in animation.stderr, animation.stderr
    for fraction in ("Cu=nan", "Cu=inf", "Cu=0.5garbage"):
        output = root / "invalid.cif"
        run("--build", "sss", "--input", source, "--frac", fraction, "--output", output, success=False)
        assert not output.exists()
    # A CIF that names its space group but lists no symmetry operations is
    # expanded (rock salt: 8 atoms); XYZ symbols in any case; the extended XYZ
    # cell survives a round trip through XYZ and CIF.
    rocksalt = root / "rocksalt.cif"
    rocksalt.write_text("data_NaCl\n_cell_length_a 5.64\n_cell_length_b 5.64\n_cell_length_c 5.64\n"
                        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                        "_space_group_IT_number 225\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
                        "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nNa1 Na 0 0 0\nCl1 Cl 0.5 0.5 0.5\n")
    expanded = root / "rocksalt.xyz"
    run("--convert", "--input", rocksalt, "--output", expanded, "--format", "xyz")
    lines = expanded.read_text().splitlines()
    assert lines[0].strip() == "8", lines[:3]
    # Orientation colours reconstructed from the geometry are not written beside it.
    assert not (root / "rocksalt.atomforge-ipf").exists(), "a reconstructed orientation sidecar was written"
    assert lines[1].startswith('Lattice="5.64 0 0 0 5.64 0 0 0 5.64"'), lines[1]
    back = root / "rocksalt_back.cif"
    run("--convert", "--input", expanded, "--output", back)
    assert "_cell_length_a 5.64" in back.read_text(), back.read_text()[:300]
    upper = root / "upper.xyz"
    upper.write_text("3\n\nCU 0 0 0\nni 1.8 0 0\nO 0 1.8 0\n")
    symbols = root / "upper_out.xyz"
    run("--convert", "--input", upper, "--output", symbols, "--format", "xyz")
    assert [line.split()[0] for line in symbols.read_text().splitlines()[2:5]] == ["Cu", "Ni", "O"], symbols.read_text()
    # Extended XYZ columns in any order (an id first, extra columns after).
    ids = root / "ids.xyz"
    ids.write_text('2\nProperties=id:I:1:species:S:1:pos:R:3:forces:R:3 Lattice="4 0 0 0 4 0 0 0 4"\n'
                   "1 Na 0 0 0 0 0 0\n2 Cl 2 2 2 0 0 0\n")
    reordered = root / "ids_out.xyz"
    run("--convert", "--input", ids, "--output", reordered, "--format", "xyz")
    assert [line.split()[0] for line in reordered.read_text().splitlines()[2:4]] == ["Na", "Cl"], reordered.read_text()
    assert reordered.read_text().splitlines()[1].startswith('Lattice="4 0 0 0 4 0 0 0 4"'), reordered.read_text()
    # POSCAR: a negative scale is the cell volume; VASP 4 names the species on
    # the title line instead of an element line.
    for name, text, expected in (
        ("vol_POSCAR", "volume\n-64.0\n1 0 0\n0 1 0\n0 0 1\nCu\n1\nDirect\n0 0 0\n", ["Cu"]),
        ("v4_POSCAR", "Na Cl\n1.0\n4 0 0\n0 4 0\n0 0 4\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n", ["Na", "Cl"]),
        ("POSCAR-pv", "labels\n1.0\n4 0 0\n0 4 0\n0 0 4\nFe_pv O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n", ["Fe", "O"]),
    ):
        poscar = root / name
        poscar.write_text(text)
        converted = root / (name + ".xyz")
        run("--convert", "--input", poscar, "--output", converted, "--format", "xyz")
        lines = converted.read_text().splitlines()
        assert lines[1].startswith('Lattice="4 0 0 0 4 0 0 0 4"'), (name, lines[1])
        assert [line.split()[0] for line in lines[2:]] == expected, (name, lines)
    # The space group sets the crystal system; fractions are exact.
    hcp = root / "hcp.xyz"
    run("--build", "bulk", "--spacegroup", "194", "--a", "3.21", "--c", "5.21", "--atom", "Mg 1/3 2/3 1/4", "--output", hcp)
    hcp_lines = hcp.read_text().splitlines()
    assert hcp_lines[1].startswith('Lattice="3.21 0 0 -1.605 2.779941546 0 0 0 5.21"'), hcp_lines[1]
    assert hcp_lines[2].split()[1:] == ["0.0000000000", "1.8532943641", "1.3025000000"], hcp_lines[2]
    run("--build", "bulk", "--system", "cubic", "--spacegroup", "194", "--atom", "Mg 0 0 0", "--output", root / "no.xyz", success=False)
    # LAMMPS data: written and read back with the same cell and atoms.
    data = root / "hcp.data"
    run("--convert", "--input", hcp, "--output", data)
    assert "1 24.305  # Mg" in data.read_text(), data.read_text()
    back = root / "hcp_back.xyz"
    run("--convert", "--input", data, "--output", back, "--format", "xyz")
    assert back.read_text().splitlines() == hcp_lines, back.read_text()
    # A disordered CIF (Cu 0.6 and Ni 0.4 on one FCC site): one atom per site,
    # the majority species, and a note saying so.
    disordered = root / "disordered.cif"
    disordered.write_text("data_x\n_cell_length_a 3.6\n_cell_length_b 3.6\n_cell_length_c 3.6\n"
                          "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                          "_symmetry_space_group_name_H-M 'F m -3 m'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
                          "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
                          "Ni1 Ni 0 0 0 0.4\nCu1 Cu 0 0 0 0.6\n")
    ordered = root / "disordered.xyz"
    note = run("--convert", "--input", disordered, "--output", ordered, "--format", "xyz")
    assert [line.split()[0] for line in ordered.read_text().splitlines()[2:]] == ["Cu"] * 4, ordered.read_text()
    assert "partial occupancy" in note.stderr, note.stderr
    # Selective dynamics: read, reported, written back; fix-atoms sets it.
    fixed = root / "fixed_POSCAR"
    fixed.write_text("sd\n1.0\n4 0 0\n0 4 0\n0 0 4\nNa Cl\n1 1\nSelective dynamics\nDirect\n0 0 0 T T T\n0.5 0.5 0.5 F F T\n")
    rewritten = root / "fixed_again.vasp"
    note = run("--convert", "--input", fixed, "--output", rewritten)
    assert "1 atom has fixed coordinates" in note.stderr, note.stderr
    lines = rewritten.read_text().splitlines()
    assert lines[7] == "Selective dynamics" and lines[9].endswith("T T T") and lines[10].endswith("F F T"), lines
    pinned = root / "pinned.vasp"
    run("--pipe", "select-element Cl | fix-atoms axes=z", "--input", rewritten, "--output", pinned)
    assert pinned.read_text().splitlines()[10].endswith("F F F"), pinned.read_text()
    # LAMMPS data: force-field labels (CA, HO) do not override the masses; a
    # six-column file that could be charge or molecular is refused.
    ff = root / "ff.data"
    ff.write_text("x\n\n2 atoms\n2 atom types\n0 4 xlo xhi\n0 4 ylo yhi\n0 4 zlo zhi\n\nMasses\n\n1 12.011 # CA\n2 1.008 # HO\n\n"
                  "Atoms # full\n\n1 1 1 0.0 1 1 1\n2 1 2 0.0 2 2 2\n")
    ff_out = root / "ff.xyz"
    run("--convert", "--input", ff, "--output", ff_out, "--format", "xyz")
    assert [line.split()[0] for line in ff_out.read_text().splitlines()[2:]] == ["C", "H"], ff_out.read_text()
    ambiguous = root / "mol.data"
    ambiguous.write_text("x\n\n2 atoms\n2 atom types\n0 4 xlo xhi\n0 4 ylo yhi\n0 4 zlo zhi\n\nMasses\n\n1 63.546\n2 58.693\n\n"
                         "Atoms\n\n1 1 2 1 1 1\n2 1 1 2 2 2\n")
    run("--convert", "--input", ambiguous, "--output", root / "mol.xyz", "--format", "xyz", success=False)
    # Two species a hair apart on one site are one site.
    straddle = root / "straddle.cif"
    straddle.write_text("data_x\n_cell_length_a 3.61\n_cell_length_b 3.61\n_cell_length_c 3.61\n"
                        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                        "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n"
                        "_atom_site_fract_z\n_atom_site_occupancy\nCu1 Cu 0.5 0.5 0.5 0.4\nNi1 Ni 0.49999 0.49999 0.49999 0.6\n")
    straddle_out = root / "straddle.xyz"
    run("--convert", "--input", straddle, "--output", straddle_out, "--format", "xyz")
    assert straddle_out.read_text().splitlines()[0].strip() == "1", straddle_out.read_text()
    # Output files named the VASP way are written as POSCAR, without an extension too.
    vasp_dir = root / "vasp.run"
    vasp_dir.mkdir()
    for name in ("POSCAR", "CONTCAR_relaxed"):
        run("--convert", "--input", source, "--output", vasp_dir / name)
        assert (vasp_dir / name).read_text().splitlines()[0] == "Cu", name
    # CIF keeps full precision: 1/6 of a 10.83 A cell stays 1.805 A.
    big = root / "big.cif"
    run("--pipe", "replicate 3 3 3", "--input", source, "--output", big)
    big_xyz = root / "big.xyz"
    run("--convert", "--input", big, "--output", big_xyz, "--format", "xyz")
    for line in big_xyz.read_text().splitlines()[2:]:
        for value in map(float, line.split()[1:]):
            assert abs(value / 1.805 - round(value / 1.805)) < 1e-9, line
    # Paths with accented names are read and written (UTF-8 on Windows too).
    accented = root / "døssier é"
    accented.mkdir()
    copied = accented / "cu ü.cif"
    run("--convert", "--input", source, "--output", copied)
    run("--convert", "--input", copied, "--output", accented / "cu ü.xyz", "--format", "xyz")
    assert (accented / "cu ü.xyz").stat().st_size > 0
    run("--build", "bulk", "--a", "nan", "--output", root / "invalid.cif", success=False)
    run("--build", "unknown", success=False)
print("CLI smoke tests passed (version, 15 help modes, 12 builders, 6 invalid inputs)")
