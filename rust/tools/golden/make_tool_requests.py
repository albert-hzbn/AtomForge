"""Requests for the science-tool golden comparison (tests/fixtures/tools/requests.json):
deterministic inputs for every ported tool, plus error cases. Run from the
repository root: python rust/tools/golden/make_tool_requests.py"""
import json
import math
import os
import random

random.seed(20261120)
R = lambda: random.uniform(-1, 1)


def fcc(n, a=3.61, jitter=0.0):
    basis = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in basis:
                    pts.append([(i + b[0]) * a + jitter * R(), (j + b[1]) * a + jitter * R(), (k + b[2]) * a + jitter * R()])
    return pts, [[n * a, 0, 0], [0, n * a, 0], [0, 0, n * a]]


def bcc(n, a=2.87):
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in [(0, 0, 0), (.5, .5, .5)]:
                    pts.append([(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a])
    return pts, [[n * a, 0, 0], [0, n * a, 0], [0, 0, n * a]]


def hcp(n, a=2.95, c=4.68):
    a1, a2 = (a, 0, 0), (-a / 2, a * math.sqrt(3) / 2, 0)
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for (u, v, w) in [(1 / 3, 2 / 3, 0.25), (2 / 3, 1 / 3, 0.75)]:
                    x = (i + u) * a1[0] + (j + v) * a2[0]
                    y = (i + u) * a1[1] + (j + v) * a2[1]
                    pts.append([x, y, (k + w) * c])
    return pts, [[n * a1[0], n * a1[1], 0], [n * a2[0], n * a2[1], 0], [0, 0, n * c]]


def diamond(n, a=5.43):
    base = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
    basis = base + [(x + .25, y + .25, z + .25) for (x, y, z) in base]
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in basis:
                    pts.append([(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a])
    return pts, [[n * a, 0, 0], [0, n * a, 0], [0, 0, n * a]]


cases = []


def add(tool, request):
    cases.append({"tool": tool, "request": request})


# ---- trajectories
frames, atoms = 7, 4
walk = [[[R() * 2, R() * 2, R() * 2] for _ in range(atoms)]]
for f in range(1, frames):
    walk.append([[walk[-1][a][k] + 0.4 * R() for k in range(3)] for a in range(atoms)])
add("msd", {"positions": walk, "timestep_fs": 2.0})
add("msd", {"positions": walk, "timestep_fs": 0.5, "remove_drift": True, "max_lag": 4})
wrapped = [[[x % 3.0 for x in atom] for atom in frame] for frame in walk]
add("msd", {"positions": wrapped, "timestep_fs": 1.0, "wrapped": True, "cell": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]})
add("msd", {"positions": wrapped, "timestep_fs": 1.0, "wrapped": True})
add("msd", {"positions": walk, "timestep_fs": 1.0, "max_lag": 9})
lag = [i * 10.0 for i in range(12)]
msd = [0.0] + [0.3 * t + 0.1 * math.sin(t) for t in lag[1:]]
add("diffusion", {"lag_fs": lag, "msd_A2": msd, "fit_range_fs": [20, 100]})
add("diffusion", {"lag_fs": lag, "msd_A2": msd, "fit_range_fs": [10, 110], "dimensions": 2})
add("diffusion", {"lag_fs": lag, "msd_A2": msd, "fit_range_fs": [10, 30], "dimensions": 4})
add("diffusion", {"lag_fs": lag, "msd_A2": msd, "fit_range_fs": [100, 300]})
vel = [[[R(), R(), R()] for _ in range(3)] for _ in range(9)]
add("vacf", {"velocities": vel, "timestep_fs": 1.0})
add("vacf", {"velocities": vel, "timestep_fs": 2.0, "normalize": True, "remove_drift": True, "max_lag": 5})
vel16 = [[[math.sin(0.7 * f + a) + 0.1 * R(), math.cos(0.3 * f * (a + 1)), R() * 0.2] for a in range(3)] for f in range(16)]
add("vibrational-spectrum", {"velocities": vel16, "timestep_fs": 1.0})
add("vibrational-spectrum", {"velocities": vel16[:15], "timestep_fs": 0.5, "window": "none", "masses": [1.0, 12.0, 16.0], "remove_mean": False})
add("vibrational-spectrum", {"velocities": vel16[:3], "timestep_fs": 0.5})
add("vibrational-spectrum", {"velocities": vel16, "timestep_fs": 0.5, "window": "blackman"})

# ---- local environments
ref, cell = fcc(2)
strain = [[1.01, 0.002, 0], [0, 0.995, 0.003], [0, 0, 1.0]]
cur = [[sum(p[i] * strain[i][j] for i in range(3)) + 0.01 * R() for j in range(3)] for p in ref]
ccell = [[sum(r[i] * strain[i][j] for i in range(3)) for j in range(3)] for r in cell]
add("local-strain", {"reference": ref, "current": cur, "cutoff_A": 3.0, "reference_cell": cell, "current_cell": ccell, "pbc": [True, True, True]})
add("local-strain", {"reference": ref[:13], "current": cur[:13], "cutoff_A": 3.0})
jit, _ = fcc(2, jitter=0.05)
add("centrosymmetry", {"positions": jit, "cutoff_A": 3.2, "cell": cell, "pbc": [True, True, True]})
add("centrosymmetry", {"positions": jit, "cutoff_A": 3.2, "neighbors": 8})
add("centrosymmetry", {"positions": jit, "cutoff_A": 3.2, "neighbors": 7})
b, bcell = bcc(3)
h, hcell = hcp(3)
add("structure-type", {"positions": jit, "cell": cell, "pbc": [True, True, True]})
add("structure-type", {"positions": b, "cell": bcell, "pbc": [True, True, True]})
add("structure-type", {"positions": h, "cell": hcell, "pbc": [True, True, True]})
add("structure-type", {"positions": jit, "cutoff_A": 3.0})
add("structure-type", {"positions": jit})
d, dcell = diamond(2)
add("diamond-structure", {"positions": d, "cell": dcell, "pbc": [True, True, True]})
dj = [[x + 0.02 * R() for x in p] for p in d]
add("diamond-structure", {"positions": dj[:40], "cutoff_A": 2.6})
add("bond-order", {"positions": jit, "cutoff_A": 3.0, "cell": cell, "pbc": [True, True, True]})
add("bond-order", {"positions": b, "cutoff_A": 2.6, "cell": bcell, "pbc": [True, True, True], "degrees": [2, 4, 6, 8, 12]})
add("bond-order", {"positions": jit, "cutoff_A": 3.0, "degrees": [6, 6]})
add("bond-order", {"positions": jit, "cutoff_A": 3.0, "degrees": [40]})
sites, scell = fcc(1)
moved = [[x + 0.3 * R() for x in p] for p in sites[1:]] + [[1.8, 1.8, 1.8], [0.1, 0.1, 3.5]]
add("wigner-seitz", {"reference_sites": sites, "positions": moved, "cell": scell, "pbc": [True, True, True]})
add("wigner-seitz", {"reference_sites": sites, "positions": moved, "cell": scell, "current_cell": [[3.7, 0, 0], [0, 3.6, 0], [0, 0, 3.61]], "pbc": [True, True, False]})
add("wigner-seitz", {"reference_sites": sites, "positions": []})
qv = [[2 * math.pi / 3.61 * h_, 2 * math.pi / 3.61 * k_, 2 * math.pi / 3.61 * l_] for (h_, k_, l_) in [(1, 1, 1), (2, 0, 0), (1, 0, 0), (2, 2, 0), (0.5, 0.3, 0.1)]]
add("structure-factor", {"positions": jit, "q_vectors": qv})
add("structure-factor", {"positions": [jit, ref], "q_vectors": qv, "weights": [29 if i % 2 else 8 for i in range(len(jit))]})
add("structure-factor", {"positions": jit, "q_vectors": qv, "weights": [0] * len(jit)})

# ---- electronic
bands = [[[-5 + b * 1.5 + 0.2 * math.sin(k + s + b) for b in range(5)] for k in range(4)] for s in range(2)]
add("band-gap", {"energies_eV": bands, "fermi_eV": -1.2})
add("band-gap", {"energies_eV": bands[0], "fermi_eV": 0.4, "tolerance_eV": 0.05})
add("band-gap", {"energies_eV": bands[0], "fermi_eV": -2.0})
add("band-gap", {"energies_eV": bands[0], "fermi_eV": 9.0})
center = [0.1, -0.05, 0.2]
kp, en = [], []
for i in range(14):
    k = [center[0] + 0.05 * R(), center[1] + 0.05 * R(), center[2] + 0.05 * R()]
    x, y, z = k[0] - center[0], k[1] - center[1], k[2] - center[2]
    kp.append(k)
    en.append(1.2 + 0.1 * x + 3.0 * x * x + 2.0 * y * y + 4.5 * z * z + 0.8 * x * y + 0.002 * R())
add("effective-mass", {"kpoints_inv_A": kp, "energies_eV": en, "center_inv_A": center})
line = [[0.01 * i, 0, 0] for i in range(12)]
add("effective-mass", {"kpoints_inv_A": line, "energies_eV": [0.5 * k[0] ** 2 for k in line], "center_inv_A": [0, 0, 0]})
dist = [0.5 * i for i in range(40)]
pot = [-8 + 9 / (1 + math.exp(-(z - 8))) + 0.001 * R() for z in dist]
add("work-function", {"distance_A": dist, "potential_eV": pot, "vacuum_range_A": [14, 19], "fermi_eV": -3.1})
add("work-function", {"distance_A": dist, "potential_eV": pot, "vacuum_range_A": [4, 12], "fermi_eV": -3.1, "max_slope_eV_per_A": 0.5})
add("work-function", {"distance_A": dist, "potential_eV": pot, "vacuum_range_A": [14, 25], "fermi_eV": -3.1})


# ---- fits
def bm(v, e0=-3.7, b0=0.88, bp=4.6, v0=11.8):
    x = (v0 / v) ** (2 / 3) - 1
    return e0 + 9 / 16 * b0 * v0 * (2 * x * x + (bp - 4) * x ** 3)


vols = [10.6 + 0.4 * i for i in range(7)]
add("equation-of-state", {"volumes_A3": vols, "energies_eV": [bm(v) + 0.0005 * R() for v in vols]})
add("equation-of-state", {"volumes_A3": vols[:4], "energies_eV": [bm(v) for v in vols[:4]]})
add("equation-of-state", {"volumes_A3": [12 + 0.2 * i for i in range(6)], "energies_eV": [bm(12 + 0.2 * i) for i in range(6)]})
C = [[250, 150, 150, 0, 0, 0], [150, 250, 150, 0, 0, 0], [150, 150, 250, 0, 0, 0], [0, 0, 0, 120, 0, 0], [0, 0, 0, 0, 120, 0], [0, 0, 0, 0, 0, 120]]
strains = [[0.01 * R() for _ in range(6)] for _ in range(12)]
stresses = [[0.3 + sum(C[i][j] * s[j] for j in range(6)) + 0.01 * R() for i in range(6)] for s in strains]
add("elastic-tensor", {"strains": strains, "stresses_GPa": stresses})
add("elastic-tensor", {"strains": strains[:6], "stresses_GPa": stresses[:6]})
unstable = [[0.3 + sum((-C[i][j] if i == j == 3 else C[i][j]) * s[j] for j in range(6)) for i in range(6)] for s in strains]
add("elastic-tensor", {"strains": strains, "stresses_GPa": unstable})

# ---- phonons
modes = [[0.0, 0.0, 0.0, 0.012 + 0.001 * q, 0.018, 0.025 + 0.002 * q] for q in range(5)]
add("phonon-dos", {"energies_eV": modes, "energy_grid_eV": [0.0005 * i for i in range(70)], "sigma_eV": 0.002})
add("phonon-dos", {"energies_eV": [[-0.002, 0.01], [0.004, 0.02]], "energy_grid_eV": [0.001 * i for i in range(30)], "weights": [1, 3]})
add("phonon-dos", {"energies_eV": modes, "energy_grid_eV": [0.1]})
add("harmonic-thermodynamics", {"energies_eV": modes, "temperatures_K": [0, 10, 100, 300, 1000]})
add("harmonic-thermodynamics", {"energies_eV": modes, "temperatures_K": [300], "weights": [1, 2, 3, 4, 5], "zero_tolerance_eV": 1e-4})
add("harmonic-thermodynamics", {"energies_eV": [[-0.01, 0.02]], "temperatures_K": [300]})

# ---- trajectory tools (extended XYZ files written here)
def write_xyz(path, frames, cell, times=None, pbc="T T T"):
    with open(path, "w") as f:
        for n, frame in enumerate(frames):
            f.write("%d\n" % len(frame))
            lattice = " ".join(repr(float(x)) for row in cell for x in row)
            extra = ' time_fs=%r' % times[n] if times else ""
            f.write('Lattice="%s" Properties=species:S:1:pos:R:3 pbc="%s"%s\n' % (lattice, pbc, extra))
            for symbol, (x, y, z) in frame:
                f.write("%s %r %r %r\n" % (symbol, x, y, z))


os.makedirs("rust/tests/fixtures/tools", exist_ok=True)
tcell = [[7.22, 0, 0], [0, 7.22, 0], [0, 0, 7.22]]
base, _ = fcc(2)
traj = []
pos = [list(q) for q in base]
for f in range(6):
    frame = []
    for i, q in enumerate(pos):
        if f:
            for k in range(3):
                q[k] += 0.15 * R() + (0.9 if i == 3 else 0.0)
        frame.append(("Cu" if i % 4 else "Ni", [q[k] % 7.22 for k in range(3)]))
    traj.append(frame)
write_xyz("rust/tests/fixtures/tools/traj.xyz", traj, tcell, times=[10.0 * f for f in range(6)])
write_xyz("rust/tests/fixtures/tools/traj_open.xyz", traj, tcell, pbc="F F F")
write_xyz("rust/tests/fixtures/tools/traj_changing.xyz", [traj[0], traj[1][:-1]], tcell)
T = "rust/tests/fixtures/tools/traj.xyz"
add("trajectory-average", {"trajectory_file": T})
add("trajectory-average", {"trajectory_file": T, "unwrap": False, "first_frame": 1, "every": 2})
add("trajectory-average", {"trajectory_file": T, "first_frame": 6})
add("trajectory-average", {"trajectory_file": "rust/tests/fixtures/tools/traj_changing.xyz"})
add("trajectory-average", {})
add("unwrap-trajectory", {"trajectory_file": T})
add("unwrap-trajectory", {"trajectory_file": T, "smooth_window": 3, "every": 2})
add("unwrap-trajectory", {"trajectory_file": "rust/tests/fixtures/tools/traj_open.xyz", "smooth_window": 4})
add("time-series", {"trajectory_file": T, "analysis": "none"})
add("time-series", {"trajectory_file": T, "analysis": "structure-type", "every": 2})
add("time-series", {"trajectory_file": T, "analysis": "centrosymmetry", "cutoff_A": 3.0, "first_frame": 2})
add("time-series", {"trajectory_file": T, "analysis": "bond-order", "property": "Coordination", "cutoff_A": 2.9})
add("frame-difference", {"trajectory_file": T})
add("frame-difference", {"trajectory_file": T, "reference_frame": 1, "frame": -2})
add("frame-difference", {"trajectory_file": "rust/tests/fixtures/tools/traj_open.xyz", "frame": 3})
add("frame-difference", {"trajectory_file": T, "analysis": "centrosymmetry", "cutoff_A": 3.0})
add("frame-difference", {"trajectory_file": T, "analysis": "structure-type", "property": "nothing"})
add("frame-difference", {"trajectory_file": T, "frame": 9})
add("frame-difference", {"trajectory_file": "rust/tests/fixtures/tools/traj_changing.xyz"})

# ---- per-atom property analyses
pa, pcell = fcc(2, jitter=0.1)
vals = [math.sin(q[0]) + 0.3 * q[2] for q in pa]
add("spatial-binning", {"positions": pa, "cell": pcell, "values": vals})
add("spatial-binning", {"positions": pa, "cell": pcell, "values": vals, "bins": [2, 3, 1], "operation": "max"})
add("spatial-binning", {"positions": pa + [[-0.2, 1, 1]], "cell": pcell, "pbc": [False, True, True], "bins": [4, 1, 2], "operation": "count"})
add("spatial-binning", {"positions": pa, "cell": pcell, "operation": "density", "bins": [3, 3, 3]})
add("spatial-binning", {"positions": pa, "cell": pcell, "operation": "median"})
add("spatial-binning", {"positions": pa, "values": vals})
add("property-histogram", {"values": vals})
add("property-histogram", {"values": vals + [float("nan")], "bins": 7, "range": [-1, 2]})
add("property-histogram", {"values": vals, "mask": [i % 3 for i in range(len(vals))], "bins": 5})
add("property-histogram", {"values": [2.5] * 5, "bins": 4})
add("property-histogram", {"values": vals, "range": [2, 1]})
add("scatter-plot", {"values_x": vals, "values_y": [v * 2 - q[1] for v, q in zip(vals, pa)]})
add("scatter-plot", {"values_x": [1, float("nan"), 3], "values_y": [2, 2, float("inf")]})
add("scatter-plot", {"values_x": [1, 1, 1], "values_y": [2, 3, 4]})
add("ambient-occlusion", {"positions": pa, "directions": 12, "resolution": 32})
add("ambient-occlusion", {"positions": pa[:9], "radii": [1.2 + 0.1 * i for i in range(9)], "directions": 7, "resolution": 40})
add("ambient-occlusion", {"positions": pa[:3], "radii": [1, 0, 1]})
add("spatial-correlation", {"positions": pa, "cell": pcell, "values_a": vals})
add("spatial-correlation", {"positions": pa, "cell": pcell, "values_a": vals, "values_b": [q[1] for q in pa], "grid_spacing_A": 0.8, "bins": 12, "r_max_A": 3.3})
add("spatial-correlation", {"positions": pa, "cell": pcell, "values_a": [float("nan")] * len(pa)})

# ---- local environments
dd, ddcell = diamond(2, a=6.36)
ddj = [[x + 0.15 * R() for x in q] for q in dd]
add("chill-plus", {"positions": dd, "cell": ddcell, "pbc": [True, True, True], "cutoff_A": 3.0})
add("chill-plus", {"positions": ddj, "cell": ddcell, "pbc": [True, True, True], "cutoff_A": 3.2})
add("chill-plus", {"positions": ddj[:30], "cutoff_A": 3.2})
add("local-entropy", {"positions": jit, "cell": cell, "pbc": [True, True, True], "cutoff_A": 4.0, "sigma_A": 0.3})
add("local-entropy", {"positions": jit, "cutoff_A": 3.5, "local_density": True, "average": False})
add("local-entropy", {"positions": jit, "cell": cell, "pbc": [True, True, True], "cutoff_A": 3.0, "average_cutoff_A": 2.8, "sigma_A": 0.25})
add("local-entropy", {"positions": jit})
gr = []
for i in range(3):
    for j in range(3):
        for (u, v) in [(1 / 3, 2 / 3), (2 / 3, 1 / 3)]:
            x = (i + u) * 2.46 + (j + v) * -1.23
            y = (j + v) * 2.46 * math.sqrt(3) / 2
            gr.append([x, y, 0.0])
grcell = [[7.38, 0, 0], [-3.69, 7.38 * math.sqrt(3) / 2, 0], [0, 0, 10]]
add("find-rings", {"positions": gr, "cell": grcell, "pbc": [True, True, False], "cutoff_A": 1.6})
add("find-rings", {"positions": d, "cell": dcell, "pbc": [True, True, True], "cutoff_A": 2.5, "max_size": 8})
add("find-rings", {"positions": dj[:40], "cutoff_A": 2.6, "max_size": 7, "max_paths": 2})
add("find-rings", {"positions": d, "cutoff_A": 2.5, "max_size": 30})
add("displacements", {"reference": ref, "current": cur})
add("displacements", {"reference": ref, "current": cur, "reference_cell": cell, "current_cell": ccell, "pbc": [True, True, True], "affine_mapping": "to_reference"})
add("displacements", {"reference": ref, "current": cur, "reference_cell": cell, "current_cell": ccell, "pbc": [True, False, True], "affine_mapping": "to_current"})
add("displacements", {"reference": ref, "current": cur, "pbc": [True, True, True]})
add("displacements", {"reference": ref, "current": cur, "reference_cell": cell, "affine_mapping": "to_current"})
add("bond-distributions", {"positions": jit, "cutoff_A": 3.0, "cell": cell, "pbc": [True, True, True], "bins": 20})
add("bond-distributions", {"positions": jit, "cutoff_A": 3.0, "types": [29 if i % 3 else 8 for i in range(len(jit))], "bins": 12})
add("bond-distributions", {"positions": jit, "cutoff_A": 3.0, "types": [1, 2]})

# ---- clusters and voids
big, bigcell = fcc(3)
mask = [1 if (q[0] < 4 or q[2] > 7) else 0 for q in big]
add("cluster-analysis", {"positions": big, "cell": bigcell, "pbc": [True, True, True], "cutoff_A": 2.6})
add("cluster-analysis", {"positions": big, "cell": bigcell, "pbc": [True, True, True], "cutoff_A": 2.6, "mask": mask})
add("cluster-analysis", {"positions": big, "cutoff_A": 2.6, "mask": [1 if (i % 7 == 0 or i % 11 == 0) else 0 for i in range(len(big))], "min_size": 2})
add("cluster-analysis", {"positions": big, "pbc": [True, False, False]})
holes = [q for q in big if not (2.5 < q[0] < 6 and 2.5 < q[1] < 6 and 2.5 < q[2] < 6)]
add("void-analysis", {"positions": holes, "cell": bigcell, "pbc": [True, True, True], "grid_spacing_A": 0.45})
add("void-analysis", {"positions": holes, "cell": bigcell, "pbc": [True, True, False], "grid_spacing_A": 0.5, "probe_radius_A": 0.5, "atomic_radius_A": 1.1})
add("void-analysis", {"positions": holes[:60], "grid_spacing_A": 0.6, "probe_radius_A": 0.3})
add("void-analysis", {"positions": holes, "cell": bigcell, "pbc": [True, True, True], "probe_radius_A": -1})

# ---- bond orders and LAMMPS regions
benz = [[1.397 * math.cos(k * math.pi / 3), 1.397 * math.sin(k * math.pi / 3), 0] for k in range(6)]
benz += [[2.481 * math.cos(k * math.pi / 3), 2.481 * math.sin(k * math.pi / 3), 0] for k in range(6)]
add("bond-orders", {"positions": benz, "types": [6] * 6 + [1] * 6})
add("bond-orders", {"positions": jit, "types": [29] * len(jit), "cell": cell, "pbc": [True, True, True], "tolerance": 1.1})
add("bond-orders", {"positions": [[0, 0, 0], [1.2, 0, 0], [2.3, 0.1, 0], [3.5, 0, 0]], "types": [6, 6, 8, 120]})
add("bond-orders", {"positions": benz, "types": [6] * 5})
regions = "region box block 0 3.5 0 7.22 INF INF\nregion ball sphere 3.6 3.6 3.6 2.0 side out # comment\nregion tube cylinder z 3.61 3.61 1.5 EDGE 5 units box"
add("lammps-regions", {"positions": jit, "cell": cell, "regions": regions, "export_mesh": False})
add("lammps-regions", {"positions": jit, "cell": cell, "regions": "region c1 cone x 3 3 1 2 0.5 6;region p1 plane 3 3 3 1 1 0;region pr prism 0 4 0 4 0 4 1 0.5 0.2 units box", "resolution": 10})
add("lammps-regions", {"positions": jit, "regions": "region half block 0 1 EDGE EDGE EDGE EDGE", "lattice_scale": 3.61, "resolution": 8})
add("lammps-regions", {"positions": jit, "regions": "fix 1 all nve"})
add("lammps-regions", {"positions": jit, "regions": "region r1 ellipsoid 0 0 0 1"})
add("lammps-regions", {"positions": jit, "regions": "region r1 sphere 0 0 0 1x"})
add("lammps-regions", {"positions": jit, "regions": "region r1 sphere 0 0 0 1 move v_x"})

# ---- surface meshes
add("surface-mesh", {"positions": holes, "cell": bigcell, "pbc": [True, True, True], "probe_radius_A": 1.5, "grid_spacing_A": 0.8})
add("surface-mesh", {"positions": jit[:20], "probe_radius_A": 2.0, "grid_spacing_A": 0.7})
add("surface-mesh", {"positions": jit, "probe_radius_A": 2.5, "grid_spacing_A": 0.9, "mask": [1 if q[2] < 3 else 0 for q in jit], "export_mesh": False})
add("surface-mesh", {"positions": jit, "cell": cell, "pbc": [True, True, False]})
add("surface-mesh", {"positions": jit, "mask": [0] * len(jit)})

# ---- Voronoi, VoroTop, CNA
add("voronoi", {"positions": jit, "cell": cell, "pbc": [True, True, True]})
add("voronoi", {"positions": jit[:14], "face_area_threshold_A2": 0.5, "relative_face_area_threshold": 0.02})
add("voronoi", {"positions": b, "cell": bcell, "pbc": [True, True, True], "edge_length_threshold_A": 0.3, "cutoff_A": 3.0})
add("voronoi", {"positions": jit, "face_area_threshold_A2": -1})
add("vorotop", {"positions": jit, "cell": cell, "pbc": [True, True, True], "samples": 3, "noise": 0.05, "export_filter": True})
add("vorotop", {"positions": b, "cell": bcell, "pbc": [True, True, True], "samples": 2, "noise": 0.08})
add("vorotop", {"positions": jit, "noise": 0.5})
add("vorotop", {"positions": jit, "cell": cell, "pbc": [True, True, True], "filter_file": "rust/tests/fixtures/tools/missing.filter"})
add("cna", {"positions": jit, "cell": cell, "pbc": [True, True, True]})
add("cna", {"positions": b, "cell": bcell, "pbc": [True, True, True]})
add("cna", {"positions": h, "cell": hcell, "pbc": [True, True, True]})
add("cna", {"positions": jit, "cell": cell, "pbc": [True, True, True], "cutoff_A": 3.0})
add("cna", {"positions": jit[:20]})
add("cna", {"positions": jit, "cutoff_A": -1})

# ---- template matching
bigj = [[x + 0.04 * R() for x in q] for q in big]
add("ptm", {"positions": bigj, "cell": bigcell, "pbc": [True, True, True]})
add("ptm", {"positions": b, "cell": bcell, "pbc": [True, True, True], "structures": ["bcc", "sc", "fcc"]})
add("ptm", {"positions": h, "cell": hcell, "pbc": [True, True, True], "rmsd_cutoff": 0})
add("ptm", {"positions": d, "cell": dcell, "pbc": [True, True, True], "structures": ["cubic-diamond", "hex-diamond"]})
add("ptm", {"positions": jit, "cell": cell, "pbc": [True, True, True], "types": [29 if (i % 4) else 79 for i in range(len(jit))]})
add("ptm", {"positions": jit, "structures": "bct"})
tw = [list(q) for q in big]
add("grain-segmentation", {"positions": bigj, "cell": bigcell, "pbc": [True, True, True], "min_grain_size": 5})
add("grain-segmentation", {"positions": b, "cell": bcell, "pbc": [True, True, True], "min_grain_size": 10, "adopt_orphans": False, "merge_threshold_deg": 2})
add("elastic-strain", {"positions": [[x * 1.01 for x in q] for q in bigj], "cell": [[x * 1.01 for x in r] for r in bigcell], "pbc": [True, True, True], "lattice_constant_A": 3.61})
add("elastic-strain", {"positions": h, "cell": hcell, "pbc": [True, True, True], "lattice_constant_A": 2.95, "lattice": "hcp", "c_over_a": 4.68 / 2.95})
add("elastic-strain", {"positions": jit, "lattice_constant_A": 3.61, "lattice": "ico"})
add("fcc-planar-faults", {"positions": bigj, "cell": bigcell, "pbc": [True, True, True]})
add("fcc-planar-faults", {"positions": h, "cell": hcell, "pbc": [True, True, True]})

# ---- symmetry
def crystal(symbols, positions, cell):
    return {"symbols": symbols, "positions": positions, "cell": cell}


cu, cucell = fcc(1)
add("symmetry", {"structure": crystal(["Cu"] * 4, cu, cucell)})
nacl = [[0, 0, 0], [0, 2.82, 2.82], [2.82, 0, 2.82], [2.82, 2.82, 0], [2.82, 0, 0], [0, 2.82, 0], [0, 0, 2.82], [2.82, 2.82, 2.82]]
add("symmetry", {"structure": crystal(["Na"] * 4 + ["Cl"] * 4, nacl, [[5.64, 0, 0], [0, 5.64, 0], [0, 0, 5.64]]), "output_cell": "primitive"})
hj, hjcell = hcp(2)
hj = [[x + 0.01 * R() for x in q] for q in hj]
add("symmetry", {"structure": crystal(["Mg"] * len(hj), hj, hjcell), "symprec_A": 0.05, "output_cell": "symmetrized"})
bb, bbcell = bcc(1)
add("symmetry", {"structure": crystal(["Fe"] * 2, bb, bbcell), "angle_tolerance_deg": 5})
sto = [[0, 0, 0], [1.9525, 1.9525, 1.9525], [1.9525, 1.9525, 0], [1.9525, 0, 1.9525], [0, 1.9525, 1.9525]]
add("symmetry", {"structure": crystal(["Sr", "Ti", "O", "O", "O"], sto, [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]])})
mono = [[0.1, 0.2, 0.3], [1.7, 2.1, 0.4], [0.9, 0.3, 2.2]]
add("symmetry", {"structure": crystal(["Si", "O", "O"], mono, [[4.1, 0, 0], [0.6, 3.7, 0], [0.3, -0.4, 5.2]])})
add("symmetry", {"structure": {"symbols": ["Cu"], "positions": [[0, 0, 0]]}})
add("symmetry", {"structure": crystal(["Cu"] * 4, cu, cucell), "output_cell": "niggli"})
add("symmetry", {"structure": crystal(["Cu"] * 4, cu, cucell), "symprec_A": 0})

# ---- diffraction
cu8 = crystal(["Cu"] * 4, cu, cucell)
add("powder-xrd", {"structure": cu8})
add("powder-xrd", {"structure": crystal(["Na"] * 4 + ["Cl"] * 4, nacl, [[5.64, 0, 0], [0, 5.64, 0], [0, 0, 5.64]]), "wavelength_A": 0.7107, "two_theta_range_deg": [5, 60], "fwhm_deg": 0.2, "debye_waller_A2": 0.5})
cu2, cu2cell = fcc(2)
add("powder-xrd", {"structure": crystal(["Cu"] * len(cu2), cu2, cu2cell), "two_theta_range_deg": [20, 100]})
add("powder-xrd", {"structure": crystal(["Mg"] * len(hj), hj, hjcell)})
add("powder-xrd", {"structure": crystal(["Si", "O", "O"], mono, [[4.1, 0, 0], [0.6, 3.7, 0], [0.3, -0.4, 5.2]]), "two_theta_range_deg": [10, 50]})
add("powder-xrd", {"structure": cu8, "two_theta_range_deg": [50, 40]})
add("powder-xrd", {"structure": crystal(["Xx"] * 4, cu, cucell)})
add("powder-xrd", {"structure": {"symbols": ["Cu"], "positions": [[0, 0, 0]]}})
add("electron-diffraction", {"structure": cu8})
add("electron-diffraction", {"structure": cu8, "zone_axis": [1, 1, 0], "voltage_kV": 300, "g_max_inv_A": 1.5})
add("electron-diffraction", {"structure": crystal(["Mg"] * len(hj), hj, hjcell), "zone_axis": [0, 0, 1]})
add("electron-diffraction", {"structure": crystal(["Sr", "Ti", "O", "O", "O"], sto, [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]]), "zone_axis": [1, 1, 1]})
add("electron-diffraction", {"structure": cu8, "zone_axis": [0, 0, 0]})
add("electron-diffraction", {"structure": cu8, "zone_axis": [1, 0]})
add("electron-diffraction", {"structure": cu8, "g_max_inv_A": 0.1})

# ---- reciprocal paths and DFT inputs
def cellxyz(symbols, fracs, cell):
    return crystal(symbols, [[sum(f[k] * cell[k][j] for k in range(3)) for j in range(3)] for f in fracs], cell)


si = cellxyz(["Si"] * 8, [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0], [.25, .25, .25], [.25, .75, .75], [.75, .25, .75], [.75, .75, .25]], [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]])
add("reciprocal-path", {"structure": si})
add("reciprocal-path", {"structure": cu8, "spacing_inv_A": 0.1})
add("reciprocal-path", {"structure": crystal(["Fe"] * 2, bb, bbcell)})
add("reciprocal-path", {"structure": crystal(["Mg"] * len(hj), hj, hjcell), "symprec_A": 0.05})
add("reciprocal-path", {"structure": crystal(["Sr", "Ti", "O", "O", "O"], sto, [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]])})
add("reciprocal-path", {"structure": crystal(["Si", "O", "O"], mono, [[4.1, 0, 0], [0.6, 3.7, 0], [0.3, -0.4, 5.2]]), "time_reversal": False})
# tetragonal body-centred (tI), orthorhombic (oP, oC, oF, oI), monoclinic (mP, mC), rhombohedral (hR), hexagonal noncentric (hP, augmented)
add("reciprocal-path", {"structure": cellxyz(["In", "In"], [[0, 0, 0], [.5, .5, .5]], [[3.25, 0, 0], [0, 3.25, 0], [0, 0, 4.95]])})
add("reciprocal-path", {"structure": cellxyz(["Sn", "Sn"], [[0, 0, 0], [.5, .5, .5]], [[5.8, 0, 0], [0, 5.8, 0], [0, 0, 3.2]])})
add("reciprocal-path", {"structure": cellxyz(["Ga", "As"], [[0, 0, 0], [.5, .4, .3]], [[3.1, 0, 0], [0, 4.2, 0], [0, 0, 5.3]])})
add("reciprocal-path", {"structure": cellxyz(["Ga", "Ga"], [[0, 0, 0], [.5, .5, 0]], [[3.1, 0, 0], [0, 4.2, 0], [0, 0, 5.3]])})
add("reciprocal-path", {"structure": cellxyz(["S"] * 4, [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], [[3.1, 0, 0], [0, 4.2, 0], [0, 0, 5.3]])})
add("reciprocal-path", {"structure": cellxyz(["U", "U"], [[0, 0, 0], [.5, .5, .5]], [[3.1, 0, 0], [0, 4.2, 0], [0, 0, 5.3]])})
add("reciprocal-path", {"structure": cellxyz(["Se", "Te"], [[0, 0, 0], [.3, .5, .2]], [[3.1, 0, 0], [0, 4.2, 0], [-0.8, 0, 5.3]])})
add("reciprocal-path", {"structure": cellxyz(["Se", "Se"], [[0, 0, 0], [.5, .5, 0]], [[3.1, 0, 0], [0, 4.2, 0], [-0.8, 0, 5.3]])})
rh = [[1.2, 0, 4.0], [-0.6, 1.0392304845413265, 4.0], [-0.6, -1.0392304845413265, 4.0]]
add("reciprocal-path", {"structure": cellxyz(["Bi"], [[0, 0, 0]], rh)})
add("reciprocal-path", {"structure": cellxyz(["Zn", "O"], [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.38]], [[3.25, 0, 0], [-1.625, 2.8145825622994254, 0], [0, 0, 5.2]]), "time_reversal": False})
add("reciprocal-path", {"structure": si, "spacing_inv_A": 0})
add("reciprocal-path", {"structure": {"symbols": ["Cu"], "positions": [[0, 0, 0]]}})
add("dft-inputs", {"structure": si})
add("dft-inputs", {"structure": crystal(["Sr", "Ti", "O", "O", "O"], sto, [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]]), "points_per_segment": 25})
add("dft-inputs", {"structure": cellxyz(["Zn", "O"], [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.38]], [[3.25, 0, 0], [-1.625, 2.8145825622994254, 0], [0, 0, 5.2]]), "time_reversal": False})
add("dft-inputs", {"structure": si, "points_per_segment": 1})
add("dft-inputs", {"structure": si, "points_per_segment": 2000})

# ---- dislocation lines
def screw(points, centre, b):
    out = []
    for (x, y, z) in points:
        theta = math.atan2(y - centre[1], x - centre[0])
        out.append([x, y, z + b * theta / (2 * math.pi)])
    return out


fp, fpcell = fcc(6)
fcentre = [fpcell[0][0] / 2 + 0.37, fpcell[1][1] / 2 + 0.41]
add("dislocation-lines", {"reference": {"symbols": ["Cu"] * len(fp), "positions": fp}, "structure": {"symbols": ["Cu"] * len(fp), "positions": screw(fp, fcentre, 3.61)}})
add("dislocation-lines", {"reference": {"symbols": ["Cu"] * len(fp), "positions": fp}, "structure": {"symbols": ["Cu"] * len(fp), "positions": screw(fp, fcentre, 3.61)}, "cluster_cutoff_A": 3.5, "min_atoms": 2, "circuit_radius_A": 6.0, "cutoff_A": 3.0})
bp, bpcell = bcc(7)
bcentre = [bpcell[0][0] / 2 + 0.31, bpcell[1][1] / 2 + 0.23]
add("dislocation-lines", {"reference": {"symbols": ["Fe"] * len(bp), "positions": bp}, "structure": {"symbols": ["Fe"] * len(bp), "positions": screw(bp, bcentre, 2.87)}})
col = []
for q in fp:
    if (q[0] - fcentre[0]) ** 2 + (q[1] - fcentre[1]) ** 2 < 9:
        col.append([x + 0.6 * R() for x in q])
    else:
        col.append(list(q))
add("dislocation-lines", {"reference": crystal(["Cu"] * len(fp), fp, fpcell), "structure": crystal(["Cu"] * len(fp), col, fpcell)})
add("dislocation-lines", {"reference": crystal(["Cu"] * len(fp), fp, fpcell), "structure": crystal(["Cu"] * len(fp), fp, fpcell)})
add("dislocation-lines", {"reference": crystal(["Cu"] * 4, cu, cucell), "structure": crystal(["Cu"] * len(fp), fp, fpcell)})
gas = [[12 * random.random(), 12 * random.random(), 12 * random.random()] for _ in range(40)]
add("dislocation-lines", {"reference": {"symbols": ["Ar"] * 40, "positions": gas}, "structure": {"symbols": ["Ar"] * 40, "positions": gas}})
add("dislocation-lines", {"reference": {"symbols": ["Ar"] * 2, "positions": [[0, 0, 0], [20, 0, 0]]}, "structure": {"symbols": ["Ar"] * 2, "positions": [[0, 0, 0], [20, 0, 0]]}})
add("dislocation-lines", {"reference": crystal(["Cu"] * len(fp), fp, fpcell), "structure": crystal(["Cu"] * len(fp), col, fpcell), "min_atoms": 0})

# ---- trajectory structure
add("trajectory-structure", {"trajectory_file": T})
add("trajectory-structure", {"trajectory_file": T, "first_frame": 1, "every": 2, "r_max_A": 3.5, "bins": 35, "cutoff_A": 2.9, "angle_bins": 36})
add("trajectory-structure", {"trajectory_file": T, "cutoff_A": 4.0, "r_max_A": 3.0, "bins": 10})
add("trajectory-structure", {"trajectory_file": "rust/tests/fixtures/tools/traj_open.xyz"})
add("trajectory-structure", {"trajectory_file": T, "first_frame": 6})
add("trajectory-structure", {"trajectory_file": T, "bins": 5})
add("trajectory-structure", {})

# ---- LAMMPS export
tri = crystal(["Si", "O", "O"], mono, [[4.1, 0, 0], [0.6, 3.7, 0], [0.3, -0.4, 5.2]])
skew = crystal(["Cu"] * 4, cu, [[3.61, 0, 0], [3.0, 3.61, 0], [-2.5, 2.4, 3.61]])
add("lammps-export", {"structure": cu8})
add("lammps-export", {"structure": tri, "task": "nvt", "temperature_K": 500, "timestep_fs": 0.5, "steps": 2000, "sample_interval": 20, "seed": 7, "calculator": {"potential": "lennard-jones", "sigma": 2.3, "epsilon": 0.0104}})
add("lammps-export", {"structure": skew, "task": "npt", "pressure_GPa": 1.5, "calculator": {"potential": "EAM", "file": "pots/Cu_u3.eam"}})
add("lammps-export", {"structure": skew, "calculator": {"potential": "eam/alloy", "file": "C:/x/CuAg.eam.alloy"}, "relax_cell": True, "fmax": 0.001})
add("lammps-export", {"structure": crystal(["Cu"] * 4, cu, cucell), "calculator": {"potential": "Finnis-Sinclair", "file": "Fe.fs", "format": "auto"}})
add("lammps-export", {"structure": si, "calculator": {"potential": "Tersoff"}})
add("lammps-export", {"structure": si, "calculator": {"potential": "stillinger-weber", "file": "my.sw"}})
mgo = crystal(["Mg"] * 4 + ["O"] * 4, nacl, [[4.21, 0, 0], [0, 4.21, 0], [0, 0, 4.21]])
add("lammps-export", {"structure": mgo, "calculator": {"potential": "Buckingham", "cutoff": 12, "charges": {"Mg": 2.0, "O": -2.0}, "pairs": [{"elements": ["Mg", "O"], "A": 1428.5, "rho": 0.2945}, {"elements": ["O", "O"], "A": 22764.0, "rho": 0.149, "C": 27.88}], "ewald_accuracy": 1e-5}})
add("lammps-export", {"structure": mgo, "calculator": {"potential": "buck", "pairs": [{"elements": ["O", "Mg"], "A": 1428.5, "rho": 0.2945}]}})
mol = {"symbols": ["O", "H", "H"], "positions": [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]}
mol2 = {"symbols": ["O", "H", "H"], "positions": [[0.1, 0, 0], [1.06, 0.1, 0], [-0.24, 0.93, 0.2]]}
add("lammps-export", {"structure": mol, "task": "neb", "final": mol2, "images": 5})
add("lammps-export", {"structure": skew, "task": "neb", "final": crystal(["Cu"] * 4, [[x + 0.1 for x in q] for q in cu], [[3.61, 0, 0], [3.0, 3.61, 0], [-2.5, 2.4, 3.61]])})
add("lammps-export", {"structure": mol, "task": "npt"})
add("lammps-export", {"structure": mol, "relax_cell": True})
add("lammps-export", {"structure": mol, "task": "neb"})
add("lammps-export", {"structure": mol, "task": "neb", "final": cu8})
add("lammps-export", {"structure": mol, "task": "md"})
add("lammps-export", {"structure": mol, "calculator": {"potential": "buck", "charges": {"O": -0.8, "H": 0.4}}})
add("lammps-export", {"structure": mgo, "calculator": {"potential": "buck", "charges": {"Mg": 2.0}}})
add("lammps-export", {"structure": crystal(["Cu"] * 4, cu, [[0, 3.61, 0], [3.61, 0, 0], [0, 0, 3.61]])})

# ---- VASP electronic structure (synthetic files)
V = "rust/tests/fixtures/tools/vasp/"
os.makedirs(V, exist_ok=True)


def band_energy(b, k, nk, spin):
    t = k / (nk - 1)
    base = [-6.0, -4.5, -2.0, 1.5, 3.0, 4.2][b]
    return base + (0.8 if b < 3 else -0.6) * math.cos(math.pi * t) + 0.05 * spin + 0.01 * R()


def write_eigenval(path, kpts, spins, nb=6):
    with open(path, "w") as f:
        f.write("    2    2    1    %d\n  0.1E+02  0.1E-08  0.1E-08  0.3E-08  0.5E-15\n  1.0E-004\n  CAR\n GaAs\n" % spins)
        f.write("      8     %d      %d\n" % (len(kpts), nb))
        for i, k in enumerate(kpts):
            f.write("\n  %.7E  %.7E  %.7E  %.7E\n" % (k[0], k[1], k[2], 1.0 / len(kpts)))
            for b in range(nb):
                es = " ".join("%12.6f" % band_energy(b, i % 5, 5, s) for s in range(spins))
                f.write("%5d %s %8.4f\n" % (b + 1, es, 1.0 if b < 3 else 0.0))


def path(points, n):
    out = []
    for a, b in points:
        for i in range(n):
            t = i / (n - 1)
            out.append([a[j] + (b[j] - a[j]) * t for j in range(3)])
    return out


G, X, W, K, L = [0, 0, 0], [0.5, 0, 0.5], [0.5, 0.25, 0.75], [0.375, 0.375, 0.75], [0.5, 0.5, 0.5]
kp = path([(G, X), (X, W)], 5)
write_eigenval(V + "EIGENVAL", kp, 1)
write_eigenval(V + "EIGENVAL.spin", kp, 2)
kpb = path([(G, X), (K, L)], 5)
write_eigenval(V + "EIGENVAL.break", kpb, 1)
with open(V + "KPOINTS", "w") as f:
    f.write("k-path\n 5\nLine-mode\nReciprocal\n0 0 0 ! \\Gamma\n0.5 0 0.5 ! X\n\n0.5 0 0.5 ! X\n0.5 0.25 0.75 ! W\n")
with open(V + "KPOINTS.break", "w") as f:
    f.write("k-path\n5\n  line\nrec\n0 0 0 GAMMA\n0.5 0 0.5 X\n\n0.375 0.375 0.75 K\n0.5 0.5 0.5 ! L\n")
with open(V + "KPOINTS.wrong", "w") as f:
    f.write("k-path\n4\nLine-mode\nReciprocal\n0 0 0 ! G\n0.5 0 0.5 ! X\n")
with open(V + "KPOINTS.mesh", "w") as f:
    f.write("mesh\n0\nGamma\n4 4 4\n")


def write_doscar(path, spins, partial):
    with open(path, "w") as f:
        f.write("    2    2    %d    0\n  0.1E+02  0.1E-08\n  1.0E-004\n  CAR\n GaAs\n" % partial)
        head = "     10.00000000    -10.00000000     21      0.45000000      1.00000000\n"
        f.write(head)
        for i in range(21):
            e = -10 + i
            vals = [max(0.0, 2.0 - abs(abs(e) - 4) / 2) + 0.1 * s for s in range(spins)]
            ints = [0.5 * i + s for s in range(spins)]
            f.write("%10.3f %s %s\n" % (e, " ".join("%.4E" % v for v in vals), " ".join("%.4E" % v for v in ints)))
        if partial:
            for a in range(2):
                f.write(head)
                for i in range(21):
                    e = -10 + i
                    cols = [0.1 * (a + 1) * (o + 1) * (1 + 0.01 * i) for o in range(3 * spins)]
                    f.write("%10.3f %s\n" % (e, " ".join("%.4E" % c for c in cols)))


write_doscar(V + "DOSCAR", 1, 1)
write_doscar(V + "DOSCAR.spin", 2, 1)
write_doscar(V + "DOSCAR.total", 1, 0)


def write_procar(path, nk, nb, phases=False):
    with open(path, "w") as f:
        f.write("PROCAR lm decomposed\n# of k-points:   %d         # of bands:   %d         # of ions:   2\n" % (nk, nb))
        for k in range(nk):
            f.write("\n k-point     %d :    0.0 0.0 0.0     weight = 0.1\n\n" % (k + 1))
            for b in range(nb):
                f.write("band     %d # energy   -5.0 # occ.  1.0\n\nion      s      py     pz     px    dxy    tot\n" % (b + 1))
                rows = []
                for ion in range(2):
                    w = [round(0.05 * ((ion + b + k + o) % 5), 3) for o in range(5)]
                    rows.append(w)
                    f.write("    %d %s %.3f\n" % (ion + 1, " ".join("%.3f" % x for x in w), sum(w)))
                f.write("tot    %s %.3f\n" % (" ".join("%.3f" % (rows[0][o] + rows[1][o]) for o in range(5)), sum(map(sum, rows))))
                if phases:
                    f.write("ion      s      py     pz     px    dxy    tot\n    1  9 9 9 9 9 9\n    2  9 9 9 9 9 9\ncharge 1 1 1 1 1 1\n")
                f.write("\n")


write_procar(V + "PROCAR", 10, 6, phases=True)
write_procar(V + "PROCAR.short", 4, 6)
with open(V + "EIGENVAL.bad", "w") as f:
    f.write("    2    2    1    3\n")
gaas = cellxyz(["Ga", "As"], [[0, 0, 0], [.25, .25, .25]], [[0, 2.83, 2.83], [2.83, 0, 2.83], [2.83, 2.83, 0]])
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "doscar_file": V + "DOSCAR", "kpoints_file": V + "KPOINTS", "structure": gaas})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL.spin", "doscar_file": V + "DOSCAR.spin", "structure": gaas, "fermi_eV": 0.3})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL.break", "kpoints_file": V + "KPOINTS.break", "fermi_eV": 0.1})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "kpoints_file": V + "KPOINTS.wrong", "doscar_file": V + "DOSCAR.total"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "PROCAR", "structure": gaas, "fermi_eV": 0.0})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "PROCAR", "projection_element": "As", "projection_orbital": "p"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "PROCAR", "structure": gaas, "projection_orbital": "dxy"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "PROCAR"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "PROCAR.short", "structure": gaas})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "kpoints_file": V + "KPOINTS.mesh"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL.bad"})
add("vasp-electronic", {"eigenval_file": V + "missing"})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "structure": {"symbols": ["Cu"], "positions": [[0, 0, 0]]}})
add("vasp-electronic", {"eigenval_file": V + "EIGENVAL", "procar_file": V + "KPOINTS"})
add("vasp-electronic", {"eigenval_file": V + "KPOINTS"})

# ---- potentials and simulations (fixtures written here)
P = "rust/tests/fixtures/tools/potentials/"
os.makedirs(P, exist_ok=True)


def table_lines(values):
    return "".join(" ".join("%.10e" % v for v in values[i:i + 5]) + "\n" for i in range(0, len(values), 5))


def eam_functions(scale_f, r0, nrho=200, drho=0.05, nr=300, dr=0.02, cut=5.5):
    F = [-scale_f * math.sqrt(i * drho) + 0.02 * (i * drho) ** 2 for i in range(nrho)]
    rho = [0.0 if i * dr >= cut else math.exp(-1.8 * (i * dr - r0)) * (1 - (i * dr / cut) ** 4) ** 2 for i in range(nr)]
    rphi = [0.0 if i * dr >= cut else (i * dr) * 0.6 * (math.exp(-2.4 * (i * dr - r0)) - 2 * math.exp(-1.2 * (i * dr - r0))) * (1 - (i * dr / cut) ** 4) ** 2 for i in range(nr)]
    return F, rho, rphi, (nrho, drho, nr, dr, cut)


F, rho, rphi, (nrho, drho, nr, dr, cut) = eam_functions(1.3, 2.55)
with open(P + "Cu_test.eam", "w") as f:
    f.write("synthetic Cu funcfl\n29 63.546 3.615 FCC\n%d %r %d %r %r\n" % (nrho, drho, nr, dr, cut))
    z = [math.sqrt(max(0.0, -v) / (27.2 * 0.529)) if i else 0.0 for i, v in enumerate(rphi)]
    f.write(table_lines(F) + table_lines(z) + table_lines(rho))
Fa, rhoa, rphia, _ = eam_functions(1.3, 2.55)
Fb, rhob, rphib, _ = eam_functions(1.1, 2.88)
with open(P + "CuAg_test.eam.alloy", "w") as f:
    f.write("synthetic\nCu Ag setfl\n\n2 Cu Ag\n%d %r %d %r %r\n" % (nrho, drho, nr, dr, cut))
    f.write("29 63.546 3.615 fcc\n" + table_lines(Fa) + table_lines(rhoa))
    f.write("47 107.8682 4.09 fcc\n" + table_lines(Fb) + table_lines(rhob))
    f.write(table_lines(rphia) + table_lines([0.5 * (x + y) for x, y in zip(rphia, rphib)]) + table_lines(rphib))
with open(P + "FeCr_test.eam.fs", "w") as f:
    f.write("synthetic\nFe Cr fs\n\n2 Fe Cr\n%d %r %d %r %r\n" % (nrho, drho, nr, dr, cut))
    f.write("26 55.845 2.86 bcc\n" + table_lines(Fa) + table_lines(rhoa) + table_lines([0.9 * x for x in rhoa]))
    f.write("24 51.996 2.88 bcc\n" + table_lines(Fb) + table_lines([0.8 * x for x in rhob]) + table_lines(rhob))
    f.write(table_lines(rphia) + table_lines([0.5 * (x + y) for x, y in zip(rphia, rphib)]) + table_lines(rphib))
with open(P + "Si_test.tersoff", "w") as f:
    f.write("# Si test\nSi Si Si 3.0 1.0 0.0 1.0039e5 16.217 -0.59825 0.78734 1.1e-6 1.7322 471.18 2.85 0.15 2.4799 1830.8\n")
with open(P + "SiC_test.tersoff", "w") as f:
    for a in ["Si", "C"]:
        for b in ["Si", "C"]:
            for c in ["Si", "C"]:
                f.write("%s %s %s 3.0 1.0 0.0 1.0039e5 16.217 -0.59825 0.78734 1.1e-6 1.7322 471.18 2.85 0.15 2.4799 1830.8 # x\n" % (a, b, c))
with open(P + "Si_test.sw", "w") as f:
    f.write("Si Si Si 2.1683 2.0951 1.80 21.0 1.20 -0.333333333333 7.049556277 0.6022245584 4.0 0.0 0.0\n")
with open(P + "bad.tersoff", "w") as f:
    f.write("Si Si Si 3.0 1.0\n")
with open(P + "badnum.sw", "w") as f:
    f.write("Si Si Si 2.1683 2.0951 1.80 21.0 1.20 -0.33x 7.049556277 0.6022245584 4.0 0.0 0.0\n")

EMT = {"potential": "EMT"}
cuj = [[x + 0.05 * R() for x in q] for q in cu]
cu4j = crystal(["Cu"] * 4, cuj, cucell)
add("relax", {"structure": cu4j, "calculator": EMT, "steps": 40})
add("relax", {"structure": cu4j, "calculator": EMT, "relax_cell": True, "pressure_GPa": 2.0, "steps": 60, "fmax": 0.005})
argon = {"symbols": ["Ar"] * 6, "positions": [[0, 0, 0], [3.8, 0.1, 0], [0, 3.7, 0.2], [0.1, 0, 3.9], [3.6, 3.7, 0.3], [-3.5, 0.2, 0.1]]}
add("relax", {"structure": argon, "calculator": {"potential": "LennardJones", "epsilon": 0.0104, "sigma": 3.4, "cutoff": 8.5}, "steps": 80})
add("relax", {"structure": argon, "calculator": {"module": "ase.calculators.lj", "attribute": "LennardJones", "kwargs": {"epsilon": 0.0104, "sigma": 3.4, "rc": 7.0}}, "steps": 10})
sij = [[x + 0.04 * R() for x in q] for q in si["positions"]]
si8 = crystal(["Si"] * 8, sij, si["cell"])
add("relax", {"structure": si8, "calculator": {"potential": "Tersoff"}, "steps": 25})
add("relax", {"structure": si8, "calculator": {"potential": "Tersoff", "file": "rust/tests/fixtures/tools/potentials/Si_test.tersoff"}, "steps": 5, "relax_cell": True})
add("relax", {"structure": si8, "calculator": {"potential": "StillingerWeber"}, "steps": 25})
add("relax", {"structure": si8, "calculator": {"potential": "sw", "file": "rust/tests/fixtures/tools/potentials/Si_test.sw"}, "steps": 5})
sic = crystal(["Si", "C"] * 4, sij, si["cell"])
add("relax", {"structure": sic, "calculator": {"potential": "Tersoff", "file": "rust/tests/fixtures/tools/potentials/SiC_test.tersoff"}, "steps": 5})
mgoj = crystal(["Mg"] * 4 + ["O"] * 4, [[x + 0.03 * R() for x in q] for q in nacl], [[4.21, 0, 0], [0, 4.21, 0], [0, 0, 4.21]])
buck = {"potential": "Buckingham", "cutoff": 6.0, "charges": {"Mg": 2.0, "O": -2.0}, "pairs": [{"elements": ["Mg", "O"], "A": 1428.5, "rho": 0.2945}, {"elements": ["O", "O"], "A": 22764.0, "rho": 0.149, "C": 27.88}], "ewald_accuracy": 1e-5}
add("relax", {"structure": mgoj, "calculator": buck, "steps": 8})
add("relax", {"structure": {"symbols": ["Mg", "O", "Mg", "O"], "positions": [[0, 0, 0], [2.0, 0.1, 0], [2.1, 2.0, 0.1], [0, 2.1, 0]]}, "calculator": buck, "steps": 8})
add("relax", {"structure": mgoj, "calculator": {"potential": "buck", "pairs": [{"elements": ["O", "Mg"], "A": 1428.5, "rho": 0.2945, "C": 2.0}], "cutoff": 7.0}, "steps": 4, "relax_cell": True})
cu8b, cu8bcell = fcc(2, a=3.615)
cu32j = crystal(["Cu"] * len(cu8b), [[x + 0.04 * R() for x in q] for q in cu8b], cu8bcell)
add("relax", {"structure": cu32j, "calculator": {"potential": "EAM", "file": "rust/tests/fixtures/tools/potentials/Cu_test.eam"}, "steps": 15})
cuag = crystal(["Cu", "Ag"] * 16, [[x + 0.04 * R() for x in q] for q in cu8b], cu8bcell)
add("relax", {"structure": cuag, "calculator": {"potential": "eam/alloy", "file": "rust/tests/fixtures/tools/potentials/CuAg_test.eam.alloy"}, "steps": 10, "relax_cell": True})
bfe, bfecell = bcc(2)
fecr = crystal(["Fe", "Cr"] * 8, [[x + 0.04 * R() for x in q] for q in bfe], bfecell)
add("relax", {"structure": fecr, "calculator": {"potential": "Finnis-Sinclair", "file": "rust/tests/fixtures/tools/potentials/FeCr_test.eam.fs"}, "steps": 10})
add("relax", {"structure": cu4j, "calculator": {"potential": "EAM", "file": "rust/tests/fixtures/tools/potentials/Cu_test.eam", "format": "funcfl"}, "steps": 0})
# potential errors
add("relax", {"structure": cu4j, "calculator": {"potential": "EAM"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "EAM", "file": "missing.eam"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "EAM", "file": "rust/tests/fixtures/tools/potentials/Cu_test.eam", "format": "xyz"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "EAM", "file": "rust/tests/fixtures/tools/potentials/CuAg_test.eam.alloy"}, "steps": 2})
add("relax", {"structure": si8, "calculator": {"potential": "Tersoff", "file": "rust/tests/fixtures/tools/potentials/bad.tersoff"}})
add("relax", {"structure": si8, "calculator": {"potential": "sw", "file": "rust/tests/fixtures/tools/potentials/badnum.sw"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "Tersoff"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "Morse"}})
add("relax", {"structure": cu4j, "calculator": {"potential": "LJ", "sigma": -1}})
add("relax", {"structure": cu4j, "calculator": {"potential": "buck"}})
add("relax", {"structure": cu4j, "calculator": "EMT"})
add("relax", {"structure": cu4j, "calculator": {"model": "EMT"}})
add("relax", {"structure": crystal(["Fe"] * 4, cu, cucell), "calculator": EMT})
add("relax", {"structure": mgoj, "calculator": {"potential": "buck", "charges": {"Mg": 2.0, "O": -1.0}}})
add("relax", {"structure": argon, "calculator": EMT, "relax_cell": True})
# dynamics
add("nvt", {"structure": cu32j, "calculator": EMT, "steps": 20, "timestep_fs": 2.0, "temperature_K": 500, "seed": 3, "sample_interval": 5})
add("nvt", {"structure": cu32j, "calculator": EMT, "steps": 12, "seed": 11, "sample_interval": 5, "production_steps": 13, "thermostat_fs": 50})
add("nvt", {"structure": argon, "calculator": {"potential": "LJ", "epsilon": 0.0104, "sigma": 3.4}, "steps": 30, "sample_interval": 3, "equilibration_fs": 10})
add("npt", {"structure": cu32j, "calculator": EMT, "steps": 20, "seed": 5, "sample_interval": 4, "pressure_GPa": 1.0, "barostat_fs": 200})
add("npt", {"structure": cu32j, "calculator": EMT, "steps": 10, "seed": 7, "sample_interval": 5, "production_steps": 10, "equilibration_fs": 4})
add("npt", {"structure": argon, "calculator": EMT, "steps": 5})
add("nvt", {"structure": {"symbols": ["Cu"], "positions": [[0, 0, 0]]}, "calculator": EMT})
add("nvt", {"structure": cu32j, "calculator": EMT, "equilibration_fs": -1})
# NEB: a vacancy hop in fcc Cu (atom 1 moves into the empty site 0)
vac = [list(q) for q in cu8b[1:]]
hop = [list(q) for q in vac]
hop[0] = list(cu8b[0])
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "final": crystal(["Cu"] * len(hop), hop, cu8bcell), "calculator_factory": EMT, "images": 5, "steps": 15, "mic": True})
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "final": crystal(["Cu"] * len(hop), hop, cu8bcell), "calculator_factory": EMT, "images": 0, "image_spacing_A": 0.6, "steps": 6, "climb": False, "relax_endpoints": True, "endpoint_fmax": 0.05})
chain = []
for i in range(5):
    t = i / 4
    frame = [("Cu", [vac[a][k] + t * (hop[a][k] - vac[a][k]) + (0.05 * math.sin(math.pi * t) if a == 0 and k == 2 else 0.0) for k in range(3)]) for a in range(len(vac))]
    chain.append(frame)
write_xyz(P + "neb_restart.xyz", chain, cu8bcell)
add("neb", {"restart_images": {"file": P + "neb_restart.xyz"}, "calculator_factory": EMT, "steps": 5})
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "calculator_factory": EMT})
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "final": cu8, "calculator_factory": EMT})
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "final": crystal(["Cu"] * len(hop), hop, fpcell), "calculator_factory": EMT})
add("neb", {"initial": crystal(["Cu"] * len(vac), vac, cu8bcell), "final": crystal(["Cu"] * len(hop), hop, cu8bcell), "calculator_factory": EMT, "images": 2})
# phonons
add("phonons", {"structure": cu8, "calculator": EMT, "supercell_min_A": 7.0, "mesh": [3, 3, 3], "spacing_inv_A": 0.2})
add("phonons", {"structure": cu8, "calculator": EMT, "symmetry_path": False, "q_points_fractional": [[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]], "supercell": [2, 2, 2], "mesh": [2, 2, 2], "temperatures_K": [100, 300]})
add("phonons", {"structure": si8, "calculator": {"potential": "Tersoff"}, "symmetry_path": False, "supercell": [1, 1, 1], "mesh": [2, 2, 2]})
add("phonons", {"structure": argon, "calculator": EMT})
add("phonons", {"structure": cu8, "calculator": EMT, "supercell": [2, 2]})
# workflows
add("elastic-constants", {"structure": crystal(["Cu"] * 4, cu, cucell), "calculator": EMT, "relax_ions": False, "strain_steps": 1})
add("elastic-constants", {"structure": cu4j, "calculator": EMT, "strain_steps": 1, "max_strain": 0.005, "fmax": 0.01, "steps": 20})
add("elastic-constants", {"structure": cu4j, "calculator": EMT, "max_strain": 0.1})
add("eos-scan", {"structure": crystal(["Cu"] * 4, cu, cucell), "calculator": EMT, "points": 7, "relax_ions": False})
add("eos-scan", {"structure": cu4j, "calculator": EMT, "points": 5, "volume_range": 0.05, "steps": 15, "fmax": 0.01})
add("eos-scan", {"structure": cu4j, "calculator": EMT, "volume_range": 0.6})
add("eos-scan", {"structure": argon, "calculator": EMT})
bulk32 = crystal(["Cu"] * len(cu8b), cu8b, cu8bcell)
add("formation-energy", {"structure": crystal(["Cu"] * len(vac), vac, cu8bcell), "bulk": bulk32, "calculator": EMT, "steps": 20})
add("formation-energy", {"structure": crystal(["Cu"] * 31 + ["Ag"], [list(q) for q in cu8b], cu8bcell), "bulk": bulk32, "calculator": EMT, "relax": False, "chemical_potentials_eV": {"Ag": -2.9}})
add("formation-energy", {"structure": crystal(["Cu"] * 31 + ["Ag"], [list(q) for q in cu8b], cu8bcell), "bulk": bulk32, "calculator": EMT, "relax": False})
add("formation-energy", {"structure": crystal(["Cu"] * len(vac), vac, cu8bcell), "bulk": bulk32, "calculator": EMT, "relax_cell": True, "steps": 10})
slab = [list(q) for q in cu8b]
for q in slab:
    if q[2] >= 3.615:
        q[0] += 0.9
add("planar-defect-energy", {"structure": crystal(["Cu"] * len(slab), slab, cu8bcell), "bulk": crystal(["Cu"] * 4, cu, [[3.615, 0, 0], [0, 3.615, 0], [0, 0, 3.615]]), "calculator": EMT, "steps": 10})
add("planar-defect-energy", {"structure": crystal(["Cu"] * len(slab), slab, cu8bcell), "bulk": crystal(["Cu"] * 4, cu, [[3.615, 0, 0], [0, 3.615, 0], [0, 0, 3.615]]), "calculator": EMT, "relax": False, "interfaces": 1, "normal_axis": 0})
add("planar-defect-energy", {"structure": crystal(["Cu"] * 3, cu[:3], cucell), "bulk": crystal(["Cu"] * 4, cu, cucell), "calculator": EMT})
add("planar-defect-energy", {"structure": crystal(["Cu"] * 4, cu, cucell), "bulk": crystal(["Cu"] * 4, cu, cucell), "calculator": EMT, "normal_axis": 3})

# ---- dispatcher errors
add("msd", {"positions": walk, "timestep": 2.0})
add("msd", {"positions": walk})
add("msdd", {})
add("band-gap", [])

with open("rust/tests/fixtures/tools/requests.json", "w") as f:
    json.dump(cases, f)
print(len(cases), "requests")
