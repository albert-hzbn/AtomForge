"""Band edges, local effective masses and vacuum-referenced work functions."""

from ._validation import array, points, positive


def band_gap(energies_eV, fermi_eV, *, tolerance_eV=1e-6):
    """Sampled band edges and smallest same-spin direct gap, relative to E_F.

    Shape is (k,band) or (spin,k,band). Detects a metal when a band crosses E_F
    on the supplied mesh, or a sampled state lies within tolerance of E_F.
    A sparse mesh can miss crossings/extrema. Smearing occupations are not
    interpreted here, and the direct gap is not an optical excitation energy.
    """
    import numpy as np
    energy = array(energies_eV, "energies_eV")
    fermi = float(array([fermi_eV], "fermi_eV", 1)[0])
    tol = positive(tolerance_eV, "tolerance_eV")
    if energy.ndim == 2:
        energy = energy[None, :, :]
    if energy.ndim != 3:
        raise ValueError("energies_eV must have shape (k,band) or (spin,k,band)")
    below, above = energy < fermi - tol, energy > fermi + tol
    if not below.any() or not above.any():
        raise ValueError("Include both occupied and unoccupied bands")
    crossing = np.any(below.any(axis=1) & above.any(axis=1))
    metal = bool(crossing or np.any(abs(energy - fermi) <= tol))
    occupied = np.where(below, energy, -np.inf)
    empty = np.where(above, energy, np.inf)
    valence = float(occupied.max()); conduction = float(empty.min())
    vertical = empty.min(axis=2) - occupied.max(axis=2)
    finite = vertical[np.isfinite(vertical)]
    return {"metal_on_sampled_mesh": metal, "gap_eV": 0.0 if metal else conduction - valence,
            "direct_gap_eV": 0.0 if metal else (float(finite.min()) if finite.size else None),
            "vbm_eV": valence, "cbm_eV": conduction,
            "vbm_indices_spin_k_band": np.argwhere(below & (abs(energy - valence) <= tol)),
            "cbm_indices_spin_k_band": np.argwhere(above & (abs(energy - conduction) <= tol))}


def effective_mass(kpoints_inv_A, energies_eV, *, center_inv_A):
    """Fit E=E0+g.dk+1/2 dk.H.dk; return signed principal masses in m_e.

    k is Cartesian wavevector in radians/A, not fractional reciprocal position
    or path length. Use a small 3D neighborhood of a band extremum, with one
    consistently tracked band. Negative masses describe valence curvature;
    positive hole masses are their negatives. This is a local quadratic model.
    """
    import numpy as np
    k = points(kpoints_inv_A, "kpoints_inv_A")
    energy = array(energies_eV, "energies_eV", 1)
    center = array(center_inv_A, "center_inv_A", 1)
    if center.shape != (3,) or energy.shape != (len(k),) or len(k) < 10:
        raise ValueError("At least ten 3D k points, matching energies and one center are required")
    x, y, z = (k - center).T
    design = np.column_stack((np.ones(len(k)), x, y, z, x*x/2, y*y/2, z*z/2, x*y, x*z, y*z))
    # Scaling protects tiny neighborhoods against misleading rank decisions.
    scale = np.linalg.norm(design, axis=0)
    if np.any(scale == 0) or np.linalg.matrix_rank(design / scale) != 10:
        raise ValueError("k points must span a full 3D quadratic fit; a line path is insufficient")
    coefficients = np.linalg.lstsq(design / scale, energy, rcond=None)[0] / scale
    hessian = np.array([[coefficients[4], coefficients[7], coefficients[8]],
                        [coefficients[7], coefficients[5], coefficients[9]],
                        [coefficients[8], coefficients[9], coefficients[6]]])
    curvature, axes = np.linalg.eigh(hessian)
    if np.any(abs(curvature) < 1e-10):
        raise ValueError("Singular band curvature: finite effective mass cannot be determined")
    # hbar^2 / m_e in eV A^2 (CODATA conversion).
    constant = 7.619964231073853
    return {"principal_masses_m_e": constant / curvature, "principal_axes_columns": axes,
            "mass_tensor_m_e": constant * np.linalg.inv(hessian), "hessian_eV_A2": hessian,
            "gradient_eV_A": coefficients[1:4], "energy_at_center_eV": float(coefficients[0]),
            "rms_fit_error_eV": float(np.sqrt(np.mean((design @ coefficients - energy) ** 2)))}


def work_function(distance_A, potential_eV, fermi_eV, *, vacuum_range_A, max_slope_eV_per_A=.01):
    """Phi=<V_vac>-E_F for a user-selected vacuum interval of a planar potential.

    Supply an electron potential-energy profile with a consistent reference,
    e.g. the appropriate VASP LOCPOT channel. A linear-fit slope and spread flag
    a poorly converged/nonflat vacuum. Calculate the two slab faces separately.
    This does not insert dipole corrections or detect vacuum automatically.
    """
    import numpy as np
    position = array(distance_A, "distance_A", 1)
    potential = array(potential_eV, "potential_eV", 1)
    region = array(vacuum_range_A, "vacuum_range_A", 1)
    fermi = float(array([fermi_eV], "fermi_eV", 1)[0])
    limit = positive(max_slope_eV_per_A, "max_slope_eV_per_A")
    if position.shape != potential.shape or np.any(np.diff(position) <= 0):
        raise ValueError("Distances must increase and match the potential")
    if region.shape != (2,) or region[0] >= region[1] or region[0] < position[0] or region[1] > position[-1]:
        raise ValueError("Vacuum interval must lie within the supplied profile")
    chosen = (position >= region[0]) & (position <= region[1])
    if chosen.sum() < 3:
        raise ValueError("At least three vacuum samples are required")
    slope, _ = np.polyfit(position[chosen], potential[chosen], 1)
    vacuum = float(np.mean(potential[chosen]))
    return {"work_function_eV": vacuum - fermi, "vacuum_level_eV": vacuum,
            "vacuum_std_eV": float(np.std(potential[chosen])), "slope_eV_per_A": float(slope),
            "flat_vacuum": bool(abs(slope) <= limit), "samples": int(chosen.sum())}
