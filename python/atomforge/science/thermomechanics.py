"""Equation of state, linear elastic response and harmonic thermal properties."""

from ._validation import array, positive


def equation_of_state(volumes_A3, energies_eV):
    """Third-order Birch-Murnaghan E(V) fit; energy and volume must share a cell basis.

    Requires at least five distinct positive volumes bracketing the fitted
    minimum. Returns bulk modulus and pressure derivative, residuals and fit
    covariance (a numerical diagnostic, not model/convergence uncertainty).
    """
    import numpy as np
    from scipy.optimize import curve_fit
    from ase.eos import birchmurnaghan
    volumes = array(volumes_A3, "volumes_A3", 1)
    energies = array(energies_eV, "energies_eV", 1)
    if volumes.shape != energies.shape or len(volumes) < 5 or np.any(volumes <= 0) or len(np.unique(volumes)) != len(volumes):
        raise ValueError("Provide at least five distinct positive volumes and matching energies")
    quadratic = np.polyfit(volumes, energies, 2)
    if quadratic[0] <= 0:
        raise ValueError("Energy-volume samples must contain a convex minimum")
    v0 = -quadratic[1] / (2 * quadratic[0])
    if not volumes.min() < v0 < volumes.max():
        raise ValueError("Volume samples must bracket the equilibrium minimum")
    initial = [np.polyval(quadratic, v0), 2 * quadratic[0] * v0, 4., v0]
    parameters, covariance = curve_fit(birchmurnaghan, volumes, energies, p0=initial,
        bounds=([-np.inf, 1e-12, -20, volumes.min()], [np.inf, np.inf, 30, volumes.max()]), maxfev=20000)
    e0, modulus, derivative, volume = parameters
    if not volumes.min() + 1e-8 < volume < volumes.max() - 1e-8:
        raise ValueError("Fitted equilibrium lies on the volume boundary; sample a wider range")
    residual = energies - birchmurnaghan(volumes, *parameters)
    return {"volume_A3": float(volume), "energy_eV": float(e0),
            "bulk_modulus_GPa": float(modulus * 160.2176634), "bulk_derivative": float(derivative),
            "residuals_eV": residual, "rms_error_eV": float(np.sqrt(np.mean(residual ** 2))),
            "parameter_covariance": covariance}


def elastic_tensor(strains, stresses_GPa):
    """Fit sigma=sigma0+C epsilon with engineering-Voigt order xx,yy,zz,yz,xz,xy.

    Strain shears are 2*epsilon_ij; stress shears are sigma_ij, with tension
    positive. Input stresses are GPa, not eV/A^3 or pressure-positive outputs.
    Requires a full-rank strain design including an intercept. Symmetrizes C
    and reports the antisymmetric mismatch and residual. Stability is the
    positive-definiteness test appropriate to small strain at zero prestress.
    """
    import numpy as np
    strain = array(strains, "strains", 2)
    stress = array(stresses_GPa, "stresses_GPa", 2)
    if strain.shape != stress.shape or strain.shape[1] != 6 or len(strain) < 7:
        raise ValueError("Provide matching (samples>=7,6) engineering strains and stresses")
    design = np.column_stack((np.ones(len(strain)), strain))
    if np.linalg.matrix_rank(design) != 7:
        raise ValueError("Strains must independently sample all six strain components")
    coefficients = np.linalg.lstsq(design, stress, rcond=None)[0]
    raw = coefficients[1:].T
    stiffness = (raw + raw.T) / 2
    eigenvalues = np.linalg.eigvalsh(stiffness)
    stable = bool(np.all(eigenvalues > 0))
    moduli = None
    if stable:
        compliance = np.linalg.inv(stiffness)
        def diagonal(m):
            return np.trace(m[:3, :3]), m[0, 1] + m[0, 2] + m[1, 2], np.trace(m[3:, 3:])
        a, b, c = diagonal(stiffness)
        bulk_v = (a + 2*b) / 9; shear_v = (a - b + 3*c) / 15
        a, b, c = diagonal(compliance)
        bulk_r = 1 / (a + 2*b); shear_r = 15 / (4*a - 4*b + 3*c)
        bulk = (bulk_v + bulk_r) / 2; shear = (shear_v + shear_r) / 2
        moduli = {"bulk_voigt_GPa": bulk_v, "bulk_reuss_GPa": bulk_r,
                  "shear_voigt_GPa": shear_v, "shear_reuss_GPa": shear_r,
                  "bulk_hill_GPa": bulk, "shear_hill_GPa": shear,
                  "young_hill_GPa": 9*bulk*shear/(3*bulk+shear),
                  "poisson_hill": (3*bulk-2*shear)/(2*(3*bulk+shear)),
                  "universal_anisotropy": 5*shear_v/shear_r + bulk_v/bulk_r - 6}
    residual = stress - (coefficients[0] + strain @ stiffness.T)
    return {"stiffness_GPa": stiffness, "prestress_GPa": coefficients[0],
            "stable_zero_prestress": stable, "stiffness_eigenvalues_GPa": eigenvalues,
            "antisymmetric_norm_GPa": float(np.linalg.norm(raw - raw.T)),
            "rms_stress_error_GPa": float(np.sqrt(np.mean(residual ** 2))), "moduli": moduli}


def _modes(energies_eV, weights):
    import numpy as np
    energies = array(energies_eV, "energies_eV", 2)
    weight = np.ones(len(energies)) if weights is None else array(weights, "weights", 1)
    if weight.shape != (len(energies),) or np.any(weight < 0) or weight.sum() <= 0:
        raise ValueError("Supply one nonnegative q-point weight per row, with positive total")
    return energies, weight / weight.sum()


def phonon_dos(energies_eV, energy_grid_eV, *, sigma_eV=.001, weights=None):
    """Gaussian harmonic DOS, normalized to number of branches per cell.

    Input is (qpoints,branches) in eV with negative energies denoting imaginary
    frequencies. Those are retained on the negative axis; they are not stable
    modes. Finite plot bounds can truncate Gaussian tails; enclosed weight is
    returned without silently renormalizing the plotted curve.
    """
    import numpy as np
    energies, weight = _modes(energies_eV, weights)
    grid = array(energy_grid_eV, "energy_grid_eV", 1)
    if len(grid) < 2 or np.any(np.diff(grid) <= 0):
        raise ValueError("Energy grid must contain at least two increasing samples")
    sigma = positive(sigma_eV, "sigma_eV")
    density = np.zeros(len(grid))
    for row, value in zip(energies, weight):
        for mode in row:
            density += value * np.exp(-.5*((grid-mode)/sigma)**2) / (sigma*np.sqrt(2*np.pi))
    return {"energy_eV": grid, "dos_per_eV": density, "total_modes": energies.shape[1],
            "enclosed_modes": float(np.trapezoid(density, grid)),
            "imaginary_weight": float(np.sum(weight[:, None] * (energies < 0)))}


def harmonic_thermodynamics(energies_eV, temperatures_K, *, weights=None, zero_tolerance_eV=1e-8):
    """Harmonic F,U,S,Cv per cell, including zero-point energy.

    Negative modes below -tolerance are rejected. Modes within +/-tolerance
    (e.g. exact Gamma acoustic zeros) are omitted and their weight reported;
    use a converged q mesh. Outputs eV/cell and eV/(cell K), not per mole.
    Electronic, configurational and anharmonic contributions are excluded.
    """
    import numpy as np
    energies, weight = _modes(energies_eV, weights)
    temperatures = array(temperatures_K, "temperatures_K", 1)
    tolerance = positive(zero_tolerance_eV, "zero_tolerance_eV")
    if np.any(temperatures < 0) or np.any(energies < -tolerance):
        raise ValueError("Nonnegative temperatures and dynamically stable phonons are required")
    selected = energies > tolerance
    mode_weights = np.broadcast_to(weight[:, None], energies.shape)[selected]
    modes = energies[selected]
    if not len(modes):
        raise ValueError("No positive-frequency modes remain")
    kb = 8.617333262145e-5
    zpe = float(np.sum(mode_weights * modes / 2))
    free, internal, entropy, capacity = [], [], [], []
    for temperature in temperatures:
        if temperature == 0:
            f = u = zpe; s = cv = 0.
        else:
            x = modes / (kb * temperature)
            decay = np.exp(-x)
            denominator = -np.expm1(-x)
            f = zpe + kb * temperature * np.sum(mode_weights * np.log(denominator))
            u = zpe + np.sum(mode_weights * modes * decay / denominator)
            s = (u-f) / temperature
            # Avoid inf*0 for very low T by omitting frozen modes (x>700).
            active = x < 700
            cv = kb * np.sum(mode_weights[active] * x[active]**2 * decay[active] / denominator[active]**2)
        free.append(f); internal.append(u); entropy.append(s); capacity.append(cv)
    return {"temperature_K": temperatures, "free_energy_eV": np.asarray(free),
            "internal_energy_eV": np.asarray(internal), "entropy_eV_per_K": np.asarray(entropy),
            "heat_capacity_eV_per_K": np.asarray(capacity), "zero_point_energy_eV": zpe,
            "omitted_zero_mode_weight": float(np.sum(weight[:, None] * ~selected))}
