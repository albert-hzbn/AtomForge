"""Trajectory statistics. Positions: Angstrom; velocities: Angstrom/fs; time: fs."""

from ._validation import array, boolean, cell_matrix, integer, minimum_image, positive


def _trajectory(values, name):
    result = array(values, name, 3)
    if result.shape[0] < 2 or result.shape[2] != 3:
        raise ValueError(name + " must have shape (at least two frames, atoms, 3)")
    return result


def mean_square_displacement(positions, timestep_fs, *, cell=None, pbc=(True, True, True),
                             wrapped=False, remove_drift=False, max_lag=None):
    """Time-origin-averaged MSD and Cartesian displacement covariance.

    Atom ordering must be stable. Unwrapped coordinates are recommended.
    With wrapped=True, unwrap successive displacements in a FIXED cell using
    its shortest image; sampling must resolve motion below the MIC ambiguity.
    Variable-cell trajectories must first be mapped to an explicit reference
    frame. remove_drift subtracts the geometric centroid at every frame.
    Covariance is <dr_a dr_b>, not a covariance about the mean displacement.
    """
    import numpy as np
    coordinates = _trajectory(positions, "positions").copy()
    wrapped = boolean(wrapped, "wrapped")
    remove_drift = boolean(remove_drift, "remove_drift")
    dt = positive(timestep_fs, "timestep_fs")
    if wrapped:
        matrix = cell_matrix(cell)
        steps = minimum_image(np.diff(coordinates, axis=0), matrix, pbc)
        coordinates[1:] = coordinates[0] + np.cumsum(steps, axis=0)
    if remove_drift:
        coordinates -= coordinates.mean(axis=1, keepdims=True)
    count = len(coordinates)
    limit = count - 1 if max_lag is None else integer(max_lag, "max_lag")
    if limit >= count:
        raise ValueError("max_lag must be smaller than frame count")
    covariance = np.zeros((limit + 1, 3, 3))
    for lag in range(1, limit + 1):
        delta = coordinates[lag:] - coordinates[:-lag]
        covariance[lag] = np.einsum("tni,tnj->ij", delta, delta) / (delta.shape[0] * delta.shape[1])
    return {"lag_fs": np.arange(limit + 1) * dt,
            "msd_A2": np.trace(covariance, axis1=1, axis2=2),
            "components_A2": np.diagonal(covariance, axis1=1, axis2=2).copy(),
            "displacement_tensor_A2": covariance,
            "origins": count - np.arange(limit + 1)}


def diffusion_coefficient(lag_fs, msd_A2, *, fit_range_fs, dimensions=3):
    """Einstein diffusion from an explicitly selected diffusive time interval.

    D = slope/(2*d); input MSD must sum ONLY the d selected directions.
    Reports slope, intercept, R-squared and OLS fit uncertainty. Overlapping
    time origins correlate MSD samples: the OLS uncertainty is NOT a physical
    confidence interval. Use independent trajectories/block estimates for that.
    A linear fit alone does not establish the diffusive regime.
    """
    import numpy as np
    from scipy.stats import linregress
    times = array(lag_fs, "lag_fs", 1)
    values = array(msd_A2, "msd_A2", 1)
    if dimensions not in (1, 2, 3) or isinstance(dimensions, bool):
        raise ValueError("dimensions must be 1, 2 or 3")
    if times.shape != values.shape or np.any(np.diff(times) <= 0) or times[0] < 0 or np.any(values < 0):
        raise ValueError("MSD must be nonnegative with matching, increasing nonnegative times")
    interval = array(fit_range_fs, "fit_range_fs", 1)
    if interval.shape != (2,) or interval[0] >= interval[1] or interval[0] < times[0] or interval[1] > times[-1]:
        raise ValueError("fit_range_fs must lie inside the sampled time interval")
    chosen = (times >= interval[0]) & (times <= interval[1])
    if chosen.sum() < 3:
        raise ValueError("At least three fit points are required")
    fit = linregress(times[chosen], values[chosen])
    if fit.slope < 0:
        raise ValueError("Negative MSD slope is not a physical diffusion estimate")
    diffusion = fit.slope / (2 * dimensions)
    return {"D_A2_per_fs": diffusion, "D_m2_per_s": diffusion * 1e-5,
            "slope_A2_per_fs": fit.slope, "intercept_A2": fit.intercept,
            "r_squared": float(fit.rvalue ** 2) if np.isfinite(fit.rvalue) else 0.0,
            "ols_slope_stderr": float(fit.stderr), "fit_points": int(chosen.sum()),
            "fit_range_fs": interval}


def velocity_autocorrelation(velocities, timestep_fs, *, max_lag=None,
                             remove_drift=False, normalize=False):
    """Unbiased time-origin average of v(t) dot v(t+lag), per atom.

    remove_drift removes each frame's geometric mean velocity, not its
    mass-weighted centre-of-mass velocity. Units before normalization: A^2/fs^2.
    """
    import numpy as np
    velocity = _trajectory(velocities, "velocities").copy()
    remove_drift = boolean(remove_drift, "remove_drift")
    normalize = boolean(normalize, "normalize")
    dt = positive(timestep_fs, "timestep_fs")
    if remove_drift:
        velocity -= velocity.mean(axis=1, keepdims=True)
    count = len(velocity)
    limit = count - 1 if max_lag is None else integer(max_lag, "max_lag")
    if limit >= count:
        raise ValueError("max_lag must be smaller than frame count")
    # FFT zero padding avoids circular correlation. Divide each lag by the
    # actual number of origins, then average atoms (Cartesian components sum).
    size = 1 << (2 * count - 1).bit_length()
    transformed = np.fft.rfft(velocity, n=size, axis=0)
    correlation = np.fft.irfft(transformed.conj() * transformed, n=size, axis=0)[:limit + 1]
    correlation = correlation.sum(axis=(1, 2)) / (velocity.shape[1] * (count - np.arange(limit + 1)))
    if normalize:
        if correlation[0] <= 0:
            raise ValueError("Cannot normalize zero velocity autocorrelation")
        correlation /= correlation[0]
    return {"lag_fs": np.arange(limit + 1) * dt, "vacf": correlation,
            "unit": "dimensionless" if normalize else "A^2/fs^2",
            "origins": count - np.arange(limit + 1)}


def vibrational_spectrum(velocities, timestep_fs, *, masses=None, window="hann", remove_mean=True):
    """One-sided velocity power spectrum, the Fourier partner of velocity correlation.

    Uses a nonnegative periodogram, with optional atomic-mass weighting, and
    returns unit-integral spectral weight per THz. This is a classical velocity
    spectrum, not quantum neutron/IR/Raman intensity. Resolution and Nyquist
    frequency are determined by duration and sampling; zero padding adds no data.
    """
    import numpy as np
    remove_mean = boolean(remove_mean, "remove_mean")
    velocity = _trajectory(velocities, "velocities").copy()
    dt = positive(timestep_fs, "timestep_fs")
    count, atoms, _ = velocity.shape
    if count < 4:
        raise ValueError("At least four velocity frames are required")
    if remove_mean:
        velocity -= velocity.mean(axis=0, keepdims=True)
    if masses is not None:
        weights = array(masses, "masses", 1)
        if weights.shape != (atoms,) or np.any(weights <= 0):
            raise ValueError("One positive mass is required per atom")
        velocity *= np.sqrt(weights)[None, :, None]
    if window not in ("hann", "none"):
        raise ValueError("window must be hann or none")
    taper = np.hanning(count) if window == "hann" else np.ones(count)
    transformed = np.fft.rfft(velocity * taper[:, None, None], axis=0)
    power = np.sum(abs(transformed) ** 2, axis=(1, 2))
    power[1:-1 if count % 2 == 0 else None] *= 2
    frequency = np.fft.rfftfreq(count, d=dt) * 1000
    area = np.trapezoid(power, frequency)
    if area <= 0:
        raise ValueError("No nonzero vibrational spectral weight")
    return {"frequency_THz": frequency, "density_per_THz": power / area,
            "resolution_THz": 1000 / (count * dt), "nyquist_THz": 500 / dt}
