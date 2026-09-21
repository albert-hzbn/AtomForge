"""Shared validation for numerical scientific tools (optional NumPy/ASE)."""

import math
from numbers import Integral


def array(value, name, ndim=None):
    import numpy as np
    result = np.asarray(value, dtype=float)
    if not result.size or not np.isfinite(result).all() or (ndim is not None and result.ndim != ndim):
        raise ValueError(name + " must be a nonempty finite array of the required dimension")
    return result


def positive(value, name):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")
    return float(value)


def integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))
    return value


def boolean(value, name):
    import numpy as np
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(name + " must be a boolean, not a string or number")
    return bool(value)


def cell_matrix(cell):
    import numpy as np
    matrix = array(cell, "cell", 2)
    if matrix.shape != (3, 3) or abs(np.linalg.det(matrix)) < 1e-12 or np.linalg.cond(matrix) > 1e10:
        raise ValueError("cell must contain three independent, well-conditioned row vectors")
    return matrix


def points(value, name="positions"):
    result = array(value, name, 2)
    if result.shape[1] != 3:
        raise ValueError(name + " must have shape (atoms, 3)")
    return result


def minimum_image(vectors, cell, pbc):
    """Use ASE's general-cell MIC, not fractional rounding in a skew cell."""
    import numpy as np
    if len(pbc) != 3 or any(not isinstance(v, (bool, np.bool_)) for v in pbc):
        raise ValueError("pbc must contain three booleans")
    if not any(pbc):
        return np.asarray(vectors, dtype=float)
    from ase.geometry import find_mic
    matrix = cell_matrix(cell)
    original = np.asarray(vectors, dtype=float)
    return find_mic(original.reshape(-1, 3), matrix, pbc=pbc)[0].reshape(original.shape)
