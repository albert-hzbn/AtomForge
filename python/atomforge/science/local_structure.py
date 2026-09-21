"""Local environment and defect analyses with explicit periodic geometry."""

from ._validation import array, cell_matrix, integer, minimum_image, points, positive


def _neighbors(positions, cutoff, cell, pbc):
    import numpy as np
    from ase.neighborlist import primitive_neighbor_list
    coordinates = points(positions)
    radius = positive(cutoff, "cutoff_A")
    if len(pbc) != 3 or any(not isinstance(v, (bool, np.bool_)) for v in pbc):
        raise ValueError("pbc must contain three booleans")
    matrix = cell_matrix(cell) if cell is not None else np.eye(3)
    if any(pbc) and cell is None:
        raise ValueError("Periodic analysis requires a cell")
    i, j, shifts, vectors = primitive_neighbor_list("ijSD", pbc, matrix, coordinates, radius,
                                                   self_interaction=False)
    if np.any(np.linalg.norm(vectors, axis=1) < 1e-10):
        raise ValueError("Coincident atoms make local environment analysis undefined")
    return coordinates, matrix, i, j, shifts, vectors


def local_strain(reference, current, cutoff_A, *, reference_cell=None, current_cell=None,
                 pbc=(False, False, False)):
    """Least-squares deformation gradient, Green-Lagrange strain and D2min.

    Atom IDs/order must correspond. Neighbors are selected in the reference.
    D2min is the SUM of residual squared distances (A^2), not a per-neighbor
    mean. Rank-deficient environments have valid=False and NaN tensors; they
    must not be interpreted as unstrained. Box deformation is included.
    """
    import numpy as np
    ref, ref_cell, centers, neighbors, shifts, vectors = _neighbors(reference, cutoff_A, reference_cell, pbc)
    now = points(current, "current")
    if now.shape != ref.shape:
        raise ValueError("Reference/current atom ordering and counts must correspond")
    cell = cell_matrix(current_cell) if current_cell is not None else ref_cell
    count = len(ref)
    deformation = np.full((count, 3, 3), np.nan)
    strain = deformation.copy()
    residual = np.full(count, np.nan)
    valid = np.zeros(count, dtype=bool)
    coordination = np.bincount(centers, minlength=count)
    for index in range(count):
        selected = centers == index
        original = vectors[selected]
        if len(original) < 3 or np.linalg.matrix_rank(original) < 3:
            continue
        displaced = now[neighbors[selected]] - now[index] + shifts[selected] @ cell
        if any(pbc):
            prediction = original @ np.linalg.solve(ref_cell, cell)
            displaced = prediction + minimum_image(displaced - prediction, cell, pbc)
        # Row-vector displacements Y = X F^T.
        transpose, _, _, _ = np.linalg.lstsq(original, displaced, rcond=None)
        deformation[index] = transpose.T
        strain[index] = (transpose @ transpose.T - np.eye(3)) / 2
        residual[index] = np.sum((displaced - original @ transpose) ** 2)
        valid[index] = True
    return {"deformation_gradient": deformation, "green_lagrange_strain": strain,
            "d2min_A2": residual, "valid": valid, "coordination": coordination}


def centrosymmetry(positions, cutoff_A, *, neighbors=12, cell=None, pbc=(False, False, False)):
    """Kelchner CSP: sum the N/2 smallest |r_i+r_j|^2 of all distinct pairs.

    Selects the nearest N neighbors within cutoff, including periodic images.
    This is the conventional smallest-pair algorithm, not disjoint optimal
    matching. Missing neighbors produce valid=False and NaN, not zero.
    """
    import numpy as np
    integer(neighbors, "neighbors", 2)
    if neighbors % 2:
        raise ValueError("neighbors must be even")
    coordinates, _, centers, _, _, vectors = _neighbors(positions, cutoff_A, cell, pbc)
    values = np.full(len(coordinates), np.nan)
    for index in range(len(coordinates)):
        bonds = vectors[centers == index]
        if len(bonds) < neighbors:
            continue
        bonds = bonds[np.argsort(np.linalg.norm(bonds, axis=1))[:neighbors]]
        a, b = np.triu_indices(neighbors, 1)
        costs = np.sum((bonds[a] + bonds[b]) ** 2, axis=1)
        values[index] = np.sort(costs)[:neighbors // 2].sum()
    return {"centrosymmetry_A2": values, "valid": np.isfinite(values)}


def bond_order(positions, cutoff_A, *, degrees=(4, 6), cell=None, pbc=(False, False, False)):
    """Local rotationally invariant Steinhardt q_l, with equal bond weights.

    q_l = sqrt(4*pi/(2*l+1) sum_m |mean_j Y_lm(r_ij)|^2).
    Reports q_l (not w_l or neighbor-averaged qbar_l). No neighbors => invalid.
    """
    import numpy as np
    from scipy.special import sph_harm_y
    orders = list(degrees)
    if not orders or len(set(orders)) != len(orders):
        raise ValueError("Supply distinct degrees")
    for order in orders:
        integer(order, "degree", 0)
        if order > 32:
            raise ValueError("degree must be <= 32")
    coordinates, _, centers, _, _, vectors = _neighbors(positions, cutoff_A, cell, pbc)
    coordination = np.bincount(centers, minlength=len(coordinates))
    result = {}
    for order in orders:
        q = np.full(len(coordinates), np.nan)
        for index in range(len(coordinates)):
            bonds = vectors[centers == index]
            if not len(bonds):
                continue
            polar = np.arccos(np.clip(bonds[:, 2] / np.linalg.norm(bonds, axis=1), -1, 1))
            azimuth = np.arctan2(bonds[:, 1], bonds[:, 0])
            harmonics = np.array([sph_harm_y(order, m, polar, azimuth).mean()
                                  for m in range(-order, order + 1)])
            q[index] = np.sqrt(4 * np.pi / (2 * order + 1) * np.sum(abs(harmonics) ** 2))
        result["q" + str(order)] = q
    return {"order": result, "coordination": coordination, "valid": coordination > 0}


def wigner_seitz(reference_sites, positions, *, cell=None, current_cell=None,
                  pbc=(False, False, False)):
    """Assign each atom to its nearest reference site; count vacancies/excess atoms.

    A site with occupancy n>1 contributes n-1 interstitial excess atoms. This
    does not identify a unique interstitial atom within that site. A supplied
    current_cell maps homogeneous box deformation back to the reference cell.
    Ties choose the lowest site index and are flagged for inspection.
    """
    import numpy as np
    sites = points(reference_sites, "reference_sites")
    atoms = np.asarray(positions, dtype=float)
    if atoms.size == 0 and atoms.shape in ((0,), (0, 3)):
        atoms = atoms.reshape(0, 3)
    else:
        atoms = points(positions)
    if current_cell is not None:
        atoms = atoms @ np.linalg.solve(cell_matrix(current_cell), cell_matrix(cell))
    assignments, distances, ties = [], [], []
    for atom in atoms:
        delta = minimum_image(sites - atom, cell, pbc)
        norms = np.linalg.norm(delta, axis=1)
        selected = int(np.argmin(norms))
        assignments.append(selected)
        distances.append(norms[selected])
        ties.append(int(np.count_nonzero(np.isclose(norms, norms[selected], rtol=1e-10, atol=1e-10))) > 1)
    occupancy = np.bincount(assignments, minlength=len(sites))
    return {"occupancy": occupancy, "site_index": np.asarray(assignments, dtype=int),
            "distance_A": np.asarray(distances), "ambiguous": np.asarray(ties, dtype=bool),
            "vacancy_sites": np.flatnonzero(occupancy == 0),
            "interstitial_sites": np.flatnonzero(occupancy > 1),
            "vacancies": int(np.sum(occupancy == 0)),
            "interstitial_excess": int(np.maximum(occupancy - 1, 0).sum())}


def static_structure_factor(positions, q_vectors, *, weights=None):
    """S(q)=<|sum_j b_j exp(i q.r_j)|^2>/sum_j |b_j|^2.

    Positions have shape (atoms,3) or (frames,atoms,3). q vectors are Cartesian
    radians/A (include 2*pi in reciprocal lattice vectors). Equal weights give
    S(0)=N. No forward-scattering subtraction, powder average or quantum effects
    are implied. Use commensurate q for wrapped periodic coordinates.
    """
    import numpy as np
    coordinates = array(positions, "positions")
    if coordinates.ndim == 2:
        coordinates = coordinates[None, :, :]
    if coordinates.ndim != 3 or coordinates.shape[-1] != 3:
        raise ValueError("positions must have shape (atoms,3) or (frames,atoms,3)")
    wavevectors = points(q_vectors, "q_vectors")
    factors = np.ones(coordinates.shape[1]) if weights is None else np.asarray(weights, dtype=complex)
    if factors.shape != (coordinates.shape[1],) or not np.isfinite(factors).all() or np.sum(abs(factors) ** 2) == 0:
        raise ValueError("weights must contain one finite scattering amplitude per atom, not all zero")
    intensity = np.zeros(len(wavevectors))
    # Bound temporary memory to one frame and a block of wave vectors.
    for frame in coordinates:
        for start in range(0, len(wavevectors), 64):
            amplitude = factors @ np.exp(1j * (frame @ wavevectors[start:start + 64].T))
            intensity[start:start + 64] += abs(amplitude) ** 2
    intensity /= len(coordinates) * np.sum(abs(factors) ** 2)
    return {"q_vectors_rad_per_A": wavevectors, "S_q": intensity}
