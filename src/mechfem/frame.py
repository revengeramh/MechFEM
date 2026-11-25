import numpy as np


def _validate_positive_finite(value, name):
    if (
        not isinstance(value, (int, float, np.number))
        or isinstance(value, (bool, np.bool_))
        or not np.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be positive and finite.")


def _validate_nodes(nodes):
    try:
        nodes = np.asarray(nodes, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid nodal coordinates.") from exc

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if len(nodes) < 2:
        raise ValueError("At least two nodes are required.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    return nodes


def _validate_elements(elements, n_nodes):
    try:
        raw = np.asarray(elements, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid element connectivity.") from exc

    if raw.ndim != 2 or raw.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(raw) == 0:
        raise ValueError("At least one element is required.")

    if not np.all(np.isfinite(raw)):
        raise ValueError("Element indices must be finite.")

    if np.any(raw != np.floor(raw)):
        raise ValueError("Element indices must be integers.")

    if np.any(raw < 0) or np.any(raw >= n_nodes):
        raise ValueError("Invalid element connectivity.")

    elements = raw.astype(int)

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    normalized = [tuple(sorted(row)) for row in elements]

    if len(set(normalized)) != len(normalized):
        raise ValueError("Duplicate elements are not allowed.")

    return elements


def _validate_fixed_dofs(fixed_dofs, n_dof):
    try:
        raw = np.asarray(fixed_dofs, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid fixed DOFs.") from exc

    if raw.ndim != 1:
        raise ValueError("fixed_dofs must be one-dimensional.")

    if not np.all(np.isfinite(raw)):
        raise ValueError("Fixed DOFs must be finite.")

    if np.any(raw != np.floor(raw)):
        raise ValueError("Fixed DOFs must contain integers.")

    fixed_dofs = raw.astype(int)

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dof):
        raise ValueError("Invalid fixed degree of freedom.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate fixed degrees of freedom.")

    return fixed_dofs


def frame_element_stiffness(E, A, I, x1, y1, x2, y2):
    """Return the global stiffness matrix of a 2D frame element."""

    _validate_positive_finite(E, "E")
    _validate_positive_finite(A, "A")
    _validate_positive_finite(I, "I")

    try:
        coordinates = np.asarray(
            [x1, y1, x2, y2],
            dtype=float
        )
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid element coordinates.") from exc

    if not np.all(np.isfinite(coordinates)):
        raise ValueError("Element coordinates must be finite.")

    x1, y1, x2, y2 = coordinates

    dx = x2 - x1
    dy = y2 - y1
    L = np.hypot(dx, dy)

    if not np.isfinite(L) or L <= 0:
        raise ValueError("Element length must be positive and finite.")

    c = dx / L
    s = dy / L

    EA_L = E * A / L
    EI = E * I

    k_local = np.array([
        [EA_L, 0, 0, -EA_L, 0, 0],
        [0, 12*EI/L**3, 6*EI/L**2,
         0, -12*EI/L**3, 6*EI/L**2],
        [0, 6*EI/L**2, 4*EI/L,
         0, -6*EI/L**2, 2*EI/L],
        [-EA_L, 0, 0, EA_L, 0, 0],
        [0, -12*EI/L**3, -6*EI/L**2,
         0, 12*EI/L**3, -6*EI/L**2],
        [0, 6*EI/L**2, 2*EI/L,
         0, -6*EI/L**2, 4*EI/L]
    ], dtype=float)

    T = np.array([
        [c, s, 0, 0, 0, 0],
        [-s, c, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0],
        [0, 0, 0, -s, c, 0],
        [0, 0, 0, 0, 0, 1]
    ], dtype=float)

    return T.T @ k_local @ T


def solve_frame(nodes, elements, E, A, I, loads, fixed_dofs):
    """Solve a linear static 2D frame structure."""

    _validate_positive_finite(E, "E")
    _validate_positive_finite(A, "A")
    _validate_positive_finite(I, "I")

    nodes = _validate_nodes(nodes)
    elements = _validate_elements(elements, len(nodes))

    n_nodes = len(nodes)
    n_dof = 3 * n_nodes

    if not isinstance(loads, dict):
        raise ValueError("Loads must be a dictionary.")

    F = np.zeros(n_dof, dtype=float)

    for dof, value in loads.items():
        if (
            not isinstance(dof, (int, np.integer))
            or isinstance(dof, (bool, np.bool_))
            or not 0 <= dof < n_dof
        ):
            raise ValueError("Invalid load degree of freedom.")

        try:
            value = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Load values must be finite numbers.") from exc

        if not np.isfinite(value):
            raise ValueError("Load values must be finite.")

        F[dof] += value

    fixed_dofs = _validate_fixed_dofs(fixed_dofs, n_dof)

    free_dofs = np.setdiff1d(
        np.arange(n_dof),
        fixed_dofs
    )

    if len(free_dofs) == 0:
        raise ValueError("No free degrees of freedom.")

    K = np.zeros((n_dof, n_dof), dtype=float)

    for n1, n2 in elements:
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        ke = frame_element_stiffness(
            E, A, I,
            x1, y1, x2, y2
        )

        dofs = [
            3*n1, 3*n1 + 1, 3*n1 + 2,
            3*n2, 3*n2 + 1, 3*n2 + 2
        ]

        K[np.ix_(dofs, dofs)] += ke

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = F[free_dofs]

    # Check stability using a diagonally scaled stiffness matrix.
    diagonal = np.diag(K_ff)

    if (
        not np.all(np.isfinite(K_ff))
        or np.any(diagonal <= 0)
    ):
        raise ValueError(
            "The structure is unstable or has insufficient constraints."
        )

    scale = np.sqrt(diagonal)
    K_scaled = K_ff / np.outer(scale, scale)

    if np.linalg.matrix_rank(K_scaled) < len(free_dofs):
        raise ValueError(
            "The structure is unstable or has insufficient constraints."
        )

    try:
        u_f = np.linalg.solve(K_ff, F_f)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "The structure is unstable or has insufficient constraints."
        ) from exc

    if not np.all(np.isfinite(u_f)):
        raise ValueError("The calculated displacements are not finite.")

    displacements = np.zeros(n_dof, dtype=float)
    displacements[free_dofs] = u_f

    reactions = K @ displacements - F

    return displacements, reactions
