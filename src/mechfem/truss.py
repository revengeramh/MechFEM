
import numpy as np


def _validate_positive_finite_scalar(value, name):
    """Validate a positive finite scalar."""
    if (
        not isinstance(value, (int, float, np.number))
        or isinstance(value, (bool, np.bool_))
        or not np.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be positive and finite.")


def _validate_nodes(nodes):
    """Validate 2D truss nodal coordinates."""
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
    """Validate connectivity without silently truncating indices."""
    try:
        raw = np.asarray(elements, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid element connectivity.") from exc

    if raw.ndim != 2 or raw.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(raw) == 0:
        raise ValueError("At least one element is required.")

    if not np.all(np.isfinite(raw)):
        raise ValueError("Element connectivity must be finite.")

    if np.any(raw != np.floor(raw)):
        raise ValueError("Element connectivity must contain integers.")

    if np.any(raw < 0) or np.any(raw >= n_nodes):
        raise ValueError("Invalid element connectivity.")

    elements = raw.astype(int)

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    normalized = [
        tuple(sorted(row))
        for row in elements
    ]

    if len(set(normalized)) != len(normalized):
        raise ValueError("Duplicate elements are not allowed.")

    return elements


def _validate_fixed_dofs(fixed_dofs, n_dofs):
    """Validate constrained global degrees of freedom."""
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

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dofs):
        raise ValueError("Invalid fixed degree of freedom.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate constrained degrees of freedom.")

    return fixed_dofs


def truss_element_stiffness(E, A, x1, y1, x2, y2):
    """
    Calculate the global stiffness matrix of a 2D truss element.

    DOF order: [u1, v1, u2, v2].
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(A, "A")

    coordinates = np.asarray([x1, y1, x2, y2], dtype=float)

    if not np.all(np.isfinite(coordinates)):
        raise ValueError("Element coordinates must be finite.")

    dx = x2 - x1
    dy = y2 - y1
    L = np.hypot(dx, dy)

    if L <= 0:
        raise ValueError("Element length must be positive.")

    c = dx / L
    s = dy / L

    k = (E * A / L) * np.array([
        [c*c, c*s, -c*c, -c*s],
        [c*s, s*s, -c*s, -s*s],
        [-c*c, -c*s, c*c, c*s],
        [-c*s, -s*s, c*s, s*s]
    ], dtype=float)

    return k


def solve_truss(nodes, elements, E, A, loads, fixed_dofs):
    """
    Solve a 2D truss structure using the direct stiffness method.

    Parameters
    ----------
    nodes : array-like
        Nodal coordinates, shape (n_nodes, 2).
    elements : array-like
        Element connectivity using zero-based node indices.
    E : float
        Young's modulus.
    A : float
        Cross-sectional area.
    loads : array-like
        Global force vector, shape (2*n_nodes,).
    fixed_dofs : array-like
        Constrained global degrees of freedom.

    Returns
    -------
    displacements : ndarray
        Global nodal displacement vector.
    reactions : ndarray
        Global support reaction vector.
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(A, "A")

    nodes = _validate_nodes(nodes)
    elements = _validate_elements(elements, len(nodes))

    try:
        loads = np.asarray(loads, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid load vector.") from exc

    n_nodes = len(nodes)
    n_dofs = 2 * n_nodes

    if loads.shape != (n_dofs,):
        raise ValueError("Loads must have length 2*n_nodes.")

    if not np.all(np.isfinite(loads)):
        raise ValueError("Loads must be finite.")

    fixed_dofs = _validate_fixed_dofs(fixed_dofs, n_dofs)

    free_dofs = np.setdiff1d(
        np.arange(n_dofs),
        fixed_dofs
    )

    if len(free_dofs) == 0:
        raise ValueError("No free degrees of freedom.")

    K = np.zeros((n_dofs, n_dofs), dtype=float)

    for n1, n2 in elements:
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        ke = truss_element_stiffness(
            E, A, x1, y1, x2, y2
        )

        dofs = [
            2*n1,
            2*n1 + 1,
            2*n2,
            2*n2 + 1
        ]

        K[np.ix_(dofs, dofs)] += ke

    Kff = K[np.ix_(free_dofs, free_dofs)]
    Ff = loads[free_dofs]

    try:
        Uf = np.linalg.solve(Kff, Ff)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "The structure is unstable or insufficiently constrained."
        ) from exc

    displacements = np.zeros(n_dofs, dtype=float)
    displacements[free_dofs] = Uf

    reactions = K @ displacements - loads

    return displacements, reactions


def calculate_truss_results(nodes, elements, E, A, displacements):
    """
    Calculate axial force and stress for each 2D truss element.

    Positive force indicates tension.
    Negative force indicates compression.
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(A, "A")

    nodes = _validate_nodes(nodes)
    elements = _validate_elements(elements, len(nodes))

    try:
        displacements = np.asarray(displacements, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Invalid displacement vector.") from exc

    if displacements.shape != (2 * len(nodes),):
        raise ValueError("Invalid displacement vector size.")

    if not np.all(np.isfinite(displacements)):
        raise ValueError("Displacements must be finite.")

    axial_forces = []
    axial_stresses = []

    for n1, n2 in elements:
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        dx = x2 - x1
        dy = y2 - y1
        L = np.hypot(dx, dy)

        if L <= 0:
            raise ValueError("Element length must be positive.")

        c = dx / L
        s = dy / L

        dofs = [
            2*n1,
            2*n1 + 1,
            2*n2,
            2*n2 + 1
        ]

        u = displacements[dofs]

        extension = np.dot(
            np.array([-c, -s, c, s]),
            u
        )

        strain = extension / L
        stress = E * strain
        force = A * stress

        axial_forces.append(force)
        axial_stresses.append(stress)

    return (
        np.asarray(axial_forces),
        np.asarray(axial_stresses)
    )


def assemble_truss_matrices(nodes, elements, E, A, density):
    """
    Assemble global stiffness and consistent mass matrices
    for a 2D truss structure.

    Parameters
    ----------
    nodes : array-like
        Nodal coordinates, shape (n_nodes, 2).
    elements : array-like
        Element connectivity, shape (n_elements, 2).
    E : float
        Young's modulus in Pa.
    A : float
        Cross-sectional area in m^2.
    density : float
        Material density in kg/m^3.

    Returns
    -------
    K : ndarray
        Global stiffness matrix.
    M : ndarray
        Global consistent mass matrix.
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(A, "A")
    _validate_positive_finite_scalar(density, "Density")

    nodes = _validate_nodes(nodes)
    elements = _validate_elements(elements, len(nodes))

    n_nodes = len(nodes)
    n_dofs = 2 * n_nodes

    K = np.zeros((n_dofs, n_dofs), dtype=float)
    M = np.zeros((n_dofs, n_dofs), dtype=float)

    identity = np.eye(2)

    for n1, n2 in elements:
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        dx = x2 - x1
        dy = y2 - y1
        L = np.hypot(dx, dy)

        if L <= 0:
            raise ValueError("Element length must be positive.")

        ke = truss_element_stiffness(
            E, A, x1, y1, x2, y2
        )

        mass_factor = density * A * L / 6.0

        me = mass_factor * np.block([
            [2.0 * identity, identity],
            [identity, 2.0 * identity]
        ])

        dofs = [
            2*n1,
            2*n1 + 1,
            2*n2,
            2*n2 + 1
        ]

        K[np.ix_(dofs, dofs)] += ke
        M[np.ix_(dofs, dofs)] += me

    return K, M
