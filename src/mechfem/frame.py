
import numpy as np


def frame_element_stiffness(E, A, I, x1, y1, x2, y2):
    """
    Calculate the global stiffness matrix of a 2D frame element.

    Degrees of freedom per node:
    [u, v, theta]

    Parameters
    ----------
    E : float
        Young's modulus.
    A : float
        Cross-sectional area.
    I : float
        Second moment of area.
    x1, y1 : float
        Coordinates of the first node.
    x2, y2 : float
        Coordinates of the second node.

    Returns
    -------
    numpy.ndarray
        6x6 global element stiffness matrix.
    """

    if E <= 0 or A <= 0 or I <= 0:
        raise ValueError("E, A, and I must be positive.")

    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx**2 + dy**2)

    if L <= 0:
        raise ValueError("Element length must be positive.")

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
    """
    Solve a linear static 2D frame structure.

    Parameters
    ----------
    nodes : array-like
        Node coordinates, shape (n_nodes, 2).
    elements : array-like
        Element connectivity, shape (n_elements, 2).
    E, A, I : float
        Material and section properties (uniform for all elements).
    loads : dict
        Global force vector indexed by global DOF.
    fixed_dofs : list
        Constrained global DOF indices.

    Returns
    -------
    displacements : numpy.ndarray
        Global displacement vector.
    reactions : numpy.ndarray
        Global reaction vector.
    """

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("elements must have shape (n_elements, 2).")

    if len(nodes) == 0 or len(elements) == 0:
        raise ValueError("At least one node and one element are required.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if E <= 0 or A <= 0 or I <= 0:
        raise ValueError("E, A, and I must be positive.")

    n_dof = 3 * len(nodes)
    K = np.zeros((n_dof, n_dof), dtype=float)
    F = np.zeros(n_dof, dtype=float)

    for dof, value in loads.items():
        if not isinstance(dof, (int, np.integer)) or not 0 <= dof < n_dof:
            raise ValueError("Invalid load degree of freedom.")
        F[dof] += float(value)

    for element in elements:
        n1, n2 = element

        if not (0 <= n1 < len(nodes) and 0 <= n2 < len(nodes)):
            raise ValueError("Invalid element connectivity.")

        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        ke = frame_element_stiffness(E, A, I, x1, y1, x2, y2)

        dofs = [
            3*n1, 3*n1+1, 3*n1+2,
            3*n2, 3*n2+1, 3*n2+2
        ]

        K[np.ix_(dofs, dofs)] += ke

    fixed_dofs = np.asarray(fixed_dofs, dtype=int)

    if fixed_dofs.ndim != 1:
        raise ValueError("fixed_dofs must be a one-dimensional list.")

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dof):
        raise ValueError("Invalid fixed degree of freedom.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate fixed degrees of freedom.")

    free_dofs = np.setdiff1d(np.arange(n_dof), fixed_dofs)

    if len(free_dofs) == 0:
        raise ValueError("No free degrees of freedom.")

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = F[free_dofs]

    try:
        u_f = np.linalg.solve(K_ff, F_f)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "The structure is unstable or has insufficient constraints."
        ) from exc

    displacements = np.zeros(n_dof)
    displacements[free_dofs] = u_f

    reactions = K @ displacements - F

    return displacements, reactions
