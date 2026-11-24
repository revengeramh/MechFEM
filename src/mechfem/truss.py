
import numpy as np


def truss_element_stiffness(E, A, x1, y1, x2, y2):
    """
    Calculate the global stiffness matrix of a 2D truss element.
    DOF order: [u1, v1, u2, v2].
    """

    if E <= 0 or A <= 0:
        raise ValueError("E and A must be positive.")

    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx**2 + dy**2)

    if L <= 0:
        raise ValueError("Element length must be positive.")

    c = dx / L
    s = dy / L

    k = (E * A / L) * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
    ])

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
    fixed_dofs : list
        Constrained global degrees of freedom.

    Returns
    -------
    displacements : ndarray
        Global nodal displacement vector.
    reactions : ndarray
        Global support reaction vector.
    """

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    loads = np.asarray(loads, dtype=float)

    if E <= 0 or A <= 0:
        raise ValueError("E and A must be positive.")

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    n_nodes = len(nodes)
    n_dofs = 2 * n_nodes

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(elements) == 0:
        raise ValueError("At least one element is required.")

    if np.any(elements < 0) or np.any(elements >= n_nodes):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    if loads.shape != (n_dofs,):
        raise ValueError("Loads must have length 2*n_nodes.")

    if not np.all(np.isfinite(loads)):
        raise ValueError("Loads must be finite.")

    fixed_dofs = np.asarray(fixed_dofs, dtype=int)

    if fixed_dofs.ndim != 1:
        raise ValueError("fixed_dofs must be a one-dimensional list.")

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dofs):
        raise ValueError("Invalid fixed degree of freedom.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate constrained degrees of freedom.")

    free_dofs = np.setdiff1d(
        np.arange(n_dofs),
        fixed_dofs
    )

    if len(free_dofs) == 0:
        raise ValueError("No free degrees of freedom.")

    K = np.zeros((n_dofs, n_dofs))

    for element in elements:
        n1, n2 = element

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

    displacements = np.zeros(n_dofs)
    displacements[free_dofs] = Uf

    reactions = K @ displacements - loads

    return displacements, reactions

def calculate_truss_results(nodes, elements, E, A, displacements):
    """
    Calculate axial force and stress for each 2D truss element.
    Positive force indicates tension; negative force indicates compression.
    """

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    displacements = np.asarray(displacements, dtype=float)

    if E <= 0 or A <= 0:
        raise ValueError("E and A must be positive.")

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(elements) == 0:
        raise ValueError("At least one element is required.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

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

        dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
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

    return np.array(axial_forces), np.array(axial_stresses)


def assemble_truss_matrices(nodes, elements, E, A, density):
    """
    Assemble global stiffness and consistent mass matrices
    for a 2D truss structure.

    Parameters
    ----------
    nodes : array-like
        Nodal coordinates in meters, shape (n_nodes, 2).
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

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)

    if E <= 0 or A <= 0 or density <= 0:
        raise ValueError("E, A, and density must be positive.")

    if not np.all(np.isfinite([E, A, density])):
        raise ValueError("Material properties must be finite.")

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(elements) == 0:
        raise ValueError("At least one element is required.")

    n_nodes = len(nodes)
    n_dofs = 2 * n_nodes

    if np.any(elements < 0) or np.any(elements >= n_nodes):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    K = np.zeros((n_dofs, n_dofs))
    M = np.zeros((n_dofs, n_dofs))

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

        # Consistent mass matrix in global coordinates
        mass_factor = density * A * L / 6.0

        identity = np.eye(2)

        me = mass_factor * np.block([
            [2.0 * identity, identity],
            [identity, 2.0 * identity]
        ])

        dofs = [
            2 * n1,
            2 * n1 + 1,
            2 * n2,
            2 * n2 + 1
        ]

        K[np.ix_(dofs, dofs)] += ke
        M[np.ix_(dofs, dofs)] += me

    return K, M
