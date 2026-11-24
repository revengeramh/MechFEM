
import numpy as np


def beam_element_stiffness(E, I, L):
    """
    Return the Euler-Bernoulli beam element stiffness matrix.
    DOF order: [v1, theta1, v2, theta2]
    """
    if E <= 0 or I <= 0 or L <= 0:
        raise ValueError("E, I, and L must be positive.")

    factor = E * I / L**3

    return factor * np.array([
        [12, 6*L, -12, 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]
    ], dtype=float)


def solve_beam(nodes, elements, E, I, loads, fixed_dofs):
    """
    Solve a 1D Euler-Bernoulli beam using FEM.

    Parameters
    ----------
    nodes : array
        Nodal coordinates.
    elements : list of tuples
        Element connectivity (node_i, node_j).
    E : float
        Young's modulus.
    I : float
        Second moment of area.
    loads : array
        Global nodal load vector.
    fixed_dofs : list
        Constrained global DOF indices.

    Returns
    -------
    displacement : ndarray
        Nodal displacement and rotation vector.
    reactions : ndarray
        Support reaction vector.
    """
    nodes = np.asarray(nodes, dtype=float)
    loads = np.asarray(loads, dtype=float)

    if E <= 0 or I <= 0:
        raise ValueError("E and I must be positive.")

    if len(nodes) < 2 or np.any(np.diff(nodes) <= 0):
        raise ValueError("Nodes must be strictly increasing.")

    n_dofs = 2 * len(nodes)

    if loads.shape != (n_dofs,):
        raise ValueError("Incorrect load vector size.")

    K = np.zeros((n_dofs, n_dofs))

    for n1, n2 in elements:
        if not (0 <= n1 < len(nodes) and 0 <= n2 < len(nodes)):
            raise ValueError("Invalid element connectivity.")

        L = nodes[n2] - nodes[n1]
        k = beam_element_stiffness(E, I, L)

        dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]

        K[np.ix_(dofs, dofs)] += k

    fixed_dofs = np.asarray(fixed_dofs, dtype=int)

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dofs):
        raise ValueError("Invalid constrained DOF.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate constrained DOFs.")

    free_dofs = np.setdiff1d(
        np.arange(n_dofs), fixed_dofs
    )

    displacement = np.zeros(n_dofs)

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = loads[free_dofs]

    if len(free_dofs) > 0:
        displacement[free_dofs] = np.linalg.solve(K_ff, F_f)

    reactions = K @ displacement - loads

    return displacement, reactions


if __name__ == "__main__":
    # Cantilever beam example
    L = 2.0
    E = 200e9
    I = 8e-6
    P = -1000.0

    nodes = [0.0, L]
    elements = [(0, 1)]

    loads = [0.0, 0.0, P, 0.0]
    fixed_dofs = [0, 1]

    displacement, reactions = solve_beam(
        nodes, elements, E, I, loads, fixed_dofs
    )

    analytical_deflection = P * L**3 / (3 * E * I)

    print("MechFEM - Beam Analysis")
    print("-----------------------")
    print(f"Tip displacement: {displacement[2]:.8e} m")
    print(f"Tip rotation:     {displacement[3]:.8e} rad")
    print(f"Support reaction: {reactions[0]:.2f} N")
    print(f"Analytical deflection: {analytical_deflection:.8e} m")
