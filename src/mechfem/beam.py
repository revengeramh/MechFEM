
import numpy as np


def _validate_positive_finite_scalar(value, name):
    """Validate a positive finite scalar."""
    if (
        not isinstance(value, (int, float, np.number))
        or not np.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be positive and finite.")


def beam_element_stiffness(E, I, L):
    """
    Return the Euler-Bernoulli beam element stiffness matrix.

    DOF order: [v1, theta1, v2, theta2]
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(I, "I")
    _validate_positive_finite_scalar(L, "L")

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

    DOF order per node: [transverse displacement, rotation].

    Returns
    -------
    displacement : ndarray
        Nodal displacement and rotation vector.

    reactions : ndarray
        Support reaction vector.
    """
    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(I, "I")

    try:
        nodes = np.asarray(nodes, dtype=float)
        loads = np.asarray(loads, dtype=float)
        raw_elements = np.asarray(elements, dtype=float)
        raw_fixed_dofs = np.asarray(fixed_dofs, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Inputs must contain valid numerical values.") from exc

    if nodes.ndim != 1 or len(nodes) < 2:
        raise ValueError("Nodes must be a one-dimensional array.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if np.any(np.diff(nodes) <= 0):
        raise ValueError("Nodes must be strictly increasing.")

    if (
        raw_elements.ndim != 2
        or raw_elements.shape[1] != 2
        or len(raw_elements) == 0
    ):
        raise ValueError(
            "Elements must be a non-empty array with shape (n_elements, 2)."
        )

    if not np.all(np.isfinite(raw_elements)):
        raise ValueError("Element connectivity must be finite.")

    if np.any(raw_elements != np.floor(raw_elements)):
        raise ValueError("Element connectivity must contain integers.")

    elements = raw_elements.astype(int)

    n_dofs = 2 * len(nodes)

    if loads.shape != (n_dofs,):
        raise ValueError("Incorrect load vector size.")

    if not np.all(np.isfinite(loads)):
        raise ValueError("Loads must contain finite values.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] >= elements[:, 1]):
        raise ValueError(
            "Each element must connect increasing node indices."
        )

    element_pairs = [tuple(row) for row in elements]

    if len(set(element_pairs)) != len(element_pairs):
        raise ValueError("Duplicate elements are not allowed.")

    K = np.zeros((n_dofs, n_dofs), dtype=float)

    for n1, n2 in elements:
        L = nodes[n2] - nodes[n1]
        k = beam_element_stiffness(E, I, L)

        dofs = [
            2 * n1,
            2 * n1 + 1,
            2 * n2,
            2 * n2 + 1
        ]

        K[np.ix_(dofs, dofs)] += k

    if raw_fixed_dofs.ndim != 1:
        raise ValueError("Fixed DOFs must be a one-dimensional array.")

    if not np.all(np.isfinite(raw_fixed_dofs)):
        raise ValueError("Fixed DOFs must contain finite values.")

    if np.any(raw_fixed_dofs != np.floor(raw_fixed_dofs)):
        raise ValueError("Fixed DOFs must contain integers.")

    fixed_dofs = raw_fixed_dofs.astype(int)

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dofs):
        raise ValueError("Invalid constrained DOF.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate constrained DOFs.")

    free_dofs = np.setdiff1d(
        np.arange(n_dofs),
        fixed_dofs
    )

    displacement = np.zeros(n_dofs, dtype=float)

    if len(free_dofs) > 0:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        F_f = loads[free_dofs]

        try:
            displacement[free_dofs] = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "The reduced stiffness matrix is singular. "
                "Check the supports, connectivity, and structural stability."
            ) from exc

    reactions = K @ displacement - loads

    return displacement, reactions


if __name__ == "__main__":
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


def assemble_beam_matrices(nodes, elements, E, I, A, density):
    """
    Assemble global stiffness and consistent mass matrices
    for an Euler-Bernoulli beam.

    DOF order per node:
    [transverse displacement, rotation]
    """
    try:
        nodes = np.asarray(nodes, dtype=float)
        raw_elements = np.asarray(elements, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "Nodes and elements must contain valid numerical values."
        ) from exc

    if nodes.ndim != 1 or len(nodes) < 2:
        raise ValueError("Nodes must be a one-dimensional array.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if np.any(np.diff(nodes) <= 0):
        raise ValueError("Nodes must be strictly increasing.")

    if raw_elements.ndim != 2 or raw_elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(raw_elements) == 0:
        raise ValueError("At least one element is required.")

    if not np.all(np.isfinite(raw_elements)):
        raise ValueError("Element connectivity must be finite.")

    if np.any(raw_elements != np.floor(raw_elements)):
        raise ValueError("Element connectivity must contain integers.")

    elements = raw_elements.astype(int)

    _validate_positive_finite_scalar(E, "E")
    _validate_positive_finite_scalar(I, "I")
    _validate_positive_finite_scalar(A, "A")
    _validate_positive_finite_scalar(density, "Density")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] >= elements[:, 1]):
        raise ValueError(
            "Each element must connect increasing node indices."
        )

    element_pairs = [tuple(row) for row in elements]

    if len(set(element_pairs)) != len(element_pairs):
        raise ValueError("Duplicate elements are not allowed.")

    n_dofs = 2 * len(nodes)

    K = np.zeros((n_dofs, n_dofs), dtype=float)
    M = np.zeros((n_dofs, n_dofs), dtype=float)

    for n1, n2 in elements:
        L = nodes[n2] - nodes[n1]

        ke = beam_element_stiffness(E, I, L)

        me = (density * A * L / 420.0) * np.array([
            [156.0, 22.0 * L, 54.0, -13.0 * L],
            [22.0 * L, 4.0 * L**2, 13.0 * L, -3.0 * L**2],
            [54.0, 13.0 * L, 156.0, -22.0 * L],
            [-13.0 * L, -3.0 * L**2, -22.0 * L, 4.0 * L**2]
        ], dtype=float)

        dofs = [
            2 * n1,
            2 * n1 + 1,
            2 * n2,
            2 * n2 + 1
        ]

        K[np.ix_(dofs, dofs)] += ke
        M[np.ix_(dofs, dofs)] += me

    return K, M
