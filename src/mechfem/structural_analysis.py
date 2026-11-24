
import numpy as np


def analyze_truss(nodes, elements, E, A, displacements):
    """
    Perform axial stress and strain analysis for a 2D truss.

    Parameters
    ----------
    nodes : array_like, shape (n_nodes, 2)
        Global nodal coordinates.
    elements : array_like, shape (n_elements, 2)
        Element connectivity (zero-based node indices).
    E : float or array_like
        Young's modulus (Pa), uniform or per element.
    A : float or array_like
        Cross-sectional area (m^2), uniform or per element.
    displacements : array_like, shape (2*n_nodes,)
        Global nodal displacement vector (m).

    Returns
    -------
    dict
        Element lengths, axial elongations, strains,
        stresses, axial forces, and strain energies.

    Sign convention
    ---------------
    Positive force/stress/strain: tension.
    Negative force/stress/strain: compression.
    """

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    displacements = np.asarray(displacements, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("nodes must have shape (n_nodes, 2).")

    if len(nodes) < 2:
        raise ValueError("At least two nodes are required.")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("elements must have shape (n_elements, 2).")

    if len(elements) == 0:
        raise ValueError("At least one element is required.")

    if displacements.shape != (2 * len(nodes),):
        raise ValueError("Invalid displacement vector size.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if not np.all(np.isfinite(displacements)):
        raise ValueError("Displacements must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    n_elements = len(elements)

    E = np.asarray(E, dtype=float)
    A = np.asarray(A, dtype=float)

    if E.ndim == 0:
        E = np.full(n_elements, float(E))

    if A.ndim == 0:
        A = np.full(n_elements, float(A))

    if E.shape != (n_elements,) or A.shape != (n_elements,):
        raise ValueError(
            "E and A must be scalars or arrays with one value per element."
        )

    if not np.all(np.isfinite(E)) or np.any(E <= 0):
        raise ValueError("Young's modulus must be positive and finite.")

    if not np.all(np.isfinite(A)) or np.any(A <= 0):
        raise ValueError("Cross-sectional areas must be positive and finite.")

    lengths = np.zeros(n_elements)
    elongations = np.zeros(n_elements)
    strains = np.zeros(n_elements)
    stresses = np.zeros(n_elements)
    axial_forces = np.zeros(n_elements)
    strain_energies = np.zeros(n_elements)

    for i, (n1, n2) in enumerate(elements):

        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        dx = x2 - x1
        dy = y2 - y1

        L = np.hypot(dx, dy)

        if L <= 0:
            raise ValueError("Element length must be positive.")

        c = dx / L
        s = dy / L

        element_dofs = np.array([
            2*n1,
            2*n1 + 1,
            2*n2,
            2*n2 + 1
        ])

        ue = displacements[element_dofs]

        # Axial deformation projected onto the element axis.
        B = np.array([-c, -s, c, s])

        delta = B @ ue

        # Engineering axial strain.
        strain = delta / L

        # Linear elastic Hooke's law.
        stress = E[i] * strain

        # Axial force.
        force = stress * A[i]

        # Elastic strain energy.
        energy = 0.5 * stress * strain * A[i] * L

        lengths[i] = L
        elongations[i] = delta
        strains[i] = strain
        stresses[i] = stress
        axial_forces[i] = force
        strain_energies[i] = energy

    return {
        "lengths": lengths,
        "elongations": elongations,
        "strains": strains,
        "stresses": stresses,
        "axial_forces": axial_forces,
        "strain_energies": strain_energies
    }
