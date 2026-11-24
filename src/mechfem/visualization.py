
import numpy as np


def plot_beam_deformation(
    nodes,
    displacements,
    scale=1.0,
    save_path=None,
    show=True
):
    """
    Plot the original and deformed shape of an Euler-Bernoulli beam.

    Parameters
    ----------
    nodes : array-like
        Beam nodal coordinates.
    displacements : array-like
        Global displacement vector [v1, theta1, ...].
    scale : float
        Deformation visualization scale factor.
    save_path : str or None
        Optional path to save the figure.
    show : bool
        Whether to display the figure.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """

    import matplotlib.pyplot as plt

    nodes = np.asarray(nodes, dtype=float)
    displacements = np.asarray(displacements, dtype=float)

    if nodes.ndim != 1 or len(nodes) < 2:
        raise ValueError("Nodes must be a 1D array with at least two entries.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if np.any(np.diff(nodes) <= 0):
        raise ValueError("Nodes must be strictly increasing.")

    if displacements.shape != (2 * len(nodes),):
        raise ValueError("Invalid displacement vector size.")

    if not np.all(np.isfinite(displacements)):
        raise ValueError("Displacements must be finite.")

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Scale must be positive and finite.")

    x_original = []
    y_original = []
    x_deformed = []
    y_deformed = []

    for i in range(len(nodes) - 1):
        x1 = nodes[i]
        x2 = nodes[i + 1]
        L = x2 - x1

        v1 = displacements[2 * i]
        theta1 = displacements[2 * i + 1]
        v2 = displacements[2 * (i + 1)]
        theta2 = displacements[2 * (i + 1) + 1]

        xi = np.linspace(0.0, 1.0, 30)

        N1 = 1 - 3*xi**2 + 2*xi**3
        N2 = L * (xi - 2*xi**2 + xi**3)
        N3 = 3*xi**2 - 2*xi**3
        N4 = L * (-xi**2 + xi**3)

        v = N1*v1 + N2*theta1 + N3*v2 + N4*theta2

        x = x1 + xi * L

        x_original.extend(x)
        y_original.extend(np.zeros_like(x))

        x_deformed.extend(x)
        y_deformed.extend(scale * v)

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(
        x_original,
        y_original,
        "--",
        label="Original beam"
    )

    ax.plot(
        x_deformed,
        y_deformed,
        label="Deformed beam"
    )

    ax.set_xlabel("Position (m)")
    ax.set_ylabel("Transverse displacement (m)")
    ax.set_title("MechFEM - Beam Deformation")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    return fig, ax

def plot_truss_axial_forces(
    nodes,
    elements,
    axial_forces,
    save_path=None,
    show=True
):
    """Plot 2D truss members colored by axial force."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    axial_forces = np.asarray(axial_forces, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(elements) == 0 or axial_forces.shape != (len(elements),):
        raise ValueError("Invalid elements or axial forces.")

    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(axial_forces)):
        raise ValueError("Inputs must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    fig, ax = plt.subplots(figsize=(8, 6))

    max_force = max(np.max(np.abs(axial_forces)), 1e-15)
    norm = Normalize(vmin=-max_force, vmax=max_force)
    cmap = plt.get_cmap("coolwarm")

    for (n1, n2), force in zip(elements, axial_forces):
        x = [nodes[n1, 0], nodes[n2, 0]]
        y = [nodes[n1, 1], nodes[n2, 1]]

        ax.plot(
            x, y,
            color=cmap(norm(force)),
            linewidth=3
        )

        midpoint = (np.mean(x), np.mean(y))
        ax.annotate(
            f"{force:.2f} N",
            midpoint,
            fontsize=8
        )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, label="Axial Force (N)")

    ax.set_title("MechFEM - Truss Axial Forces")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.axis("equal")
    ax.grid(True)

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    return fig, ax

def plot_truss_stresses(
    nodes,
    elements,
    axial_stresses,
    save_path=None,
    show=True
):
    """Plot 2D truss members colored by axial stress."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    axial_stresses = np.asarray(axial_stresses, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if len(elements) == 0 or axial_stresses.shape != (len(elements),):
        raise ValueError("Invalid elements or axial stresses.")

    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(axial_stresses)):
        raise ValueError("Inputs must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if np.any(elements[:, 0] == elements[:, 1]):
        raise ValueError("An element cannot connect a node to itself.")

    fig, ax = plt.subplots(figsize=(8, 6))

    max_stress = max(np.max(np.abs(axial_stresses)), 1e-15)
    norm = Normalize(vmin=-max_stress, vmax=max_stress)
    cmap = plt.get_cmap("coolwarm")

    for (n1, n2), stress in zip(elements, axial_stresses):
        x = [nodes[n1, 0], nodes[n2, 0]]
        y = [nodes[n1, 1], nodes[n2, 1]]

        ax.plot(
            x, y,
            color=cmap(norm(stress)),
            linewidth=3
        )

        midpoint = (np.mean(x), np.mean(y))
        ax.annotate(
            f"{stress / 1e6:.2f} MPa",
            midpoint,
            fontsize=8
        )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, label="Axial Stress (Pa)")

    ax.set_title("MechFEM - Truss Axial Stresses")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.axis("equal")
    ax.grid(True)

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    return fig, ax


def plot_frame_deformation(
    nodes,
    elements,
    displacements,
    scale=1.0,
    save_path=None,
    show=True
):
    """
    Plot the original and deformed 2D frame configuration.

    Each node has three DOFs: [ux, uy, theta].
    """

    import numpy as np
    import matplotlib.pyplot as plt

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    displacements = np.asarray(displacements, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("elements must have shape (n_elements, 2).")

    if displacements.shape != (3 * len(nodes),):
        raise ValueError("Invalid displacement vector size.")

    if not np.all(np.isfinite(displacements)):
        raise ValueError("Displacements must be finite.")

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Scale must be a positive finite number.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    if save_path is not None and not isinstance(save_path, (str, bytes)):
        raise ValueError("save_path must be a string or None.")

    deformed = nodes + scale * displacements.reshape(-1, 3)[:, :2]

    fig, ax = plt.subplots(figsize=(9, 6))

    for n1, n2 in elements:
        ax.plot(
            nodes[[n1, n2], 0],
            nodes[[n1, n2], 1],
            "k--",
            linewidth=1.5
        )

        ax.plot(
            deformed[[n1, n2], 0],
            deformed[[n1, n2], 1],
            "r-",
            linewidth=2
        )

    ax.plot([], [], "k--", label="Original frame")
    ax.plot([], [], "r-", label="Deformed frame")

    ax.set_xlabel("X coordinate (m)")
    ax.set_ylabel("Y coordinate (m)")
    ax.set_title("2D Frame Deformation")
    ax.axis("equal")
    ax.grid(True)
    ax.legend()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig, ax
