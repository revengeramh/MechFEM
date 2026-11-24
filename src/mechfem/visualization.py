
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
