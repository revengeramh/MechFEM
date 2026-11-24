
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from mechfem.visualization import plot_beam_deformation


def test_beam_deformation_plot():
    nodes = np.array([0.0, 1.0, 2.0])

    displacements = np.array([
        0.0, 0.0,
        -0.01, -0.02,
        -0.03, 0.0
    ])

    fig, ax = plot_beam_deformation(
        nodes,
        displacements,
        scale=10.0,
        show=False
    )

    assert fig is not None
    assert ax is not None

    assert len(ax.lines) == 2
    assert ax.get_title() == "MechFEM - Beam Deformation"

    plt.close(fig)


def test_invalid_displacement_vector():
    nodes = np.array([0.0, 1.0])

    displacements = np.array([0.0, 0.0])

    try:
        plot_beam_deformation(
            nodes,
            displacements,
            show=False
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError")


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
