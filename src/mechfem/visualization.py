
import numpy as np
import matplotlib.pyplot as plt


def plot_beam_deformation(
    nodes,
    elements,
    displacements,
    scale=1.0,
    save_path=None,
    show=True
):
    """Plot Euler-Bernoulli beam deformation."""

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    displacements = np.asarray(displacements, dtype=float)

    if nodes.ndim != 1 or len(nodes) < 2:
        raise ValueError("Invalid beam nodes.")

    if displacements.shape != (2 * len(nodes),):
        raise ValueError("Invalid displacement vector.")

    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(displacements)):
        raise ValueError("Input values must be finite.")

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Scale must be positive and finite.")

    fig, ax = plt.subplots(figsize=(9, 5))

    for n1, n2 in elements:
        if not (0 <= n1 < len(nodes) and 0 <= n2 < len(nodes)):
            plt.close(fig)
            raise ValueError("Invalid element connectivity.")

        x1, x2 = nodes[n1], nodes[n2]
        v1, theta1 = displacements[2*n1:2*n1+2]
        v2, theta2 = displacements[2*n2:2*n2+2]

        L = x2 - x1

        if L <= 0:
            plt.close(fig)
            raise ValueError("Element length must be positive.")

        xi = np.linspace(0, 1, 30)

        N1 = 1 - 3*xi**2 + 2*xi**3
        N2 = L*(xi - 2*xi**2 + xi**3)
        N3 = 3*xi**2 - 2*xi**3
        N4 = L*(-xi**2 + xi**3)

        v = N1*v1 + N2*theta1 + N3*v2 + N4*theta2

        x = np.linspace(x1, x2, 30)

        ax.plot([x1, x2], [0, 0], "k--")
        ax.plot(x, scale*v, "r-")

    ax.set_xlabel("X coordinate (m)")
    ax.set_ylabel("Deflection (m)")
    ax.set_title("Beam Deformation")
    ax.grid(True)
    ax.legend(["Original beam", "Deformed beam"])

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

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
    """Plot axial force distribution in a 2D truss."""

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    axial_forces = np.asarray(axial_forces, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if axial_forces.shape != (len(elements),):
        raise ValueError("Invalid axial force vector.")

    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(axial_forces)):
        raise ValueError("Input values must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (n1, n2) in enumerate(elements):
        if n1 == n2:
            plt.close(fig)
            raise ValueError("Invalid zero-length element.")

        x = nodes[[n1, n2], 0]
        y = nodes[[n1, n2], 1]

        ax.plot(x, y, "o-", linewidth=3)

        midpoint_x = np.mean(x)
        midpoint_y = np.mean(y)

        ax.text(
            midpoint_x,
            midpoint_y,
            f"{axial_forces[i]:.2f} N",
            fontsize=9
        )

    ax.set_xlabel("X coordinate (m)")
    ax.set_ylabel("Y coordinate (m)")
    ax.set_title("Truss Axial Force Distribution")
    ax.axis("equal")
    ax.grid(True)

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

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
    """Plot axial stress distribution in a 2D truss."""

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    axial_stresses = np.asarray(axial_stresses, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if axial_stresses.shape != (len(elements),):
        raise ValueError("Invalid stress vector.")

    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(axial_stresses)):
        raise ValueError("Input values must be finite.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (n1, n2) in enumerate(elements):
        if n1 == n2:
            plt.close(fig)
            raise ValueError("Invalid zero-length element.")

        x = nodes[[n1, n2], 0]
        y = nodes[[n1, n2], 1]

        ax.plot(x, y, "o-", linewidth=3)

        ax.text(
            np.mean(x),
            np.mean(y),
            f"{axial_stresses[i]/1e6:.2f} MPa",
            fontsize=9
        )

    ax.set_xlabel("X coordinate (m)")
    ax.set_ylabel("Y coordinate (m)")
    ax.set_title("Truss Stress Distribution")
    ax.axis("equal")
    ax.grid(True)

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig, ax


def plot_frame_deformation(
    nodes,
    elements,
    displacements,
    scale=1.0,
    save_path=None,
    show=True,
    points_per_element=30
):
    """
    Plot original and deformed 2D frame using
    Euler-Bernoulli beam Hermite interpolation.

    Each node has three DOFs:
    [ux, uy, theta]
    """

    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    displacements = np.asarray(displacements, dtype=float)

    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("Nodes must have shape (n_nodes, 2).")

    if elements.ndim != 2 or elements.shape[1] != 2:
        raise ValueError("Elements must have shape (n_elements, 2).")

    if displacements.shape != (3 * len(nodes),):
        raise ValueError("Invalid displacement vector size.")

    if len(elements) == 0:
        raise ValueError("At least one element is required.")

    if not np.all(np.isfinite(nodes)):
        raise ValueError("Node coordinates must be finite.")

    if not np.all(np.isfinite(displacements)):
        raise ValueError("Displacements must be finite.")

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Scale must be positive and finite.")

    if not isinstance(points_per_element, (int, np.integer)):
        raise ValueError("points_per_element must be an integer.")

    if points_per_element < 2:
        raise ValueError("At least two points per element are required.")

    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("Invalid element connectivity.")

    fig, ax = plt.subplots(figsize=(9, 6))

    xi = np.linspace(0.0, 1.0, points_per_element)

    for n1, n2 in elements:

        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        dx = x2 - x1
        dy = y2 - y1

        L = np.hypot(dx, dy)

        if L <= 0:
            plt.close(fig)
            raise ValueError("Element length must be positive.")

        c = dx / L
        s = dy / L

        dofs = [
            3*n1, 3*n1+1, 3*n1+2,
            3*n2, 3*n2+1, 3*n2+2
        ]

        ue = displacements[dofs]

        # Transform global displacements to local coordinates.
        u1 = c*ue[0] + s*ue[1]
        v1 = -s*ue[0] + c*ue[1]
        theta1 = ue[2]

        u2 = c*ue[3] + s*ue[4]
        v2 = -s*ue[3] + c*ue[4]
        theta2 = ue[5]

        x_local = xi * L

        # Cubic Hermite interpolation.
        N1 = 1 - 3*xi**2 + 2*xi**3
        N2 = L*(xi - 2*xi**2 + xi**3)
        N3 = 3*xi**2 - 2*xi**3
        N4 = L*(-xi**2 + xi**3)

        v_local = (
            N1*v1 + N2*theta1 +
            N3*v2 + N4*theta2
        )

        # Linear interpolation of axial displacement.
        u_local = (1-xi)*u1 + xi*u2

        original_x = x1 + c*x_local
        original_y = y1 + s*x_local

        deformed_x = (
            original_x + scale*(c*u_local - s*v_local)
        )

        deformed_y = (
            original_y + scale*(s*u_local + c*v_local)
        )

        ax.plot(
            [x1, x2],
            [y1, y2],
            "k--",
            linewidth=1.5
        )

        ax.plot(
            deformed_x,
            deformed_y,
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
