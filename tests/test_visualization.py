
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from mechfem.visualization import (
    plot_beam_deformation,
    plot_truss_axial_forces,
    plot_truss_stresses,
)


def test_beam_deformation_plot():
    nodes = np.array([0.0, 1.0, 2.0])
    displacements = np.array([
        0.0, 0.0,
        -0.01, -0.02,
        -0.03, 0.0
    ])

    fig, ax = plot_beam_deformation(
        nodes, displacements, scale=10.0, show=False
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
        plot_beam_deformation(nodes, displacements, show=False)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError")


def test_truss_axial_forces_plot():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 1.0]
    ])
    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2]
    ])
    forces = np.array([1000.0, -500.0, 250.0])

    fig, ax = plot_truss_axial_forces(
        nodes, elements, forces, show=False
    )

    assert fig is not None
    assert len(ax.lines) == 3
    assert ax.get_title() == "MechFEM - Truss Axial Forces"
    plt.close(fig)


def test_truss_stresses_plot():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 1.0]
    ])
    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2]
    ])
    stresses = np.array([10e6, -5e6, 2.5e6])

    fig, ax = plot_truss_stresses(
        nodes, elements, stresses, show=False
    )

    assert fig is not None
    assert len(ax.lines) == 3
    assert ax.get_title() == "MechFEM - Truss Axial Stresses"
    plt.close(fig)
