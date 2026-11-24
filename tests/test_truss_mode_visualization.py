import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from mechfem.visualization import plot_truss_mode_shape


def test_truss_mode_shape_plot():
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

    mode_shape = np.array([
        0.0, 0.0,
        0.01, 0.02,
        -0.01, 0.03
    ])

    fig, ax = plot_truss_mode_shape(
        nodes,
        elements,
        mode_shape,
        scale=10.0,
        show=False
    )

    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "MechFEM - Truss Mode Shape"
    assert len(ax.lines) == 8

    plt.close(fig)


def test_invalid_mode_shape():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0]
    ])

    elements = np.array([[0, 1]])
    mode_shape = np.array([0.0, 0.0])

    with pytest.raises(ValueError):
        plot_truss_mode_shape(
            nodes, elements, mode_shape, show=False
        )


def test_invalid_mode_shape_scale():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0]
    ])

    elements = np.array([[0, 1]])
    mode_shape = np.zeros(4)

    with pytest.raises(ValueError):
        plot_truss_mode_shape(
            nodes,
            elements,
            mode_shape,
            scale=-1.0,
            show=False
        )


def test_invalid_mode_shape_connectivity():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0]
    ])

    elements = np.array([[0, 2]])
    mode_shape = np.zeros(4)

    with pytest.raises(ValueError):
        plot_truss_mode_shape(
            nodes, elements, mode_shape, show=False
        )
