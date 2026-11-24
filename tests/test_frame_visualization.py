
import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")

from mechfem.visualization import plot_frame_deformation


def test_frame_deformation_plot():
    nodes = np.array([
        [0.0, 0.0],
        [0.0, 3.0],
        [4.0, 3.0],
        [4.0, 0.0]
    ])

    elements = np.array([
        [0, 1],
        [1, 2],
        [2, 3]
    ])

    displacements = np.array([
        0.0, 0.0, 0.0,
        0.01, 0.02, 0.001,
        0.03, 0.01, 0.002,
        0.0, 0.0, 0.0
    ])

    fig, ax = plot_frame_deformation(
        nodes,
        elements,
        displacements,
        scale=10,
        show=False
    )

    assert fig is not None
    assert ax is not None
    assert len(ax.lines) == 9

    import matplotlib.pyplot as plt
    plt.close(fig)


def test_invalid_frame_displacement_size():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        plot_frame_deformation(
            nodes,
            elements,
            np.array([0.0, 0.0]),
            show=False
        )


def test_invalid_frame_scale():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        plot_frame_deformation(
            nodes,
            elements,
            np.zeros(6),
            scale=-1,
            show=False
        )
