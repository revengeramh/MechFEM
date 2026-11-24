
import numpy as np
import pytest

from mechfem.beam import assemble_beam_matrices
from mechfem.modal import analyze_beam_modal


def test_beam_modal_single_element_analytical():
    E = 200e9
    I = 8e-6
    A = 0.01
    density = 7850.0
    L = 2.0

    nodes = np.array([0.0, L])
    elements = np.array([[0, 1]])
    fixed_dofs = [0, 1]

    K, M = assemble_beam_matrices(
        nodes, elements, E, I, A, density
    )

    Kff = K[2:4, 2:4]
    Mff = M[2:4, 2:4]

    # Independent analytical solution of the
    # 2x2 generalized eigenvalue characteristic equation.
    a = np.linalg.det(Mff)

    b = -(
        Kff[0, 0] * Mff[1, 1]
        + Kff[1, 1] * Mff[0, 0]
        - Kff[0, 1] * Mff[1, 0]
        - Kff[1, 0] * Mff[0, 1]
    )

    c = np.linalg.det(Kff)

    expected_eigenvalues = np.sort(
        np.roots([a, b, c]).real
    )

    expected_frequencies = (
        np.sqrt(expected_eigenvalues) / (2.0 * np.pi)
    )

    results = analyze_beam_modal(
        nodes,
        elements,
        E,
        I,
        A,
        density,
        fixed_dofs,
    )

    assert np.allclose(
        results["natural_frequencies"],
        expected_frequencies,
    )

    assert np.allclose(
        results["eigenvalues"],
        expected_eigenvalues,
    )


def test_beam_modal_mode_shapes_respect_boundary_conditions():
    nodes = np.array([0.0, 1.0, 2.0])
    elements = np.array([
        [0, 1],
        [1, 2],
    ])

    results = analyze_beam_modal(
        nodes,
        elements,
        200e9,
        8e-6,
        0.01,
        7850.0,
        [0, 1],
    )

    mode_shapes = results["mode_shapes"]

    assert mode_shapes.shape[0] == 6
    assert np.allclose(mode_shapes[0, :], 0.0)
    assert np.allclose(mode_shapes[1, :], 0.0)


def test_beam_modal_num_modes():
    nodes = np.array([0.0, 1.0, 2.0])
    elements = np.array([
        [0, 1],
        [1, 2],
    ])

    results = analyze_beam_modal(
        nodes,
        elements,
        200e9,
        8e-6,
        0.01,
        7850.0,
        [0, 1],
        num_modes=1,
    )

    assert len(results["natural_frequencies"]) == 1
    assert results["mode_shapes"].shape == (6, 1)


def test_beam_modal_invalid_num_modes():
    nodes = np.array([0.0, 1.0])
    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        analyze_beam_modal(
            nodes,
            elements,
            200e9,
            8e-6,
            0.01,
            7850.0,
            [0, 1],
            num_modes=10,
        )
