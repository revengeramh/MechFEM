
import numpy as np
import pytest

from mechfem.beam import (
    beam_element_stiffness,
    assemble_beam_matrices,
)


def test_beam_matrices_single_element():
    E = 200e9
    I = 8e-6
    A = 0.01
    density = 7850.0
    L = 2.0

    nodes = np.array([0.0, L])
    elements = np.array([[0, 1]])

    K, M = assemble_beam_matrices(
        nodes, elements, E, I, A, density
    )

    expected_K = beam_element_stiffness(E, I, L)

    expected_M = (density * A * L / 420.0) * np.array([
        [156, 22 * L, 54, -13 * L],
        [22 * L, 4 * L**2, 13 * L, -3 * L**2],
        [54, 13 * L, 156, -22 * L],
        [-13 * L, -3 * L**2, -22 * L, 4 * L**2]
    ])

    assert K.shape == (4, 4)
    assert M.shape == (4, 4)
    assert np.allclose(K, expected_K)
    assert np.allclose(M, expected_M)


def test_beam_matrices_are_symmetric():
    nodes = np.array([0.0, 1.0, 2.0])
    elements = np.array([[0, 1], [1, 2]])

    K, M = assemble_beam_matrices(
        nodes, elements, 200e9, 8e-6, 0.01, 7850
    )

    assert np.allclose(K, K.T)
    assert np.allclose(M, M.T)


def test_beam_mass_matrix_positive_definite():
    nodes = np.array([0.0, 2.0])
    elements = np.array([[0, 1]])

    _, M = assemble_beam_matrices(
        nodes, elements, 200e9, 8e-6, 0.01, 7850
    )

    assert np.all(np.linalg.eigvalsh(M) > 0)


def test_beam_matrices_invalid_material():
    nodes = np.array([0.0, 1.0])
    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        assemble_beam_matrices(
            nodes, elements, 200e9, 8e-6, 0.01, 0
        )


def test_beam_matrices_invalid_connectivity():
    nodes = np.array([0.0, 1.0])
    elements = np.array([[0, 2]])

    with pytest.raises(ValueError):
        assemble_beam_matrices(
            nodes, elements, 200e9, 8e-6, 0.01, 7850
        )


def test_beam_matrices_empty_elements():
    nodes = np.array([0.0, 1.0])
    elements = np.empty((0, 2), dtype=int)

    with pytest.raises(ValueError):
        assemble_beam_matrices(
            nodes, elements, 200e9, 8e-6, 0.01, 7850
        )
