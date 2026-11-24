
import numpy as np
import pytest

from mechfem.truss import assemble_truss_matrices


def test_horizontal_truss_matrix_shapes():
    nodes = [[0.0, 0.0], [2.0, 0.0]]
    elements = [[0, 1]]

    K, M = assemble_truss_matrices(
        nodes, elements, 200e9, 0.01, 7850
    )

    assert K.shape == (4, 4)
    assert M.shape == (4, 4)


def test_stiffness_matrix_symmetry():
    nodes = [[0.0, 0.0], [2.0, 0.0]]
    elements = [[0, 1]]

    K, M = assemble_truss_matrices(
        nodes, elements, 200e9, 0.01, 7850
    )

    assert np.allclose(K, K.T)
    assert np.allclose(M, M.T)


def test_horizontal_truss_axial_stiffness():
    E = 200e9
    A = 0.01
    L = 2.0

    nodes = [[0.0, 0.0], [L, 0.0]]
    elements = [[0, 1]]

    K, M = assemble_truss_matrices(
        nodes, elements, E, A, 7850
    )

    expected = E * A / L

    assert np.isclose(K[0, 0], expected)
    assert np.isclose(K[0, 2], -expected)


def test_consistent_mass_matrix():
    density = 7850
    area = 0.01
    length = 2.0

    nodes = [[0.0, 0.0], [length, 0.0]]
    elements = [[0, 1]]

    K, M = assemble_truss_matrices(
        nodes, elements, 200e9, area, density
    )

    factor = density * area * length / 6.0

    assert np.isclose(M[0, 0], 2.0 * factor)
    assert np.isclose(M[0, 2], factor)
    assert np.isclose(M[1, 1], 2.0 * factor)
    assert np.isclose(M[1, 3], factor)


def test_diagonal_truss_matrices_are_symmetric():
    nodes = [[0.0, 0.0], [3.0, 4.0]]
    elements = [[0, 1]]

    K, M = assemble_truss_matrices(
        nodes, elements, 200e9, 0.01, 7850
    )

    assert np.allclose(K, K.T)
    assert np.allclose(M, M.T)


def test_multiple_elements_are_assembled():
    nodes = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0]
    ]

    elements = [[0, 1], [1, 2]]

    K, M = assemble_truss_matrices(
        nodes, elements, 200e9, 0.01, 7850
    )

    assert K.shape == (6, 6)
    assert M.shape == (6, 6)

    assert np.allclose(K, K.T)
    assert np.allclose(M, M.T)


@pytest.mark.parametrize(
    "E, A, density",
    [
        (0, 0.01, 7850),
        (200e9, 0, 7850),
        (200e9, 0.01, 0),
        (-200e9, 0.01, 7850),
        (200e9, -0.01, 7850),
        (200e9, 0.01, -7850),
    ]
)
def test_invalid_material_properties(E, A, density):
    nodes = [[0.0, 0.0], [1.0, 0.0]]
    elements = [[0, 1]]

    with pytest.raises(ValueError):
        assemble_truss_matrices(
            nodes, elements, E, A, density
        )


def test_invalid_node_shape():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            [[0.0, 0.0, 0.0]],
            [[0, 0]],
            200e9,
            0.01,
            7850
        )


def test_invalid_connectivity():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            [[0.0, 0.0], [1.0, 0.0]],
            [[0, 2]],
            200e9,
            0.01,
            7850
        )


def test_zero_length_element():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            [[0.0, 0.0], [0.0, 0.0]],
            [[0, 1]],
            200e9,
            0.01,
            7850
        )


def test_empty_elements():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            [[0.0, 0.0], [1.0, 0.0]],
            [],
            200e9,
            0.01,
            7850
        )
