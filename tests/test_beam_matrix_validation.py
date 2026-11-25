
import numpy as np
import pytest

from mechfem.beam import assemble_beam_matrices


E = 200e9
I = 8e-6
A = 0.01
DENSITY = 7850.0


def test_rejects_fractional_element_indices():
    with pytest.raises(ValueError, match="integers"):
        assemble_beam_matrices(
            [0.0, 1.0, 2.0],
            [[0, 1], [1.5, 2]],
            E, I, A, DENSITY
        )


def test_rejects_non_finite_element_indices():
    with pytest.raises(ValueError):
        assemble_beam_matrices(
            [0.0, 1.0],
            [[0, np.inf]],
            E, I, A, DENSITY
        )


def test_rejects_duplicate_elements():
    with pytest.raises(ValueError, match="Duplicate"):
        assemble_beam_matrices(
            [0.0, 1.0],
            [[0, 1], [0, 1]],
            E, I, A, DENSITY
        )


def test_rejects_reversed_connectivity():
    with pytest.raises(ValueError, match="increasing"):
        assemble_beam_matrices(
            [0.0, 1.0],
            [[1, 0]],
            E, I, A, DENSITY
        )


def test_rejects_non_finite_nodes():
    with pytest.raises(ValueError):
        assemble_beam_matrices(
            [0.0, np.inf],
            [[0, 1]],
            E, I, A, DENSITY
        )


def test_rejects_invalid_density():
    with pytest.raises(ValueError):
        assemble_beam_matrices(
            [0.0, 1.0],
            [[0, 1]],
            E, I, A, 0.0
        )


def test_global_matrices_are_symmetric():
    K, M = assemble_beam_matrices(
        [0.0, 1.0, 2.0],
        [[0, 1], [1, 2]],
        E, I, A, DENSITY
    )

    assert np.allclose(K, K.T)
    assert np.allclose(M, M.T)


def test_global_mass_matrix_is_positive_definite():
    _, M = assemble_beam_matrices(
        [0.0, 1.0, 2.0],
        [[0, 1], [1, 2]],
        E, I, A, DENSITY
    )

    assert np.all(np.linalg.eigvalsh(M) > 0)


def test_global_matrices_have_correct_dimensions():
    K, M = assemble_beam_matrices(
        [0.0, 1.0, 2.0],
        [[0, 1], [1, 2]],
        E, I, A, DENSITY
    )

    assert K.shape == (6, 6)
    assert M.shape == (6, 6)
