
import numpy as np
import pytest

from mechfem.modal import modal_analysis


def test_modal_analysis_basic():
    K = np.diag([100.0, 800.0])
    M = np.diag([1.0, 2.0])

    results = modal_analysis(K, M, [])

    expected = np.array([10.0, 20.0]) / (2.0 * np.pi)

    assert np.allclose(
        results["natural_frequencies"],
        expected
    )


def test_modal_analysis_with_fixed_dof():
    K = np.diag([999.0, 100.0, 800.0])
    M = np.diag([1.0, 1.0, 2.0])

    results = modal_analysis(K, M, [0])

    assert results["mode_shapes"].shape == (3, 2)
    assert np.allclose(results["mode_shapes"][0, :], 0.0)

    expected = np.array([10.0, 20.0]) / (2.0 * np.pi)

    assert np.allclose(
        results["natural_frequencies"],
        expected
    )


def test_modal_analysis_num_modes():
    K = np.diag([100.0, 800.0])
    M = np.diag([1.0, 2.0])

    results = modal_analysis(K, M, [], num_modes=1)

    assert len(results["natural_frequencies"]) == 1
    assert results["mode_shapes"].shape == (2, 1)


def test_mode_shapes_are_mass_orthonormal():
    K = np.diag([100.0, 800.0])
    M = np.diag([1.0, 2.0])

    results = modal_analysis(K, M, [])

    phi = results["mode_shapes"]

    assert np.allclose(
        phi.T @ M @ phi,
        np.eye(2)
    )


def test_modal_analysis_invalid_stiffness_shape():
    K = np.array([[1.0, 2.0]])
    M = np.eye(2)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [])


def test_modal_analysis_invalid_mass_shape():
    K = np.eye(2)
    M = np.eye(3)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [])


def test_modal_analysis_non_symmetric_matrix():
    K = np.array([
        [100.0, 10.0],
        [0.0, 200.0]
    ])
    M = np.eye(2)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [])


def test_modal_analysis_non_positive_mass_matrix():
    K = np.eye(2)
    M = np.diag([1.0, 0.0])

    with pytest.raises(ValueError):
        modal_analysis(K, M, [])


def test_modal_analysis_unstable_structure():
    K = np.diag([100.0, -10.0])
    M = np.eye(2)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [])


def test_modal_analysis_invalid_num_modes():
    K = np.eye(2)
    M = np.eye(2)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [], num_modes=3)


def test_modal_analysis_invalid_fixed_dof():
    K = np.eye(2)
    M = np.eye(2)

    with pytest.raises(ValueError):
        modal_analysis(K, M, [5])
