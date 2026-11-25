
import numpy as np

from mechfem.modal import modal_analysis


def test_generalized_eigenvalue_equation():
    K = np.array([
        [6.0, -2.0],
        [-2.0, 4.0]
    ])

    M = np.array([
        [2.0, 0.2],
        [0.2, 1.0]
    ])

    results = modal_analysis(K, M, [])

    phi = results["mode_shapes"]
    eigenvalues = results["eigenvalues"]

    for i in range(len(eigenvalues)):
        residual = K @ phi[:, i] - (
            eigenvalues[i] * M @ phi[:, i]
        )

        assert np.linalg.norm(residual) < 1e-10


def test_natural_frequencies_are_sorted():
    K = np.array([
        [6.0, -2.0],
        [-2.0, 4.0]
    ])

    M = np.array([
        [2.0, 0.2],
        [0.2, 1.0]
    ])

    results = modal_analysis(K, M, [])

    frequencies = results["natural_frequencies"]

    assert np.all(np.diff(frequencies) > 0)


def test_mode_shapes_are_stiffness_orthogonal():
    K = np.array([
        [6.0, -2.0],
        [-2.0, 4.0]
    ])

    M = np.array([
        [2.0, 0.2],
        [0.2, 1.0]
    ])

    results = modal_analysis(K, M, [])

    phi = results["mode_shapes"]

    assert np.allclose(
        phi.T @ K @ phi,
        np.diag(results["eigenvalues"]),
        atol=1e-10
    )


def test_fixed_dofs_remain_zero_in_mode_shapes():
    K = np.diag([100.0, 200.0, 300.0])
    M = np.diag([1.0, 2.0, 3.0])

    results = modal_analysis(K, M, [0])

    assert np.allclose(
        results["mode_shapes"][0, :],
        0.0
    )


def test_angular_and_natural_frequencies_are_consistent():
    K = np.diag([100.0, 400.0])
    M = np.eye(2)

    results = modal_analysis(K, M, [])

    assert np.allclose(
        results["natural_frequencies"],
        results["angular_frequencies"] / (2.0 * np.pi)
    )
