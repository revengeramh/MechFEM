
import numpy as np
import pytest

from mechfem.modal import analyze_truss_modal


def test_single_element_truss_natural_frequency():
    E = 200e9
    A = 0.01
    density = 7850
    L = 2.0

    nodes = [
        [0.0, 0.0],
        [L, 0.0]
    ]

    elements = [[0, 1]]

    # Fix both DOFs at the first node and
    # the transverse DOF at the second node.
    fixed_dofs = [0, 1, 3]

    results = analyze_truss_modal(
        nodes,
        elements,
        E,
        A,
        density,
        fixed_dofs
    )

    expected_frequency = (
        1.0 / (2.0 * L)
        * np.sqrt(E / density)
    )

    assert len(results["natural_frequencies"]) == 1

    assert np.isclose(
        results["natural_frequencies"][0],
        expected_frequency,
        rtol=1e-10
    )


def test_truss_mode_shape_has_zero_fixed_dofs():
    nodes = [
        [0.0, 0.0],
        [2.0, 0.0]
    ]

    elements = [[0, 1]]
    fixed_dofs = [0, 1, 3]

    results = analyze_truss_modal(
        nodes,
        elements,
        200e9,
        0.01,
        7850,
        fixed_dofs
    )

    mode_shapes = results["mode_shapes"]

    assert np.allclose(
        mode_shapes[fixed_dofs, :],
        0.0
    )


def test_truss_modal_analysis_returns_finite_values():
    nodes = [
        [0.0, 0.0],
        [2.0, 0.0]
    ]

    elements = [[0, 1]]

    results = analyze_truss_modal(
        nodes,
        elements,
        200e9,
        0.01,
        7850,
        [0, 1, 3]
    )

    assert np.all(
        np.isfinite(results["natural_frequencies"])
    )

    assert np.all(
        np.isfinite(results["mode_shapes"])
    )


def test_truss_modal_analysis_num_modes():
    nodes = [
        [0.0, 0.0],
        [2.0, 0.0]
    ]

    elements = [[0, 1]]

    results = analyze_truss_modal(
        nodes,
        elements,
        200e9,
        0.01,
        7850,
        [0, 1, 3],
        num_modes=1
    )

    assert results["mode_shapes"].shape == (4, 1)
    assert len(results["natural_frequencies"]) == 1


def test_truss_modal_analysis_invalid_material():
    nodes = [
        [0.0, 0.0],
        [2.0, 0.0]
    ]

    with pytest.raises(ValueError):
        analyze_truss_modal(
            nodes,
            [[0, 1]],
            200e9,
            0.01,
            0,
            [0, 1, 3]
        )
