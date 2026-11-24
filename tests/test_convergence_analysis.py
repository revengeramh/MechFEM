
import numpy as np
import pytest

from mechfem.convergence import beam_mesh_convergence


def test_beam_mesh_convergence_returns_expected_keys():
    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[2, 4, 8, 16],
        analytical_value=-1.0
    )

    expected_keys = {
        "element_counts",
        "displacements",
        "absolute_errors",
        "relative_errors",
        "convergence_rates",
    }

    assert set(result.keys()) == expected_keys


def test_mesh_counts_are_sorted():
    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[8, 2, 4]
    )

    assert np.array_equal(
        result["element_counts"],
        [2, 4, 8]
    )


def test_displacements_are_finite():
    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[2, 4, 8]
    )

    assert np.all(np.isfinite(result["displacements"]))


def test_error_decreases_with_mesh_refinement():
    length = 4.0
    E = 200e9
    I = 8e-6
    load = -1000.0

    analytical = load * length**3 / (48 * E * I)

    result = beam_mesh_convergence(
        length=length,
        E=E,
        I=I,
        load=load,
        element_counts=[2, 4, 8, 16],
        analytical_value=analytical
    )

    errors = result["absolute_errors"]

    assert errors[-1] <= errors[0]
    assert np.all(np.isfinite(errors))


def test_relative_errors_are_nonnegative():
    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[2, 4, 8],
        analytical_value=-1.0
    )

    assert np.all(result["relative_errors"] >= 0)


def test_convergence_rates_shape():
    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[2, 4, 8, 16],
        analytical_value=-1.0
    )

    rates = result["convergence_rates"]

    assert len(rates) == 4
    assert np.isnan(rates[0])


@pytest.mark.parametrize("length", [0, -1, np.nan, np.inf])
def test_invalid_length(length):
    with pytest.raises(ValueError):
        beam_mesh_convergence(
            length=length,
            E=200e9,
            I=8e-6,
            load=-1000.0,
            element_counts=[2, 4]
        )


@pytest.mark.parametrize("E", [0, -200e9, np.nan, np.inf])
def test_invalid_young_modulus(E):
    with pytest.raises(ValueError):
        beam_mesh_convergence(
            length=4.0,
            E=E,
            I=8e-6,
            load=-1000.0,
            element_counts=[2, 4]
        )


@pytest.mark.parametrize("I", [0, -8e-6, np.nan, np.inf])
def test_invalid_moment_of_inertia(I):
    with pytest.raises(ValueError):
        beam_mesh_convergence(
            length=4.0,
            E=200e9,
            I=I,
            load=-1000.0,
            element_counts=[2, 4]
        )


def test_invalid_element_counts():
    with pytest.raises(ValueError):
        beam_mesh_convergence(
            length=4.0,
            E=200e9,
            I=8e-6,
            load=-1000.0,
            element_counts=[2]
        )


def test_element_counts_must_be_unique():
    with pytest.raises(ValueError):
        beam_mesh_convergence(
            length=4.0,
            E=200e9,
            I=8e-6,
            load=-1000.0,
            element_counts=[2, 2, 4]
        )


def test_analytical_errors_are_returned():
    analytical = -0.01

    result = beam_mesh_convergence(
        length=4.0,
        E=200e9,
        I=8e-6,
        load=-1000.0,
        element_counts=[2, 4, 8],
        analytical_value=analytical
    )

    expected = np.abs(
        result["displacements"] - analytical
    )

    assert np.allclose(
        result["absolute_errors"],
        expected
    )
