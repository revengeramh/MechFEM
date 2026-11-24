
import numpy as np
import pytest

from mechfem.beam import beam_element_stiffness, solve_beam


def test_beam_element_stiffness_is_symmetric():
    k = beam_element_stiffness(E=200e9, I=8e-6, L=2.0)

    assert k.shape == (4, 4)
    np.testing.assert_allclose(k, k.T)


def test_cantilever_tip_load_matches_analytical_solution():
    length = 2.0
    E = 200e9
    I = 8e-6
    force = -1000.0

    displacement, reactions = solve_beam(
        nodes=[0.0, length],
        elements=[(0, 1)],
        E=E,
        I=I,
        loads=[0.0, 0.0, force, 0.0],
        fixed_dofs=[0, 1],
    )

    expected_deflection = force * length**3 / (3 * E * I)
    expected_rotation = force * length**2 / (2 * E * I)

    np.testing.assert_allclose(
        displacement[2], expected_deflection
    )

    np.testing.assert_allclose(
        displacement[3], expected_rotation
    )

    np.testing.assert_allclose(reactions[0], -force)
    np.testing.assert_allclose(reactions[1], -force * length)


def test_invalid_beam_properties_raise_value_error():
    with pytest.raises(ValueError):
        beam_element_stiffness(E=0, I=1, L=1)

