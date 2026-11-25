
import numpy as np
import pytest

from mechfem.beam import solve_beam, beam_element_stiffness


E = 200e9
I = 8e-6


def test_beam_element_stiffness_is_symmetric():
    k = beam_element_stiffness(E, I, 2.0)
    assert np.allclose(k, k.T)


def test_beam_element_stiffness_rejects_invalid_values():
    with pytest.raises(ValueError):
        beam_element_stiffness(-E, I, 2.0)


def test_solve_beam_rejects_non_finite_nodes():
    with pytest.raises(ValueError):
        solve_beam(
            [0.0, np.inf],
            [(0, 1)],
            E,
            I,
            [0.0, 0.0, 0.0, 0.0],
            [0, 1],
        )


def test_solve_beam_rejects_duplicate_elements():
    with pytest.raises(ValueError):
        solve_beam(
            [0.0, 1.0, 2.0],
            [(0, 1), (0, 1)],
            E,
            I,
            [0.0] * 6,
            [0, 1],
        )


def test_solve_beam_rejects_invalid_connectivity():
    with pytest.raises(ValueError):
        solve_beam(
            [0.0, 1.0],
            [(1, 0)],
            E,
            I,
            [0.0] * 4,
            [0, 1],
        )


def test_solve_beam_rejects_duplicate_fixed_dofs():
    with pytest.raises(ValueError):
        solve_beam(
            [0.0, 1.0],
            [(0, 1)],
            E,
            I,
            [0.0] * 4,
            [0, 0],
        )


def test_solve_beam_detects_singular_system():
    with pytest.raises(ValueError, match="singular"):
        solve_beam(
            [0.0, 1.0],
            [(0, 1)],
            E,
            I,
            [0.0] * 4,
            [],
        )


def test_solve_beam_cantilever_equilibrium():
    load = -1000.0

    displacement, reactions = solve_beam(
        [0.0, 2.0],
        [(0, 1)],
        E,
        I,
        [0.0, 0.0, load, 0.0],
        [0, 1],
    )

    assert np.isclose(reactions[0], -load)
    assert np.isclose(reactions[1], -load * 2.0)
    assert np.all(np.isfinite(displacement))
