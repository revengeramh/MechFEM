import numpy as np
import pytest

from mechfem.truss import (
    solve_truss,
    calculate_truss_results,
    assemble_truss_matrices,
    truss_element_stiffness,
)


NODES = [[0.0, 0.0], [1.0, 0.0]]
ELEMENTS = [[0, 1]]
E = 200e9
A = 0.01


@pytest.mark.parametrize(
    "nodes",
    [
        [[0.0, 0.0]],
        [[0.0, 0.0, 1.0], [1.0, 0.0, 2.0]],
        [[0.0, 0.0], [np.nan, 0.0]],
        [[0.0, 0.0], [np.inf, 0.0]],
    ],
)
def test_invalid_nodes(nodes):
    with pytest.raises(ValueError):
        solve_truss(nodes, ELEMENTS, E, A, np.zeros(4), [0, 1])


@pytest.mark.parametrize(
    "elements",
    [
        [[0.0, 1.5]],
        [[0, 2]],
        [[0, 0]],
        [[0, 1], [1, 0]],
        [],
    ],
)
def test_invalid_elements(elements):
    with pytest.raises(ValueError):
        solve_truss(NODES, elements, E, A, np.zeros(4), [0, 1])


@pytest.mark.parametrize(
    "loads",
    [
        [0.0, 0.0],
        [0.0, 0.0, 0.0, np.nan],
        [0.0, 0.0, 0.0, np.inf],
    ],
)
def test_invalid_load_vectors(loads):
    with pytest.raises(ValueError):
        solve_truss(NODES, ELEMENTS, E, A, loads, [0, 1])


@pytest.mark.parametrize(
    "fixed_dofs",
    [
        [0.5, 1],
        [-1, 1],
        [0, 4],
        [0, 0],
        [0, np.nan],
    ],
)
def test_invalid_fixed_dofs(fixed_dofs):
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(4),
            fixed_dofs,
        )


def test_unstable_truss_raises_error():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(4),
            [0],
        )


def test_invalid_displacement_vector():
    with pytest.raises(ValueError):
        calculate_truss_results(
            NODES,
            ELEMENTS,
            E,
            A,
            [0.0, 0.0, np.nan, 0.0],
        )


@pytest.mark.parametrize(
    "density",
    [0, -1, np.nan, np.inf],
)
def test_invalid_density(density):
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            NODES,
            ELEMENTS,
            E,
            A,
            density,
        )


def test_truss_element_rejects_nonfinite_coordinates():
    with pytest.raises(ValueError):
        truss_element_stiffness(
            E, A, 0.0, 0.0, np.inf, 0.0
        )


def test_truss_element_rejects_zero_length():
    with pytest.raises(ValueError):
        truss_element_stiffness(
            E, A, 0.0, 0.0, 0.0, 0.0
        )
