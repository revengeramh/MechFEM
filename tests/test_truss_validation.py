import numpy as np
import pytest

from mechfem.truss import (
    solve_truss,
    calculate_truss_results,
    assemble_truss_matrices,
)


NODES = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
]

ELEMENTS = [
    [0, 1],
    [0, 2],
    [1, 2],
]

E = 200e9
A = 0.01


def test_reject_fractional_element_indices():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            [[0, 1.5]],
            E,
            A,
            np.zeros(6),
            [0, 1, 3],
        )


def test_reject_out_of_range_element_indices():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            NODES,
            [[0, 5]],
            E,
            A,
            7850,
        )


def test_reject_self_connected_element():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            NODES,
            [[1, 1]],
            E,
            A,
            7850,
        )


def test_reject_duplicate_elements():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            NODES,
            [[0, 1], [1, 0]],
            E,
            A,
            7850,
        )


def test_reject_non_finite_node_coordinates():
    invalid_nodes = [
        [0.0, 0.0],
        [np.nan, 1.0],
    ]

    with pytest.raises(ValueError):
        assemble_truss_matrices(
            invalid_nodes,
            [[0, 1]],
            E,
            A,
            7850,
        )


def test_reject_invalid_node_shape():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            [[0.0, 0.0, 1.0], [1.0, 0.0, 2.0]],
            [[0, 1]],
            E,
            A,
            7850,
        )


def test_reject_non_finite_loads():
    loads = np.zeros(6)
    loads[2] = np.nan

    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            loads,
            [0, 1, 3],
        )


def test_reject_wrong_load_vector_size():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(5),
            [0, 1, 3],
        )


def test_reject_fractional_fixed_dofs():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(6),
            [0, 1, 3.5],
        )


def test_reject_out_of_range_fixed_dofs():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(6),
            [0, 1, 6],
        )


def test_reject_duplicate_fixed_dofs():
    with pytest.raises(ValueError):
        solve_truss(
            NODES,
            ELEMENTS,
            E,
            A,
            np.zeros(6),
            [0, 1, 1, 3],
        )


def test_reject_non_finite_displacements():
    displacements = np.zeros(6)
    displacements[2] = np.inf

    with pytest.raises(ValueError):
        calculate_truss_results(
            NODES,
            ELEMENTS,
            E,
            A,
            displacements,
        )


def test_reject_invalid_density():
    with pytest.raises(ValueError):
        assemble_truss_matrices(
            NODES,
            ELEMENTS,
            E,
            A,
            0,
        )


def test_reject_unstable_structure():
    nodes = [
        [0.0, 0.0],
        [1.0, 0.0],
    ]

    elements = [[0, 1]]
    loads = np.zeros(4)

    with pytest.raises(ValueError):
        solve_truss(
            nodes,
            elements,
            E,
            A,
            loads,
            [0, 1],
        )
