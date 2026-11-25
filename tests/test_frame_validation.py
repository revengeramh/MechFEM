import numpy as np
import pytest

from mechfem.frame import (
    frame_element_stiffness,
    solve_frame,
)


NODES = [
    [0.0, 0.0],
    [2.0, 0.0],
]

ELEMENTS = [[0, 1]]

E = 200e9
A = 0.01
I = 8e-6


@pytest.mark.parametrize(
    "elements",
    [
        [[0, 1.5]],
        [[0, 2]],
        [[0, 0]],
        [[0, 1], [1, 0]],
        [],
    ],
)
def test_invalid_elements(elements):
    with pytest.raises(ValueError):
        solve_frame(
            NODES, elements, E, A, I, {}, [0, 1, 2]
        )


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
        solve_frame(
            nodes, ELEMENTS, E, A, I, {}, [0, 1, 2]
        )


@pytest.mark.parametrize(
    "property_value",
    [0, -1, np.nan, np.inf],
)
def test_invalid_material_properties(property_value):
    with pytest.raises(ValueError):
        solve_frame(
            NODES,
            ELEMENTS,
            property_value,
            A,
            I,
            {},
            [0, 1, 2],
        )


@pytest.mark.parametrize(
    "loads",
    [
        [],
        [1000.0],
        {6: np.nan},
        {6: np.inf},
        {12: 1000.0},
        {1.5: 1000.0},
        {"1": 1000.0},
    ],
)
def test_invalid_loads(loads):
    with pytest.raises(ValueError):
        solve_frame(
            NODES,
            ELEMENTS,
            E,
            A,
            I,
            loads,
            [0, 1, 2],
        )


@pytest.mark.parametrize(
    "fixed_dofs",
    [
        [0.5, 1, 2],
        [-1, 1, 2],
        [0, 1, 6],
        [0, 1, 1],
        [0, np.nan, 2],
    ],
)
def test_invalid_fixed_dofs(fixed_dofs):
    with pytest.raises(ValueError):
        solve_frame(
            NODES,
            ELEMENTS,
            E,
            A,
            I,
            {},
            fixed_dofs,
        )


def test_unstable_frame_raises_error():
    with pytest.raises(ValueError):
        solve_frame(
            NODES,
            ELEMENTS,
            E,
            A,
            I,
            {},
            [0, 1],
        )


def test_frame_element_rejects_zero_length():
    with pytest.raises(ValueError):
        frame_element_stiffness(
            E, A, I,
            0.0, 0.0,
            0.0, 0.0,
        )


def test_frame_element_rejects_nonfinite_coordinates():
    with pytest.raises(ValueError):
        frame_element_stiffness(
            E, A, I,
            0.0, 0.0,
            np.inf, 0.0,
        )


def test_valid_frame_solution():
    loads = {3: 1000.0}

    displacements, reactions = solve_frame(
        NODES,
        ELEMENTS,
        E,
        A,
        I,
        loads,
        [0, 1, 2],
    )

    assert displacements.shape == (6,)
    assert reactions.shape == (6,)
    assert np.all(np.isfinite(displacements))
    assert np.all(np.isfinite(reactions))
