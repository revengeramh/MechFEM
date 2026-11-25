import numpy as np

from mechfem.frame import solve_frame


def test_frame_cantilever_tip_load_matches_analytical_solution():
    E = 200e9
    A = 0.01
    I = 8e-6
    L = 2.0
    P = 1000.0

    nodes = [
        [0.0, 0.0],
        [L, 0.0],
    ]

    elements = [[0, 1]]

    loads = {
        4: -P
    }

    fixed_dofs = [0, 1, 2]

    displacements, reactions = solve_frame(
        nodes,
        elements,
        E,
        A,
        I,
        loads,
        fixed_dofs,
    )

    expected_deflection = -P * L**3 / (3 * E * I)
    expected_rotation = -P * L**2 / (2 * E * I)

    assert np.isclose(
        displacements[4],
        expected_deflection,
        rtol=1e-9,
    )

    assert np.isclose(
        displacements[5],
        expected_rotation,
        rtol=1e-9,
    )

    assert np.isclose(reactions[1], P, rtol=1e-9)
    assert np.isclose(reactions[2], P * L, rtol=1e-9)


def test_frame_axial_displacement_matches_analytical_solution():
    E = 200e9
    A = 0.01
    I = 8e-6
    L = 2.0
    P = 10000.0

    nodes = [
        [0.0, 0.0],
        [L, 0.0],
    ]

    elements = [[0, 1]]

    loads = {
        3: P
    }

    fixed_dofs = [0, 1, 2]

    displacements, reactions = solve_frame(
        nodes,
        elements,
        E,
        A,
        I,
        loads,
        fixed_dofs,
    )

    expected_displacement = P * L / (E * A)

    assert np.isclose(
        displacements[3],
        expected_displacement,
        rtol=1e-9,
    )

    assert np.isclose(reactions[0], -P, rtol=1e-9)


def test_frame_global_equilibrium():
    E = 200e9
    A = 0.01
    I = 8e-6

    nodes = [
        [0.0, 0.0],
        [2.0, 0.0],
    ]

    elements = [[0, 1]]

    loads = {
        3: 5000.0,
        4: -1000.0,
        5: 200.0,
    }

    fixed_dofs = [0, 1, 2]

    displacements, reactions = solve_frame(
        nodes,
        elements,
        E,
        A,
        I,
        loads,
        fixed_dofs,
    )

    total_forces = reactions + np.array([
        0.0, 0.0, 0.0,
        5000.0, -1000.0, 200.0,
    ])

    assert np.allclose(
        total_forces[[0, 1]],
        [0.0, 0.0],
        atol=1e-8,
    )

    total_moment = (
        reactions[2]
        + 200.0
        + 2.0 * (-1000.0)
    )

    assert np.isclose(total_moment, 0.0, atol=1e-8)
