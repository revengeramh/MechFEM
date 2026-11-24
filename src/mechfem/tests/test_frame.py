
import numpy as np
import pytest

from mechfem.frame import (
    frame_element_stiffness,
    solve_frame
)


def test_frame_element_stiffness_symmetry():
    k = frame_element_stiffness(
        E=200e9,
        A=0.01,
        I=8e-6,
        x1=0,
        y1=0,
        x2=4,
        y2=0
    )

    assert k.shape == (6, 6)
    assert np.allclose(k, k.T)


def test_frame_element_stiffness_positive_semidefinite():
    k = frame_element_stiffness(
        200e9, 0.01, 8e-6,
        0, 0, 4, 0
    )

    eigenvalues = np.linalg.eigvalsh(k)

    assert np.all(eigenvalues >= -1e-5)


def test_frame_element_invalid_length():
    with pytest.raises(ValueError):
        frame_element_stiffness(
            200e9, 0.01, 8e-6,
            0, 0, 0, 0
        )


def test_cantilever_frame_axial_displacement():
    E = 200e9
    A = 0.01
    I = 8e-6
    L = 4.0
    P = 10000.0

    nodes = [
        [0, 0],
        [L, 0]
    ]

    elements = [[0, 1]]

    loads = {3: P}
    fixed_dofs = [0, 1, 2]

    u, reactions = solve_frame(
        nodes, elements, E, A, I,
        loads, fixed_dofs
    )

    expected = P * L / (E * A)

    assert u[3] == pytest.approx(expected, rel=1e-8)
    assert reactions[0] == pytest.approx(-P)


def test_cantilever_frame_transverse_displacement():
    E = 200e9
    A = 0.01
    I = 8e-6
    L = 4.0
    P = 1000.0

    nodes = [
        [0, 0],
        [L, 0]
    ]

    elements = [[0, 1]]

    loads = {4: P}
    fixed_dofs = [0, 1, 2]

    u, reactions = solve_frame(
        nodes, elements, E, A, I,
        loads, fixed_dofs
    )

    expected = P * L**3 / (3 * E * I)

    assert u[4] == pytest.approx(expected, rel=1e-8)
    assert reactions[1] == pytest.approx(-P)


def test_frame_invalid_material_properties():
    with pytest.raises(ValueError):
        frame_element_stiffness(
            0, 0.01, 8e-6,
            0, 0, 4, 0
        )


def test_frame_invalid_connectivity():
    with pytest.raises(ValueError):
        solve_frame(
            nodes=[[0, 0], [4, 0]],
            elements=[[0, 2]],
            E=200e9,
            A=0.01,
            I=8e-6,
            loads={},
            fixed_dofs=[0, 1, 2]
        )
