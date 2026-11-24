
import numpy as np
import pytest

from mechfem.truss import (
    truss_element_stiffness,
    solve_truss,
    calculate_truss_results
) 


def test_truss_stiffness_symmetry():
    k = truss_element_stiffness(
        E=200e9,
        A=0.01,
        x1=0,
        y1=0,
        x2=3,
        y2=4
    )

    assert k.shape == (4, 4)
    assert np.allclose(k, k.T)


def test_horizontal_truss():
    k = truss_element_stiffness(
        E=200e9,
        A=0.01,
        x1=0,
        y1=0,
        x2=2,
        y2=0
    )

    assert np.isclose(k[0, 0], 1e9)
    assert np.isclose(k[2, 2], 1e9)
    assert np.isclose(k[0, 2], -1e9)


def test_invalid_material_properties():
    with pytest.raises(ValueError):
        truss_element_stiffness(
            E=0,
            A=0.01,
            x1=0,
            y1=0,
            x2=1,
            y2=0
        )


def test_zero_length_element():
    with pytest.raises(ValueError):
        truss_element_stiffness(
            E=200e9,
            A=0.01,
            x1=0,
            y1=0,
            x2=0,
            y2=0
        )

def test_triangular_truss_solution():
    nodes = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]
    ]

    elements = [
        [0, 1],
        [0, 2],
        [1, 2]
    ]

    E = 200e9
    A = 0.01

    loads = np.zeros(6)
    loads[5] = -1000.0

    fixed_dofs = [0, 1, 3]

    displacements, reactions = solve_truss(
        nodes,
        elements,
        E,
        A,
        loads,
        fixed_dofs
    )

    assert displacements.shape == (6,)
    assert reactions.shape == (6,)

    assert np.all(np.isfinite(displacements))
    assert np.all(np.isfinite(reactions))

    assert np.allclose(
        reactions[fixed_dofs],
        [0.0, 1000.0, 0.0],
        atol=1e-6
    )

    assert np.isclose(
        np.sum(reactions[1::2] + loads[1::2]),
        0.0,
        atol=1e-6
    )

def test_truss_axial_force_and_stress():
    E = 200e9
    A = 0.01
    L = 2.0

    nodes = [
        [0.0, 0.0],
        [L, 0.0]
    ]

    elements = [[0, 1]]

    # 1 mm extension at the second node
    displacements = np.array([
        0.0,
        0.0,
        0.001,
        0.0
    ])

    forces, stresses = calculate_truss_results(
        nodes,
        elements,
        E,
        A,
        displacements
    )

    expected_stress = 100e6
    expected_force = 1e6

    assert np.isclose(stresses[0], expected_stress)
    assert np.isclose(forces[0], expected_force)
    
def test_truss_example():
    import runpy
    from pathlib import Path

    example_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "truss_example.py"
    )

    example = runpy.run_path(str(example_path))
    example["main"]( )
