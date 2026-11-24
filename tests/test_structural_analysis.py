
import numpy as np
import pytest

from mechfem.structural_analysis import analyze_truss


def test_single_element_tension():
    nodes = np.array([
        [0.0, 0.0],
        [2.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    E = 200e9
    A = 0.01

    displacements = np.array([
        0.0, 0.0,
        0.001, 0.0
    ])

    result = analyze_truss(
        nodes, elements, E, A, displacements
    )

    assert result["lengths"][0] == pytest.approx(2.0)
    assert result["elongations"][0] == pytest.approx(0.001)
    assert result["strains"][0] == pytest.approx(0.0005)
    assert result["stresses"][0] == pytest.approx(100e6)
    assert result["axial_forces"][0] == pytest.approx(1e6)
    assert result["strain_energies"][0] == pytest.approx(500.0)


def test_single_element_compression():
    nodes = np.array([
        [0.0, 0.0],
        [2.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    result = analyze_truss(
        nodes,
        elements,
        200e9,
        0.01,
        np.array([0.0, 0.0, -0.001, 0.0])
    )

    assert result["strains"][0] == pytest.approx(-0.0005)
    assert result["stresses"][0] == pytest.approx(-100e6)
    assert result["axial_forces"][0] == pytest.approx(-1e6)


def test_diagonal_element():
    nodes = np.array([
        [0.0, 0.0],
        [3.0, 4.0]
    ])

    elements = np.array([[0, 1]])

    # 0.005 m displacement along the element axis.
    displacement = 0.005 / 5.0

    displacements = np.array([
        0.0, 0.0,
        3.0 * displacement,
        4.0 * displacement
    ])

    result = analyze_truss(
        nodes, elements, 100e9, 0.002, displacements
    )

    assert result["lengths"][0] == pytest.approx(5.0)
    assert result["elongations"][0] == pytest.approx(0.005)
    assert result["strains"][0] == pytest.approx(0.001)
    assert result["stresses"][0] == pytest.approx(100e6)


def test_element_specific_material_properties():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [3.0, 0.0]
    ])

    elements = np.array([
        [0, 1],
        [1, 2]
    ])

    displacements = np.array([
        0.0, 0.0,
        0.001, 0.0,
        0.003, 0.0
    ])

    result = analyze_truss(
        nodes,
        elements,
        np.array([100e9, 200e9]),
        np.array([0.01, 0.02]),
        displacements
    )

    assert result["stresses"][0] == pytest.approx(100e6)
    assert result["stresses"][1] == pytest.approx(200e6)


@pytest.mark.parametrize("E,A", [
    (0, 0.01),
    (-200e9, 0.01),
    (200e9, 0),
    (200e9, -0.01)
])
def test_invalid_material_properties(E, A):
    nodes = np.array([[0.0, 0.0], [1.0, 0.0]])
    elements = np.array([[0, 1]])
    displacements = np.zeros(4)

    with pytest.raises(ValueError):
        analyze_truss(nodes, elements, E, A, displacements)


def test_invalid_displacement_size():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0]])
    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        analyze_truss(nodes, elements, 200e9, 0.01, np.zeros(3))


def test_invalid_connectivity():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0]])
    elements = np.array([[0, 2]])

    with pytest.raises(ValueError):
        analyze_truss(nodes, elements, 200e9, 0.01, np.zeros(4))


def test_zero_length_element():
    nodes = np.array([
        [0.0, 0.0],
        [0.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    with pytest.raises(ValueError):
        analyze_truss(nodes, elements, 200e9, 0.01, np.zeros(4))


def test_rigid_body_translation():
    nodes = np.array([
        [0.0, 0.0],
        [2.0, 0.0]
    ])

    elements = np.array([[0, 1]])

    displacements = np.array([
        0.01, 0.02,
        0.01, 0.02
    ])

    result = analyze_truss(
        nodes, elements, 200e9, 0.01, displacements
    )

    assert result["strains"][0] == pytest.approx(0.0)
    assert result["stresses"][0] == pytest.approx(0.0)
    assert result["axial_forces"][0] == pytest.approx(0.0)
