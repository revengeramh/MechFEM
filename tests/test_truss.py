
import numpy as np
import pytest

from mechfem.truss import truss_element_stiffness


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
