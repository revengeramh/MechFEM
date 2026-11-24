
import numpy as np


def truss_element_stiffness(E, A, x1, y1, x2, y2):
    """
    Calculate the global stiffness matrix of a 2D truss element.
    DOF order: [u1, v1, u2, v2].
    """

    if E <= 0 or A <= 0:
        raise ValueError("E and A must be positive.")

    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx**2 + dy**2)

    if L <= 0:
        raise ValueError("Element length must be positive.")

    c = dx / L
    s = dy / L

    k = (E * A / L) * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
    ])

    return k
