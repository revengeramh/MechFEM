
import numpy as np

from mechfem.beam import solve_beam


def test_simply_supported_beam_mesh_convergence():
    E = 200e9
    I = 8e-6
    L = 2.0
    q = -1000.0

    analytical_deflection = (
        5 * q * L**4 / (384 * E * I)
    )

    errors = []

    for n_elements in [2, 4, 8]:
        nodes = np.linspace(0.0, L, n_elements + 1)

        elements = [
            (i, i + 1)
            for i in range(n_elements)
        ]

        loads = np.zeros(2 * len(nodes))

        h = L / n_elements

        # Approximate uniform load using nodal forces
        loads[0] = q * h / 2
        loads[-2] = q * h / 2

        for i in range(1, n_elements):
            loads[2 * i] = q * h

        # Simply supported beam
        fixed_dofs = [0, 2 * n_elements]

        displacements, reactions = solve_beam(
            nodes,
            elements,
            E,
            I,
            loads,
            fixed_dofs
        )

        midpoint_index = n_elements // 2
        fem_deflection = displacements[2 * midpoint_index]

        relative_error = abs(
            (fem_deflection - analytical_deflection)
            / analytical_deflection
        )

        errors.append(relative_error)

    # Error must decrease as the mesh is refined
    assert errors[1] < errors[0]
    assert errors[2] < errors[1]

    # Final relative error must be below 2 percent
    assert errors[2] < 0.02

