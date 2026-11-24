
import numpy as np

from mechfem.truss import (
    solve_truss,
    calculate_truss_results
)


def main():
    E = 200e9
    A = 0.01

    nodes = np.array([
        [0.0, 0.0],
        [4.0, 0.0],
        [2.0, 3.0]
    ])

    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2]
    ])

    loads = np.zeros(6)
    loads[5] = -10000.0

    fixed_dofs = [0, 1, 3]

    displacements, reactions = solve_truss(
        nodes,
        elements,
        E,
        A,
        loads,
        fixed_dofs
    )

    forces, stresses = calculate_truss_results(
        nodes,
        elements,
        E,
        A,
        displacements
    )

    print("Nodal Displacements:")
    print(displacements.reshape(-1, 2))

    print("\nSupport Reactions:")
    print(reactions.reshape(-1, 2))

    print("\nAxial Forces:")
    print(forces)

    print("\nAxial Stresses:")
    print(stresses)


if __name__ == "__main__":
    main()
