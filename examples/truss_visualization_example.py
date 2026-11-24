
import numpy as np
import matplotlib.pyplot as plt

from mechfem.truss import solve_truss, calculate_truss_results
from mechfem.visualization import (
    plot_truss_axial_forces,
    plot_truss_stresses,
)


def main():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 0.8],
    ])

    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2],
    ])

    E = 200e9
    A = 0.001

    loads = np.zeros(6)
    loads[5] = -10000.0

    fixed_dofs = [0, 1, 3]

    displacements, reactions = solve_truss(
        nodes, elements, E, A, loads, fixed_dofs
    )

    forces, stresses = calculate_truss_results(
        nodes, elements, E, A, displacements
    )

    print("Displacements:")
    print(displacements)

    print("Axial forces (N):")
    print(forces)

    print("Axial stresses (Pa):")
    print(stresses)

    plot_truss_axial_forces(
        nodes, elements, forces,
        save_path="truss_forces.png",
        show=False
    )

    plot_truss_stresses(
        nodes, elements, stresses,
        save_path="truss_stresses.png",
        show=False
    )

    plt.show()


if __name__ == "__main__":
    main()
