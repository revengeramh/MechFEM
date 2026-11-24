
import numpy as np

from mechfem.frame import solve_frame


def main():
    # Material and section properties
    E = 200e9
    A = 0.01
    I = 8e-6

    # Two-dimensional portal frame
    # Units: meters
    nodes = np.array([
        [0.0, 0.0],
        [0.0, 3.0],
        [4.0, 3.0],
        [4.0, 0.0]
    ])

    elements = np.array([
        [0, 1],
        [1, 2],
        [2, 3]
    ])

    # Global degrees of freedom:
    # Node 0: 0, 1, 2
    # Node 1: 3, 4, 5
    # Node 2: 6, 7, 8
    # Node 3: 9, 10, 11

    # Horizontal load at node 2
    loads = {
        6: 10000.0
    }

    # Fixed supports at nodes 0 and 3
    fixed_dofs = [
        0, 1, 2,
        9, 10, 11
    ]

    displacements, reactions = solve_frame(
        nodes,
        elements,
        E,
        A,
        I,
        loads,
        fixed_dofs
    )

    print("2D Portal Frame Analysis")
    print("=" * 40)

    for i in range(len(nodes)):
        ux = displacements[3*i]
        uy = displacements[3*i + 1]
        theta = displacements[3*i + 2]

        print(f"\nNode {i}:")
        print(f"Horizontal displacement: {ux:.8e} m")
        print(f"Vertical displacement:   {uy:.8e} m")
        print(f"Rotation:                 {theta:.8e} rad")

    print("\nSupport Reactions")
    print("=" * 40)

    for dof, value in enumerate(reactions):
        if dof in fixed_dofs:
            print(f"DOF {dof}: {value:.4f} N or N.m")


if __name__ == "__main__":
    main()
