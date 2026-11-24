
import numpy as np

from mechfem import solve_truss, analyze_truss


def main():
    # Geometry: 2D triangular truss
    nodes = np.array([
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 1.5]
    ])

    elements = np.array([
        [0, 1],
        [0, 2],
        [1, 2]
    ])

    # Material and cross-section
    E = 200e9       # Young's modulus (Pa)
    A = 0.01        # Area (m^2)

    # External loads: [Fx0, Fy0, Fx1, Fy1, ...]
    loads = np.zeros(2 * len(nodes))
    loads[2 * 2 + 1] = -10000.0

    # Fixed degrees of freedom
    fixed_dofs = [
        0, 1,       # Node 0
        3           # Vertical DOF of node 1
    ]

    # Solve structural displacements
    displacements = solve_truss(
        nodes,
        elements,
        E,
        A,
        loads,
        fixed_dofs
    )

    # Calculate element-level structural results
    results = analyze_truss(
        nodes,
        elements,
        E,
        A,
        displacements
    )

    print("\nMECHFEM - STRUCTURAL ANALYSIS")
    print("=" * 55)

    print("\nNODAL DISPLACEMENTS")
    for i in range(len(nodes)):
        ux = displacements[2*i]
        uy = displacements[2*i + 1]

        print(
            f"Node {i}: "
            f"Ux = {ux:.6e} m, "
            f"Uy = {uy:.6e} m"
        )

    print("\nELEMENT RESULTS")
    print("-" * 55)

    for i in range(len(elements)):
        print(f"\nElement {i + 1}: {elements[i]}")
        print(f"Length: {results['lengths'][i]:.4f} m")
        print(
            f"Axial deformation: "
            f"{results['elongations'][i]:.6e} m"
        )
        print(
            f"Axial strain: "
            f"{results['strains'][i]:.6e}"
        )
        print(
            f"Axial stress: "
            f"{results['stresses'][i] / 1e6:.4f} MPa"
        )
        print(
            f"Axial force: "
            f"{results['axial_forces'][i] / 1e3:.4f} kN"
        )
        print(
            f"Strain energy: "
            f"{results['strain_energies'][i]:.6f} J"
        )

    total_energy = np.sum(results["strain_energies"])

    print("\n" + "=" * 55)
    print(f"Total strain energy: {total_energy:.6f} J")
    print("Analysis completed successfully.")


if __name__ == "__main__":
    main()
