
import numpy as np

from mechfem import (
    analyze_truss_modal,
    plot_truss_mode_shape,
)


def main():
    E = 200e9
    A = 0.01
    density = 7850

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

    fixed_dofs = [0, 1, 3]

    results = analyze_truss_modal(
        nodes,
        elements,
        E,
        A,
        density,
        fixed_dofs
    )

    print("Natural Frequencies (Hz):")
    print(results["natural_frequencies"])

    print("\nFirst Mode Shape:")
    first_mode = results["mode_shapes"][:, 0]
    print(first_mode)

    plot_truss_mode_shape(
        nodes,
        elements,
        first_mode,
        scale=100.0,
        show=True
    )


if __name__ == "__main__":
    main()
