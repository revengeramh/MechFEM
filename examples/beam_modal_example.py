
import numpy as np

from mechfem import analyze_beam_modal


def main():
    E = 200e9
    I = 8e-6
    A = 0.01
    density = 7850.0
    length = 2.0
    number_of_elements = 10

    nodes = np.linspace(
        0.0,
        length,
        number_of_elements + 1
    )

    elements = np.array([
        [i, i + 1]
        for i in range(number_of_elements)
    ])

    fixed_dofs = [0, 1]

    results = analyze_beam_modal(
        nodes,
        elements,
        E,
        I,
        A,
        density,
        fixed_dofs,
        num_modes=5,
    )

    print("MechFEM - Beam Modal Analysis")
    print("--------------------------------")

    print("\nNatural Frequencies (Hz):")

    for i, frequency in enumerate(
        results["natural_frequencies"],
        start=1
    ):
        print(f"Mode {i}: {frequency:.6f} Hz")

    print("\nFirst Mode Shape:")
    print(results["mode_shapes"][:, 0])


if __name__ == "__main__":
    main()
