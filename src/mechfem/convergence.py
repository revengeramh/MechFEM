
import numpy as np

from .beam import solve_beam


def beam_mesh_convergence(
    length,
    E,
    I,
    load,
    element_counts,
    fixed_dofs=None,
    target_dof=None,
    analytical_value=None,
):
    if not np.isfinite(length) or length <= 0:
        raise ValueError("length must be positive and finite.")

    if not np.isfinite(E) or E <= 0:
        raise ValueError("E must be positive and finite.")

    if not np.isfinite(I) or I <= 0:
        raise ValueError("I must be positive and finite.")

    if not np.isfinite(load):
        raise ValueError("load must be finite.")

    element_counts = np.asarray(element_counts)

    if element_counts.ndim != 1 or len(element_counts) < 2:
        raise ValueError("At least two mesh sizes are required.")

    if not np.issubdtype(element_counts.dtype, np.integer):
        raise ValueError("Element counts must be integers.")

    if np.any(element_counts < 2):
        raise ValueError("Each mesh must contain at least two elements.")

    if len(np.unique(element_counts)) != len(element_counts):
        raise ValueError("Element counts must be unique.")

    element_counts = np.sort(element_counts)

    if analytical_value is not None and not np.isfinite(analytical_value):
        raise ValueError("analytical_value must be finite.")

    displacements = []
    errors = []
    relative_errors = []
    rates = []

    for n_elements in element_counts:
        nodes = np.linspace(0.0, length, n_elements + 1)

        elements = np.column_stack((
            np.arange(n_elements),
            np.arange(1, n_elements + 1)
        ))

        loads = np.zeros(2 * len(nodes))

        mid_node = int(np.argmin(np.abs(nodes - length / 2)))
        loads[2 * mid_node] = load

        if fixed_dofs is None:
            constraints = [0, 2 * n_elements]
        else:
            constraints = list(fixed_dofs)

        if target_dof is None:
            dof = 2 * mid_node
        else:
            dof = int(target_dof)

        if dof < 0 or dof >= 2 * len(nodes):
            raise ValueError("target_dof is outside the mesh DOF range.")

        u, reactions = solve_beam(
            nodes,
            elements,
            E,
            I,
            loads,
            constraints
        )

        value = float(u[dof])
        displacements.append(value)

        if analytical_value is not None:
            error = abs(value - analytical_value)
            errors.append(error)

            if analytical_value != 0:
                relative_errors.append(error / abs(analytical_value))
            else:
                relative_errors.append(np.nan)

    if analytical_value is not None:
        for i in range(len(errors)):
            if i == 0:
                rates.append(np.nan)
            elif errors[i] > 0 and errors[i - 1] > 0:
                ratio = element_counts[i] / element_counts[i - 1]
                rate = np.log(errors[i - 1] / errors[i]) / np.log(ratio)
                rates.append(float(rate))
            else:
                rates.append(np.nan)

    return {
        "element_counts": element_counts,
        "displacements": np.asarray(displacements),
        "absolute_errors": np.asarray(errors),
        "relative_errors": np.asarray(relative_errors),
        "convergence_rates ": np.asarray(rates),
    }
