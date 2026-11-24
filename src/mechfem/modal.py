
import numpy as np


def modal_analysis(
    stiffness_matrix,
    mass_matrix,
    fixed_dofs,
    num_modes=None,
):
    """
    Perform modal analysis using the generalized eigenvalue problem.

    K * phi = omega^2 * M * phi

    Returns
    -------
    dict
        natural_frequencies : Frequencies in Hz.
        angular_frequencies : Frequencies in rad/s.
        mode_shapes : Full displacement mode shapes.
        eigenvalues : Squared angular frequencies.
    """

    K = np.asarray(stiffness_matrix, dtype=float)
    M = np.asarray(mass_matrix, dtype=float)

    if K.ndim != 2 or K.shape[0] != K.shape[1]:
        raise ValueError("Stiffness matrix must be square.")

    if M.ndim != 2 or M.shape != K.shape:
        raise ValueError(
            "Mass matrix must have the same square shape as K."
        )

    if not np.all(np.isfinite(K)) or not np.all(np.isfinite(M)):
        raise ValueError("Matrices must contain only finite values.")

    if not np.allclose(K, K.T):
        raise ValueError("Stiffness matrix must be symmetric.")

    if not np.allclose(M, M.T):
        raise ValueError("Mass matrix must be symmetric.")

    n_dofs = K.shape[0]

    fixed_dofs = np.asarray(fixed_dofs)

    if fixed_dofs.ndim != 1:
        raise ValueError("fixed_dofs must be one-dimensional.")

    if fixed_dofs.size == 0:
        fixed_dofs = np.array([], dtype=int)

    elif not np.issubdtype(fixed_dofs.dtype, np.integer):
        raise ValueError("fixed_dofs must contain integers.")

    if np.any(fixed_dofs < 0) or np.any(fixed_dofs >= n_dofs):
        raise ValueError("Invalid fixed degree of freedom.")

    if len(np.unique(fixed_dofs)) != len(fixed_dofs):
        raise ValueError("Duplicate constrained degrees of freedom.")

    free_dofs = np.setdiff1d(
        np.arange(n_dofs),
        fixed_dofs
    )

    if len(free_dofs) == 0:
        raise ValueError("No free degrees of freedom.")

    if num_modes is not None:
        if (
            not isinstance(num_modes, (int, np.integer))
            or isinstance(num_modes, (bool, np.bool_))
            or num_modes < 1
            or num_modes > len(free_dofs)
        ):
            raise ValueError("Invalid number of requested modes.")

    Kff = K[np.ix_(free_dofs, free_dofs)]
    Mff = M[np.ix_(free_dofs, free_dofs)]

    try:
        L = np.linalg.cholesky(Mff)

        left = np.linalg.solve(L, Kff)

        transformed_K = np.linalg.solve(
            L,
            left.T
        ).T

        transformed_K = (
            transformed_K + transformed_K.T
        ) / 2.0

        eigenvalues, eigenvectors = np.linalg.eigh(
            transformed_K
        )

    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Mass matrix must be positive definite and "
            "the eigenvalue problem must be solvable."
        ) from exc

    if np.any(eigenvalues <= 0):
        raise ValueError(
            "The constrained structure has non-positive "
            "eigenvalues. Check stability and boundary conditions."
        )

    if num_modes is not None:
        eigenvalues = eigenvalues[:num_modes]
        eigenvectors = eigenvectors[:, :num_modes]

    mode_shapes_free = np.linalg.solve(
        L.T,
        eigenvectors
    )

    mode_shapes = np.zeros(
        (n_dofs, len(eigenvalues))
    )

    mode_shapes[free_dofs, :] = mode_shapes_free

    angular_frequencies = np.sqrt(eigenvalues)

    natural_frequencies = (
        angular_frequencies / (2.0 * np.pi)
    )

    return {
        "natural_frequencies": natural_frequencies,
        "angular_frequencies": angular_frequencies,
        "eigenvalues": eigenvalues,
        "mode_shapes": mode_shapes,
    }
    
def analyze_truss_modal(
    nodes,
    elements,
    E,
    A,
    density,
    fixed_dofs,
    num_modes=None,
):
    """
    Perform modal analysis of a 2D truss structure.

    Returns natural frequencies and mode shapes.
    """

    from .truss import assemble_truss_matrices

    K, M = assemble_truss_matrices(
        nodes,
        elements,
        E,
        A,
        density
    )

    return modal_analysis(
        K,
        M,
        fixed_dofs,
        num_modes=num_modes
    )

