
# MechFEM - Numerical Verification and Validation

## 1. Introduction

This document describes the numerical verification and validation strategy of the MechFEM finite element framework.

The purpose of these tests is to evaluate the mathematical consistency, numerical behavior, and implementation reliability of the developed finite element modules.

The current validation activities include:

- Euler-Bernoulli beam analysis.
- Two-dimensional truss analysis.
- Structural post-processing.
- Mesh convergence analysis.
- Modal analysis.
- Numerical matrix verification.

## 2. Beam Element Validation

The Euler-Bernoulli beam element stiffness matrix is checked for symmetry.

A cantilever beam subjected to a transverse tip load is compared against the analytical Euler-Bernoulli beam solution.

The analytical tip displacement is:

\[
\delta = \frac{FL^3}{3EI}
\]

The analytical tip rotation is:

\[
\theta = \frac{FL^2}{2EI}
\]

The numerical results are compared with these analytical expressions.

## 3. Two-Dimensional Truss Validation

The two-dimensional truss element stiffness matrix is checked for symmetry.

A triangular truss model is used to verify:

- Nodal displacement calculations.
- Support reaction calculations.
- Structural equilibrium.
- Axial force calculations.
- Normal stress calculations.

The element length is calculated as:

\[
L=\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}
\]

The axial deformation of an element is:

\[
\delta =
c(u_2-u_1)+s(v_2-v_1)
\]

The axial strain is:

\[
\varepsilon=\frac{\delta}{L}
\]

The corresponding linear elastic stress is:

\[
\sigma=E\varepsilon
\]

## 4. Analytical Displacement Validation

For a one-dimensional axial bar subjected to an axial force, the analytical displacement is:

\[
\delta=\frac{FL}{AE}
\]

Where:

- \(F\): Applied axial force.
- \(L\): Bar length.
- \(A\): Cross-sectional area.
- \(E\): Young's modulus.

The numerical displacement is compared against the analytical solution.

## 5. Structural Analysis Validation

The structural post-processing module is evaluated using simple axial member configurations.

The validation includes:

- Element length calculation.
- Axial elongation calculation.
- Axial strain calculation.
- Normal stress calculation.
- Axial force calculation.
- Element strain energy calculation.

For a linearly elastic axial member, the strain energy is:

\[
U=\frac{1}{2}N\delta
\]

The numerical results are checked against the corresponding analytical relationships.

## 6. Mesh Convergence Validation

Mesh convergence analysis is evaluated using multiple beam discretizations.

The numerical displacement results are compared against an analytical reference solution.

The absolute error is:

\[
e_{\mathrm{abs}}=
|u_h-u_{\mathrm{ref}}|
\]

The relative error is:

\[
e_{\mathrm{rel}}=
\frac{|u_h-u_{\mathrm{ref}}|}
{|u_{\mathrm{ref}}|}
\]

The validation checks include:

- Correct mesh refinement ordering.
- Finite displacement values.
- Non-negative error values.
- Expected error reduction with mesh refinement.
- Convergence rate calculation.

## 7. Modal Analysis Validation

Modal analysis is verified using a generalized eigenvalue problem:

\[
\mathbf{K}\boldsymbol{\phi}
=
\omega^2\mathbf{M}\boldsymbol{\phi}
\]

The numerical implementation is evaluated through:

- Natural frequency calculation.
- Angular frequency calculation.
- Mode shape reconstruction.
- Boundary condition enforcement.
- Positive eigenvalue verification.
- Mass matrix symmetry and positive definiteness.
- Finite numerical results.

### 7.1 Analytical Single-Element Truss Benchmark

A single horizontal truss element is considered.

The first node is fully constrained, and the transverse degree of freedom at the second node is constrained.

The remaining free degree of freedom represents axial vibration.

For a truss element with consistent mass formulation, the reduced stiffness and mass matrices are:

\[
K_{ff}=\frac{EA}{L}
\]

\[
M_{ff}=\frac{\rho AL}{3}
\]

The angular natural frequency is therefore:

\[
\omega=
\sqrt{\frac{K_{ff}}{M_{ff}}}
\]

Simplifying:

\[
\omega=
\frac{1}{L}\sqrt{\frac{3E}{\rho}}
\]

The analytical natural frequency in Hertz is:

\[
f=
\frac{1}{2\pi L}
\sqrt{\frac{3E}{\rho}}
\]

The numerical natural frequency calculated by MechFEM is compared against this analytical result.

This benchmark verifies the implementation of the consistent mass matrix and the generalized eigenvalue solver for the specified single-element configuration.

## 8. Automated Testing

The project uses pytest for automated unit and numerical verification.

GitHub Actions automatically executes the test suite when changes are pushed to the repository.

The automated tests cover:

- Element stiffness matrix properties.
- Analytical displacement benchmarks.
- Structural stability checks.
- Material and geometry validation.
- Structural post-processing.
- Mesh convergence behavior.
- Modal analysis.
- Visualization functionality.
- Example script execution.

Passing automated tests confirms that the tested cases satisfy their specified assertions. It does not establish accuracy for every possible structural configuration.

## 9. Current Validation Scope

The current validation framework covers selected benchmark cases for:

- Euler-Bernoulli beam elements.
- Two-dimensional truss elements.
- Two-dimensional frame elements.
- Linear static structural analysis.
- Structural post-processing.
- Beam mesh convergence.
- Undamped modal analysis.
- Consistent truss mass matrix assembly.
- Numerical visualization.

Additional benchmark problems and independent comparisons with established FEM software are planned for future development.

## 10. Future Validation Activities

Future validation work may include:

- Three-dimensional structural benchmarks.
- Plane stress and plane strain benchmark problems.
- Nonlinear material behavior.
- Geometric nonlinear analysis.
- Dynamic response verification.
- Additional modal benchmark structures.
- Comparison with established commercial and open-source FEM software.

---

**Project:** MechFEM

**Repository:** https://github.com/revengeramh/MechFEM
