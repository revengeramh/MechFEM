
# Theoretical Foundation of MechFEM

## 1. Introduction

MechFEM is based on the Finite Element Method (FEM), a numerical technique used to approximate solutions of engineering problems governed by differential equations.

The main objective of FEM is to transform a continuous engineering system into a finite number of discrete elements connected through nodes.

This document describes the theoretical foundations of the current beam, two-dimensional truss, frame, structural analysis, mesh convergence, and modal analysis implementations.

---

## 2. Fundamental Concepts of FEM

The finite element procedure generally consists of:

1. Discretization of the physical domain.
2. Selection of element interpolation functions.
3. Derivation of element stiffness matrices.
4. Assembly of the global stiffness matrix.
5. Application of boundary conditions.
6. Solution of the resulting algebraic equations.
7. Calculation of engineering quantities such as displacements, reactions, forces, and stresses.

The general static equilibrium equation is:

\[
\mathbf{K}\mathbf{u}=\mathbf{F}
\]

Where:

- \(\mathbf{K}\): Global stiffness matrix.
- \(\mathbf{u}\): Global displacement vector.
- \(\mathbf{F}\): Global external force vector.

---

## 3. Euler-Bernoulli Beam Element

### 3.1 Governing Assumptions

The Euler-Bernoulli beam theory is based on the following assumptions:

- The beam material is linearly elastic.
- Deflections and rotations are small.
- Plane cross-sections remain plane and perpendicular to the neutral axis.
- Shear deformation is neglected.
- The beam has constant material and geometric properties within each element.

### 3.2 Degrees of Freedom

Each two-node beam element has four degrees of freedom:

\[
\mathbf{u}_e =
\begin{bmatrix}
v_1 & \theta_1 & v_2 & \theta_2
\end{bmatrix}^{T}
\]

Where:

- \(v_1,v_2\): Transverse nodal displacements.
- \(\theta_1,\theta_2\): Nodal rotations.

### 3.3 Element Stiffness Matrix

For a beam element with length \(L\), Young's modulus \(E\), and second moment of area \(I\), the local stiffness matrix is:

\[
\mathbf{k}_e =
\frac{EI}{L^3}
\begin{bmatrix}
12 & 6L & -12 & 6L\\
6L & 4L^2 & -6L & 2L^2\\
-12 & -6L & 12 & -6L\\
6L & 2L^2 & -6L & 4L^2
\end{bmatrix}
\]

### 3.4 Global Assembly

The global stiffness matrix is constructed by assembling the element stiffness matrices according to the connectivity of the mesh.

The resulting system is:

\[
\mathbf{K}\mathbf{u}=\mathbf{F}
\]

After applying prescribed displacement boundary conditions, the unknown nodal displacements can be calculated.

### 3.5 Bending Moment and Curvature

The Euler-Bernoulli bending relation is:

\[
M(x)=EI\frac{d^2v}{dx^2}
\]

The exact sign depends on the adopted sign convention.

---

## 4. Two-Dimensional Truss Element

### 4.1 Basic Assumptions

A two-dimensional truss element is characterized by:

- Two nodes.
- Axial deformation only.
- Linear elastic material behavior.
- Small displacements and rotations.
- Loads applied at the nodes.
- Pin-connected members.

Each node has two translational degrees of freedom.

### 4.2 Element Geometry

For two nodes with coordinates:

\[
(x_1,y_1),\qquad (x_2,y_2)
\]

The element length is:

\[
L=\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}
\]

The direction cosines are:

\[
c=\frac{x_2-x_1}{L}
\]

\[
s=\frac{y_2-y_1}{L}
\]

### 4.3 Element Stiffness Matrix

The global-coordinate stiffness matrix of a two-dimensional truss element is:

\[
\mathbf{k}_e=\frac{EA}{L}
\begin{bmatrix}
c^2 & cs & -c^2 & -cs\\
cs & s^2 & -cs & -s^2\\
-c^2 & -cs & c^2 & cs\\
-cs & -s^2 & cs & s^2
\end{bmatrix}
\]

Where:

- \(E\): Young's modulus.
- \(A\): Cross-sectional area.
- \(L\): Element length.
- \(c,s\): Direction cosines.

### 4.4 Axial Deformation

The axial deformation of the element is:

\[
\delta =
c(u_2-u_1)+s(v_2-v_1)
\]

The axial strain is:

\[
\varepsilon=\frac{\delta}{L}
\]

Using Hooke's law:

\[
\sigma=E\varepsilon
\]

The axial force is:

\[
N=A\sigma
\]

Positive axial force represents tension, while negative axial force represents compression under the adopted convention.

---

## 5. Two-Dimensional Frame Element

### 5.1 Basic Assumptions

The two-dimensional frame element combines axial and bending behavior.

The current linear frame formulation assumes:

- Linear elastic material behavior.
- Small displacements and rotations.
- Constant material and geometric properties within each element.
- Two-dimensional structural behavior.
- Axial and bending deformations.

### 5.2 Degrees of Freedom

Each frame node has three degrees of freedom:

\[
\mathbf{u}_i =
\begin{bmatrix}
u_i & v_i & \theta_i
\end{bmatrix}^{T}
\]

Where:

- \(u_i\): Horizontal displacement.
- \(v_i\): Vertical displacement.
- \(\theta_i\): Rotation.

A two-node frame element therefore has six degrees of freedom.

### 5.3 Local Element Stiffness Matrix

For a frame element of length \(L\), Young's modulus \(E\), cross-sectional area \(A\), and second moment of area \(I\), the local stiffness matrix is:

\[
\mathbf{k}_e =
\begin{bmatrix}
\frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0\\
0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2}\\
0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L}\\
-\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0\\
0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2}\\
0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L}
\end{bmatrix}
\]

### 5.4 Global Coordinate Transformation

For an element inclined at an angle \(\theta\), the direction cosines are:

\[
c=\cos\theta,\qquad s=\sin\theta
\]

The global stiffness matrix is obtained using the transformation matrix:

\[
\mathbf{k}_g=\mathbf{T}^{T}\mathbf{k}_e\mathbf{T}
\]

The global element matrices are assembled according to the structural connectivity.

---

## 6. Structural Analysis

### 6.1 Element Deformation

For a two-dimensional truss element, the axial deformation is calculated from the nodal displacement vector projected onto the element axis.

\[
\delta =
c(u_2-u_1)+s(v_2-v_1)
\]

### 6.2 Strain and Stress

The axial strain is:

\[
\varepsilon=\frac{\delta}{L}
\]

For a linear elastic material:

\[
\sigma=E\varepsilon
\]

### 6.3 Axial Force

The axial force is:

\[
N=A\sigma
\]

### 6.4 Element Strain Energy

For a linearly elastic truss element:

\[
U_e=\frac{1}{2}N\delta
\]

Equivalently:

\[
U_e=\frac{1}{2}EA L\varepsilon^2
\]

These expressions assume linear elastic axial behavior.

---

## 7. Mesh Convergence Analysis

### 7.1 Introduction

Mesh convergence analysis evaluates how the numerical solution changes as the finite element mesh is refined.

A numerical solution should approach a stable value as the element size decreases, provided that the formulation and model are appropriate.

### 7.2 Displacement Comparison

Let \(u_h\) represent a numerical displacement obtained using a mesh with characteristic element size \(h\), and let \(u_{\mathrm{ref}}\) represent an analytical or reference solution.

The absolute error is:

\[
e_{\mathrm{abs}}=|u_h-u_{\mathrm{ref}}|
\]

The relative error is:

\[
e_{\mathrm{rel}}=
\frac{|u_h-u_{\mathrm{ref}}|}
{|u_{\mathrm{ref}}|}
\]

This relative error is defined when the reference displacement is nonzero.

### 7.3 Convergence Rate

For two successive mesh refinements, an approximate convergence rate can be estimated from:

\[
p=
\frac{\ln(e_1/e_2)}
{\ln(h_1/h_2)}
\]

Where:

- \(e_1,e_2\): Errors corresponding to two mesh sizes.
- \(h_1,h_2\): Characteristic element sizes.

The convergence rate depends on the element formulation, solution smoothness, and refinement strategy.

### 7.4 MechFEM Implementation

The current beam mesh convergence functionality provides:

- Multiple mesh sizes.
- Displacement calculations.
- Absolute error evaluation.
- Relative error evaluation.
- Approximate convergence rate estimation.
- Comparison with an analytical reference solution.

---

## 8. Modal Analysis

### 8.1 Introduction

Modal analysis is used to determine the natural frequencies and corresponding mode shapes of a structural system.

For an undamped linear elastic structure undergoing free vibration, the governing equation is:

\[
\mathbf{M}\ddot{\mathbf{u}}+\mathbf{K}\mathbf{u}=0
\]

Where:

- \(\mathbf{M}\): Global mass matrix.
- \(\mathbf{K}\): Global stiffness matrix.
- \(\mathbf{u}\): Displacement vector.
- \(\ddot{\mathbf{u}}\): Acceleration vector.

### 8.2 Harmonic Motion Assumption

Assuming harmonic motion:

\[
\mathbf{u}(t)=\boldsymbol{\phi}\sin(\omega t)
\]

Substitution into the governing equation gives:

\[
(\mathbf{K}-\omega^2\mathbf{M})\boldsymbol{\phi}=0
\]

This leads to the generalized eigenvalue problem:

\[
\mathbf{K}\boldsymbol{\phi}
=
\lambda\mathbf{M}\boldsymbol{\phi}
\]

Where:

\[
\lambda=\omega^2
\]

The eigenvalues determine the natural frequencies, while the eigenvectors represent the corresponding mode shapes.

### 8.3 Natural Frequencies

The angular natural frequency is:

\[
\omega_i=\sqrt{\lambda_i}
\]

The natural frequency in Hertz is:

\[
f_i=\frac{\omega_i}{2\pi}
\]

### 8.4 Consistent Mass Matrix for 2D Truss Elements

For a two-node truss element with uniform density, cross-sectional area, and length, the consistent mass matrix in global coordinates is:

\[
\mathbf{M}_e=
\frac{\rho A L}{6}
\begin{bmatrix}
2&0&1&0\\
0&2&0&1\\
1&0&2&0\\
0&1&0&2
\end{bmatrix}
\]

Where:

- \(\rho\): Material density.
- \(A\): Cross-sectional area.
- \(L\): Element length.

The global mass matrix is assembled from the individual element mass matrices using the structural connectivity.

### 8.5 Boundary Conditions in Modal Analysis

Prescribed degrees of freedom are removed from the eigenvalue problem.

The reduced system is:

\[
\mathbf{K}_{ff}\boldsymbol{\phi}_f
=
\lambda\mathbf{M}_{ff}\boldsymbol{\phi}_f
\]

After solving, the complete mode shapes are reconstructed by inserting zeros at the constrained degrees of freedom.

### 8.6 Modal Analysis Implementation in MechFEM

The MechFEM modal analysis module provides:

- Generalized eigenvalue problem solution.
- Natural frequency calculation in Hz.
- Angular frequency calculation in rad/s.
- Mode shape extraction.
- Boundary condition handling.
- Optional selection of the number of modes.
- Consistent mass matrix assembly for 2D truss structures.

The implementation uses NumPy numerical routines and supports linear undamped modal analysis.

---

## 9. Boundary Conditions

Boundary conditions define the constraints imposed on the structural model.

For example, a fixed displacement condition is:

\[
u_i=0
\]

The global equilibrium equations are modified to account for prescribed degrees of freedom before solving for the unknown displacements.

Support reactions are subsequently recovered from the original equilibrium system:

\[
\mathbf{R}=\mathbf{K}\mathbf{u}-\mathbf{F}
\]

---

## 10. Numerical Solution

For a stable and properly constrained linear elastic model, the reduced system can be written as:

\[
\mathbf{K}_{ff}\mathbf{u}_f=\mathbf{F}_f
\]

Therefore:

\[
\mathbf{u}_f=\mathbf{K}_{ff}^{-1}\mathbf{F}_f
\]

In numerical implementations, direct linear system solvers are generally preferred over explicitly computing the inverse matrix.

---

## 11. Verification and Validation

Numerical results should be evaluated through:

- Comparison with analytical solutions.
- Equilibrium checks.
- Mesh refinement studies.
- Boundary condition verification.
- Dimensional consistency.
- Unit testing of individual element formulations.
- Modal frequency verification.

A successful automated test indicates that the tested conditions passed. It does not independently establish accuracy for every possible engineering configuration.

---

## 12. Current Scope and Future Development

The current theoretical scope of MechFEM includes:

- Linear Euler-Bernoulli beam elements.
- Two-dimensional linear truss elements.
- Two-dimensional linear frame elements.
- Linear elastic material behavior.
- Static structural analysis.
- Nodal displacement and reaction calculations.
- Axial force and stress recovery.
- Structural result visualization.
- Beam mesh convergence analysis.
- Undamped modal analysis for supported structural systems.
- Consistent mass matrix formulation for 2D truss elements.

Potential future developments include:

- Two-dimensional continuum elements.
- Plane stress and plane strain formulations.
- Plate and shell elements.
- Three-dimensional structural analysis.
- Nonlinear material models.
- Geometric nonlinearity.
- Advanced dynamic finite element analysis.
- Advanced numerical solvers.
- Engineering optimization.

---

## 13. References

1. Cook, R. D., Malkus, D. S., Plesha, M. E., and Witt, R. J., *Concepts and Applications of Finite Element Analysis*.
2. Bathe, K. J., *Finite Element Procedures*.
3. Logan, D. L., *A First Course in the Finite Element Method*.
4. Hutton, D. V., *Fundamentals of Finite Element Analysis*.
5. Clough, R. W., and Penzien, J., *Dynamics of Structures*.

---

**Project:** MechFEM  
**Repository:** https://github.com/revengeramh/MechFEM
