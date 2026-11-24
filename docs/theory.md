
# Theoretical Foundation of MechFEM

## 1. Introduction

MechFEM is based on the Finite Element Method (FEM), a numerical technique used to approximate solutions of engineering problems governed by differential equations.

The main objective of FEM is to transform a continuous engineering system into a finite number of discrete elements connected through nodes.

This document describes the theoretical foundations of the current beam and two-dimensional truss implementations.

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

The general equilibrium equation is:

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

## 5. Boundary Conditions

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

## 6. Numerical Solution

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

## 7. Verification and Validation

Numerical results should be evaluated through:

- Comparison with analytical solutions.
- Equilibrium checks.
- Mesh refinement studies.
- Boundary condition verification.
- Dimensional consistency.
- Unit testing of individual element formulations.

A successful automated test indicates that the tested conditions passed. It does not independently establish accuracy for every possible engineering configuration.

---

## 8. Current Scope and Future Development

The current theoretical scope of MechFEM includes:

- Linear Euler-Bernoulli beam elements.
- Two-dimensional linear truss elements.
- Linear elastic material behavior.
- Static structural analysis.
- Nodal displacement and reaction calculations.
- Axial force and stress recovery.
- Structural result visualization.

Potential future developments include:

- Two-dimensional continuum elements.
- Plane stress and plane strain formulations.
- Plate and shell elements.
- Three-dimensional structural analysis.
- Nonlinear material models.
- Geometric nonlinearity.
- Dynamic finite element analysis.
- Advanced numerical solvers.

---

## 9. References

1. Cook, R. D., Malkus, D. S., Plesha, M. E., and Witt, R. J., *Concepts and Applications of Finite Element Analysis*.
2. Bathe, K. J., *Finite Element Procedures*.
3. Logan, D. L., *A First Course in the Finite Element Method*.
4. Hutton, D. V., *Fundamentals of Finite Element Analysis*.

---

**Project:** MechFEM  
**Repository:** https://github.com/revengeramh/MechFEM
