
# MechFEM

**An Open-Source Python Framework for Computational Mechanics and Finite Element Analysis**

MechFEM is an open-source Python framework designed for computational mechanics, structural analysis, and Finite Element Method (FEM) applications.

The project focuses on transparent numerical implementations, modular architecture, reproducible simulations, and engineering verification.

## Current Features

### 1. Beam Analysis
- Euler-Bernoulli beam element formulation.
- Global stiffness matrix assembly.
- Boundary condition application.
- Nodal displacement calculation.
- Support reaction calculation.
- Analytical verification of beam responses.

### 2. Truss Analysis
- Two-dimensional truss element formulation.
- Global stiffness matrix assembly.
- Nodal displacement calculation.
- Support reaction calculation.
- Axial force and normal stress calculation.

### 3. Frame Analysis
- Two-dimensional frame element formulation.
- Six-degree-of-freedom element stiffness matrix.
- Global stiffness matrix assembly.
- Linear static structural analysis.
- Nodal displacement and reaction calculation.
- Frame deformation visualization.

### 4. Structural Analysis
- Element length calculation.
- Axial deformation analysis.
- Strain and stress calculation.
- Axial force calculation.
- Element strain energy calculation.
- Support for element-specific material properties.

### 5. Mesh Convergence Analysis
- Beam mesh refinement studies.
- Displacement comparison between mesh sizes.
- Absolute error calculation.
- Relative error calculation.
- Convergence rate estimation.
- Analytical reference comparison.

### 6. Modal Analysis
- Generalized eigenvalue problem solution.
- Natural frequency calculation.
- Angular frequency calculation.
- Mode shape extraction.
- Boundary condition support.
- Consistent mass matrix assembly for 2D truss structures.
- Optional extraction of a specified number of modes.

### 7. Visualization
- Beam deformation plots.
- Truss structural visualization.
- Axial force distribution.
- Axial stress distribution.
- Frame deformation plots.
- Truss modal shape visualization.
- Original and deformed mode shape comparison.
- Matplotlib-based visualization.

### 8. Automated Testing
- Unit testing with pytest.
- Analytical benchmark verification.
- Structural analysis tests.
- Modal analysis tests.
- Visualization tests.
- Mesh convergence tests.
- Continuous integration using GitHub Actions.

## Documentation

- [Theoretical Foundation and FEM Formulations](docs/theory.md)
- [Validation Documentation](VALIDATION.md)

## Technology Stack

- Python
- NumPy
- Matplotlib
- pytest
- GitHub Actions

## Installation

Clone the repository:

```bash
git clone https://github.com/revengeramh/MechFEM.git
cd MechFEM
```

Install the project:

```bash
pip install -e ".[test,visualization]"
```

## Usage

### Beam Analysis

```python
from mechfem import solve_beam
```

### Truss Analysis

```python
from mechfem import (
    solve_truss,
    calculate_truss_results
)
```

### Frame Analysis

```python
from mechfem import solve_frame
```

### Structural Analysis

```python
from mechfem import analyze_truss
```

### Mesh Convergence

```python
from mechfem.convergence import beam_mesh_convergence
```

### Modal Analysis

```python
from mechfem import (
    modal_analysis,
    analyze_truss_modal,
    assemble_truss_matrices
)
```

### Modal Visualization

```python
from mechfem import (
    analyze_truss_modal,
    plot_truss_mode_shape
)

results = analyze_truss_modal(
    nodes,
    elements,
    E,
    A,
    density,
    fixed_dofs
)

first_mode = results["mode_shapes"][:, 0]

plot_truss_mode_shape(
    nodes,
    elements,
    first_mode,
    scale=100.0
)
```

A complete working example is available at:

`examples/truss_modal_visualization_example.py`

Example scripts are available in the `examples/` directory.

## Project Structure

```text
MechFEM/
├── src/
│   └── mechfem/
│       ├── __init__.py
│       ├── beam.py
│       ├── truss.py
│       ├── frame.py
│       ├── visualization.py
│       ├── structural_analysis.py
│       ├── convergence.py
│       └── modal.py
├── tests/
│   ├── test_beam.py
│   ├── test_truss.py
│   ├── test_frame.py
│   ├── test_visualization.py
│   ├── test_frame_visualization.py
│   ├── test_structural_analysis.py
│   ├── test_convergence.py
│   ├── test_convergence_analysis.py
│   ├── test_modal.py
│   ├── test_truss_modal.py
│   ├── test_truss_modal_analysis.py
│   ├── test_truss_modal_example.py
│   └── test_truss_modal_visualization_example.py
├── examples/
│   ├── truss_modal_example.py
│   └── truss_modal_visualization_example.py
├── docs/
│   └── theory.md
├── .github/
│   └── workflows/
│       └── tests.yml
├── pyproject.toml
├── VALIDATION.md
└── README.md
```

## Running Tests

Run all tests:

```bash
pytest -v
```

GitHub Actions automatically runs the test suite when changes are pushed to the repository.

## Verification and Validation

MechFEM uses analytical solutions, numerical benchmarks, and automated tests to evaluate its implementations.

Validation activities include:

- Analytical benchmark comparisons.
- Structural equilibrium checks.
- Mesh refinement studies.
- Modal frequency verification.
- Numerical consistency checks.
- Automated regression testing.

Passing automated tests does not guarantee accuracy for every possible engineering configuration. Independent verification is required for engineering applications.

## Research and Development

MechFEM is being developed as a research-oriented computational mechanics framework.

Future development directions include:

- Advanced beam and plate elements.
- Computational solid mechanics.
- Material constitutive modeling.
- Dynamic finite element analysis.
- Nonlinear structural analysis.
- Mesh generation and post-processing.
- Advanced numerical solvers.
- Engineering optimization.

## Development Status

**Status: Active Development**

The project currently includes beam, truss, and frame analysis modules, structural post-processing, mesh convergence analysis, modal analysis, visualization tools, and automated tests.

Additional capabilities are planned as the framework evolves.

## Contribution

Contributions, suggestions, bug reports, and research collaborations are welcome.

## License

This project is distributed under the MIT License.

---

**Maintainer:** [revengeramh](https://github.com/revengeramh)

**Repository:** https://github.com/revengeramh/MechFEM
