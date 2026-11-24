
# MechFEM

**An Open-Source Python Framework for Computational Mechanics and Finite Element Analysis**

MechFEM is an open-source Python project focused on developing modular and accessible computational tools for mechanical engineering, structural mechanics, and Finite Element Analysis (FEA).

The project aims to provide transparent numerical implementations, reproducible simulations, and engineering analysis tools with an emphasis on mathematical formulation, verification, and extensibility.

## Current Features

### 1. Beam Analysis

- Euler-Bernoulli beam element formulation.
- Global stiffness matrix assembly.
- Boundary condition application.
- Nodal displacement calculation.
- Support reaction calculation.
- Beam deformation visualization.

### 2. Truss Analysis

- Two-dimensional truss element stiffness formulation.
- Global structural stiffness matrix assembly.
- Nodal displacement calculation.
- Support reaction calculation.
- Axial force calculation.
- Normal stress calculation.

### 3. Visualization

MechFEM includes visualization tools based on Matplotlib:

- Beam deformation plots.
- 2D truss structural visualization.
- Axial force distribution plots.
- Axial stress distribution plots.
- Color-coded structural results.
- Saving visualization outputs as image files.

### 4. Automated Testing

- Unit tests for beam analysis.
- Unit tests for truss analysis.
- Visualization tests.
- Analytical verification tests.
- Mesh convergence testing.
- Automated testing through GitHub Actions.

## Documentation

The theoretical foundations and mathematical formulations of the implemented finite elements are documented in the following file:

- [Theoretical Foundation and FEM Formulations](docs/theory.md)

This documentation covers:

- Fundamental concepts of the Finite Element Method.
- Euler-Bernoulli beam element formulation.
- Two-dimensional truss element formulation.
- Element stiffness matrices.
- Boundary conditions.
- Numerical solution procedures.
- Verification and validation principles.
- Current project scope and future development directions.

## Technology Stack

- Python
- NumPy
- SciPy
- Matplotlib
- pytest

## Installation

Clone the repository:

```bash
git clone https://github.com/revengeramh/MechFEM.git
cd MechFEM
```

Install the project with testing and visualization dependencies:

```bash
pip install -e ".[test,visualization]"
```

## Usage

### Beam Analysis

Import the beam solver:

```python
from mechfem.beam import solve_beam
```

The solver supports structural analysis using nodes, elements, material properties, loads, and boundary conditions.

### Truss Analysis

Import the truss solver:

```python
from mechfem.truss import solve_truss
```

Additional functions are available for calculating axial forces and stresses:

```python
from mechfem.truss import calculate_truss_results
```

### Visualization

Import the visualization functions:

```python
from mechfem.visualization import (
    plot_beam_deformation,
    plot_truss_axial_forces,
    plot_truss_stresses
)
```

Example scripts are provided in the `examples/` directory.

Run the truss visualization example:

```bash
python examples/truss_visualization_example.py
```

## Project Structure

```text
MechFEM/
├── src/
│   └── mechfem/
│       ├── __init__.py
│       ├── beam.py
│       ├── truss.py
│       └── visualization.py
├── tests/
│   ├── test_beam.py
│   ├── test_truss.py
│   ├── test_convergence.py
│   └── test_visualization.py
├── examples/
│   ├── __init__.py
│   ├── truss_example.py
│   └── truss_visualization_example.py
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

Run the test suite:

```bash
pytest -v
```

The project also uses GitHub Actions for automated testing.

## Verification and Validation

MechFEM uses automated tests and analytical comparisons to evaluate its numerical implementations.

The validation documentation is available here:

[VALIDATION.md](VALIDATION.md)

Validation activities include:

- Analytical benchmark comparisons.
- Structural equilibrium checks.
- Mesh refinement studies.
- Automated unit testing.

Passing automated tests does not guarantee accuracy for every possible engineering configuration. Further verification is required as new elements and capabilities are introduced.

## Research and Development

MechFEM is being developed as a research-oriented computational mechanics project.

Future development directions include:

- One-dimensional and two-dimensional finite element formulations.
- Frame structural analysis.
- Advanced beam and plate elements.
- Computational solid mechanics.
- Material modeling.
- Numerical methods and solver improvements.
- Mesh generation and post-processing.
- Advanced engineering visualization.
- Dynamic and nonlinear finite element analysis.

## Development Status

**Status: Initial Development**

The beam and 2D truss analysis modules, visualization tools, example scripts, theoretical documentation, and automated tests are currently implemented.

Features listed under future development are planned and may not yet be available.

## Contribution

Contributions, suggestions, bug reports, and research collaborations are welcome.

## License

This project is distributed under the MIT License.

---

**Maintainer:** [revengeramh](https://github.com/revengeramh)

**Repository:** https://github.com/revengeramh/MechFEM
