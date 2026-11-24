import runpy
from pathlib import Path


def test_truss_modal_example():
    example_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "truss_modal_example.py"
    )

    example = runpy.run_path(str(example_path))

    example["main"]()
