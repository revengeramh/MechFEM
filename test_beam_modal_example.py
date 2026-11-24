import runpy
from pathlib import Path


def test_beam_modal_example_runs():
    example_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "beam_modal_example.py"
    )

    example = runpy.run_path(str(example_path))
    example["main"]()
