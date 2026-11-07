"""Setup and command-line tests."""

from pathlib import Path

from ortam_core.cli import main


def test_list(capsys):
    assert main(["list"]) == 0
    assert "groundwater.three-point" in capsys.readouterr().out


def test_module_runs_standalone(capsys):
    example = Path(__file__).parents[2] / "modules/groundwater/examples/three_point.yaml"
    assert main(["run", "groundwater.three-point", str(example)]) == 0
    assert "Flow direction" in capsys.readouterr().out
