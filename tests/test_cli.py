import json

from shift.cli import main


def test_info_command(capsys) -> None:
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert "SHIFT" in out
    assert "thresholds are frozen from uncertainty" in out


def test_doctor_command(capsys) -> None:
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "Python:" in out


def test_poiseuille_benchmark_command_emits_machine_readable_result(capsys) -> None:
    assert main(["benchmark", "poiseuille"]) == 0
    payload = json.loads(capsys.readouterr().out)
    distances = payload["distances"]
    assert distances["d_tau"] == 0.0
    assert distances["d_q_invariant"] == 0.0
    assert distances["d_lu"] is None
    assert distances["d_m"] > 0.07
