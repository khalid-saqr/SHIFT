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
